"""
PSD-RMS Visualisation
=====================

Plots derived from the ``psd_rms`` workflow step: time-series, clock plots,
hour-maps, grid-maps, and daily stacks.  Ported and adapted from the
SeismoRMS / SeismoSocialDistancing project
(https://github.com/ThomasLecocq/SeismoRMS).

CLI usage::

    msnoise qc plot_psd_rms BE.UCC..HHZ --type timeseries --band 1.0-10.0
    msnoise qc plot_psd_rms BE.UCC..HHZ --type clockplot  --timezone Europe/Brussels
    msnoise qc plot_psd_rms BE.UCC..HHZ --type hourmap
    msnoise qc plot_psd_rms BE.UCC..HHZ --type gridmap
    msnoise qc plot_psd_rms BE.UCC..HHZ --type dailyplot

Notebook usage::

    from msnoise.core.db import connect
    from msnoise.results import MSNoiseResult
    from msnoise.plots.psd_rms import load_rms, plot_timeseries, plot_clockplot

    db = connect()
    r = MSNoiseResult.from_ids(db, psd=1, psd_rms=1)
    data = load_rms(r, ["BE.UCC..HHZ"])
    plot_timeseries(data, band="1.0-10.0")

"""

from __future__ import annotations

import datetime
import textwrap
from dataclasses import dataclass
from typing import Optional

import matplotlib
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

matplotlib.rcParams["pdf.fonttype"] = 42  # editable text in Illustrator/Inkscape

__all__ = [
    "DayNightWindow",
    "load_rms",
    "plot_timeseries",
    "plot_clockplot",
    "plot_hourmap",
    "plot_gridmap",
    "plot_dailyplot",
    "main",
    # public helpers for custom notebooks
    "localize_tz_and_reindex",
    "pivot_for_hourmap",
    "stack_wday_time",
    "radial_hours",
    "clock24_plot_commons",
]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_DAYS = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
_WRAPPER = textwrap.TextWrapper(width=15, break_long_words=False)


# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Day / Night window
# ---------------------------------------------------------------------------

@dataclass
class DayNightWindow:
    """Define the local "daytime" and "nighttime" hour windows.

    Used by all plot functions to:

    * draw the daytime median overlay on time-series plots,
    * shade night hours on clockplots, gridmaps, and dailyplots,
    * shade night sectors on hourmaps.

    Hours are expressed as floats in 0..24 in **local** time (after timezone
    conversion).  The window wraps midnight correctly, so day_start=22 with
    day_end=6 defines a "daytime" from 22:00 to 06:00 (e.g. a city station
    where the quiet hours are at night).

    Parameters
    ----------
    day_start:
        Start of daytime in decimal hours (default ``7.0`` = 07:00 local).
    day_end:
        End of daytime in decimal hours (default ``19.0`` = 19:00 local).
    day_color:
        Matplotlib colour for daytime shading (default light yellow).
    night_color:
        Matplotlib colour for nighttime shading (default light blue).
    day_alpha:
        Alpha for daytime shading (default ``0.15``).
    night_alpha:
        Alpha for nighttime shading (default ``0.12``).

    Examples
    --------
    >>> # Standard business-hours window
    >>> dnw = DayNightWindow(day_start=8, day_end=18)

    >>> # Seismological quiet window: anthropogenic noise lowest 22:00-05:00
    >>> dnw = DayNightWindow(day_start=5, day_end=22)

    >>> # Marine station: no meaningful day/night -- disable shading
    >>> dnw = DayNightWindow(day_start=0, day_end=24)
    """

    day_start: float = 7.0
    day_end: float = 19.0
    day_color: str = "#fffacd"    # lemon chiffon
    night_color: str = "#cce5ff"  # pale blue
    day_alpha: float = 0.15
    night_alpha: float = 0.12

    def __post_init__(self) -> None:
        if not (0.0 <= self.day_start < 24.0):
            raise ValueError(
                f"day_start must be in [0, 24), got {self.day_start}"
            )
        if not (0.0 < self.day_end <= 24.0):
            raise ValueError(
                f"day_end must be in (0, 24], got {self.day_end}"
            )

    @property
    def day_label(self) -> str:
        """Human-readable label, e.g. '7-19 h'."""
        return f"{self.day_start:.4g}–{self.day_end:.4g} h"

    def shade_cartesian_hour(self, ax) -> None:
        """Shade day / night horizontal bands on a Cartesian hour-of-day y-axis.

        Assumes the y-axis represents hours 0-24 (gridmap).
        """
        if self.day_start < self.day_end:
            if self.day_start > 0:
                ax.axhspan(0, self.day_start,
                           color=self.night_color, alpha=self.night_alpha,
                           zorder=-10, linewidth=0)
            if self.day_end < 24:
                ax.axhspan(self.day_end, 24,
                           color=self.night_color, alpha=self.night_alpha,
                           zorder=-10, linewidth=0)
            ax.axhspan(self.day_start, self.day_end,
                       color=self.day_color, alpha=self.day_alpha,
                       zorder=-10, linewidth=0)
        else:
            # Wraps midnight: e.g. day_start=22, day_end=6
            ax.axhspan(0, self.day_end,
                       color=self.day_color, alpha=self.day_alpha,
                       zorder=-10, linewidth=0)
            ax.axhspan(self.day_end, self.day_start,
                       color=self.night_color, alpha=self.night_alpha,
                       zorder=-10, linewidth=0)
            ax.axhspan(self.day_start, 24,
                       color=self.day_color, alpha=self.day_alpha,
                       zorder=-10, linewidth=0)

    def shade_cartesian_time(
        self,
        ax,
        x_min: datetime.datetime,
        x_max: datetime.datetime,
    ) -> None:
        """Shade night windows as vertical spans on a calendar time x-axis.

        Used by plot_timeseries to shade each nightly period.
        """
        day = x_min.date()
        while day <= x_max.date():
            midnight = datetime.datetime.combine(day, datetime.time(0))
            if self.day_start < self.day_end:
                # Two night blocks per calendar day
                if self.day_start > 0:
                    ax.axvspan(midnight,
                               midnight + datetime.timedelta(hours=self.day_start),
                               color=self.night_color, alpha=self.night_alpha,
                               zorder=-10, linewidth=0)
                if self.day_end < 24:
                    ax.axvspan(
                               midnight + datetime.timedelta(hours=self.day_end),
                               midnight + datetime.timedelta(hours=24),
                               color=self.night_color, alpha=self.night_alpha,
                               zorder=-10, linewidth=0)
            else:
                # Wraps midnight: night = [day_end, day_start] within the day
                ax.axvspan(midnight,
                           midnight + datetime.timedelta(hours=self.day_end),
                           color=self.night_color, alpha=self.night_alpha,
                           zorder=-10, linewidth=0)
                ax.axvspan(midnight + datetime.timedelta(hours=self.day_start),
                           midnight + datetime.timedelta(hours=24),
                           color=self.night_color, alpha=self.night_alpha,
                           zorder=-10, linewidth=0)
            day += datetime.timedelta(days=1)

    def shade_polar(self, ax, rmin: float, rmax: float) -> None:
        """Shade night sectors on a polar (clockplot / hourmap) axis.

        Parameters
        ----------
        ax:
            Polar Axes with clockwise direction and North at top
            (as set by clock24_plot_commons).
        rmin, rmax:
            Radial extent of the shaded wedge.
        """
        def _wedge(h_start: float, h_end: float, color: str, alpha: float) -> None:
            a0 = 2 * np.pi * h_start / 24.0
            a1 = 2 * np.pi * h_end / 24.0
            theta = np.linspace(a0, a1, 60)
            ax.fill_between(theta,
                            np.full_like(theta, rmin),
                            np.full_like(theta, rmax),
                            color=color, alpha=alpha, zorder=-5, linewidth=0)

        if self.day_start < self.day_end:
            if self.day_start > 0:
                _wedge(0, self.day_start, self.night_color, self.night_alpha)
            if self.day_end < 24:
                _wedge(self.day_end, 24, self.night_color, self.night_alpha)
        else:
            # Wraps midnight
            _wedge(self.day_end, self.day_start, self.night_color, self.night_alpha)


#: Default window: 07:00-19:00 daytime (mid-latitude urban stations)
DEFAULT_DAY_NIGHT = DayNightWindow()

# ---------------------------------------------------------------------------
# Public helper functions (importable for custom notebooks)
# ---------------------------------------------------------------------------

def localize_tz_and_reindex(
    df: pd.DataFrame,
    freq: str = "15min",
    time_zone: str = "UTC",
) -> pd.DataFrame:
    """Convert a UTC-indexed DataFrame to *time_zone* and resample to *freq*.

    Parameters
    ----------
    df:
        DataFrame with a timezone-naive UTC :class:`~pandas.DatetimeIndex`.
    freq:
        Pandas offset string for resampling (default ``"15min"``).
    time_zone:
        Target timezone string (default ``"UTC"``).

    Returns
    -------
    :class:`~pandas.DataFrame`
    """
    out = (
        df.copy()
        .tz_localize("UTC")
        .dropna()
        .tz_convert(time_zone)
        .tz_localize(None)
        .resample(freq)
        .mean()
    )
    return out if isinstance(out, pd.DataFrame) else out.to_frame()


def pivot_for_hourmap(
    data: pd.DataFrame,
    columns: str = "angles",
) -> pd.DataFrame:
    """Pivot a single-column RMS DataFrame into a (days x time) matrix.

    Parameters
    ----------
    data:
        DataFrame with a :class:`~pandas.DatetimeIndex` and a **single** RMS
        column (already scaled to the desired display unit).
    columns:
        ``"angles"`` (default) converts the time axis to radians 0-2pi for
        polar pcolormesh; ``"hours"`` leaves them as float hours 0-24 for the
        Cartesian grid-map.

    Returns
    -------
    :class:`~pandas.DataFrame`
        Rows = elapsed integer days since the first sample; columns = hour or
        angle values.
    """
    band = data.columns[0]
    df = data[[band]].copy()
    df["day"] = [d.year * 365 + d.dayofyear for d in df.index]
    df["time"] = [d.hour + d.minute / 60.0 for d in df.index]
    pivot = df.pivot(index="day", columns="time", values=band)
    pivot.index = (pivot.index - pivot.index[0]).astype(float)
    if columns == "angles":
        pivot.columns = 2 * np.pi * pivot.columns / 24.0
    return pivot


def stack_wday_time(
    df: pd.DataFrame,
    scale: float = 1.0,
) -> pd.DataFrame:
    """Median-stack a single-column RMS DataFrame into a (hour x weekday) table.

    Parameters
    ----------
    df:
        DataFrame with a :class:`~pandas.DatetimeIndex` and a single column.
        The column name is ignored — only the first column is used.
    scale:
        Multiply all values by this factor (default ``1.0``).

    Returns
    -------
    :class:`~pandas.DataFrame`
        Index = float hour-of-day; columns = weekday names ordered
        Monday to Sunday.
    """
    tmp = df.iloc[:, :1].copy()
    tmp.columns = ["rms"]
    tmp["wday"] = [t.day_name() for t in tmp.index]
    tmp["hour"] = [t.hour + t.minute / 60.0 for t in tmp.index]
    pivot = (
        tmp.groupby(["hour", "wday"])["rms"]
        .median()
        .unstack("wday")
        .reindex(columns=_DAYS)
    )
    return pivot * scale


def radial_hours(N: int) -> np.ndarray:
    """Return *N* evenly-spaced angles for a 24-hour clock plot (closed polygon).

    The last element equals the first so that ``DataFrame.plot()`` closes the
    polygon around the clock face.
    """
    hours = np.deg2rad(np.linspace(0, 360, N - 1, endpoint=False))
    return np.append(hours, hours[0])


def clock24_plot_commons(ax, unit: str = "nm") -> None:
    """Apply standard 24-hour polar axis decoration in-place."""
    from matplotlib.ticker import FixedLocator, FuncFormatter
    ax.set_xticks(np.linspace(0, 2 * np.pi, 24, endpoint=False))
    ax.set_xticklabels([f"{h:d} h" for h in range(24)], fontsize=8)
    # Fix radial ticks before labelling to avoid FixedLocator warning
    rticks = ax.get_yticks()
    ax.yaxis.set_major_locator(FixedLocator(rticks))
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda v, _: f"{v:.2g} {unit}")
    )
    ax.yaxis.set_tick_params(labelsize=7)
    ax.set_rlabel_position(0)
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2.0)
    plt.xlabel("Hour (local time)", fontsize=10)
    plt.grid(True)


# ---------------------------------------------------------------------------
# Private drawing helpers
# ---------------------------------------------------------------------------

def _draw_timeseries_annotations(ax, annotations: dict) -> None:
    """Draw vertical lines + legend entries for event annotations."""
    from obspy import UTCDateTime
    for i, (date_str, label) in enumerate(annotations.items()):
        ax.axvline(
            UTCDateTime(date_str).datetime,
            color=f"C{i}",
            linewidth=2,
            linestyle=["-", "--", "-.", ":"][i % 4],
            path_effects=[pe.withStroke(linewidth=4, foreground="k")],
            zorder=-9,
            label="\n".join(_WRAPPER.wrap(label)),
        )


def _draw_polar_annotations(
    ax, annotations: dict, origin: datetime.datetime
) -> None:
    """Draw radial spokes + markers on a polar hour-map for event annotations."""
    from obspy import UTCDateTime
    rticks = []
    for i, (date_str, label) in enumerate(annotations.items()):
        t = UTCDateTime(date_str).datetime
        r = (t - origin).days
        if r <= 0:
            continue
        angle = (t.hour / 24.0 + t.minute / 60.0 / 24.0) * 2 * np.pi
        rticks.append(r)
        ax.plot(
            angle, r, "o",
            label="\n".join(_WRAPPER.wrap(label)),
            color=f"C{i}",
            path_effects=[
                pe.withStroke(linewidth=5, foreground="w"),
                pe.withStroke(linewidth=3, foreground="k"),
            ],
        )
    if rticks:
        ax.set_rticks(rticks)


def _clockpanel(
    ax, s_rms: pd.DataFrame, rmax: float, unit: str, title: str
) -> None:
    """Draw one polar clock panel from a single-column 'rms' DataFrame.

    Builds the weekday/hour median stack, closes the polygon, and plots it.
    """
    stacked = stack_wday_time(s_rms, scale=1.0)
    # Close the polygon: append first row at the end with a new numeric index
    closed = pd.concat([stacked, stacked.iloc[[0]]], ignore_index=True)
    closed.index = radial_hours(len(closed))
    closed.plot(ax=ax)
    ax.set_title(title, fontsize=12)
    clock24_plot_commons(ax, unit=unit)
    ax.set_rmax(rmax)
    ax.set_rmin(0)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_rms(
    result,
    seed_ids: Optional[list[str]] = None,
    *,
    psd_id: int = 1,
    psd_rms_id: int = 1,
) -> dict[str, pd.DataFrame]:
    """Load PSD-RMS results via an :class:`~msnoise.results.MSNoiseResult`.

    This is the preferred entry point for notebooks and scripts.  It delegates
    all path resolution to ``MSNoiseResult.get_psd_rms()``, which correctly
    handles the ``psd_N/psd_rms_N/_output/`` lineage without hard-coding step
    names or output folders.

    Parameters
    ----------
    result:
        An :class:`~msnoise.results.MSNoiseResult` that contains ``psd_rms``
        in its lineage, obtained via::

            from msnoise.results import MSNoiseResult
            r = MSNoiseResult.from_ids(db, psd=1, psd_rms=1)

    seed_ids:
        Optional list of SEED IDs to load (e.g. ``["BE.UCC..HHZ"]``).
        When ``None`` (default) every station found on disk is returned.

    Returns
    -------
    dict
        Mapping *seed_id* -> :class:`~pandas.DataFrame` (columns = band labels,
        index = :class:`~pandas.DatetimeIndex` in UTC).  Stations with no
        data on disk are silently omitted.

    Examples
    --------
    **Preferred -- MSNoiseResult:**

    >>> from msnoise.core.db import connect
    >>> from msnoise.results import MSNoiseResult
    >>> from msnoise.plots.psd_rms import load_rms
    >>> db = connect()
    >>> r = MSNoiseResult.from_ids(db, psd=1, psd_rms=1)
    >>> data = load_rms(r)                   # all stations auto-discovered
    >>> data = load_rms(r, ["BE.UCC..HHZ"])  # single station

"""

    if seed_ids is not None:
        out: dict[str, pd.DataFrame] = {}
        for sid in seed_ids:
            df = result.get_psd_rms(seed_id=sid, format="dataframe")
            if df is not None and not (isinstance(df, pd.DataFrame) and df.empty):
                out[sid] = df
        return out

    # No seed_ids -> auto-discover all stations on disk
    raw = result.get_psd_rms(format="dataframe")
    return {sid: df for sid, df in raw.items() if df is not None}


# ---------------------------------------------------------------------------
# Plot functions
# ---------------------------------------------------------------------------

def plot_timeseries(
    data: dict[str, pd.DataFrame],
    band: Optional[str] = None,
    scale: float = 1e9,
    unit: str = "nm",
    annotations: Optional[dict] = None,
    time_zone: str = "UTC",
    resample_freq: str = "30min",
    day_night: Optional[DayNightWindow] = None,
    logo: Optional[str] = None,
    show: bool = True,
    outfile: Optional[str] = None,
) -> plt.Figure:
    """Time-series plot of RMS amplitudes for one or multiple stations.

    Parameters
    ----------
    data:
        dict ``{seed_id: DataFrame}`` as returned by :func:`load_rms`.
    band:
        Frequency band label (column name in the DataFrame, e.g. ``"1.0-10.0"``).
        Defaults to the first available band.
    scale:
        Multiply amplitudes before display (default ``1e9`` -> nm for DISP).
    unit:
        Y-axis unit label (default ``"nm"``).
    annotations:
        Optional dict ``{date_string: label}`` for vertical event lines.
    time_zone:
        Time zone for x-axis localisation (default ``"UTC"``).
    resample_freq:
        Resampling frequency string (default ``"30min"``).
    day_night:
        :class:`DayNightWindow` defining local daytime/nighttime hours.
        Controls: (1) the daytime median overlay time window and label,
        and (2) night shading applied as vertical spans.
        Defaults to :data:`DEFAULT_DAY_NIGHT` (07:00-19:00).
        Pass ``DayNightWindow(day_start=0, day_end=24)`` to disable shading.
    logo:
        URL or path to a logo image to embed bottom-left (optional).
    show:
        Call ``plt.show()`` when done.
    outfile:
        Save figure to this path.

    Returns
    -------
    :class:`matplotlib.figure.Figure`
    """
    if not data:
        raise ValueError("data dict is empty")

    dnw = day_night if day_night is not None else DEFAULT_DAY_NIGHT
    # Build "HH:MM" strings for between_time
    h0 = int(dnw.day_start)
    m0 = int(round((dnw.day_start - h0) * 60))
    h1 = int(dnw.day_end) % 24
    m1 = int(round((dnw.day_end - int(dnw.day_end)) * 60))
    t_day_start = f"{h0:02d}:{m0:02d}"
    t_day_end   = f"{h1:02d}:{m1:02d}" if h1 > 0 else "23:59"

    fig, ax = plt.subplots(figsize=(12, 5))
    if logo is not None:
        try:
            img = plt.imread(logo)
            fig.figimage(img, 40, 40, alpha=0.4, zorder=1)
        except Exception:
            pass

    all_series: list[pd.DataFrame] = []
    for sid, df in data.items():
        b = band or df.columns[0]
        s = localize_tz_and_reindex(df[[b]], resample_freq, time_zone) * scale
        all_series.append(s)
        ax.plot(s.index, s.iloc[:, 0], label=sid)
        # Daytime median overlay (user-defined window)
        daytime = s.between_time(t_day_start, t_day_end).resample("1D").median()
        midday = (dnw.day_start + dnw.day_end) / 2.0
        daytime.index = daytime.index + pd.Timedelta(hours=midday)
        ax.plot(daytime.index, daytime.iloc[:, 0], "--",
                label=f"{sid} ({dnw.day_label} median)", linewidth=1.2)

    # Night shading (vertical spans)
    all_idx = pd.concat(all_series).index
    if len(all_idx):
        dnw.shade_cartesian_time(ax, all_idx.min().to_pydatetime(),
                                 all_idx.max().to_pydatetime())

    if annotations:
        _draw_timeseries_annotations(ax, annotations)

    all_vals = pd.concat(all_series).iloc[:, 0].dropna()
    if len(all_vals):
        ax.set_ylim(0, float(all_vals.quantile(0.95)) * 1.5)

    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.3g}"))
    ax.set_ylabel(f"Amplitude ({unit})")
    b_label = band or next(iter(data.values())).columns[0]
    ax.set_title(f"Seismic Noise RMS -- band [{b_label}] Hz")
    ax.grid(True, zorder=-1)
    ax.set_axisbelow(True)
    fig.autofmt_xdate()
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.tight_layout()

    if outfile:
        fig.savefig(outfile, bbox_inches="tight", facecolor="w")
    if show:
        plt.show()
    return fig


def plot_clockplot(
    data: dict[str, pd.DataFrame],
    band: Optional[str] = None,
    scale: float = 1e9,
    unit: str = "nm",
    split_date: Optional[str] = None,
    time_zone: str = "UTC",
    resample_freq: str = "30min",
    day_night: Optional[DayNightWindow] = None,
    show: bool = True,
    outfile: Optional[str] = None,
) -> plt.Figure:
    """Polar 24-hour clock plot of median RMS per weekday / hour.

    When *split_date* is provided two panels are drawn (before / after),
    useful for comparing noise regimes across an event or intervention.

    Parameters
    ----------
    data:
        dict ``{seed_id: DataFrame}`` as returned by :func:`load_rms`.
        Only the first station is used; for multi-station comparisons build
        a side-by-side figure using :func:`stack_wday_time` directly (see the
        multi-station notebook example).
    band:
        Frequency band label.
    scale:
        Amplitude scale factor (default ``1e9`` -> nm).
    unit:
        Unit label for radial tick labels.
    split_date:
        ISO date string (``"YYYY-MM-DD"``) to split the data.  When ``None``
        a single panel is drawn covering the full time range.
    time_zone:
        Local time zone for hour-of-day grouping (default ``"UTC"``).
    resample_freq:
        Resampling frequency string (default ``"30min"``).
    day_night:
        :class:`DayNightWindow` for shading night sectors on the polar axes.
        Defaults to :data:`DEFAULT_DAY_NIGHT` (07:00-19:00).
    show:
        Call ``plt.show()``.
    outfile:
        Save figure to this path.

    Returns
    -------
    :class:`matplotlib.figure.Figure`
    """
    if not data:
        raise ValueError("data dict is empty")

    sid = next(iter(data))
    df = data[sid]
    b = band or df.columns[0]

    # Localise, resample, scale into a single-column 'rms' DataFrame
    s_full = localize_tz_and_reindex(df[[b]], resample_freq, time_zone) * scale
    s_full.columns = ["rms"]

    # Derive rmax from the full series so both panels share the same scale
    rmax = float(s_full["rms"].quantile(0.95)) * 1.5

    dnw = day_night if day_night is not None else DEFAULT_DAY_NIGHT
    n_panels = 2 if split_date else 1
    fig = plt.figure(figsize=(6 * n_panels, 6))

    if split_date:
        ax1 = fig.add_subplot(121, polar=True)
        _clockpanel(ax1, s_full.loc[:split_date], rmax, unit,
                    f"Before {split_date}")
        dnw.shade_polar(ax1, rmin=0, rmax=rmax)

        ax2 = fig.add_subplot(122, polar=True, sharey=ax1)
        post = s_full.loc[split_date:]
        if len(post):
            _clockpanel(ax2, post, rmax, unit, f"After {split_date}")
        else:
            ax2.set_title(f"After {split_date} -- no data", fontsize=12)
        dnw.shade_polar(ax2, rmin=0, rmax=rmax)
    else:
        ax1 = fig.add_subplot(111, polar=True)
        _clockpanel(ax1, s_full, rmax, unit, f"{sid}  [{b}] Hz")
        dnw.shade_polar(ax1, rmin=0, rmax=rmax)

    plt.suptitle(f"Day/Hour Median Noise -- {sid}  [{b}] Hz", fontsize=14)
    plt.subplots_adjust(top=0.85)

    if outfile:
        fig.savefig(outfile, bbox_inches="tight", facecolor="w")
    if show:
        plt.show()
    return fig


def plot_hourmap(
    data: dict[str, pd.DataFrame],
    band: Optional[str] = None,
    scale: float = 1e9,
    unit: str = "nm",
    annotations: Optional[dict] = None,
    time_zone: str = "UTC",
    resample_freq: str = "30min",
    day_night: Optional[DayNightWindow] = None,
    show: bool = True,
    outfile: Optional[str] = None,
) -> plt.Figure:
    """Polar hour-map: radial axis = elapsed days, angular axis = hour of day.

    Equivalent to SeismoRMS ``clockmap`` / ``hourmap``.  Best for datasets
    up to ~1-2 years; beyond that the polar readability degrades and
    :func:`plot_gridmap` is preferable.

    Parameters
    ----------
    data:
        dict ``{seed_id: DataFrame}`` as returned by :func:`load_rms`.
    band:
        Frequency band label.
    scale / unit:
        Amplitude scale and unit label.
    annotations:
        Optional dict ``{date_str: label}`` to mark events as radial markers.
    time_zone:
        Local time zone.
    resample_freq:
        Resampling frequency.
    day_night:
        :class:`DayNightWindow` for shading night sectors on the polar axes.
        Defaults to :data:`DEFAULT_DAY_NIGHT` (07:00-19:00).
    show / outfile:
        Display / save options.

    Returns
    -------
    :class:`matplotlib.figure.Figure`
    """
    if not data:
        raise ValueError("data dict is empty")

    sid = next(iter(data))
    df = data[sid]
    b = band or df.columns[0]
    s = localize_tz_and_reindex(df[[b]], resample_freq, time_zone) * scale
    origin_time = s.index[0]
    origin_text = origin_time.strftime("%Y-%m-%d")

    vmin = float(s.iloc[:, 0].quantile(0.01))
    vmax = float(s.iloc[:, 0].quantile(0.95))
    pivoted = pivot_for_hourmap(s, columns="angles")

    fig = plt.figure(figsize=(7, 9))
    ax = fig.add_subplot(111, projection="polar")
    ax.set_xticks(np.linspace(0, np.pi * 2 * 23 / 24, 24))
    ax.set_xticklabels([f"{h:d} h" for h in range(24)])
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    X = np.append(pivoted.columns.values, 2 * np.pi)
    Y = np.append(pivoted.index.values, pivoted.index.values[-1] + 1)
    pm = ax.pcolormesh(X, Y, pivoted.values,
                       vmin=vmin, vmax=vmax, rasterized=True, antialiased=True)
    cb = plt.colorbar(pm, ax=ax, orientation="horizontal", shrink=0.8)
    cb.ax.set_xlabel(f"Amplitude ({unit})")
    ax.set_rorigin(max(Y) / -4)
    ax.text(np.pi, max(Y) / -4, origin_text, ha="center", va="center")
    ax.set_rmax(max(Y))
    ax.grid(color="w")

    # Night sector shading on the polar axes
    dnw = day_night if day_night is not None else DEFAULT_DAY_NIGHT
    dnw.shade_polar(ax, rmin=0, rmax=max(Y))

    if annotations:
        _draw_polar_annotations(ax, annotations, origin_time.to_pydatetime())
        plt.legend(loc="lower left", bbox_to_anchor=(0.0, -0.2),
                   ncol=2, borderaxespad=0, frameon=False)

    ax.set_title(
        f"Seismic Noise Hour-Map -- {sid}  [{b}] Hz\n{origin_text}", pad=15
    )

    if outfile:
        fig.savefig(outfile, bbox_inches="tight", facecolor="w")
    if show:
        plt.show()
    return fig


def plot_gridmap(
    data: dict[str, pd.DataFrame],
    band: Optional[str] = None,
    scale: float = 1e9,
    unit: str = "nm",
    annotations: Optional[dict] = None,
    time_zone: str = "UTC",
    resample_freq: str = "30min",
    day_night: Optional[DayNightWindow] = None,
    show: bool = True,
    outfile: Optional[str] = None,
) -> plt.Figure:
    """Cartesian grid-map: x = calendar date, y = hour of day.

    Equivalent to SeismoRMS ``gridmap``.  Handles long time series (years)
    better than :func:`plot_hourmap`.

    Parameters
    ----------
    data:
        dict ``{seed_id: DataFrame}`` as returned by :func:`load_rms`.
    band:
        Frequency band label.
    scale / unit:
        Amplitude scale and unit label.
    annotations:
        Optional dict ``{date_str: label}`` for event markers.
    time_zone:
        Local time zone.
    resample_freq:
        Resampling frequency.
    day_night:
        :class:`DayNightWindow` for shading night bands on the hour-of-day
        y-axis.  Defaults to :data:`DEFAULT_DAY_NIGHT` (07:00-19:00).
    show / outfile:
        Display / save options.

    Returns
    -------
    :class:`matplotlib.figure.Figure`
    """
    if not data:
        raise ValueError("data dict is empty")

    sid = next(iter(data))
    df = data[sid]
    b = band or df.columns[0]
    s = localize_tz_and_reindex(df[[b]], resample_freq, time_zone) * scale
    origin_text = s.index[0].strftime("%Y-%m-%d")

    vmin = float(s.iloc[:, 0].quantile(0.01))
    vmax = float(s.iloc[:, 0].quantile(0.95))
    pivoted = pivot_for_hourmap(s, columns="hours")

    fig, ax = plt.subplots(figsize=(16, 5))
    X = pd.date_range(origin_text, periods=len(pivoted) + 1).to_pydatetime()
    Y = np.append(pivoted.columns.values, 24.0)
    pm = ax.pcolormesh(X, Y, pivoted.values.T,
                       vmin=vmin, vmax=vmax, rasterized=True, antialiased=True)
    plt.colorbar(pm, ax=ax, shrink=0.7, pad=0.01).set_label(f"Amplitude ({unit})")
    ax.set_xticks(pd.date_range(X[0], X[-1], freq="W-MON").to_pydatetime())
    ax.set_yticks(np.arange(25))
    ax.set_yticklabels([f"{h:d} h" for h in range(25)])
    ax.grid(True, which="both", c="k")

    # Night band shading (horizontal bands on the hour-of-day y-axis)
    dnw = day_night if day_night is not None else DEFAULT_DAY_NIGHT
    dnw.shade_cartesian_hour(ax)

    if annotations:
        from obspy import UTCDateTime
        for i, (date_str, label) in enumerate(annotations.items()):
            t = UTCDateTime(date_str).datetime
            ax.plot(t, t.hour + t.minute / 60.0, "o",
                    label="\n".join(_WRAPPER.wrap(label)),
                    color=f"C{i}",
                    path_effects=[
                        pe.withStroke(linewidth=5, foreground="w"),
                        pe.withStroke(linewidth=3, foreground="k"),
                    ])
        plt.legend(loc="lower left", bbox_to_anchor=(0.0, -0.25),
                   ncol=2, borderaxespad=0, frameon=False)

    ax.set_title(f"Seismic Noise Grid-Map -- {sid}  [{b}] Hz")
    fig.autofmt_xdate()
    plt.tight_layout()

    if outfile:
        fig.savefig(outfile, bbox_inches="tight", facecolor="w")
    if show:
        plt.show()
    return fig


def plot_dailyplot(
    data: dict[str, pd.DataFrame],
    band: Optional[str] = None,
    scale: float = 1e9,
    unit: str = "nm",
    split_date: Optional[str] = None,
    time_zone: str = "UTC",
    resample_freq: str = "30min",
    day_night: Optional[DayNightWindow] = None,
    show: bool = True,
    outfile: Optional[str] = None,
) -> plt.Figure:
    """Hour-of-day stacked plot, one curve per weekday.

    Solid curves cover the full dataset (or the period before *split_date*
    when provided); dashed curves represent the period after *split_date*.

    Parameters
    ----------
    data:
        dict ``{seed_id: DataFrame}`` as returned by :func:`load_rms`.
    band:
        Frequency band label.
    scale / unit:
        Amplitude scale and unit label.
    split_date:
        Optional date string to draw a pre/post comparison (solid vs dashed).
    time_zone:
        Local time zone for hour-of-day grouping.
    resample_freq:
        Resampling frequency.
    day_night:
        :class:`DayNightWindow` for shading night bands on the hour-of-day
        x-axis.  Defaults to :data:`DEFAULT_DAY_NIGHT` (07:00-19:00).
    show / outfile:
        Display / save options.

    Returns
    -------
    :class:`matplotlib.figure.Figure`
    """
    if not data:
        raise ValueError("data dict is empty")

    sid = next(iter(data))
    df = data[sid]
    b = band or df.columns[0]

    s = localize_tz_and_reindex(df[[b]], resample_freq, time_zone) * scale
    s.columns = ["rms"]

    fig, ax = plt.subplots(figsize=(14, 6))
    cmap = plt.get_cmap("tab20")

    pre = s if split_date is None else s.loc[:split_date]
    stack_wday_time(pre).plot(ax=ax, cmap=cmap, linewidth=1.5)

    if split_date:
        post = s.loc[split_date:]
        if len(post):
            stack_wday_time(post).plot(ax=ax, ls="--", cmap=cmap,
                                       legend=False, linewidth=1.5)

    # Night band shading (vertical spans on the hour x-axis)
    dnw = day_night if day_night is not None else DEFAULT_DAY_NIGHT
    if dnw.day_start < dnw.day_end:
        if dnw.day_start > 0:
            ax.axvspan(0, dnw.day_start,
                       color=dnw.night_color, alpha=dnw.night_alpha,
                       zorder=-10, linewidth=0)
        if dnw.day_end < 24:
            ax.axvspan(dnw.day_end, 24,
                       color=dnw.night_color, alpha=dnw.night_alpha,
                       zorder=-10, linewidth=0)
        ax.axvspan(dnw.day_start, dnw.day_end,
                   color=dnw.day_color, alpha=dnw.day_alpha,
                   zorder=-10, linewidth=0)
    else:
        # Wraps midnight
        ax.axvspan(0, dnw.day_end,
                   color=dnw.day_color, alpha=dnw.day_alpha,
                   zorder=-10, linewidth=0)
        ax.axvspan(dnw.day_end, dnw.day_start,
                   color=dnw.night_color, alpha=dnw.night_alpha,
                   zorder=-10, linewidth=0)
        ax.axvspan(dnw.day_start, 24,
                   color=dnw.day_color, alpha=dnw.day_alpha,
                   zorder=-10, linewidth=0)

    ax.set_title(f"Daily Noise Levels -- {sid}  [{b}] Hz")
    ax.set_xlabel("Hour of day (local time)")
    ax.set_ylabel(f"Amplitude ({unit})")
    ax.set_xlim(0, 24)
    ax.set_ylim(0, float(s["rms"].quantile(0.95)) * 1.5)
    ax.grid(True)
    plt.tight_layout()

    if outfile:
        fig.savefig(outfile, bbox_inches="tight", facecolor="w")
    if show:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
# Unified entry point (CLI + script)
# ---------------------------------------------------------------------------

def main(
    seed_ids: list[str],
    psd_id: int = 1,
    psd_rms_id: int = 1,
    plot_type: str = "timeseries",
    band: Optional[str] = None,
    scale: float = 1e9,
    unit: str = "nm",
    time_zone: str = "UTC",
    day_start: float = 7.0,
    day_end: float = 19.0,
    split_date: Optional[str] = None,
    annotations: Optional[dict] = None,
    outfile: Optional[str] = None,
    show: bool = True,
    loglevel: str = "INFO",
) -> None:
    """CLI/script entry point: load RMS data and dispatch to plot function(s).

    Constructs an :class:`~msnoise.results.MSNoiseResult` from the database,
    loads data via :func:`load_rms`, then dispatches to the requested plot
    function(s).

    Parameters
    ----------
    seed_ids:
        List of SEED IDs (e.g. ``["BE.UCC..HHZ"]``).  Empty list -> all
        stations found on disk for this lineage.
    psd_id:
        Config-set number for the ``psd`` step (default ``1``).
    psd_rms_id:
        Config-set number for the ``psd_rms`` step (default ``1``).
    plot_type:
        One of ``"timeseries"``, ``"clockplot"``, ``"hourmap"``,
        ``"gridmap"``, ``"dailyplot"``, or ``"all"``.
    band:
        Frequency band column label (e.g. ``"1.0-10.0"``).
    scale / unit:
        Amplitude scale and label.
    time_zone:
        Time zone string for local-time plots.
    day_start:
        Start of daytime in decimal hours (default ``7.0``).
    day_end:
        End of daytime in decimal hours (default ``19.0``).
    split_date:
        Date string for before/after split (clockplot, dailyplot).
    annotations:
        Dict ``{date_str: label}`` for event markers (timeseries, hourmap,
        gridmap).
    outfile:
        Base filename for saved figures; a ``_<type>`` suffix is inserted
        before the extension for each type when ``plot_type="all"``.
    show:
        Call ``plt.show()``.
    loglevel:
        Logging level string.
    """
    from ..core.db import connect, get_logger

    logger = get_logger("msnoise.plots.psd_rms", loglevel)
    db = connect()
    from ..results import MSNoiseResult
    result = MSNoiseResult.from_ids(db, psd=psd_id, psd_rms=psd_rms_id)

    data = load_rms(result, seed_ids if seed_ids else None)
    db.close()

    if not data:
        logger.error(
            "No RMS data found for psd=%d psd_rms=%d %s -- "
            "have you run `msnoise qc compute_psd_rms`?",
            psd_id, psd_rms_id, seed_ids or "(all stations)",
        )
        return

    dnw = DayNightWindow(day_start=day_start, day_end=day_end)

    def _outfile(suffix: str) -> Optional[str]:
        if not outfile:
            return None
        base, _, ext = outfile.rpartition(".")
        return f"{base}_{suffix}.{ext}" if ext else f"{outfile}_{suffix}.png"

    _dispatch = {
        "timeseries": lambda: plot_timeseries(
            data, band=band, scale=scale, unit=unit,
            annotations=annotations, time_zone=time_zone,
            day_night=dnw, show=show, outfile=_outfile("timeseries"),
        ),
        "clockplot": lambda: plot_clockplot(
            data, band=band, scale=scale, unit=unit,
            split_date=split_date, time_zone=time_zone,
            day_night=dnw, show=show, outfile=_outfile("clockplot"),
        ),
        "hourmap": lambda: plot_hourmap(
            data, band=band, scale=scale, unit=unit,
            annotations=annotations, time_zone=time_zone,
            day_night=dnw, show=show, outfile=_outfile("hourmap"),
        ),
        "gridmap": lambda: plot_gridmap(
            data, band=band, scale=scale, unit=unit,
            annotations=annotations, time_zone=time_zone,
            day_night=dnw, show=show, outfile=_outfile("gridmap"),
        ),
        "dailyplot": lambda: plot_dailyplot(
            data, band=band, scale=scale, unit=unit,
            split_date=split_date, time_zone=time_zone,
            day_night=dnw, show=show, outfile=_outfile("dailyplot"),
        ),
    }

    types = (
        list(_dispatch.keys())
        if plot_type in ("all", "*")
        else [plot_type]
    )

    for pt in types:
        if pt not in _dispatch:
            logger.warning("Unknown plot type '%s' -- skipping", pt)
            continue
        logger.info("Plotting %s", pt)
        _dispatch[pt]()
