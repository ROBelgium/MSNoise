"""
This plot shows the final output of MSNoise using the wavelet.


Example:

``msnoise cc dvv plot wct``

"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import matplotlib.dates as mdates
import pandas as pd
from ..api import *
from datetime import datetime, timedelta

def plot_dvv_heatmap(data_type, dvv_df, pair, rolling, start, end, low, high, logger, mincoh=0.5):
    # Extracting relevant data from dvv_df
    dvv_df = dvv_df.loc[start:end]
    if dvv_df is None or dvv_df.empty:
        logger.error(f"No data available for {pair} between {start} and {end}. Exiting function.")
        return None, None
    rolling_window = int(rolling)
    
    dvv_freq = dvv_df['dvv']
    coh_freq = dvv_df['coh']

    dvv_freq = dvv_freq.rolling(window=rolling_window, min_periods=1).mean()
    coh_freq = coh_freq.rolling(window=rolling_window, min_periods=1).mean()

    fig, ax = plt.subplots(figsize=(16, 10))   
    # Scatter plot of dv/v data
    #norm1 = plt.Normalize(vmin=np.min(dvv_freq.T), vmax=np.max(dvv_freq.T))

    if data_type == 'dvv':
        low_per = np.nanpercentile(dvv_freq, 1)
        high_per = np.nanpercentile(dvv_freq, 99)
        
        ax.pcolormesh(np.asarray(dvv_freq.index), np.asarray(dvv_freq.columns), dvv_freq.T,
            cmap=mpl.cm.seismic, edgecolors='none', vmin=low_per, vmax=high_per)
        save_name = f"{pair[0]}_{low}_{high}_Hz_m{rolling_window}_dvv_heatmap"
        color_bar_label = 'dv/v (%)'
        
    elif data_type == 'coh':
        ax.pcolormesh(np.asarray(coh_freq.index), np.asarray(coh_freq.columns), coh_freq.T,
                cmap='RdYlGn', edgecolors='none', vmin=mincoh, vmax=1)
        save_name = f"{pair[0]}_{low}_{high}_Hz_m{rolling_window}_coh_heatmap"
        color_bar_label = 'Coherence value'
    else:
        logger.error("Unknown data type: %s, write 'dvv' or 'coh'? " % data_type)
        return None, None

    #if current_config.get('plot_event', False):
    #    plot_events(ax, current_config['event_list'], start, end)
                
    ax.set_xlim(pd.to_datetime(start), pd.to_datetime(end))
        
    ax.set_ylabel('Frequency (Hz)', fontsize=18)
    ax.set_title(save_name, fontsize=22)

    cbar1 = plt.colorbar(ax.collections[0], ax=ax, pad=0.02)
    cbar1.set_label(color_bar_label, fontsize=18)
    #norm1 = Normalize(vmin=0, vmax=1)

    ax.tick_params(axis='both', which='both', labelsize=16, width=2, length=5)
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    fig.autofmt_xdate()
    
    # Adjust the layout if necessary
    fig.subplots_adjust(right=0.85)
    fig.tight_layout()
    return fig, save_name

def plot_dvv_scatter(dvv_df, pair, rolling, start, end, ranges, logger):
    # Extracting relevant data from dvv_df
    dvv_df = dvv_df.loc[start:end]
    if dvv_df is None or dvv_df.empty:
        logger.error(f"No data available for {pair} between {start} and {end}. Exiting function.")
        return None, None
    rolling_window = int(rolling)

    color = ['Blues', 'Reds','Greens','Greys'] #'Purples'
    color2 = ['blue', 'red', 'green', 'grey']  # Colors for different frequency ranges
    freq_names = []
    legend_handles = []
    ranges_list  = [list(map(float, r.strip().strip('[]').split(','))) for r in ranges.split('], [')]

    fig, ax = plt.subplots(figsize=(16, 10))
    # Loop through the frequency ranges specified in current_config
    for i, freqrange in enumerate(ranges_list):#[[0.5, 1.0], [1.0, 2.0], [2.0, 4.0]]):#, [0.5, 2.0]]):
        freq_name = f"{freqrange[0]}-{freqrange[1]} Hz"
        freq_names.append(freq_name)

        freqs = np.asarray(dvv_df['dvv'].columns)
        filtered_freqs = freqs[(freqs >= freqrange[0]) & (freqs <= freqrange[1])].tolist()
        dvv_freq = dvv_df['dvv'][filtered_freqs]        
        coh_freq = dvv_df['coh'][filtered_freqs]

        dvv_freq = dvv_freq.rolling(window=rolling_window, min_periods=1).mean()
        coh_freq = coh_freq.rolling(window=rolling_window, min_periods=1).mean()

        # Scatter plot of dv/v data
        norm1 = plt.Normalize(vmin=0, vmax=1)
        ax.scatter([0,1], [0,1], c=[0,1], cmap=color[-1])

        sc = ax.scatter(dvv_freq.index, dvv_freq.mean(axis=1), c=coh_freq.mean(axis=1), cmap=color[i], norm=norm1, label=freq_name)
        legend_handles.append(Line2D([0], [0], marker='o', color=color2[i], markerfacecolor=color2[i], markersize=10, label=freq_name))

    #if current_config.get('plot_event', False):
    #    plot_events(ax, current_config['event_list'], start, end)
                
    ax.set_xlim(pd.to_datetime(start), pd.to_datetime(end))
    #if current_config.get('same_dvv_scale', False):
    #    ax.set_ylim(current_config['dvv_min'], current_config['dvv_max'])
        
    ax.set_ylabel('dv/v (%)', fontsize=18)
    ax.set_title(f"{pair[0]} dv/v scatter plot", fontsize=22)

    legend1 = ax.legend(handles=legend_handles, fontsize=22, loc='upper left')
    ax.add_artist(legend1)

    cbar1 = plt.colorbar(ax.collections[0], ax=ax, pad=0.02)
    cbar1.set_label('Coherence value \n(darkness of the point)', fontsize=18)
    norm1 = Normalize(vmin=0, vmax=1)

    ax.tick_params(axis='both', which='both', labelsize=16, width=2, length=5)
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    fig.autofmt_xdate()
    
    # Adjust the layout if necessary
    fig.subplots_adjust(right=0.85)
    fig.tight_layout()
    
    return fig, f"{pair[0]}_{ np.min(ranges_list)}_{np.max(ranges_list)}_Hz_m{rolling_window}_wctscatter"
    

def save_figure(fig, filename, logger, mov_stack, components,filterid, visualize, plot_all_period=False, start=None, end=None, outfile=None):
    fig_path = os.path.join('Figures' if plot_all_period else 'Figures/Zooms')
    create_folder(fig_path, logger)
    mov_stack= mov_stack[0]
    if start and end:
        filename = f'{filename}_{str(start)[:10]}_{str(end)[:10]}'
    filepath = os.path.join(fig_path, f'{filename}.png')

    if outfile:
        if outfile.startswith("?"):
            if len(mov_stack) == 1:
                outfile = outfile.replace('?', '%s-f%i-m%s_%s-%s' % (components,
                                                                   filterid,
                                                                   mov_stack[0],
                                                                   mov_stack[1],
                                                                   visualize))
            else:
                outfile = outfile.replace('?', '%s-f%i-%s' % (components,
                                                               filterid,
                                                               visualize))
        filepath = "wct " + outfile
        logger.info("output to: %s" % outfile)

    fig.savefig(filepath, dpi=300, bbox_inches='tight', transparent=True)
    
def create_folder(folder_path, logger):
    try:
        os.makedirs(folder_path)
        logger.info(f"Folder '{folder_path}' created successfully.")
    except FileExistsError:
        pass

def xr_get_wct_pair(pair, components, filterid, mov_stack, logger):
    fn = os.path.join("WCT", "%02i" % filterid,
                      "%s_%s" % (mov_stack[0], mov_stack[1]),
                      "%s" % components, "%s.nc" % (pair))
    if not os.path.isfile(fn):
        logger.error("FILE DOES NOT EXIST: %s, skipping" % fn)
    
    data = xr_create_or_open(fn, name="WCT")
    if data is None:
        logger.error(f"Empty file for pair {pair}.")
        return None
    data = data.to_dataframe().unstack(level='frequency')
    return data

def xr_get_wct(components, filterid, mov_stack, logger):
    fn = os.path.join("WCT", "%02i" % filterid,
                      "%s_%s" % (mov_stack[0], mov_stack[1]),
                      "%s" % components, "*.nc" )
    matching_files = glob.glob(fn)
    if not matching_files:
       logger.error(f"No files found matching pattern: {fn}")

    all_wct = []
    for fil in matching_files:
        data = xr_create_or_open(fil, name="WCT")
        data = data.to_dataframe().unstack(level='frequency')
        all_wct.append(data)
    combined = pd.concat(all_wct, axis=0)
    dvv = combined.groupby(combined.index).mean()
    return dvv

def validate_and_adjust_date(date_string, end_date, logger):
    try:
        start_date = datetime.strptime(date_string, '%Y-%m-%d')
    except ValueError:
        try:
            days_delta = int(date_string)
            start_date = datetime.strptime(end_date, '%Y-%m-%d') + timedelta(days=days_delta)
        except ValueError:
            logger.error(f"Invalid start string: {date_string}")
            return None
    return start_date

def main(mov_stackid=None, components='ZZ', filterid=1,
        pairs=[], showALL=False, start="1970-01-01", end="2100-01-01", visualize='dvv', ranges="[0.5, 1.0], [1.0, 2.0], [2.0, 4.0]", show=True,outfile=None, loglevel="INFO"):
    logger = get_logger('msnoise.cc_dvv_plot_dvv', loglevel,
                        with_pid=True)
    db = connect()
    params = get_params(db)
    mincoh = params.dtt_mincoh

    # Check start and end dates
    if start == "1970-01-01":
         start= params.startdate
    else:
        start = validate_and_adjust_date(start, end, logger)
    if end == "2100-01-01":
        end = params.enddate

    # TODO clearer  mov_stackid to additionnal rolling
    if mov_stackid and mov_stackid != "": #if mov_stackid given
        try:
            mov_stack = params.mov_stack[mov_stackid - 1]
            if mov_stack in params.mov_stack:  # Check if mov_stack is in params.mov_stack
                mov_stacks = [mov_stack, ]
                rolling = 1
            else:
                rolling = mov_stackid  # Assign  mov_stack to rolling
        except:
            mov_stack = params.mov_stack[0]
            if mov_stack in params.mov_stack:  # Check if mov_stack is in params.mov_stack new format
                mov_stacks = [mov_stack, ]
                rolling = 1  # Keeping the mov_stack result
            else:
                rolling = mov_stack  # Assign mov_stack to rolling
    else:
        mov_stacks = params.mov_stack
        rolling = int(params.mov_stack[0][0][0])

    if components.count(","):
        components = components.split(",")
    else:
        components = [components, ]

    filter = get_filters(db, ref=filterid)
    low = float(filter.low)
    high = float(filter.high)

    for i, mov_stack in enumerate(mov_stacks):
        for comps in components:
            # Get the data
            if not pairs:
                dvv = xr_get_wct(comps, filterid, mov_stack, logger)
                pairs = ["all stations",]
            else:
                try:
                    dvv = xr_get_wct_pair(pairs, comps, filterid, mov_stack, logger)    
                except FileNotFoundError as fullpath:
                    logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                    continue
            # Plotting
            if dvv is None:
                logger.error(f"No data available for {pairs}. Skipping plot.")
                continue
            if visualize == 'dvv':
                fig, savename = plot_dvv_heatmap('dvv', dvv, pairs, rolling, start, end, low, high, logger, mincoh)
            elif visualize == 'coh':
                fig, savename = plot_dvv_heatmap('coh', dvv, pairs, rolling, start, end, low, high, logger, mincoh)
            elif visualize == 'curve':
                fig, savename = plot_dvv_scatter(dvv, pairs, rolling, start, end, ranges, logger)
            else:
                logger.error("PLOT TYPE DOES NOT EXIST: %s" % visualize)
            # Save and show the figure
            if fig is not None :
                save_figure(fig, savename, logger, mov_stacks, comps, filterid, visualize, plot_all_period=False, start=start, end=end, outfile=outfile)
                if show:
                    plt.show()
            else:
                logger.error("Figure was not created. Skipping save.")

if __name__ == "__main__":
    main()

