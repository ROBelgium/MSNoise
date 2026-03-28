"""MSNoise configuration management, parameter merging, and plot filename helpers."""
import traceback

from .db import get_logger
from ..msnoise_table_def import Config


_STEP_ABBREVS = {
    "preprocess":     "pre",
    "cc":             "cc",
    "filter":         "f",
    "stack":          "stk",
    "refstack":       "ref",
    "mwcs":           "mwcs",
    "mwcs_dtt":       "dtt",
    "mwcs_dtt_dvv":   "dvv",
    "stretching":     "str",
    "stretching_dvv": "sdvv",
    "wavelet":        "wct",
    "wavelet_dtt":    "wdtt",
    "wavelet_dtt_dvv": "wdvv",
    "psd":            "psd",
    "psd_rms":        "rms",
}

def get_config(session, name=None, isbool=False, plugin=None, category='global', set_number=None):
    """Get the value of one or all config bits from the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type name: str
    :param name: The name of the config bit to get. If omitted, a dictionary
        with all config items will be returned
    :type isbool: bool
    :param isbool: if True, returns True/False for config `name`. Defaults to
        False
    :type plugin: str
    :param plugin: if provided, gives the name of the Plugin config to use. E.g.
        if "Amazing" is provided, MSNoise will try to load the "AmazingConfig"
        entry point. See :doc:`plugins` for details.

    :rtype: str, bool or dict
    :returns: the value for `name` or a dict of all config values
    """
    if plugin:
        from importlib.metadata import entry_points
        for ep in entry_points(group='msnoise.plugins.table_def'):
            if ep.name.replace("Config", "") == plugin:
                table = ep.load()
    else:
        table = Config
    if name:
        query = session.query(table).filter(table.name == name)
        if hasattr(table, 'category') and category is not None:
            query = query.filter(table.category == category)
            if category == 'global':
                set_number = 1
        if hasattr(table, 'set_number') and set_number is not None:
            query = query.filter(table.set_number == set_number)
        config = query.first()
        if config is not None:
            if isbool:
                if config.value in [True, 'True', 'true', 'Y', 'y', '1', 1]:
                    config = True
                else:
                    config = False
            else:
                config = config.value
        else:
            config = ''
    else:
        config = {}
        configs = session.query(Config).all()
        for c in configs:
            config[c.name] = c.value
    return config



def update_config(session, name, value, plugin=None, category='global', set_number=None):
    """Update one config bit in the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :type name: str
    :param name: The name of the config bit to set.

    :type value: str
    :param value: The value of parameter `name`. Can also be NULL if you don't
        want to use this particular parameter.

    :type plugin: str
    :param plugin: if provided, gives the name of the Plugin config to use. E.g.
        if "Amazing" is provided, MSNoise will try to load the "AmazingConfig"
        entry point. See :doc:`plugins` for details.

    :type category: str
    :param category: The config category (default 'global'). Use e.g. 'stack',
        'filter', 'mwcs' to target a specific config set.

    :type set_number: int or None
    :param set_number: The config set number within the category (default None
        for global config). Use e.g. 1 for stack_1.

    """
    if plugin:
        from importlib.metadata import entry_points
        for ep in entry_points(group='msnoise.plugins.table_def'):
            if ep.name.replace("Config", "") == plugin:
                table = ep.load()
    else:
        table = Config
    query = session.query(table).filter(table.name == name)
    if hasattr(table, 'category') and category is not None:
        query = query.filter(table.category == category)
    if hasattr(table, 'set_number') and set_number is not None:
        query = query.filter(table.set_number == set_number)
    config = query.first()
    if config is None:
        return
    if "NULL" in str(value):
        config.value = None
    else:
        config.value = value
    session.commit()
    return



def create_config_set(session, set_name):
    """
    Create a configuration set for a workflow step.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type set_name: str
    :param set_name: The name of the workflow step (e.g., 'mwcs', 'mwcs_dtt', etc.)

    :rtype: int or None
    :returns: The set_number if set was created successfully, None otherwise
    """
    import os
    import csv
    from sqlalchemy import func

    # Define the config file path
    config_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config', f'config_{set_name}.csv')
    if not os.path.exists(config_file):
        return None

    # Find the next available set_number for this category
    max_set_number = session.query(func.max(Config.set_number)).filter(
        Config.category == set_name
    ).scalar()
    next_set_number = (max_set_number + 1) if max_set_number is not None else 1

    # Read the CSV file and add config entries
    with open(config_file, 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            name = row['name']
            default_value = row['default']
            definition = row.get('definition', '')
            param_type = row.get('type', 'str')
            possible_values = row.get('possible_values', '')

            # Create new config entry with set_number
            config = Config(
                name=name,
                category=set_name,
                set_number=next_set_number,
                value=default_value,
                param_type=param_type,
                default_value=default_value,
                description=definition,
                possible_values=possible_values,
                used_in=f"[{set_name}]",
                used=True
            )
            session.add(config)

    session.commit()
    return next_set_number



def delete_config_set(session, set_name, set_number):
    """
    Delete a configuration set for a workflow step.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type set_name: str
    :param set_name: The name of the workflow step (e.g., 'mwcs', 'mwcs_dtt', etc.)
    :type set_number: int
    :param set_number: The set number to delete
    :rtype: bool
    :returns: True if set was deleted successfully, False otherwise
    """
    try:
        from sqlalchemy import func
        from ..msnoise_table_def import Config

        # Validate input parameters
        if not isinstance(set_number, int) or set_number < 1:
            raise ValueError("set_number must be a positive integer")

        if not set_name or not isinstance(set_name, str):
            raise ValueError("set_name must be a non-empty string")

        # Check if the config set exists
        config_count = session.query(func.count(Config.ref)).filter(
            Config.category == set_name,
            Config.set_number == set_number
        ).scalar()

        if config_count == 0:
            return False  # Set doesn't exist

        # Delete all config entries for this set
        deleted_count = session.query(Config).filter(
            Config.category == set_name,
            Config.set_number == set_number
        ).delete()

        session.commit()

        get_logger("msnoise", "INFO").info(f"Deleted config set '{set_name}' set_number {set_number} "
                          f"({deleted_count} entries)")

        return True

    except Exception as e:
        session.rollback()
        traceback.print_exc()
        get_logger("msnoise").error(f"Failed to delete config set '{set_name}' set_number {set_number}: {str(e)}")
        return False



def list_config_sets(session, set_name=None):
    """
    List all configuration sets, optionally filtered by category.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object
    :type set_name: str or None
    :param set_name: Optional category filter (e.g., 'mwcs', 'mwcs_dtt', etc.)
    :rtype: list
    :returns: List of tuples (category, set_number, entry_count)
    """
    try:
        from sqlalchemy import func
        from ..msnoise_table_def import Config

        query = session.query(
            Config.category,
            Config.set_number,
            func.count(Config.ref).label('entry_count')
        ).filter(
            Config.set_number.isnot(None)  # Exclude global configs
        )

        if set_name:
            query = query.filter(Config.category == set_name)

        results = query.group_by(
            Config.category,
            Config.set_number
        ).order_by(
            Config.category,
            Config.set_number
        ).all()

        return [(row.category, row.set_number, row.entry_count) for row in results]

    except Exception as e:
        traceback.print_exc()
        get_logger("msnoise").error(f"Failed to list config sets: {str(e)}")

        return []



def get_config_set_details(session, set_name, set_number, format="list"):
    """
    Get details of a specific configuration set.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object
    :type set_name: str
    :param set_name: The category name
    :type set_number: int
    :param set_number: The set number
    :rtype: list
    :returns: List of config entries in the set
    """
    try:
        from ..msnoise_table_def import Config

        configs = session.query(Config).filter(
            Config.category == set_name,
            Config.set_number == set_number
        ).order_by(Config.name).all()

        if format == "list":
            return [{'name': c.name, 'value': c.value, 'description': c.description}
             for c in configs]
        elif format == "AttribDict":
            from obspy.core import AttribDict
            import pydoc
            attrib_dict = AttribDict()
            for c in configs:
                if c.value is None or c.value == '':
                    if c.param_type == 'str':
                        attrib_dict[c.name] = ''
                    # else: skip — can't convert empty string to int/float/eval
                    continue
                itemtype = pydoc.locate(c.param_type)
                if itemtype is bool:
                    if c.value in [True, 'True', 'true', 'Y', 'y', '1', 1]:
                        attrib_dict[c.name] = True
                    else:
                        attrib_dict[c.name] = False
                else:
                    attrib_dict[c.name] = itemtype(c.value)
            return attrib_dict


    except Exception as e:
        from obspy.core import AttribDict
        get_logger("msnoise", loglevel="ERROR").error(f"Failed to get config set details: {str(e)}")
        return AttribDict() if format == "AttribDict" else []



def get_config_categories_definition():
    """Get the standard configuration categories with display names, order, and indent level.

    Each entry is ``(category_key, display_name, level)`` where *level* is the depth
    relative to ``global`` (0).  Used by the config-sets admin page for visual
    indentation that mirrors the workflow graph.
    """
    return [
        # Tree order: children immediately follow their parent
        ('global',      'Global Parameters',  0),
        ('preprocess',  'Preprocessing',       1),
        ('cc',          'Cross-Correlation',   2),
        ('filter',      'Filters',             3),
        ('stack',       'Moving Stacks',       4),
        ('refstack',    'Reference Stacks',    5),
        ('mwcs',          'MWCS',                    6),
        ('mwcs_dtt',      'MWCS dt/t',               7),
        ('mwcs_dtt_dvv',  'MWCS dv/v Aggregate',     8),
        ('stretching',    'Stretching',               6),
        ('stretching_dvv','Stretching dv/v Aggregate',7),
        ('wavelet',       'Wavelet',                  6),
        ('wavelet_dtt',   'Wavelet dt/t',             7),
        ('wavelet_dtt_dvv',   'WCT dv/v Aggregate',       8),
        ('psd',         'PSD',                 1),
        ('psd_rms',     'PSD RMS',             2),
    ]



def get_config_sets_organized(session):
    """Get configuration sets organized by category in the standard order"""
    from sqlalchemy import func
    from .. import msnoise_table_def as schema

    # Get category definitions
    category_order = get_config_categories_definition()

    # Get all config sets grouped by category and set_number
    all_sets = session.query(
        schema.Config.category,
        schema.Config.set_number,
        func.count(schema.Config.ref).label('param_count')
    ).group_by(
        schema.Config.category,
        schema.Config.set_number
    ).all()

    # Migrate legacy NULL set_number rows (written by pre-refactor code
    # where global config had no set_number).  Fold them into set_number=1
    # so they are accessible via the normal workflow path.  Skip if a
    # set_number=1 row already exists for that category (no overwrite).
    null_categories = {cat for cat, sn, _ in all_sets if sn is None}
    existing_1 = {cat for cat, sn, _ in all_sets if sn == 1}
    for null_cat in null_categories - existing_1:
        session.query(schema.Config).filter(
            schema.Config.category == null_cat,
            schema.Config.set_number == None,  # noqa: E711
        ).update({"set_number": 1}, synchronize_session=False)
    if null_categories - existing_1:
        session.commit()
        # Re-query after migration
        all_sets = session.query(
            schema.Config.category,
            schema.Config.set_number,
            func.count(schema.Config.ref).label('param_count')
        ).group_by(
            schema.Config.category,
            schema.Config.set_number
        ).all()

    # Group sets by category, skipping any remaining NULL rows (e.g.
    # where set_number=1 already existed so migration was skipped).
    sets_by_category = {}
    for category, set_number, param_count in all_sets:
        if set_number is None:
            continue
        if category not in sets_by_category:
            sets_by_category[category] = []
        sets_by_category[category].append({
            'category': category,
            'set_number': set_number,
            'param_count': param_count
        })

    # Sort sets within each category by set_number
    for category in sets_by_category:
        sets_by_category[category].sort(key=lambda x: x['set_number'])

    # Create ordered list of categories with their sets
    ordered_categories = []
    for category_key, display_name, level in category_order:
        if category_key in sets_by_category:
            ordered_categories.append({
                'category': category_key,
                'display_name': display_name,
                'level': level,
                'sets': sets_by_category[category_key]
            })

    # Add any additional categories not in the predefined order
    for category in sets_by_category:
        if category not in [cat[0] for cat in category_order]:
            ordered_categories.append({
                'category': category,
                'display_name': category.title(),
                'level': 0,
                'sets': sets_by_category[category]
            })

    return ordered_categories


# ── Parameters ─────────────────────────────────────────────


def get_params(session):
    """Get global config parameters as a single-layer :class:`~msnoise.params.LayeredParams`.

    Access via ``params.global_.<key>`` e.g. ``params.global_.hpc``.
    For full pipeline params (per-step config included), use
    :func:`get_merged_params_for_lineage` instead.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :rtype: :class:`~msnoise.params.LayeredParams`
    """
    from obspy.core.util.attribdict import AttribDict
    from ..default import default
    from ..params import LayeredParams
    s = session
    global_attrib = AttribDict()
    for name in default.keys():
        itemtype = default[name].type
        if itemtype is bool:
            global_attrib[name] = get_config(s, name, isbool=True)
        else:
            global_attrib[name] = itemtype(get_config(s, name))
    p = LayeredParams()
    p._add_layer("global", global_attrib)
    return p



def lineage_to_plot_tag(lineage_names):
    """Convert a list of step names to a compact plot filename tag.

    Each step is abbreviated using :data:`_STEP_ABBREVS` and the set number is
    appended.  Steps are joined with ``-``.

    Examples::

        lineage_to_plot_tag(["preprocess_1","cc_1","filter_1","stack_1"])
        # → "pre1-cc1-f1-stk1"

        lineage_to_plot_tag(["preprocess_1","cc_1","filter_1","stack_1",
                              "refstack_1","mwcs_1","mwcs_dtt_1","mwcs_dtt_dvv_1"])
        # → "pre1-cc1-f1-stk1-ref1-mwcs1-dtt1-dvv1"

    :param lineage_names: List of step name strings from
        :attr:`MSNoiseResult.lineage_names`.
    :returns: Hyphen-joined abbreviated tag string.
    """
    parts = []
    for step in lineage_names:
        tail = step.rsplit("_", 1)
        if len(tail) == 2 and tail[1].isdigit():
            base, num = tail[0], tail[1]
        else:
            base, num = step, ""
        abbrev = _STEP_ABBREVS.get(base, base)
        parts.append(f"{abbrev}{num}")
    return "-".join(parts)



def build_plot_outfile(outfile, plot_name, lineage_names, *,
                       pair=None, components=None,
                       mov_stack=None, extra=None):
    """Resolve the output filename for a plot.

    When *outfile* starts with ``"?"``, replaces ``"?"`` with a standardised
    auto-generated tag and prepends *plot_name*, producing a filename of the
    form::

        <plot_name>__<lineage>[__<pair>][__<comp>][__m<ms>][__<extra>].<ext>

    If *outfile* does not start with ``"?"``, it is returned unchanged (the
    caller supplied an explicit path).

    If *outfile* is ``None`` or empty, returns ``None``.

    :param outfile: Raw ``outfile`` argument passed to the plot function
        (e.g. ``"?.png"``).
    :param plot_name: Short plot identifier, e.g. ``"ccftime"`` or
        ``"dvv_mwcs"``.  Must not contain ``__``.
    :param lineage_names: List of step name strings (e.g.
        ``["preprocess_1", "cc_1", "filter_1", "stack_1"]``) from
        :attr:`MSNoiseResult.lineage_names`.
    :param pair: Optional station pair string (``"NET.STA.LOC-NET.STA.LOC"``).
        Dots are kept; colons are replaced with ``-``.
    :param components: Optional component string, e.g. ``"ZZ"``.
    :param mov_stack: Optional moving-stack tuple or string,
        e.g. ``("1D", "1D")`` → ``"m1D-1D"``.
    :param extra: Optional extra qualifier string appended last.
    :returns: Resolved filename string, or *outfile* unchanged.
    """
    if not outfile:
        return outfile
    if not outfile.startswith("?"):
        return outfile

    # Build tag components
    lin_tag = lineage_to_plot_tag(lineage_names)
    parts = [lin_tag]

    if pair:
        pair_clean = str(pair).replace(":", "-").replace(" ", "")
        parts.append(pair_clean)
    if components:
        parts.append(str(components))
    if mov_stack is not None:
        if isinstance(mov_stack, (tuple, list)) and len(mov_stack) == 2:
            parts.append(f"m{mov_stack[0]}-{mov_stack[1]}")
        else:
            parts.append(f"m{mov_stack}")
    if extra:
        parts.append(str(extra))

    tag = "__".join(parts)
    # Replace the leading '?' — preserve any suffix the user appended (e.g. '.png')
    resolved = outfile.replace("?", tag, 1)
    return f"{plot_name}__{resolved}"



def _load_step_config(db, step):
    # step must contain config_category + config_set_number (or equivalent)
    return get_config_set_details(
        db,
        step.category,
        step.set_number,
        format="AttribDict",
    )



def get_merged_params_for_lineage(db, orig_params, step_params, lineage):
    """Build a :class:`~msnoise.params.LayeredParams` for the given lineage.

    Returns ``(lineage, lineage_names, LayeredParams)``.
    """
    from ..params import _build_layered_params

    lineage = [s for s in lineage if s.category not in {"global"}]
    lineage_names = [s.step_name for s in lineage]
    lineage_cfgs = [_load_step_config(db, s) for s in lineage]

    # orig_params is now a single-layer LayeredParams; unwrap the global AttribDict
    global_attrib = (
        orig_params.global_
        if hasattr(orig_params, "_layers")
        else orig_params  # backward compat if raw AttribDict passed
    )
    params = _build_layered_params(
        global_attrib=global_attrib,
        lineage_steps=lineage,
        lineage_names=lineage_names,
        step_configs=lineage_cfgs,
    )
    return lineage, lineage_names, params
