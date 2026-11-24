"""MSNoise configuration management, parameter merging, and plot filename helpers."""

__all__ = [
    "build_plot_outfile",
    "create_config_set",
    "create_project_from_yaml",
    "export_project_to_yaml",
    "delete_config_set",
    "get_config",
    "get_config_categories_definition",
    "get_config_set_details",
    "get_config_sets_organized",
    "get_merged_params_for_lineage",
    "get_params",
    "lineage_to_plot_tag",
    "list_config_sets",
    "parse_config_key",
    "update_config",
]

import traceback

from .db import get_logger
from ..msnoise_table_def import Config


def _get_step_abbrevs():
    """Return category→abbreviation map from workflow (plugin-aware)."""
    from .workflow import get_step_abbrevs
    return get_step_abbrevs()


def parse_config_key(key):
    """Parse a dot-notation config key into ``(category, set_number, name)``.

    Accepted forms::

        output_folder            →  ("global", 1, "output_folder")   # global shorthand
        cc.cc_sampling_rate      →  ("cc",     1, "cc_sampling_rate")
        cc.2.cc_sampling_rate    →  ("cc",     2, "cc_sampling_rate")
        mwcs.2.mwcs_wlen         →  ("mwcs",   2, "mwcs_wlen")

    :raises ValueError: for malformed keys (too many parts, non-integer set
        number, or bare name that is not a global parameter).
    """
    parts = key.split(".")
    if len(parts) == 1:
        # Bare name → global shorthand (e.g. "output_folder")
        return ("global", 1, parts[0])
    elif len(parts) == 2:
        # category.name  → set_number=1
        category, name = parts
        return (category, 1, name)
    elif len(parts) == 3:
        # category.set_number.name
        category, set_str, name = parts
        try:
            set_number = int(set_str)
        except ValueError:
            raise ValueError(
                f"Invalid set number {set_str!r} in key {key!r}. "
                "Use the form 'category.N.name' where N is an integer."
            )
        return (category, set_number, name)
    else:
        raise ValueError(
            f"Too many parts in config key {key!r}. "
            "Expected 'name', 'category.name', or 'category.N.name'."
        )


def _cast_config_value(name, value_str, param_type):
    """Validate and cast *value_str* to the type declared for this parameter.

    :param name: Parameter name (for error messages).
    :param value_str: Raw string value from the CLI.
    :param param_type: One of ``'str'``, ``'int'``, ``'float'``, ``'bool'``, ``'eval'``.
    :returns: The validated string to store in the DB (always stored as string).
    :raises ValueError: if the value is not valid for *param_type*.
    """
    if param_type == "bool":
        normalised = value_str.strip().upper()
        if normalised in ("Y", "YES", "TRUE", "1"):
            return "Y"
        elif normalised in ("N", "NO", "FALSE", "0"):
            return "N"
        else:
            raise ValueError(
                f"'{name}' expects a boolean (Y/N). Got {value_str!r}. "
                "Use Y, N, True, False, 1, or 0."
            )
    elif param_type == "int":
        try:
            int(value_str)
        except ValueError:
            raise ValueError(
                f"'{name}' expects an integer. Got {value_str!r}."
            )
    elif param_type == "float":
        try:
            float(value_str)
        except ValueError:
            raise ValueError(
                f"'{name}' expects a float. Got {value_str!r}."
            )
    # 'str' and 'eval' accept any string — no further validation
    return value_str


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
    return config



def create_config_set(session, set_name, set_number_hint=None):
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

    # If caller hints a set_number and rows already exist for it, reuse — don't duplicate.
    if set_number_hint is not None:
        from sqlalchemy import func as _f
        existing_count = session.query(_f.count(Config.ref)).filter(
            Config.category == set_name,
            Config.set_number == set_number_hint,
        ).scalar()
        if existing_count:
            return set_number_hint

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



def create_project_from_yaml(session, yaml_path):
    """Seed a freshly-initialised database from a **project YAML** file.

    The project YAML uses ``category_N`` keys so that multiple config sets of
    the same category are unambiguous.  Each entry may declare an ``after``
    field (string or list of strings) naming the ``category_N`` steps that
    must precede it; one :class:`~msnoise.msnoise_table_def.WorkflowLink` row
    is created per ``after`` entry.  Steps whose ``after`` target is absent
    from the YAML produce a warning instead of an error, so partial YAMLs and
    incremental DB builds are both supported.

    Example::

        msnoise_project_version: 1

        global_1:
          startdate: "2013-04-01"
          enddate:   "2014-10-31"

        preprocess_1:
          after: global_1
          cc_sampling_rate: 20.0

        cc_1:
          after: preprocess_1
          cc_type_single_station_AC: PCC
          whitening: "N"

        filter_1:
          after: cc_1          # single parent
          freqmin: 1.0
          freqmax: 2.0
          AC: "Y"

        filter_2:
          after: cc_1          # fan-out from the same cc step
          freqmin: 0.5
          freqmax: 1.0
          AC: "Y"

        stack_1:
          after: [filter_1, filter_2]
          mov_stack: "(('2D','1D'))"

        refstack_1:
          after: [filter_1, filter_2]   # sibling of stack, not child
          ref_begin: "2013-04-01"
          ref_end:   "2014-10-31"

        mwcs_1:
          after: [stack_1, refstack_1]   # list: requires both parents
          freqmin: 1.0
          freqmax: 2.0

    :param session: SQLAlchemy session (from :func:`~msnoise.core.db.connect`).
    :param yaml_path: Path to a project YAML file.
    :returns: ``(created_steps, warnings)`` — list of ``category_N`` strings
        created and list of warning strings for missing ``after`` targets.
    :raises ValueError: if the YAML is missing ``msnoise_project_version`` or
        a key is not in ``category_N`` format.
    """
    import re
    import logging
    import yaml
    from ..core.workflow import create_workflow_step, create_workflow_link

    logger = logging.getLogger(__name__)

    with open(yaml_path, encoding="utf-8") as fh:
        doc = yaml.safe_load(fh)

    if doc.get("msnoise_project_version") != 1:
        raise ValueError(
            f"{yaml_path!r} is not a project YAML "
            f"(missing or wrong msnoise_project_version). "
            f"Note: MSNoiseParams.to_yaml() produces a per-lineage YAML; "
            f"use a project YAML with 'category_N' keys instead."
        )

    _SKIP = {"msnoise_project_version", "generated", "description", "data_sources", "stations"}
    _KEY_RE = re.compile(r"^(?P<category>.+)_(?P<set_number>\d+)$")

    # Parse entries, preserving YAML order (insertion order in py3.7+ dicts).
    entries = []       # [(yaml_key, category, set_number_hint, after_list, overrides)]
    for key, body in doc.items():
        if key in _SKIP:
            continue
        m = _KEY_RE.match(key)
        if not m:
            raise ValueError(
                f"Project YAML key {key!r} is not in 'category_N' format."
            )
        body = body or {}
        after_raw = body.pop("after", None)
        after_list = (
            [after_raw] if isinstance(after_raw, str)
            else list(after_raw) if after_raw
            else []
        )
        entries.append((key, m.group("category"), int(m.group("set_number")),
                        after_list, body))

    # Pass 1: create config sets and workflow steps, collect step objects.
    step_by_yaml_key = {}   # yaml_key → WorkflowStep
    created_steps = []
    warnings = []

    for yaml_key, category, _set_number_hint, _after, overrides in entries:
        set_number = create_config_set(session, category, set_number_hint=_set_number_hint)
        if set_number is None:
            w = f"No config CSV for category {category!r} (key {yaml_key!r}) — skipped."
            warnings.append(w)
            logger.warning(w)
            continue
        for key, value in overrides.items():
            if isinstance(value, list):
                if value and isinstance(value[0], (list, tuple)):
                    value = str(tuple(tuple(v) for v in value))
                else:
                    value = ",".join(str(v) for v in value)
            result = update_config(session, key, str(value),
                                   category=category, set_number=set_number)
            if result is None:
                w = f"{yaml_key!r}: unknown key {key!r} — ignored (typo?)."
                warnings.append(w)
                logger.warning(w)
        step = create_workflow_step(
            session, f"{category}_{set_number}", category, set_number,
            description=f"Created from {yaml_path}",
        )
        step_by_yaml_key[yaml_key] = step
        created_steps.append(f"{category}_{set_number}")

    # Pass 2: create links from after declarations.
    for yaml_key, _cat, _sn, after_list, _overrides in entries:
        if yaml_key not in step_by_yaml_key:
            continue   # step was skipped (no CSV)
        to_step = step_by_yaml_key[yaml_key]
        for parent_key in after_list:
            if parent_key not in step_by_yaml_key:
                w = (f"Link {parent_key!r} → {yaml_key!r}: "
                     f"{parent_key!r} not found in this YAML — "
                     f"add the link manually via the admin UI.")
                warnings.append(w)
                logger.warning(w)
                continue
            from_step = step_by_yaml_key[parent_key]
            create_workflow_link(session, from_step.step_id, to_step.step_id)

    # Pass 3: data sources
    # The first declared data source overwrites the installer default
    # (ref=1), which is what populate/scan_archive use when stations have no
    # explicit data_source_id. Further entries are added as new sources.
    created_ds_ids = []
    from ..core.stations import (add_data_source, get_data_source,
                                 update_data_source)
    for i, ds_dict in enumerate(doc.get("data_sources", []) or []):
        fields = dict(
            name=ds_dict["name"],
            uri=ds_dict.get("uri", ""),
            data_structure=ds_dict.get("data_structure", "SDS"),
            auth_env=ds_dict.get("auth_env", "MSNOISE"),
        )
        if i == 0 and get_data_source(session, id=1) is not None:
            ds = update_data_source(session, 1, **fields)
        else:
            ds = add_data_source(session, **fields)
        created_ds_ids.append(ds.ref)

    # Pass 4: stations — endpoint fetch (level=response recommended so that
    # import_stationxml saves responses to response_path automatically) or
    # explicit list (same format as MSNoiseResult bundle provenance export,
    # so bundle provenance can be copy-pasted directly into a project YAML).
    stations_section = doc.get("stations")
    if stations_section:
        from ..core.stations import import_stationxml, update_station
        if isinstance(stations_section, dict):
            endpoint = stations_section.get("station_endpoint")
            if endpoint:
                try:
                    created_s, updated_s, _ = import_stationxml(session, endpoint)
                    logger.info(
                        f"Stations via endpoint: {created_s} created, "
                        f"{updated_s} updated."
                    )
                except Exception as e:
                    w = f"station_endpoint fetch failed: {e}"
                    warnings.append(w)
                    logger.warning(w)
        elif isinstance(stations_section, list):
            for s in stations_section:
                update_station(
                    session,
                    net=s["net"], sta=s["sta"],
                    X=s["X"], Y=s["Y"],
                    altitude=s.get("altitude", 0.0),
                    coordinates=s.get("coordinates", "DEG"),
                    used=s.get("used", True),
                )

    # Pass 4b: if YAML declared datasources, assign the first one to all
    # stations whose data_source_id is still NULL (installer default).
    if created_ds_ids:
        from ..msnoise_table_def import Station
        target_ds_id = created_ds_ids[0]
        stations = session.query(Station).filter(Station.data_source_id == None).all()
        for sta in stations:
            sta.data_source_id = target_ds_id
        session.commit()
        if len(created_ds_ids) > 1:
            w = (f"Multiple data_sources declared; assigned datasource id={target_ds_id} "
                 f"to all stations. Update manually for per-station assignment.")
            warnings.append(w)
            logger.warning(w)

    return created_steps, warnings


def export_project_to_yaml(session, yaml_path, only_non_defaults=True):
    """Snapshot the current project DB state as a **project YAML** file.

    The inverse of :func:`create_project_from_yaml`.  Reads all active
    :class:`~msnoise.msnoise_table_def.WorkflowStep` rows, their config sets,
    the :class:`~msnoise.msnoise_table_def.WorkflowLink` topology, all
    :class:`~msnoise.msnoise_table_def.DataSource` rows, and all
    :class:`~msnoise.msnoise_table_def.Station` rows, and writes a
    ``msnoise_project_version: 1`` YAML that can be re-applied with
    ``msnoise db init --from-yaml``.

    :param session: SQLAlchemy session.
    :param yaml_path: Destination path for the YAML file.
    :param only_non_defaults: When ``True`` (default) only keys whose value
        differs from the CSV default are written, keeping the file minimal.
        Set to ``False`` to write every key explicitly.
    :returns: Path string of the written file.
    """
    import os
    import csv
    import datetime
    import yaml
    from ..core.workflow import get_workflow_steps, get_workflow_links, get_workflow_order
    from ..msnoise_table_def import DataSource, Station

    # ── load CSV defaults once ────────────────────────────────────────────
    _config_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "config")

    def _csv_defaults(category):
        path = os.path.join(_config_dir, f"config_{category}.csv")
        if not os.path.exists(path):
            return {}
        with open(path, newline="", encoding="utf-8") as fh:
            return {row["name"]: row["default"] for row in csv.DictReader(fh)}

    # ── build step_name → parent step_names map from WorkflowLinks ───────
    steps = get_workflow_steps(session)
    links = get_workflow_links(session)

    step_by_id = {s.step_id: s for s in steps}
    parents_of = {s.step_name: [] for s in steps}   # step_name → [parent step_names]
    for lnk in links:
        child = step_by_id.get(lnk.to_step_id)
        parent = step_by_id.get(lnk.from_step_id)
        if child and parent:
            parents_of[child.step_name].append(parent.step_name)

    # ── sort steps by canonical workflow order ────────────────────────────
    _order = get_workflow_order()
    def _sort_key(s):
        try:
            return (_order.index(s.category), s.set_number)
        except ValueError:
            return (len(_order), s.set_number)

    steps_sorted = sorted(steps, key=_sort_key)

    # ── build YAML document (ordered dict preserves insertion order) ──────
    doc = {
        "msnoise_project_version": 1,
        "generated": datetime.datetime.now(datetime.timezone.utc).isoformat(
            timespec="seconds") + "Z",
    }

    for step in steps_sorted:
        entry = {}

        # after
        pnames = sorted(parents_of.get(step.step_name, []))
        if len(pnames) == 1:
            entry["after"] = pnames[0]
        elif len(pnames) > 1:
            entry["after"] = pnames

        # config values
        defaults = _csv_defaults(step.category)
        rows = get_config_set_details(session, step.category, step.set_number,
                                      format="list")
        for row in rows:
            name, value = row["name"], row["value"]
            if only_non_defaults and value == defaults.get(name):
                continue
            entry[name] = value

        doc[step.step_name] = entry

    # ── data sources ──────────────────────────────────────────────────────
    data_sources = session.query(DataSource).all()
    if data_sources:
        doc["data_sources"] = [
            {
                "name":           ds.name,
                "uri":            ds.uri or "",
                "data_structure": ds.data_structure or "SDS",
                "auth_env":       ds.auth_env or "MSNOISE",
            }
            for ds in data_sources
        ]

    # ── stations ──────────────────────────────────────────────────────────
    stations = session.query(Station).order_by(Station.net, Station.sta).all()
    if stations:
        doc["stations"] = [
            {
                "net":            s.net,
                "sta":            s.sta,
                "X":              s.X,
                "Y":              s.Y,
                "altitude":       s.altitude,
                "coordinates":    s.coordinates or "DEG",
                "data_source_id": 0,   # index into data_sources list above
            }
            for s in stations
        ]

    with open(yaml_path, "w", encoding="utf-8") as fh:
        yaml.dump(doc, fh, default_flow_style=False, sort_keys=False,
                  allow_unicode=True)

    return yaml_path


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
    """Return ordered display metadata for all workflow categories.

    Each entry is a ``(category_key, display_name, level)`` tuple where *level*
    is the DAG depth (global=0, preprocess/psd=1 — both branch off global, …).

    Derived from :func:`~msnoise.core.workflow.get_category_display_info` so
    plugin-added categories appear automatically.  The list is ordered for UI
    rendering (CC branch first, PSD branch at end).
    """
    from .workflow import get_category_display_info
    return [
        (d["category"], d["display_name"], d["level"])
        for d in get_category_display_info()
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

    # Group sets by category
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
    """Build a single-layer :class:`~msnoise.params.MSNoiseParams` for global config.

    Queries the ``Config`` table directly for all ``(category='global',
    set_number=1)`` rows and casts each value to its declared ``param_type``.
    Does **not** depend on ``default.py``.

    For full pipeline params (per-step config included), use
    :func:`get_merged_params_for_lineage` instead.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :rtype: :class:`~msnoise.params.MSNoiseParams`
    """
    from obspy.core.util.attribdict import AttribDict
    from ..params import MSNoiseParams

    _BOOL_TRUE = {"y", "yes", "true", "1"}
    _TYPE_MAP  = {"int": int, "float": float, "bool": None, "str": str, "eval": str}

    rows = (
        session.query(Config)
        .filter(Config.category == "global")
        .filter(Config.set_number == 1)
        .all()
    )

    global_attrib = AttribDict()
    for row in rows:
        raw   = row.value if row.value is not None else ""
        ptype = (row.param_type or "str").lower()
        if ptype == "bool":
            global_attrib[row.name] = raw.strip().lower() in _BOOL_TRUE
        elif ptype in _TYPE_MAP and _TYPE_MAP[ptype] is not None:
            try:
                global_attrib[row.name] = _TYPE_MAP[ptype](raw)
            except (ValueError, TypeError):
                global_attrib[row.name] = raw
        else:
            global_attrib[row.name] = raw

    p = MSNoiseParams()
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
        abbrev = _get_step_abbrevs().get(base, base)
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
    """Build a :class:`~msnoise.params.MSNoiseParams` for the given lineage.

    Returns ``(lineage, lineage_names, MSNoiseParams)``.
    """
    from ..params import _build_msnoise_params

    lineage = [s for s in lineage if s.category not in {"global"}]
    lineage_names = [s.step_name for s in lineage]
    lineage_cfgs = [_load_step_config(db, s) for s in lineage]

    # orig_params is now a single-layer MSNoiseParams; unwrap the global AttribDict
    global_attrib = (
        orig_params.global_
        if hasattr(orig_params, "_layers")
        else orig_params  # backward compat if raw AttribDict passed
    )
    params = _build_msnoise_params(
        global_attrib=global_attrib,
        lineage_steps=lineage,
        lineage_names=lineage_names,
        step_configs=lineage_cfgs,
    )
    return lineage, lineage_names, params
