.. include:: configs.hrst

******************************
Extending MSNoise with Plugins
******************************

.. versionadded:: 1.4
.. versionchanged:: 2.0 New entry point groups for workflow topology, config sets,
   and job types. Stable import surface via :mod:`msnoise.plugins`.

MSNoise supports plugins — external Python packages that extend the default
``preprocess → CC → stack → mwcs/stretching/wavelet → dv/v`` workflow with
new processing steps, new job types, or new web-admin pages.

.. contents::
    :local:


Stable import surface
---------------------

Plugin code should import from :mod:`msnoise.plugins` rather than from
``msnoise.core.*`` directly.  This protects plugins from internal restructuring:

.. code-block:: python

    from msnoise.plugins import (
        connect,
        get_config, update_config,
        get_params,
        get_workflow_chains, get_workflow_order,
        get_next_lineage_batch, is_next_job_for_step,
        massive_update_job, propagate_downstream,
        xr_get_ccf, xr_save_ccf,
        get_stations, get_station_pairs,
    )

.. automodule:: msnoise.plugins
   :members:


pyproject.toml entry points
----------------------------

Declare all entry points in your ``pyproject.toml``:

.. code-block:: toml

    [project.entry-points."msnoise.plugins.commands"]
    myplugin = "myplugin.cli:cli"

    [project.entry-points."msnoise.plugins.jobtypes"]
    myplugin = "myplugin.jobs:get_jobtypes"

    [project.entry-points."msnoise.plugins.workflow_chains"]
    myplugin = "myplugin.workflow:get_chains"

    [project.entry-points."msnoise.plugins.workflow_order"]
    myplugin = "myplugin.workflow:get_order"

    [project.entry-points."msnoise.plugins.table_def"]
    myplugin = "myplugin.tabledef:MyPluginConfig"

Alternatively, using the legacy ``setup.py`` format:

.. code-block:: python

    entry_points = {
        "msnoise.plugins.commands":       ["myplugin = myplugin.cli:cli"],
        "msnoise.plugins.jobtypes":       ["myplugin = myplugin.jobs:get_jobtypes"],
        "msnoise.plugins.workflow_chains":["myplugin = myplugin.workflow:get_chains"],
        "msnoise.plugins.workflow_order": ["myplugin = myplugin.workflow:get_order"],
        "msnoise.plugins.table_def":      ["myplugin = myplugin.tabledef:MyPluginConfig"],
    }


Entry point groups
------------------

msnoise.plugins.commands
~~~~~~~~~~~~~~~~~~~~~~~~

Adds a :class:`click.Group` (or :class:`click.Command`) to the ``msnoise``
CLI.  Once the plugin is declared in the project's ``plugins`` config key and
the package is installed, the command appears under ``msnoise plugin``:

.. code-block:: python

    # myplugin/cli.py
    import click

    @click.group()
    def cli():
        """My Plugin commands."""

    @cli.command()
    def compute():
        """Run my custom computation."""
        from .compute import main
        main()

.. code-block:: sh

    $ msnoise plugin myplugin compute


msnoise.plugins.jobtypes
~~~~~~~~~~~~~~~~~~~~~~~~

Registers custom job types that are created during ``msnoise new_jobs``.
The callable must return a list of dicts with ``name`` and ``after`` keys:

.. code-block:: python

    # myplugin/jobs.py
    def get_jobtypes():
        return [
            {"name": "MYPLUGIN", "after": "scan_archive"},
        ]

Supported ``after`` values: ``"scan_archive"``, ``"new_files"``.


msnoise.plugins.workflow_chains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adds new categories to the workflow DAG.  The callable must return a dict
with the same schema as :data:`msnoise.core.workflow._BUILTIN_WORKFLOW_CHAINS`.
Once registered, the category is automatically valid in ``msnoise new_jobs
--after``, in job propagation, ``MSNoiseResult``, and
``msnoise utils run_workflow``:

.. code-block:: python

    # myplugin/workflow.py
    def get_chains():
        return {
            "my_dtt": {
                "next_steps":    ["my_dtt_dvv"],
                "is_entry_point": False,
                "is_terminal":   False,
                "abbrev":        "mdtt",
                "display_name":  "My dt/t",
                "level":         7,
                "cli":           ["plugin", "myplugin", "compute_dtt"],  # ← required for run_workflow
            },
            "my_dtt_dvv": {
                "next_steps":    [],
                "is_entry_point": False,
                "is_terminal":   True,
                "abbrev":        "mdvv",
                "display_name":  "My dv/v Aggregate",
                "level":         8,
                "cli":           ["plugin", "myplugin", "compute_dvv"],  # ← required for run_workflow
            },
        }

.. note::

    The ``"cli"`` key maps the category to the ``msnoise`` sub-command token
    list used by ``msnoise utils run_workflow``.  Pass-through categories
    (no worker, no jobs) should use ``"cli": None``.  If omitted, the category
    is silently skipped by ``run_workflow``.


msnoise.plugins.workflow_order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Appends new category names to the canonical display order used by the
admin UI and ``msnoise config list``:

.. code-block:: python

    # myplugin/workflow.py
    def get_order():
        return ["my_dtt", "my_dtt_dvv"]


msnoise.plugins.table_def
~~~~~~~~~~~~~~~~~~~~~~~~~~

Exposes a plugin's SQLAlchemy config table class so that
:func:`~msnoise.core.config.get_config` and
:func:`~msnoise.core.config.update_config` can access it via the
``plugin=`` argument:

.. code-block:: python

    # myplugin/tabledef.py
    from sqlalchemy import Column, String
    from sqlalchemy.orm import declarative_base

    Base = declarative_base()

    class MyPluginConfig(Base):
        __tablename__ = "myplugin-config"
        name  = Column(String(255), primary_key=True)
        value = Column(String(255))

    def get_table_def():
        return MyPluginConfig

Usage from Python:

.. code-block:: python

    from msnoise.plugins import connect, get_config
    db = connect()
    val = get_config(db, "myparam", plugin="MyPlugin")


Declaring the plugin in a project
----------------------------------

After installing the plugin package, declare its **package name** in the
project's ``plugins`` global configuration key (comma-separated, no spaces):

.. code-block:: sh

    $ msnoise config set plugins myplugin
    # or for multiple plugins:
    $ msnoise config set plugins myplugin,anotherplugin

The plugin will then appear in ``msnoise plugin`` and its workflow categories
will be available in ``msnoise new_jobs --after``.


Minimal plugin structure
-------------------------

.. code-block:: text

    msnoise-myplugin/
    ├── pyproject.toml
    └── myplugin/
        ├── __init__.py
        ├── cli.py        ← click commands
        ├── compute.py    ← worker main()
        ├── jobs.py       ← get_jobtypes()
        ├── tabledef.py   ← SQLAlchemy config table
        └── workflow.py   ← get_chains(), get_order()

