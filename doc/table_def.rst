Table Definitions
==================

.. module:: msnoise
   :no-index:

The following SQLAlchemy ORM classes represent the MSNoise database schema.
Column attributes are described in each class docstring.

Job
---

.. autoclass:: msnoise.msnoise_table_def.Job
   :members: config_category, config_set_number, lineage, step_name

Station
-------

.. autoclass:: msnoise.msnoise_table_def.Station
   :members: chans, locs

Config
------

.. autoclass:: msnoise.msnoise_table_def.Config
   :members: get_config_dict, get_global_config, get_workflow_config, get_typed_value, is_global_config, is_workflow_config, validate_value

DataAvailability
----------------

.. autoclass:: msnoise.msnoise_table_def.DataAvailability

DataSource
----------

.. autoclass:: msnoise.msnoise_table_def.DataSource

Lineage
-------

.. autoclass:: msnoise.msnoise_table_def.Lineage

WorkflowStep
------------

.. autoclass:: msnoise.msnoise_table_def.WorkflowStep
   :members: get_config_params

WorkflowLink
------------

.. autoclass:: msnoise.msnoise_table_def.WorkflowLink

Filter
------

.. autoclass:: msnoise.msnoise_table_def.Filter
