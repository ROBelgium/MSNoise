.. include:: ../configs.hrst

.. _populate:

Managing Stations
=================

In MSNoise 2.x, stations and their metadata are managed via the
:ref:`webadmin` and the ``msnoise utils import-stationxml`` command.

Importing station metadata
--------------------------

The preferred way to add stations is to import a StationXML file or
query an FDSN station web service:

.. code-block:: sh

    # From a local inventory file
    msnoise utils import-stationxml inventory.xml

    # From an FDSN URL
    msnoise utils import-stationxml https://eida.ethz.ch/fdsnws/station/1/query?...

The inventory is saved to the project's ``|global.response_path|`` directory
automatically, making instrument responses available to the preprocessing step.

Assigning data sources
-----------------------

Each station can be assigned to a ``DataSource`` (local SDS archive or
remote FDSN/EIDA service).  See :ref:`workflow_concepts` for the full
:ref:`concepts_configsets` explanation.

.. code-block:: sh

    msnoise admin   # → Data Sources → Stations
