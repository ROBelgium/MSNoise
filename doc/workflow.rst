.. _workflow:

Workflow
==========



.. image:: .static/Figure01_workflow_summary_cmyk.png


The general workflow of MSNoise is presented above. The left part of the workflow is executed only once, while the right part can be executed on a scheduled basis (cron job).

This section only presents the "installation" and configuration of MSNoise (read "the first startup of MSNoise"), not the installation of the required software, which is described in :ref:`installation`.


Installer
----------

.. automodule:: msnoise.s000installer



Configurator
-------------

.. automodule:: msnoise.s001configurator



Populate Station Table
-----------------------

.. automodule:: msnoise.s002populate_station_table



Scan Archive
---------------
.. automodule:: msnoise.s01scan_archive



New Jobs
-----------

.. automodule:: msnoise.s02new_jobs



Compute Cross-Correlations
----------------------------

.. automodule:: msnoise.s03compute_cc



Stack
--------
.. automodule:: msnoise.s04stack

	

Compute MWCS
-----------------

.. automodule::  msnoise.s05compute_mwcs




Compute dt/t
---------------

.. automodule::  msnoise.s06compute_dtt

