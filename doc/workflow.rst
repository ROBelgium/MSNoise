.. _workflow:

Workflow
==========



.. image:: .static/Figure01_workflow_summary_cmyk.png


The general workflow of MSNoise is presented above. The left part of the workflow is executed only once, while the right part can be executed on a scheduled basis (cron job).

This section only presents the "installation" and configuration of MSNoise (read "the first startup of MSNoise"), not the installation of the required software, which is described in :ref:`installation`.


Installer
----------

.. note:: TODO

....

Configurator
-------------

.. automodule:: MSNoise.s001configurator

....

Populate Station Table
-----------------------

.. automodule:: MSNoise.s002populate_station_table

....

Scan Archive
---------------
.. automodule:: MSNoise.s01scan_archive
	:members:

....

New Jobs
-----------

.. automodule:: MSNoise.s02new_jobs

....

Compute Cross-Correlations
----------------------------

.. automodule:: MSNoise.s03compute_cc

....

Stack
--------
.. automodule:: MSNoise.s04stack
	:members:
	
....	

Compute MWCS
-----------------

.. automodule::  MSNoise.s05compute_mwcs


....

Compute dt/t
---------------

.. automodule::  MSNoise.s06compute_dtt

....

Plot dt/t
----------