Plotting
========
MSNoise comes with some default plotting tools.

All plotting commands accept the ``--outfile`` argument. If provided, the
figure will be saved to the disk. Names can be explicit, or tell the code to
generate the filename automatically (using the `?` question mark), for example:

.. code-block:: sh

    # automatic naming, save to PNG
    msnoise plot dvv -o ?.png

    # automatic naming, save to PDF
    msnoise plot dvv -o ?.pdf

    # explicit naming, save to JPG
    msnoise plot dvv -o mydvv.jpg

.. contents::
    :local:


Data Availability Plot
----------------------

.. automodule:: msnoise.plots.data_availability


Station Map
-----------

.. automodule:: msnoise.plots.station_map


Interferogram Plot
------------------

.. automodule:: msnoise.plots.interferogram


CCF vs Time
-----------

.. automodule:: msnoise.plots.ccftime


MWCS Plot
---------

.. automodule:: msnoise.plots.mwcs


Distance Plot
-------------

.. automodule:: msnoise.plots.distance


dv/v Plot
---------

.. automodule:: msnoise.plots.dvv

dt/t Plot
---------

.. automodule:: msnoise.plots.dtt
