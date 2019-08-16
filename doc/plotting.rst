********
Plotting
********

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

Customizing Plots
-----------------

All plots commands can be overriden using a `-c` agument *in front of the
plot command* !!

Examples:

* ``msnoise -c plot distance``
* ``msnoise -c plot ccftime YA.UV02 YA.UV06 -m 5``
* etc.

To make this work, one has to copy the plot script from the msnoise install
directory to the project directory (where your db.ini file is located, then
edit it to one's desires. The first thing to edit in the code is the import of
the :doc:`../api`:

``from ..api import *``

to

``from msnoise.api import *``

and it should work.

.. versionadded:: 1.4


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


CCF's spectrum vs Time
----------------------

.. automodule:: msnoise.plots.spectime


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
