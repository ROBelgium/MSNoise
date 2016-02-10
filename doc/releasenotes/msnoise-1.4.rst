.. include:: ../configs.hrst

MSNoise 1.4
=============

Release date: XX XXXXXX 2016


Release type: major


Release notes:

.. contents::
    :local:

Introduction
------------
X months after the last bugfix release ( :doc:`msnoise-1.3.1`), and less than a
year after the last major release (:doc:`msnoise-1.3`) we are proud to announce
the new :doc:`msnoise-1.4`. It is a **major** release, with a massive amount of
work since the last release: in `GitHub numbers
<https://github.com/ROBelgium/MSNoise/graphs/contributors?from=2015-04-1&to=2016-03-01&type=c>`_
, it's over XXX commits and about XXX new lines of code and documentation added
! MSNoise 1.4 introduces two **major new features** : a new web-based admin
interface and the support for plugins ! + PWS !!!!


TESTS

This version has benefited from outputs/ideas/pull requests/questions from
several users:

***


Thanks to all for using MSNoise, and please, let us know why/how you use it
(and please cite it!)!

To date, we found/are aware of 11 publications using MSNoise ! That's the best
validation of our project ever ! See the full list on the
`MSNoise website <http://www.msnoise.org/they-cite-msnoise/>`_.


*Thomas Lecocq & Corentin Caudron*


new:
* config --set option

~~~~

PS: if you use MSNoise for your research and prepare publications, please
consider citing it:

**Lecocq, T., C. Caudron, et F. Brenguier (2014)**, MSNoise, a Python Package
for Monitoring Seismic Velocity Changes Using Ambient Seismic Noise,
*Seismological Research Letters*, 85(3), 715‑726, doi:10.1785/0220130073.


Web-based Admin Interface
-------------------------

For this release, we have replaced the Configurator by a more intuitive web-
based configuration interface. All fields present in the Configurator are
present, and *more* !



Pros:

* easier to customise
* modern
* less bandwidth when working remotely
* removes dependency for traits/traitui
* allows to customise "Views" to provide more information
* allows the validation of fields before saving to database
* will allow interactive plotting in the future

Cons:

* adds dependency to flask & flask-admin


Plugin support
--------------


MSNoise support plugin ...


Phase Weighted Stack
--------------------

...

Command Line changes
--------------------


* ``msnoise admin``: new command to start the web interface
* ``msnoise config``: accepts a ``--set name=value`` option, to rapidly
  change a configuration parameter in the database.
* The ``msnoise`` command accepts a ``-c`` option that triggers the "custom"
  mode, currently only for plots. See below.

All commands are now documented: :doc:`../clickhelp/msnoise`.


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


New plots
---------



* plot dtt



Math updates & bugfixes
-----------------------
Some improvements to the maths have been done for MSNoise 1.4:

* should we add Aurélien's less-agressive whiten ?


Performance improvements ------------------------ Improvements in terms of
performances have also been done for MSNoise 1.4:

* ``keep_all``: if set to ``Y`` (=True) in the config, all CCF are now stored in
* a single HDF5 file, which makes it much nicer to backup/transfer/delete. --
* XXX not used actually !!! ``compute_cc``: if only ZZ components are to be
* computed, the whitened windows are pre-computed, which makes the process
* faster. This could lead to memory issues if the job contains a lot of
* stations, a lot of filters are configured and a large number of windows. --
* Not sure it is a good idea !!!





Upgrading an existing project to MSNoise 1.4
--------------------------------------------

Some users will want to keep their current project without recomputing
everything. This requires adding a few configuration parameters to the database

Running the following command will take care of the upgrade from 1.3 to 1.4:

.. code-block:: sh

    msnoise upgrade_db


A note on parallel processing
-----------------------------
Although the ``msnoise`` command accepts the
``-t INTEGER`` argument to launch a number of threads in parallel, it currently
only works with ``scan_archive``: ``msnoise -t 4 scan_archive`` will run the
scan on four folders in parallel. For the other steps, one has still to run
multiple commands in a console. This should change in the future.


A final note about development pace and choices
-----------------------------------------------

* MSNoise team is

  * 1 developper (Thomas)
  * 1 dedicated debugger (Corentin)
  * less than 5 really *active* users, providing feedback and/or lines of codes
    (Esteban, Raphaël, Aurélien, Carmelo, ...)

* All software engineering ideas are coming from too infrequent beerstormings
  between Thomas & others
* The web-interface and the plugin support were developed during Thomas'
  holidays

If you need help, please ask your question on the mailing list. Don't be afraid
to ask. If you have ideas, please share them. If you develop codes to supplement
MSNoise, please share them, even if very small, even if you don't master gitHub.
If you have complaints, post them too, but remember that the package you are
using has been coded by 1 person, and that it's not his full time job. So
MSNoise is provided "as-is", carefully written and tested, but there will be
bugs, issues, incompatibility with certain python installations, OS or module
versions. If you **want or need** developments made and you can afford it,
contact Thomas via email directly, you can contract with the ROB for
paid-developments. If the developments you want are within the focus of the
developers' research, then a collaboration, i.e. resulting in a co-authored
peer reviewed publication, can be another option.
