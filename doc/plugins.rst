.. include:: ../configs.hrst

Extending MSNoise with Plugins
==============================

.. versionadded:: 1.4

Starting with :doc:`releasenotes/msnoise-1.4`, MSNoise supports Plugins, this
means the default workflow "from archive to dv/v" can be branched at any step!


What is a Plugin and how to declare it in MSNoise
-------------------------------------------------

A plugin is a python package, properly structured, that can be imported from
msnoise, i.e. it has to be "installed" like any other python package.


After installing a plugin, its **package name** must be declared in the
``plugins`` parameter in the configuration. This must be done **PER PROJECT**.
This configuration field supports a list of plugins, separated by a simple comma
(!no space), e.g. ``msnoise_amazing,msnoise_plugin101``.


Once configured in a project, the plugin should appear when calling the
``msnoise plugin`` command:

.. code-block:: sh

    $ msnoise plugin

      Usage: msnoise-script.py plugin [OPTIONS] COMMAND [ARGS]...

      Runs a command in a named plugin

    Options:
      --help  Show this message and exit.

    Commands:
      amazing  Example Amazing Plugin for MSNoise


Plugin minimal structure
------------------------

A plugin is a python package, so its minimal structure is:

.. code-block:: sh

    msnoise-amazingplugin
    ├── __init__.py
    ├── setup.py
    └── msnoise_amazingplugin
        ├── __init__.py
        └── plugin_definition.py


The ``setup.py`` declares where the plugin actually hooks into MSNoise:

.. code-block:: python

    from setuptools import setup, find_packages

    setup(
        name='msnoise_amazing',
        version='0.1a',
        packages=find_packages(),
        include_package_data=True,
        install_requires=['msnoise',
                          'obspy'],
        entry_points = {
            'msnoise.plugins.commands': [
                'amazing = msnoise_amazing.plugin_definition:amazing',
                ],
            },
        author = "Thomas Lecocq & MSNoise dev team",
        author_email = "Thomas.Lecocq@seismology.be",
        description = "An example plugin",
        license = "EUPL-1.1",
        url = "http://www.msnoise.org",
        keywords="amazing seismology"
    )

The most important line of this file is the one declaring the ``amazing`` entry
point in ``msnoise.plugins.commands`` and linking it to the plugin's
``plugin_definition.py`` file.

The content of ``plugin_definition.py`` must then provide at least one
:class:`click.Command`, or more commonly, one :class:`click.Group` and
many :class:`click.Command`.

.. code-block:: python

    import click

    @click.group()
    def amazing():
        """Example Amazing Plugin for MSNoise"""
        pass

    @click.command()
    def sayhi():
        """A Very Polite Command"""
        print("Hi")

    amazing.add_command(sayhi)

This way, once properly installed and activated (declared in the ``plugins``
config), the plugin will be callable from msnoise:

.. code-block:: sh

    $ msnoise plugin amazing

      Usage: msnoise-script.py plugin amazing [OPTIONS] COMMAND [ARGS]...

      Example Amazing Plugin for MSNoise

    Options:
      --help  Show this message and exit.

    Commands:
      sayhi  A Very Polite Command

and its command too:

.. code-block:: sh

    $ msnoise plugin amazing sayhi

    Hi

Amazing, isn't it ?

Declaring Job Types - Hooking
-----------------------------

Plugin-based job types are defined by providing a ``register_job_types`` method
in ``plugin_definition.py``. A new job type is defined with two parameters:

* ``name``: the actual job name (acronym style) used all over (example: CC2, TEST)
* ``after``: when is this job added to the database.

Current supported "after" are:

* ``new_files``: will be created when running the ``new_jobs`` command and will
  create a job with those parameters (nf is a new file identified in the
  scan_archive procedure). In this specific case, the ``pair`` field of the job
  will only be NET.STA, not a "pair". A job will only be inserted if the station
  is "Used" in the configuration.

  .. code-block:: python

        all_jobs.append({"day": current_date, "pair": "%s.%s"%(nf.net,nf.sta),
                         "jobtype": jobtype, "flag": "T",
                         "lastmod": datetime.datetime.utcnow()})

* ``scan_archive``: will be created when running the ``new_jobs`` command, in
  parallel to ``CC`` jobs. This is, for example, useful when one wants to
  compute relative amplitude ratios between station pairs. In this case, the
  ``pair`` field of the job is set to the pair name.

* ``refstack``: will be created when running the ``stack`` command and when
  a new REF stack needed to be calculated. This is, for example, useful when one
  wants to work on the REF stacks using a Ambient Seismic Noise Tomography code.

Plugin's Job Types are first declared in ``setup.py`` (in Entry Points):

.. code-block:: python

    'msnoise.plugins.jobtypes': [
    'register = msnoise_amazing.plugin_definition:register_job_types',
    ],

.. code-block:: python

    def register_job_types():
        jobtypes = []
        jobtypes.append( {"name":"AMAZ1", "after":"new_files"} )
        return jobtypes


Then, adding a ``compute`` command to the ``plugin_definition.py``:

.. code-block:: python

    @click.command()
    def compute():
        """Compute an Amazing Value"""
        from .compute import main()
        main()

    amazing.add_command(compute)

and creating a  ``compute.py`` file in the plugin folder:

.. code-block:: python

    import os
    from obspy.core import UTCDateTime, read

    from msnoise.api import connect, is_next_job, get_next_job, \
        get_data_availability, get_config, update_job

    def main():
        db = connect()
        while is_next_job(db, jobtype='AMAZ1'):
            jobs = get_next_job(db, jobtype='AMAZ1')
            for job in jobs:
                net, sta = job.pair.split('.')
                gd = UTCDateTime(job.day).datetime
                print("Processing %s.%s for day %s"%(net,sta, job.day))
                files = get_data_availability(
                        db, net=net, sta=sta, starttime=gd, endtime=gd,
                        comp="Z")
                for file in files:
                    fn = os.path.join(file.path, file.file)
                    st = read(fn, starttime=UTCDateTime(job.day), endtime=UTCDateTime(job.day)+86400)
                    print(st)

Aaaand:

.. code-block:: sh

    $ msnoise plugin amazing compute

    Processing YA.UV05 for day 2010-09-01
    1 Trace(s) in Stream:
    YA.UV05.00.HHZ | 2010-09-01T00:00:00.000000Z - 2010-09-01T23:59:59.990000Z | 100.0 Hz, 8640000 samples
    Processing YA.UV06 for day 2010-09-01
    1 Trace(s) in Stream:
    YA.UV06.00.HHZ | 2010-09-01T00:00:00.000000Z - 2010-09-01T23:59:59.990000Z | 100.0 Hz, 8640000 samples
    Processing YA.UV10 for day 2010-09-01
    1 Trace(s) in Stream:
    YA.UV10.00.HHZ | 2010-09-01T00:00:00.000000Z - 2010-09-01T23:59:59.990000Z | 100.0 Hz, 8640000 samples

Provided you have reset the DataAvailability rows with a "M" or "N" flag so that
when you ran ``new_jobs`` it actually inserted the ``AMAZ1`` jobs !

Because job-based stuff always requires a lot of trial-and-error, rememeber that
the ``msnoise reset`` command is your best friend. In this example, we would
need to ``msnoise reset AMAZ1`` to reset "I"n Progress jobs, or
``msnoise reset AMAZ1 --all`` to reset all ``AMAZ1`` jobs to "T"o Do.

.. note::

    * Currently, not all MSNoise workflow steps use the ``is_next_job`` -
      ``get_next_job`` logic, but it'll be the case for MSNoise 1.5

    * Only three hooks are currently present, of course, more will be added in
      in the future.
