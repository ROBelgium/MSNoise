.. _installation:

************
Installation
************

MSNoise requires Python ≥ 3.10 and a database backend (SQLite for quick
starts, PostgreSQL or MariaDB for production use).  This page covers the
recommended path — conda-forge + MSNoiseDB — and optional notes for
advanced users who prefer to manage their own database server.

.. contents::
   :local:
   :depth: 2


Step 1 — Install MSNoise
=========================

We recommend `Miniconda <https://docs.anaconda.com/free/miniconda/>`_ and a
dedicated conda environment.  All dependencies are available on
**conda-forge**.

Create the environment
-----------------------

Save the following as ``environment.yml`` — this is the exact file used by
MSNoise's own GitHub CI:

.. code-block:: yaml

    name: msnoise
    channels:
      - conda-forge
    dependencies:
      - numpy
      - scipy
      - obspy>=1.4.1
      - sqlalchemy
      - sqlalchemy-utils
      - flask
      - flask-admin<2
      - flask-wtf
      - flask-babel
      - markdown
      - wtforms<3.2
      - logbook
      - pytables
      - click
      - click-plugins
      - pandas
      - folium
      - pymysql
      - xarray
      - netCDF4
      - pooch
      - pip
      - pip:
          - git+https://github.com/regeirk/pycwt

Then create and activate it:

.. code-block:: sh

    conda env create -f environment.yml
    conda activate msnoise

Install MSNoise
---------------

The latest stable release from conda-forge:

.. code-block:: sh

    conda install -c conda-forge msnoise

Or directly from GitHub (development version, not recommended for production):

.. code-block:: sh

    pip install git+https://github.com/ROBelgium/MSNoise.git

Verify the installation:

.. code-block:: sh

    msnoise --version
    msnoise utils bugreport -s -m


Step 2 — Set up a database with MSNoiseDB
==========================================

Setting up MySQL or PostgreSQL from scratch has always been a friction point
for new users.  **MSNoiseDB** eliminates that: it bundles a self-contained,
user-run PostgreSQL server that requires no system privileges, no passwords,
and no server administration.

.. note::

   MSNoiseDB is the recommended database backend for all new MSNoise
   projects.  SQLite remains available for quick local tests (see
   :ref:`sqlite_option` below) but has limitations when running many
   parallel workers.

Install MSNoiseDB
-----------------

.. code-block:: sh

    pip install msnoisedb

or via conda-forge (once available):

.. code-block:: sh

    conda install -c conda-forge msnoisedb

See `github.com/ROBelgium/msnoise-db <https://github.com/ROBelgium/msnoise-db>`_
for the full documentation.

Start the database server
--------------------------

.. important::

   MSNoiseDB creates its database files **in the current working directory**
   the first time it starts.  Always run it from the **same dedicated folder**
   every time.  If you start it from a different directory, it will create a
   new, empty database there instead.

Choose or create a permanent home for MSNoiseDB — separate from your MSNoise
projects:

.. code-block:: sh

    mkdir ~/msnoisedb
    cd ~/msnoisedb
    msnoisedb start

The first run initialises the PostgreSQL cluster inside ``~/msnoisedb/``.
Subsequent runs just start the server.  The command prints the connection
details, for example:

.. code-block:: text

    MSNoiseDB started.
    Host : localhost
    Port : 5050

Keep this terminal open (or run it in the background) while MSNoise is
processing.  To stop the server:

.. code-block:: sh

    cd ~/msnoisedb
    msnoisedb stop

Create a new MSNoise project
-----------------------------

In a **separate** folder (your project folder, not the MSNoiseDB folder):

.. code-block:: sh

    mkdir ~/my_msnoise_project
    cd ~/my_msnoise_project
    msnoise db init

When prompted, select **postgresql** and use the host and port printed by
MSNoiseDB (``localhost:5050`` by default).  MSNoise will create a new
database automatically.  You are ready to go!

.. code-block:: sh

    msnoise admin   # opens the web configurator at http://localhost:5000


Step 3 — Verify
================

Run the fast smoke-test suite from your project folder to confirm everything
is wired up correctly:

.. code-block:: sh

    msnoise utils test --fast

All 32 tests should pass in under 30 seconds.  If any fail, run:

.. code-block:: sh

    msnoise utils bugreport -a

and check the output for missing dependencies.

Proceed to the :ref:`workflow_000` section to start your first MSNoise run!


.. _sqlite_option:

Advanced — SQLite (quick local tests only)
==========================================

SQLite requires no server setup and is useful for exploring MSNoise on a
laptop or running the test suite without MSNoiseDB:

.. code-block:: sh

    mkdir ~/msnoise_sqlite_test
    cd ~/msnoise_sqlite_test
    msnoise db init   # choose sqlite when prompted

.. warning::

   SQLite does not support concurrent writes.  Running more than one worker
   (``msnoise -t 2 cc compute``) on an SQLite project will cause database
   lock errors.  Use PostgreSQL via MSNoiseDB for any real processing.


Advanced — Self-managed MariaDB / PostgreSQL
=============================================

If your institution already runs a database server, or you need multi-user
access to a shared MSNoise database, you can connect directly.

MariaDB / MySQL
---------------

1. Create a database and user on your server:

   .. code-block:: sql

       CREATE DATABASE msnoise;
       CREATE USER 'msnoise'@'localhost' IDENTIFIED BY 'secret';
       GRANT ALL PRIVILEGES ON msnoise.* TO 'msnoise'@'localhost';
       FLUSH PRIVILEGES;

2. Run ``msnoise db init`` and choose **mysql** with your server's host,
   port, user and password.

PostgreSQL
----------

.. code-block:: sh

    createdb msnoise

Then run ``msnoise db init`` and choose **postgresql**.

See :ref:`aboutdbandperformances` for guidance on when a dedicated server
pays off over MSNoiseDB.


Building this documentation
============================

.. code-block:: sh

    conda install -c conda-forge sphinx sphinx-rtd-theme numpydoc sphinx-gallery
    cd doc/
    make html

The built documentation appears in ``doc/_build/html/``.

To include the gallery examples, define:

.. code-block:: sh

    export MSNOISE_DOC=/path/to/example/data


.. _Miniconda: https://docs.anaconda.com/free/miniconda/
