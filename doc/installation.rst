.. _installation:

******************************
Installation
******************************

MSNoise is a python package that uses a database (sqlite or MySQL) for storing
station and files metadata together with jobs. When installed, it provides a top
level command ``msnoise`` in the console.

Note that MSNoise is always tested against the latest release versions of the main packages, so older installations that are not maintained/updated regularly (years) could encounter issues. Please make sure you have the latest version of Numpy and Scipy (and MKL), as performance gets better and better (especially since Anaconda Inc. released its fast MKL implementations for all users, in the conda-forge channel).

To run MSNoise, you need:

* A recent version of Python (3.11 recommended). We suggest using Miniconda_
  and creating a fresh environment for running msnoise.
  MSNoise is tested "continuously" on GitHub for the last 3 most recent python version, and the three OS.


* Database: MariaDB or Postgresql: if you want to use a database,
  you need to install and configure a dabatase Server beforehand.
  This is not needed for sqlite. Read :ref:`aboutdbandperformances` for
  more information. We recommend using a database server when the number of
  stations/jobs is large, as parallelism of the workflow will allow you to
  process more data at once, compared with using sqlite.


Full Installation
-----------------

1. Download and install Miniconda_ for your machine, make sure Miniconda's Python
   is the default python for your user


2. Execute the following command to install the missing packages:
   
   .. code-block:: sh

        conda install -c conda-forge flask-admin flask-wtf markdown folium pymysql logbook pandas pytables pip xarray
        conda install -c conda-forge obspy
        conda install -c conda-forge msnoise


3. Prepare a portable database server:

   * The following instructions are for MariaDB similar to https://www.mariadb.education/install-portable, but can be replicated easily for postgresql.
   * Download the zip/tarball version of MariaDB portable (MariaDBs_)
   * Extract the zip/tarball in a folder and open a console into the folder:

   * Navigate to the bin/ directory (Windows) or scripts/ directory (Linux)
   * execute the ``mariadb-install-db.exe`` (Windows) or ``mariadb-install-db`` (Linux)
   * test the server by running the ``bin/mysqld --console``  command, you should see the server starting.
   * Keep the server running for now (later you'll use CTRL-C to kill the server)


4. Create a database:

   * In a new (!!) console (so keep the server running in the other console)
   * Go to the ``bin/`` directory and execute ``mysqladmin -u root  flush-privileges password "SECRET"`` where SECRET is a password (not very important to make it secure here)
   * test the connection with ``mysql -u root -p`` where you should be prompted for the password
   * the prompt should look like ``MariaDB [(none)]>`` now.
   * execute the command ``CREATE DATABASE msnoise;`` (the ; semicolumn is important)
   * List the databases with ``SHOW DATABASES;`` which should show the native MariaDB dbs, and our msnoise db.
   * From here, we can continue using the root user, or create a msnoise user, but this is not essential (see https://phoenixnap.com/kb/how-to-create-mariadb-user-grant-privileges) for instructions.
   * Close the connection by executing the ``quit;`` command.


5. Check which required packages you are still missing by executing the
   ``msnoise utils bugreport`` command. (See :ref:`testing`)


6. To be sure all is running OK, one could start the ``msnoise utils test`` command.
   This will start the standard MSNoise test suite, which should end with a
   "Ran xx tests in yy seconds : OK".


7. Proceed to the :ref:`Workflow` description to start MSNoise!

Done !


Database Structure - Tables
----------------------------
MSNoise will create the tables automatically upon running the installer script
(see :ref:`Workflow`).


Building this documentation
---------------------------

To build this documentation, some modules are required:

.. code-block:: sh

    conda install -c conda-forge "sphinx<6" sphinx_bootstrap_theme numpydoc sphinx-gallery
    pip install "sphinx_rtd_theme>1"
    pip install pillow==9.0.0


Then, this should simply work:

.. code-block:: sh

    make html
    
it will create a .build folder containing the documentation.

You can also build the doc to Latex and then use your favorite Latex-to-PDF
tool.


Using the development version
-----------------------------

This is not recommended, but users willing to test the latest development
(hopefully stable) version of MSNoise can:

.. code-block:: sh

    pip uninstall msnoise
    pip install http://msnoise.org/master.zip

Please note this version most probably uses the very latest version of every
package: Release versions of `numpy`, `scipy`, etc obtained from conda-forge
and "master" version of `obspy`. The development version (master) of obspy can
be installed from github: (warning regular Windows users, you might not be able to build the obspy package)

.. code-block:: sh

    pip uninstall obspy
    pip install https://github.com/obspy/obspy/archive/master.zip

If you are using the master version, please use the issue tracker of github to
communicate about bugs and not the mailing list, preferably used for Releases.


.. _obspy: http://www.obspy.org
.. _Miniconda: https://docs.anaconda.com/free/miniconda/#latest-miniconda-installer-links
.. _MariaDBs: https://mariadb.org/download/?t=mariadb&p=mariadb&r=10.11.6&os=windows&cpu=x86_64&pkg=zip&m=serverion