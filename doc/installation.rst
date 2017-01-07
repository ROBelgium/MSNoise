.. _installation:


Installation
------------

MSNoise is a set of Python codes that use a database (sqlite or MySQL) and
the `find` command. When installed, it provides a top level command ``msnoise``
in the console.

.. note:: MSNoise version 1.4 introduces new config parameters in the database.
    For each project, to upgrade your database, remember to run
    ``msnoise upgrade_db`` in the project folder.


To run MSNoise, you need:

* A recent version of Python (3.x recommended). We suggest using Anaconda_
  with extra modules (those modules are already distributed with Anaconda_):

  * setuptools
  * numpy
  * scipy
  * pandas
  * matplotlib
  * statsmodels
  * sqlalchemy
  * click
  * flask (new in 1.4)
  * pymysql (new in 1.4)

* Not shipped with Anaconda_:

  * obspy
  * flask-admin (new in 1.4)
  * multiprocessing_logging (new in 1.4)
  * markdown (new in 1.4)
  * folium (new in 1.4)


* MySQL: if you want to use MySQL, you need to install and configure a
  MySQL Server beforehand. This is not needed for sqlite.
  Read :ref:`aboutdbandperformances` for more information.
  We recommend using MySQL.

* The `find` command: present by default on linux and available with gow_
  on Windows.


Python and Packages Installation
--------------------------------

If you don't know which Python distribution to use and even if your system comes
with a python distribution, we suggest installing Anaconda_, as it comes with most of the
above-mentionned tools (those with [*]), and provides the easy_install tool
to install the remaining ones.

From now on, we suppose you installed Anaconda_, here are the instructions for installing
the remaining packages. If you don't use Anaconda, all the packages are available through `pip` or `conda`.

To know which packages you are missing, use the bug_reporter script!


Full Installation
-----------------

1. Download and install Anaconda_ for your machine, make sure Anaconda's Python is the default python for your user

2. Execute the following command to install the missing packages:
   
   .. code-block:: sh
    
        pip install flask-admin
        conda install -c obspy obspy

3. Install a MySQL server and phpMyAdmin:

   * WINDOWS: Download and install EasyPHP_ (their "Dev" version is OK - read their installation page carefully, you need to install some Microsoft Visual C redistribuable manually too). Then start EasyPHP.
   * LINUX:

        .. code-block:: sh

            sudo apt-get install mysql-server mysql-client phpmyadmin

4. Create a privileged user and a database:
      
   * Connect to your local host: http://localhost/phpmyadmin (or http://127.0.0.1/phpmyadmin)
   * Click on "Privileges" and create a new user, with all privileges (Select all). Ideally, create user "msnoise" with password "msnoise".

5. WINDOWS ONLY: Install gow_ and make sure its /bin directory is in the PATH (Control Panel -> Environment Variables -> PATH)

6. Install MSNoise:

   .. code-block:: sh

        pip install msnoise

7. Check which required packages you are still missing by executing the ``msnoise bugreport`` command. (See :ref:`testing`)

8. Proceed to the :ref:`Workflow` description to start MSNoise!


Done !




MySQL Server
------------
.. warning:: MySQL is not compulsory, one *can* work only using sqlite database. See :ref:`aboutdbandperformances`. for more info.

MSNoise requires a database in order to store waveform metadata, configuration bits and jobs.
If you choose to use MySQL, a running MySQL server must be available, either on the network or on localhost and have a privileged user and a database.

Windows
~~~~~~~
The simplest option to install a MySQL server on your machine is to install EasyPHP_, a small AMP (Apache, MySQL, PHP) server.

Linux
~~~~~

If you don't have a MySQL server on the network, you need to install one locally on your computer.
MySQL is usually prepackaged for every distribution, on Ubuntu/Debian you should:

.. code-block:: sh

    sudo apt-get install mysql-server mysql-client

We recommend to install phpmyadmin too, as it is a handy tool to edit the database directly

.. code-block:: sh

    sudo apt-get install phpmyadmin

This will also install apache2 and php, needed to run phpmyadmin. Once installed, it should be available through http://localhost/phpmyadmin.


Database Structure - Tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MSNoise will create the tables automatically upon running the installer script (see :ref:`Workflow`).


Building this documentation
---------------------------

To build this documentation, some modules are required:

.. code-block:: sh

    pip install sphinx
    pip install sphinx_bootstrap_theme
    
Then, this should simply work:

.. code-block:: sh

    make html
    
it will create a .build folder containing the documentation.

You can also build the doc to Latex and then use your favorite Latex-to-PDF tool.

.. _gow: https://github.com/downloads/bmatzelle/gow/Gow-0.7.0.exe
.. _EasyPHP: http://www.easyphp.org/
.. _obspy: http://www.obspy.org
.. _Anaconda: http://www.continuum.io/downloads