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
  * flask
  * pymysql
  * wtforms

* Not shipped with Anaconda_:

  * obspy
  * flask-admin
  * multiprocessing_logging
  * markdown
  * folium


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
above-mentionned tools (those with [*]), and provides the `pip` or `conda` tools
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
        conda install -c conda-forge obspy

3. Install a MySQL server and MySQL Workbench:

   Download and install MySQL Community Server (MySQLs_ ) and MySQL Workbench (MySQLw_ ) ; On Windows one can also use the MySQL installer (MySQLi_ ).

   On Linux, the MySQL server can also be installed using the following command:

        .. code-block:: sh

            sudo apt-get install mysql-server

4. Create a privileged user and a database:

   * Start MySQL Workbench and connect to the local database
   * Click on "Privileges" and create a new user, with all privileges (Select all). Ideally, create user "msnoise" with password "msnoise".

5. WINDOWS ONLY: Install gow_ and make sure its /bin directory is in the beginning of the System PATH (Control Panel -> Environment Variables -> Machine -> PATH).
   After that, the find command should be available in cmd.exe. If this fails, you will have to provide the full path to find.exe in the configuration (``find_command``).

6. Install the latest release version of MSNoise:

   .. code-block:: sh

        pip install msnoise

   Power user could install the development version too, but it is not recommended.

7. Check which required packages you are still missing by executing the ``msnoise bugreport`` command. (See :ref:`testing`)

8. To be sure all is running OK, one could start the ``msnoise test`` command in an empty directory.
   This will start the standard MSNoise test suite, which should end with a "Ran xx tests in yy seconds : OK".

8. Proceed to the :ref:`Workflow` description to start MSNoise!


Done !

MySQL Server and Workbench
--------------------------

Using the MySQL Server and Workbench is fairly easy and lots of tutorials are available online as text or videos.

Once both are installed, start Workbench and you should see the local MySQL server automatically identified:

.. image:: .static/workbench_1.png

And by clicking on "Local Instance ..." another tab should open, connected to the local database.

Create a msnoise user
~~~~~~~~~~~~~~~~~~~~~

Select "Users and Privileges" in the left sidebar, then "Add Account". Define the username and the password (msnoise:msnoise could do, although "weak"):

.. image:: .static/workbench_2.png

Then, under "Administrative Roles", grant this user the *DBA* mode (user can perform all tasks on the database server) and click "Apply".

.. image:: .static/workbench_3.png

Create an empty database
~~~~~~~~~~~~~~~~~~~~~~~~

Each "project" needs a database. That is, if one has two different volcanoes and wants to run MSNoise the two datasets, one needs to create two empty databases.

Click on the "Create new schema" button in the taskbar:

.. image:: .static/workbench_4.png

and give the database a name (for example msnoise; or msnoise_project1, or project1, or else, you choose) ; and click "Apply":

.. image:: .static/workbench_5.png

and click "Apply" again and it should state all is OK:

.. image:: .static/workbench_6.png

.. image:: .static/workbench_7.png

When done, the database we created is present in the left sidebar:

.. image:: .static/workbench_8.png

And you're ready to start your first project: :ref:`Workflow`.


When moving your project to a larger server, HPC or else, just add the
connection to this server in Workbench and you're good to go with the very
same interface/tool !


Database Structure - Tables
----------------------------
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
.. _MySQLi: https://dev.mysql.com/downloads/installer
.. _MySQLs: https://dev.mysql.com/downloads/mysql
.. _MySQLw: https://dev.mysql.com/downloads/workbench
