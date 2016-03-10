.. _installation:


Installation
============


.. contents::
    :local:


Introduction
------------

MSNoise is a set of Python codes that use a database (sqlite or MySQL) and
the `find` command. When installed, it provides a top level command ``msnoise``
in the console.

.. note:: MSNoise version 1.4 introduces new config parameters in the database.
    For each project, to upgrade your database, remember to run
    ``msnoise upgrade_db`` in the project folder.


To run MSNoise, you need:

* A recent version of Python (2.7.x recommended). We suggest using Anaconda_
  with extra modules ([+] modules are already distributed with Anaconda_):

  * setuptools [+]
  * numpy [+]
  * scipy [+]
  * pandas [+]
  * matplotlib [+]
  * statsmodels [+]
  * sqlalchemy [+]
  * flask [+] (new in 1.4)
  * click [+]
  * pymysql [+] (new in 1.4)
  * scikits.samplerate
  * obspy
  * flask-admin (new in 1.4)
  * multiprocessing_logging (new in 1.4)
  * markdown (new in 1.4)


* MySQL: if you want to use MySQL, you need to install and configure a
  :ref:`mysql` beforehand. This is not needed for sqlite.
  Read :ref:`aboutdbandperformances` for more information.
  We recommend using MySQL.

* The `find` command: present by default on linux and available with gnufind_
  on Windows.


Quick Start - Windows
----------------------

1. Download and install Anaconda_ for your machine, make sure Anaconda's Python is the default python for your user

2. Execute the following command to install the missing packages:
   
   .. code-block:: sh
    
        pip install flask-admin
        conda install -c obspy obspy

3. Download and install scikits.samplerate from here: http://www.lfd.uci.edu/~gohlke/pythonlibs/#scikits.samplerate
   see examples below (filenames might be different):

   .. code-block:: sh

        pip install scikits.samplerate-0.3.3-cp27-none-win_amd64.whl

4. Install a MySQL server: download and install EasyPHP_ (their "Dev" version is OK - read their installation page carefully, you need to install some Microsoft Visual C redistribuable manually too).

5. Start EasyPHP and create a privileged user and a database:
      
   * Connect to your local host: http://localhost/phpmyadmin (or http://127.0.0.1/phpmyadmin)
   * Click on "Privileges" and create a new user, with all privileges (Select all). Ideally, create user "msnoise" with password "msnoise".

6. Install gnufind_ and make sure its /bin directory is in the PATH (Control Panel -> Environment Variables -> PATH)

7. Install MSNoise:

   .. code-block:: sh

        pip install msnoise

8. Check which required packages you are still missing by executing the ``msnoise bugreport`` command. (See :ref:`testing`)

9. Proceed to the :ref:`Workflow` description to start MSNoise!


Done !


Quick Start - Linux
-------------------

1. Download and install Anaconda_ for your machine, make sure Anaconda's Python is the default python for your user

2. Execute the following commands to install the missing packages:
   
   .. code-block:: sh
    
        pip install flask-admin
        conda install -c obspy obspy
 
   .. code-block:: sh
        
        sudo apt-get install libsamplerate0 libsamplerate0-dev
        pip install scikits.samplerate
    
   If this fails, follow those instructions: :ref:`samplerate`.

3. Install a MySQL server and phpMyAdmin:
   
   .. code-block:: sh
    
        sudo apt-get install mysql-server mysql-client phpmyadmin


5. Create a privileged user and a database:
 
   * Connect to your local host: http://localhost/phpmyadmin (or http://127.0.0.1/phpmyadmin)
   * Click on "Privileges" and create a new user, with all privileges (Select all). Ideally, create user "msnoise" with password "msnoise".

6. Install MSNoise:

   .. code-block:: sh

        pip install msnoise

7. Check which required packages you are still missing by executing the ``msnoise bugreport`` command. (See :ref:`testing`)

8. Proceed to the :ref:`Workflow` description to start MSNoise!

Done !

UltraFast Start on Linux
------------------------
If one starts with a vanilla fresh Linux install (e.g. on a new virtual machine)
, the install can be eased with an installer script we have prepared. Indeed,
to run de tests on TravisCI, we had to prepare a pre-install script. This is
only valid for linux x86_64 (Debian or Ubuntu):

.. code-block:: sh

    wget https://raw.githubusercontent.com/ROBelgium/MSNoise/master/misc/install_debian_x86_64.sh
    chmod +x install_debian_x86_64.sh
    ./install_debian_x86_64.sh
    pip install msnoise

Done !


Python and Packages Installation
--------------------------------

If you don't know which Python distribution to use and even if your system comes
with a python distribution, we suggest installing Anaconda_, as it comes with most of the
above-mentionned tools (those with [*]), and provides the easy_install tool
to install the remaining ones.

From now on, we suppose you installed Anaconda_, here are the instructions for installing
the remaining packages. If you don't use Anaconda, all the packages are available through 'easy_install'.
Windows users are recommended to check the prebuilt binaries when advised.

To know which packages you are missing, use the bug_reporter script (see :ref:`troubleshooting`) !

.. _samplerate:

scikits.samplerate
~~~~~~~~~~~~~~~~~~

.. warning:: Scikits.samplerate isn't easy to install on Python3. If not
    successful, one has to use the "Decimate" resampling method.

https://pypi.python.org/pypi/scikits.samplerate is a wrapper to the Secret Rabbit Code (aka libsamplerate) (de Castro Lopo, 2013)

Windows
+++++++

Download and install the right version from here: http://www.lfd.uci.edu/~gohlke/pythonlibs/#scikits.samplerate

Linux
+++++

You first need to install the SRC library:

.. code-block:: sh

    sudo apt-get install libsamplerate0 libsamplerate0-dev

This python package will probably be the most tricky to install. If you are lucky, you can just

.. code-block:: sh

    pip install scikits.samplerate

On my Ubuntu 12.04, this results in an error because the SRC library path is not found. The reason is that the setup searches SRC in /usr/lib and not in /usr/lib/x86_64-linux-gnu where the library is actually present. To install, you need to download the archive from pypi and edit some configuration file:

.. code-block:: sh

    wget https://pypi.python.org/packages/source/s/scikits.samplerate/scikits.samplerate-0.3.3.tar.gz#md5=96c8d8ba3aa95a9db15994f78792efb4
    tar -xvf scikits.samplerate-0.3.3.tar.gz
    cd scikits.samplerate-0.3.3

then edit the site.cfg example file and insert the following lines:

.. code-block:: sh

    [samplerate]
    library_dirs=/usr/lib/x86_64-linux-gnu
    include_dirs=/usr/include

To know where the SRC library is on you machine:

.. code-block:: sh

    sudo dpkg -L libsamplerate0
    sudo dpkg -L libsamplerate0-dev

then, build and install:

.. code-block:: sh

    python setup.py build
    python setup.py install


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

    easy_install sphinx
    easy_install sphinx_bootstrap_theme
    
Then, this should simply work:

.. code-block:: sh

    make html
    
it will create a .build folder containing the documentation.

You can also build the doc to Latex and then use your favorite Latex-to-PDF tool.

.. _gnufind: http://sourceforge.net/projects/getgnuwin32/files/
.. _EasyPHP: http://www.easyphp.org/
.. _obspy: http://www.obspy.org
.. _Anaconda: http://www.continuum.io/downloads