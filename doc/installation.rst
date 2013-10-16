.. _installation:


Installation
=============

Introduction
------------

MSNoise is a set of Python codes that use a database (sqlite or MySQL) and
the `find` command. 

To run MSNoise, you need:

*   A recent version of Python (2.7.x recommended). We suggest using Anaconda_ with extra modules ([+] modules are already distributed with Anaconda_):
    
    * numpy [+]
    * scipy [+]
    * pandas [+]
    * matplotlib [+]
    * statsmodels [+]
    * SqlAlchemy [+]
    * Enthought Tool Suite
    * scikits.samplerate
    * obspy

*   MySQL: if you want to use MySQL, you need to install and configure a :ref:`mysql` beforehand. This is not needed for sqlite.
    Read :ref:`aboutdbandperformances` for more information. We recommend using MySQL.

*   The `find` command: present by default on linux and available with gnufind_ on Windows.


Fast way
---------

1. Download and install Anaconda_ for your machine, make sure Anaconda's Python is the default python for your user
2. Execute the following command:
   
   .. code-block:: sh
    
        easy_install traitsui traits obspy
        easy_install scikits.samplerate

   Note: on Windows the latter will probably not work! See :ref:`samplerate`.

3. Install a MySQL server:
   
   a) Linux:
      
      .. code-block:: sh
        
            sudo apt-get install mysql-server mysql-client
   
   b) Windows: Install EasyPHP_.

   c) Create a privileged user and a database:
      
      TODO

4. Windows-only: install gnufind_ and make sure its /bin directory is in the PATH (Control Panel -> Environment Variables -> PATH)

5. Download MSNoise from `GitHub <https://github.com/ROBelgium/MSNoise>`_.

6. Proceed to the :ref:`Workflow` description to start MSNoise!

Done !

....

Python Installation
-------------------

If you don't know which Python distribution to use and even if your system comes
with a python distribution, we suggest installing Anaconda_, as it comes with most of the
above-mentionned tools (those with [*]), and provides the easy_install tool
to install the remaining ones.

Required Packages (vs Anaconda)
-----------------------------------

We suppose you installed Anaconda_, here are the instructions for installing the remaining packages. If you don't use Anaconda, all the packages are available through 'easy_install'.
Windows users are recommended to check the prebuilt binaries when advised.

Obspy
~~~~~

http://www.obspy.org (Beyreuther et al., 2010; Megies et al., 2011)

.. code-block:: sh

	easy_install obspy

Enthought Tools Suite
~~~~~~~~~~~~~~~~~~~~~

Most of the suite should be present, one just needs to install the traitsui package and its dependencies (traits, pyface, 
), which easy_install will do for you:

.. code-block:: sh

	easy_install traitsui

.. _samplerate:

scikits.samplerate
~~~~~~~~~~~~~~~~~~
https://pypi.python.org/pypi/scikits.samplerate is a wrapper to the Secret Rabbit Code (aka libsamplerate) (de Castro Lopo, 2013)

Windows
++++++++

Download and install the right version from here: http://www.lfd.uci.edu/~gohlke/pythonlibs/#scikits.samplerate

Linux
+++++++

You first need to install the SRC library:

.. code-block:: sh

	sudo apt-get install libsamplerate0 libsamplerate0-dev

This python package will probably be the most tricky to install. If you are lucky, you can just

.. code-block:: sh

	easy_install scikits.samplerate

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


SQLAlchemy
~~~~~~~~~~
Windows
++++++++
Download and install the right version from here: http://www.lfd.uci.edu/~gohlke/pythonlibs/#sqlalchemy


Linux:
+++++++

.. code-block:: sh

	easy_install sqlalchemy

.. _mysql:

MySQL Server
-------------
.. warning:: MySQL is not compulsory, one *can* work only using sqlite database. See :ref:`aboutdbandperformances`. for more info.
MSNoise requires a database in order to store waveform metadata, configuration bits and jobs.
If you choose to use MySQL, a running MySQL server must be available, either on the network or on localhost and have a privileged user and a database.

Windows
~~~~~~~~~~
The simplest option to install a MySQL server on your machine is to install EasyPHP_, a small AMP (Apache, MySQL, PHP) server.

Linux
~~~~~~~~~~

If you don't have a MySQL server on the network, you need to install one locally on your computer.
MySQL is usually prepackaged for every distribution, on Ubuntu/Debian you should:

.. code-block:: sh

	sudo apt-get install mysql-server mysql-client

We recommend to install phpmyadmin too, as it is a handy tool to edit the database directly

.. code-block:: sh

	sudo apt-get install phpmyadmin

This will also install apache2 and php, needed to run phpmyadmin. Once installed, it should be available through http://localhost/phpmyadmin.


Database Structure - Tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MSNoise will create the tables automatically upon running the installer script.


Find command for Windows
-------------------------

One has to install gnufind_ in order to be able to search for recently modified files in the data archive.


Testing the Dependencies
-------------------------

Once installed, you should be able to import the python packages in a python console. 
For testing purposes, it might be a good idea to install IP[y]thon as the latest
version provides the fantastic notebook, which will allow you to test the different
functions/calls of MSNoise interactively.


Building this documentation
-----------------------------

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