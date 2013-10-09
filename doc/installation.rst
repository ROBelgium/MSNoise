.. _installation:

Installation
=============

In order to run MSNoise, you need to have Python 2.7 installed on a Windows/Linux/MacOSX machine, with the following packages:

*    numpy and scipy: http://www.scipy.org (Jones et al., 2001; Oliphant, 2006) [*]
*    obspy: http://www.obspy.org (Beyreuther et al., 2010; Megies et al., 2011)
*    scikits.samplerate: https://pypi.python.org/pypi/scikits.samplerate which is a wrapper to the Secret Rabbit Code (aka libsamplerate) (de Castro Lopo, 2013)
*    Enthought Tool Suite: http://code.enthought.com
*    Matplotlib: http://www.matplotlib.org (Hunter, 2007) [*]
*    Basemap: http://matplotlib.org/basemap [*]
*    Pandas: http://pandas.pydata.org (McKinney, 2012) [*]
*    statsmodels http://statsmodels.sourceforge.net [*]
*    MySQL-python https://pypi.python.org/pypi/MySQL-python


Python Installation
-------------------

If you don't know which Python distribution to use and even if your system comes
with a python distribution, we suggest installing Anaconda 
(http://continuum.io/downloads), as it comes with most of the
above-mentionned tools (those with [*]), and provides the easy_install tool
to install the remaining ones.

Required Packages (vs Anaconda)
-----------------------------------

We suppose you installed Anaconda, here are the instructions for installing the remaining packages. If you don't use Anaconda, all the packages are available through 'easy_install'.
Windows users are recommended to check the prebuilt binaries when advised.

Obspy
~~~~~
.. code-block:: sh

	easy_install obspy

Enthought Tools Suite
~~~~~~~~~~~~~~~~~~~~~

Most of the suite should be present, one just needs to install the traitsui package and its dependencies (traits, pyface, 
), which easy_install will do for you:

.. code-block:: sh

	easy_install traitsui


scikits.samplerate
~~~~~~~~~~~~~~~~~~

.. note:: Windows users: http://www.lfd.uci.edu/~gohlke/pythonlibs/#scikits.samplerate

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

.. note:: Windows users: http://www.lfd.uci.edu/~gohlke/pythonlibs/#sqlalchemy

.. code-block:: sh

	easy_install sqlalchemy

MySQL-python
~~~~~~~~~~~~
.. warning:: MySQL-python was needed in the original MSNoise release. The current version (currently developped in the sqlvolution branch in GitHub) 
	uses SQLAlchemy to connect to a MySQL server.

.. note:: Windows users: http://www.lfd.uci.edu/~gohlke/pythonlibs/#mysql-python
	
Install the MySQL server first ! (see below)

you need to install the dependencies of the module and then the module itself:

.. code-block:: sh

	sudo apt-get build-dep python-mysqldb
	sudo apt-get install libmysqlclient-dev
	easy_install mysql-python

MySQL Database
-------------------
.. warning:: MySQL is not compulsory, one *can* work only using sqlite database. See :ref:`aboutdbandperformances`. for more info.

MSNoise requires MySQL database in order to store waveform metadata, configuration bits and jobs. A running MySQL database must be available, either on the network or on localhost, and you must somehow have the rights
to create a new database.

If you don't have a MySQL server on the network, you need to install one locally on your computer. MySQL is usually prepackaged for every distribution, on Ubuntu/Debian you should:

.. code-block:: sh

	sudo apt-get install mysql-server

We recommend to install phpmyadmin too, as it is a handy tool to edit the database directly

.. code-block:: sh

	sudo apt-get install phpmyadmin

This will also install apach2 and php, needed to run phpmyadmin. Once installed, it should be available through http://localhost/phpmyadmin.



Database Structure - Tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. warning:: In the sqlvolution branch, those tables are created on the fly !

MSNoise will create the tables automatically upon running the installer script. If you prefer asking your network manager to install it, provide him the following SQL:

Windows - Prebuilt Binaries
---------------------------

When working on Windows, you can either install the packages 'from source', as above, or install them using pre-compiled binaries prepared by Christoph Gohlke.

Find command for Windows
~~~~~~~~~~~~~~~~~~~~~~~~

One has to install gnufind in order to be able to search for recently modified files in the data archive.

.. note::
    link here

MySQL  Apache  PhpMyAdmin
~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest option to install a MySQL server on your machine is to install EasyPHP, a small AMP (Apache, MySQL, PHP) server.

Testing the Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~

Once installed, you should be able to import the python packages in a python console. 
For testing purposes, it might be a good idea to install IP[y]thon as the latest
version provides the fantastic notebook, which will allow you to test the different
functions/calls of MSNoise interactively.


