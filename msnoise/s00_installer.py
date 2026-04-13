"""
Initialize a MSNoise Project
=============================

This script initialises a new MSNoise project in the **current directory**.
It creates the ``db.ini`` connection file and sets up all database tables with
their default configuration values.

Run it once per project with:

.. code-block:: sh

    msnoise db init

.. important::

   Always run ``msnoise db init`` (and all subsequent ``msnoise`` commands)
   from inside the **project folder** — MSNoise looks for ``db.ini`` there.


Prerequisites — MSNoiseDB (recommended)
-----------------------------------------

The easiest way to get a database running is **MSNoiseDB**, a self-contained
PostgreSQL server that requires no system privileges or server administration.

Install it once:

.. code-block:: sh

    pip install msnoisedb

Then, from its **own dedicated folder** (not your project folder):

.. code-block:: sh

    mkdir ~/msnoisedb
    cd ~/msnoisedb
    msnoisedb start

.. important::

   MSNoiseDB stores the PostgreSQL cluster files in the directory where it is
   started.  **Always start it from the same folder.**  Starting from a
   different directory creates a new, empty cluster there instead.

The server prints its connection details on startup:

.. code-block:: text

    MSNoiseDB started.
    Host : localhost
    Port : 5099

Keep this terminal open (or background the process) while MSNoise is running.
To stop the server:

.. code-block:: sh

    cd ~/msnoisedb
    msnoisedb stop


Create a database for this project
-------------------------------------

Each project needs its own PostgreSQL database.  Create one with:

.. code-block:: sh

    msnoisedb create-db msnoise_myproject

Replace ``msnoise_myproject`` with any name you like.  Verify it exists:

.. code-block:: sh

    msnoisedb list-db

This step is needed **once per project** — the database persists inside
``~/msnoisedb/`` across server restarts.


Initialise the project
------------------------

Create and enter a fresh project folder, then run the initialiser:

.. code-block:: sh

    mkdir ~/my_msnoise_project
    cd ~/my_msnoise_project
    msnoise db init

The script asks a few questions.  For a typical MSNoiseDB setup:

.. code-block:: text

    What database technology do you want to use?
     [1] sqlite
     [2] mysql
     [3] postgresql
    Choice: 3
    Hostname [localhost:5099] (hostname or hostname:port, e.g. localhost:5099 for MSNoiseDB): localhost:5099
    Database name [msnoise]: msnoise_myproject
    Username [msnoise]: msnoise
    Password: msnoise

MSNoiseDB creates the ``msnoise`` user with password ``msnoise`` by default —
this is intentional for a local, single-user convenience server.

.. warning::

   The credentials are saved in plain text in ``db.ini`` in the project
   folder.  Keep this file private if your server is network-accessible.


Quick-start summary
--------------------

.. code-block:: sh

    # 1. Start MSNoiseDB (from its own folder, every time)
    cd ~/msnoisedb && msnoisedb start

    # 2. Create a database for this project (once only)
    msnoisedb create-db msnoise_myproject

    # 3. Create the project folder and initialise
    mkdir ~/my_msnoise_project && cd ~/my_msnoise_project
    msnoise db init
    # → postgresql | localhost:5099 | msnoise_myproject | msnoise / msnoise

    # 4. Open the configurator
    msnoise admin


Alternative backends
---------------------

**SQLite** (no server, single worker only):

When prompted, choose ``sqlite``.  MSNoise creates a ``.sqlite`` file in the
current folder — no ``msnoisedb create-db`` step needed.

.. warning::

   SQLite does not support concurrent writes.  Running more than one worker
   (``msnoise -t 2 cc compute``) will cause database lock errors.  Use
   PostgreSQL via MSNoiseDB for any real processing.

**Self-managed MySQL / MariaDB or PostgreSQL**:

Create an empty database on your server first, then provide the host, database
name, user and password when prompted.  See :ref:`aboutdbandperformances` for
guidance on when a dedicated server is worth the extra setup.


After initialisation
----------------------

Launch the web configurator to set project parameters before processing:

.. code-block:: sh

    msnoise admin

This opens ``http://localhost:5000`` where you can configure dates, archive
paths, filters, preprocessing parameters, and moving-stack windows.

See :doc:`/workflow_000/001_msnoise_admin` for a tour of the admin interface,
and :ref:`workflow_000` to continue with the next steps.


What the initialiser does
--------------------------

1. Writes ``db.ini`` with the database connection string.
2. Creates all MSNoise tables (``Config``, ``Station``, ``Job``,
   ``WorkflowStep``, ``WorkflowLink``, …) in the database.
3. Populates ``Config`` with default values from the built-in CSV files.
4. Creates the default workflow steps and links (the processing DAG).
5. Creates the default ``DataSource`` record (local SDS archive, id=1).

Running ``msnoise db init`` again in the same folder updates ``db.ini`` but
leaves existing tables and data intact.
"""

import logging
import os
import sys
from getpass import getpass
from sqlalchemy import create_engine, text
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import sessionmaker
from sqlalchemy_utils import database_exists, create_database
import pickle

from .msnoise_table_def import declare_tables


if sys.version_info[0] >= 3:
    raw_input = input


DEFAULT_INPUTS = {
        'sqlite_filename': 'msnoise.sqlite',
        'table_prefix': '',
        'mysql_host': 'localhost',
        'pg_host': 'localhost:5099',   # MSNoiseDB default
        'mysql_db': 'msnoise',
        'mysql_user': 'msnoise',
        }



def create_database_inifile(tech, hostname, database, username, password,
                            prefix=""):
    """Creates the db.ini file based on supplied parameters.

    :type tech: int
    :param tech: The database technology used: 1=sqlite 2=mysql 3=postgresql
    :type hostname: string
    :param hostname: The hostname of the server (if tech>=2) or the name of the
        sqlite file if tech=1)
    :type database: string
    :param database: The database name
    :type username: string
    :param username: The user name
    :type password: string
    :param prefix: The prefix to use for all tables
    :type prefix: string
    :param password: The password of `user`

    :return: None
    """
    f = open(os.path.join(os.getcwd(), 'db.ini'), 'wb')
    pickle.dump([tech, hostname, database, username, password, prefix], f,
                 protocol=2)
    f.close()


def create_indices(session, prefix):
    """
    Create indices in the database, using the provided sqlalchemy session.
    """
    if prefix:
        prefix += '_'
    try:
        session.execute(text("CREATE UNIQUE INDEX job_index ON %sjobs (day, pair, "
                             "jobtype)" % prefix))
        session.commit()
    except Exception:
        logging.info("It looks like the v1.5 'job_index' is already in the DB")
        session.rollback()

    try:
        session.execute(text("CREATE INDEX job_index2 ON %sjobs (jobtype, flag)"
                             % prefix))
        session.commit()
    except Exception:
        logging.info("It looks like the v1.6 'job_index2' is already in the DB")
        session.rollback()

    try:
        session.execute(text("CREATE UNIQUE INDEX da_index ON %sdata_availability ("
                             "path, file, net, sta, comp)" % prefix))
        session.commit()
    except Exception:
        logging.info("It looks like the v1.5 'da_index' is already in the DB")
        session.rollback()


def main(tech=None, hostname=None, username=None, password=None,
         database="msnoise", filename=None, prefix=None):
    """
    Create the db.ini file and create database.

    Interactively ask database type and connection information to the user,
    write them to the db.ini file and create the database tables.
    Information input by the user must also be allowed to be passed as function
    argument for the automatic tests.
    """
    if tech is None:
        database = None
        print("Welcome to MSNoise")
        print()
        print("What database technology do you want to use?")
        while not tech:
            if tech == 0:
                print('Sorry, your choice is invalid. Please enter 1 or 2.')
            print(" [1] sqlite")
            print(" [2] mysql")
            print(" [3] postgresql")
            try:
                tech = int(raw_input('Choice: '))
            except ValueError:
                tech = 0
            else:
                if tech not in (1, 2, 3):
                    tech = 0

        def ask(prompt, default, input_func=raw_input):
            return input_func(prompt.format(default)) or default

        if tech == 1:
            if filename is None:
                filename = ask('Filename: [{}]: ',
                               DEFAULT_INPUTS['sqlite_filename'])
            if prefix is None:
                prefix = ask('Table prefix: []: ', '')
            database = None
            username = None
            password = None
        elif tech == 2:
            if hostname is None:
                hostname = ask('Hostname [{}]: ', DEFAULT_INPUTS['mysql_host'])
            if database is None:
                database = ask('Database name [{}]: ', DEFAULT_INPUTS['mysql_db'])
            if username is None:
                username = ask('Username [{}]: ', DEFAULT_INPUTS['mysql_user'])
            if password is None:
                password = ''
                while not password:
                    password = ask('Password: ', '', getpass)
                    if not password:
                        print('Sorry, you must define a password.')
            if prefix is None:
                prefix = ask('Table prefix [{}]: ',
                             DEFAULT_INPUTS['table_prefix'])
        else:
            if hostname is None:
                hostname = ask('Hostname [{}] (hostname or hostname:port, e.g. localhost:5099 for MSNoiseDB): ',
                               DEFAULT_INPUTS['pg_host'])
            if database is None:
                database = ask('Database name [{}]: ', DEFAULT_INPUTS['mysql_db'])
            if username is None:
                username = ask('Username [{}]: ', DEFAULT_INPUTS['mysql_user'])
            if password is None:
                password = ''
                while not password:
                    password = ask('Password: ', '', getpass)
                    if not password:
                        print('Sorry, you must define a password.')
            if prefix is None:
                prefix = ask('Table prefix [{}]: ',
                             DEFAULT_INPUTS['table_prefix'])
    else:
        tech = int(tech)
        prefix = ""

    if tech == 1:
        if not filename:
            filename = "db.sqlite"
        engine = create_engine('sqlite:///%s' % filename, echo=False)
        database = None
        username = None
        password = None
        hostname = filename
    elif tech == 2:
        engine = create_engine('mysql+pymysql://%s:%s@%s/%s'
                               % (username, password, hostname, database),
                               echo=False)
    elif tech == 3:
        engine = create_engine('postgresql+psycopg2://%s:%s@%s/%s'
                               % (username, password, hostname, database),
                               echo=False)



    if not database_exists(engine.url):
        logging.info("Database does not exist. Creating it right away!")
        create_database(engine.url)

    create_database_inifile(tech, hostname, database, username, password,
                            prefix)

    schema = declare_tables(prefix)
    schema.Base.metadata.create_all(engine)

    Session = sessionmaker(bind=engine)
    session = Session()

    # Add default configuration values to the database
    from .core.config import create_config_set


    try:
        create_config_set(session, "global")
    except IntegrityError:
        session.rollback()
        logging.error("The database seems to already exist and is not empty, "
                      "cannot continue")
        return 1

    # Create the default DataSource (id=1, local SDS archive)
    DataSource = schema.DataSource
    default_source = DataSource(
        name="local",
        uri="",
        data_structure="SDS",
        auth_env="MSNOISE",
    )
    session.add(default_source)
    session.commit()
    logging.info("Created default DataSource 'local' (id=1)")

    create_indices(session, prefix)

    session.close()
    print("Installation Done! - Go to Configuration Step!")
    return 0
