"""
This console script is responsible for asking questions about the database
connection, to create the db.ini file and to create the tables in the database.

Questions are:

* What database technology do you want to use?

  * sqlite: this will create a file in the current folder and use it as DB

  * mysql: this will connect to a local or remote mysql server, additional
    information is then required:

    - hostname: of the mysql server, defaults to 127.0.0.1
    - database: must already exist on `hostname`
    - username: as registered in the privileged users of the mysql server
    - password: his password
    - prefix: useful when users have only access to a single database. Similar
      to the way wordpress handles prefixes. The tables will be named
      ``%prefix%_config`` (etc) instead of ``config``, for example.

The SQLite choice will create a xxx.sqlite file in the current (project) folder,
while, for MySQL, one has to create an empty database first on the mysql server,
see :ref:`how to do this <emptydb>` .

To run this script:

.. include:: ../clickhelp/msnoise-db-init.rst


.. warning:: The credentials will be saved in a flat text file in the current
    directory. It's not very safe, but until now we haven't thought of another
    solution.
"""

import argparse
import logging
import os
import sys
from getpass import getpass
from sqlalchemy import create_engine
from sqlalchemy.exc import *
from sqlalchemy.orm import sessionmaker
try:
    import cPickle
except ImportError:
    import pickle as cPickle

from .default import default
from .msnoise_table_def import declare_tables


if sys.version_info[0] >= 3:
    raw_input = input


DEFAULT_INPUTS = {
        'sqlite_filename': 'msnoise.sqlite',
        'table_prefix': '',
        'mysql_host': 'localhost',
        'mysql_db': 'msnoise',
        'mysql_user': 'msnoise',
        }



def create_database_inifile(tech, hostname, database, username, password,
                            prefix=""):
    """Creates the db.ini file based on supplied parameters.

    :type tech: int
    :param tech: The database technology used: 1=sqlite 2=mysql
    :type hostname: string
    :param hostname: The hostname of the server (if tech=2) or the name of the
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
    cPickle.dump([tech, hostname, database, username, password, prefix], f,
                 protocol=2)
    f.close()


def create_indices(session, prefix):
    """
    Create indices in the database, using the provided sqlalchemy session.
    """
    # TODO move those calls to a def and call it from install / msnoise.scripts
    if prefix:
        prefix += '_'
    try:
        session.execute("CREATE UNIQUE INDEX job_index ON %sjobs (day, pair, "
                        "jobtype)" % prefix)
        session.commit()
    except:
        logging.info("It looks like the v1.5 'job_index' is already in the DB")
        session.rollback()

    try:
        session.execute("CREATE INDEX job_index2 ON %sjobs (jobtype, flag)"
                        % prefix)
        session.commit()
    except:
        logging.info("It looks like the v1.6 'job_index2' is already in the DB")
        session.rollback()

    try:
        session.execute("CREATE UNIQUE INDEX da_index ON %sdata_availability ("
                        "path, file, net, sta, comp)" % prefix)
        session.commit()
    except:
        logging.info("It looks like the v1.5 'da_index' is already in the DB")
        session.rollback()


def main(tech=None, hostname=None, username=None, password=None,
         database="msnoise", filename="msnoise.sqlite", prefix=None):
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
            try:
                tech = int(raw_input('Choice: '))
            except ValueError:
                tech = 0
            else:
                if tech not in (1, 2):
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
        else:
            if hostname is None:
                hostname = ask('Server: [{}]: ', DEFAULT_INPUTS['mysql_host'])
            if database is None:
                database = ask('Database: [{}]: ', DEFAULT_INPUTS['mysql_db'])
            if username is None:
                username = ask('Username: [{}]: ', DEFAULT_INPUTS['mysql_user'])
            if password is None:
                password = ''
                while not password:
                    password = ask('Password (not shown as you type): ',
                                   '', getpass)
                    if not password:
                        print('Sorry, you must define a password.')
            if prefix is None:
                prefix = ask('Table prefix: [{}]: ',
                             DEFAULT_INPUTS['table_prefix'])
    else:
        tech = int(tech)
        prefix = ""

    if tech == 1:
        engine = create_engine('sqlite:///%s' % filename, echo=False)
        database = None
        username = None
        password = None
        hostname = filename
    else:
        engine = create_engine('mysql+pymysql://%s:%s@%s/%s'
                               % (username, password, hostname, database),
                               echo=False)

    create_database_inifile(tech, hostname, database, username, password,
                            prefix)

    schema = declare_tables(prefix)
    schema.Base.metadata.create_all(engine)

    Session = sessionmaker(bind=engine)
    session = Session()

    # Add default configuration values to the database
    for name in default.keys():
        session.add(schema.Config(name=name, value=default[name][1]))

    try:
        session.commit()
    except IntegrityError:
        logging.error("The database seems to already exist and is not empty, "
                      "cannot continue")
        return 1

    create_indices(session, prefix)

    session.close()
    print("Installation Done! - Go to Configuration Step!")
    return 0
