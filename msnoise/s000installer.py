"""
This console script is responsible asking questions about the database
connection, to create the db.ini file in order to store the answers and
to create the tables in the database.

Questions are:

* What database technology do you want to use?

  * sqlite: this will create a file in the current folder and use it as DB

  * mysql: this will connect to a local or remote mysql server, additional
    information is then required:
 
    - hostname: of the mysql server, defaults to 127.0.0.1  
    - database: must already exist on `hostname`
    - username: as registered in the privileged users of the mysql server
    - password: his password

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
import sys
import logging
from getpass import getpass

from sqlalchemy import create_engine
from sqlalchemy.exc import *
from sqlalchemy.orm import sessionmaker


from .default import default


if sys.version_info[0] >= 3:
    raw_input = input

import os
try:
    import cPickle
except:
    import pickle as cPickle

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


def main(tech=None, hostname="localhost", username="msnoise",
         password="msnoise", database="msnoise", filename="msnoise.sqlite",
         prefix=""):
    global install_mode
    install_mode = True
    if tech is None:
        print("Welcome to MSNoise")
        print()
        print("What database technology do you want to use?")
        print(" [1] sqlite")
        print(" [2] mysql")
        tech = int(raw_input('Choice:'))

        if tech == 1:
            a = raw_input('Filename: [msnoise.sqlite]: ')
            hostname = a if len(a) != 0 else "msnoise.sqlite"

            a = raw_input('Table prefix: []: ')
            prefix = a if len(a) != 0 else ""

            database = None
            username = None
            password = None
        else:
            a = raw_input('Server: [127.0.0.1]: ')
            hostname = a if len(a) != 0 else "127.0.0.1"
            a = raw_input('Database: [msnoise]: ')
            database = a if len(a) != 0 else "msnoise"

            a = raw_input('Username: [msnoise]: ')
            username = a if len(a) != 0 else "msnoise"
            a = getpass('Password: [msnoise]: ')
            password = a if len(a) != 0 else "msnoise"

            a = raw_input('Table prefix: []: ')
            prefix = a if len(a) != 0 else ""
    else:
        tech = int(tech)

    if tech == 1:
        engine = create_engine('sqlite:///%s' % filename, echo=False)
        database = None
        username = None
        password = None
        hostname = filename
    else:
        engine = create_engine('mysql+pymysql://%s:%s@%s/%s' % (username,
                                                                password,
                                                                hostname,
                                                                database),
                               echo=False)

    # from .api import create_database_inifile
    create_database_inifile(tech, hostname, database, username, password,
                            prefix)

    from .msnoise_table_def import PrefixerBase, Config
    PrefixerBase._the_prefix = prefix
    # create tables
    PrefixerBase.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()
    
    for name in default.keys():
        session.add(Config(name=name, value=default[name][1]))
    
    try:
        session.commit()
    except IntegrityError:
        print("The database seems to already exist and is not empty, cannot"
              " continue")
        return "Integrity Error - DB already exists"
    if prefix != "":
        prefix = prefix + "_"
    # TODO move those calls to a def and call it from install / msnoise.scripts
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

    session.close()
    msg = "Installation Done! - Go to Configuration Step!"
    return msg
