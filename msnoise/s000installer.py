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

.. include:: ../clickhelp/msnoise-install.rst


.. warning:: The credentials will be saved in a flat text file in the current
    directory. It's not very safe, but until now we haven't thought of another
    solution.
"""
import argparse
import sys
from getpass import getpass

from sqlalchemy import create_engine
from sqlalchemy.exc import *
from sqlalchemy.orm import sessionmaker

from .api import create_database_inifile
from .default import default
from .msnoise_table_def import *

if sys.version_info[0] >= 3:
    raw_input = input


def main(tech=None, hostname="localhost", username="msnoise",
         password="msnoise", database="msnoise", filename="msnoise.sqlite"):
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
    
    create_database_inifile(tech, hostname, database, username, password)

    # create tables
    Base.metadata.create_all(engine)
    
    Session = sessionmaker(bind=engine)
    session = Session()
    
    for name in default.keys():
        session.add(Config(name=name, value=default[name][-1]))
    
    try:
        session.commit()
    except IntegrityError:
        print("The database seems to already exist and is not empty, cannot"
              " continue")
        return "Integrity Error - DB already exists"
    session.close()
    msg = "Installation Done! - Go to Configuration Step!"
    return msg
