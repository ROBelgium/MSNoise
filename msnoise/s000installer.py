"""
This console scripts is responsible asking questions about the database
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

To run this script:

.. include:: clickhelp/msnoise-install.rst


.. warning:: The credentials will be saved in a flat text file in the current
    directory. It's not very safe, but until now we haven't thought of another
    solution.
"""
import argparse
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import *
from getpass import getpass
from .msnoise_table_def import *
from .default import default
from .api import create_database_inifile



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
        engine = create_engine('sqlite:///%s'%filename, echo=False)
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
    
    configs = []
    for name in default.keys():
        session.add(Config(name=name,value=default[name][-1]))
    
    # session.add_all(configs)
    try:
        session.commit()
    except IntegrityError:
        print("The database seems to already exist and is not empty, cannot continue")
        return("Integrity Error - DB already exist")
    session.close()
    msg = "Installation Done! - Go to Configuration Step!"
    return msg
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Creates the database connection file (db.ini)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tech', type=int, default=None,
                        help='DB technology to use [1] sqlite or [2] mysql.')
    parser.add_argument('-f', '--filename', type=str, default='msnoise.sqlite',
                        help='File name of the SQLite DB')
    parser.add_argument('-d', '--hostname', type=str, default='localhost',
                        help='Hostname for the MySQL DB')
    parser.add_argument('-u', '--username', type=str, default='msnoise',
                        help='Username for the MySQL DB')
    parser.add_argument('-p', '--password', type=str, default='msnoise',
                        help='Password for the MySQL DB')
    parser.add_argument('-b', '--database', type=str, default='msnoise',
                        help='Database name for the MySQL DB')
    args = parser.parse_args()
    main(args.tech, args.hostname, args.username, args.password, args.database, args.filename)
    
