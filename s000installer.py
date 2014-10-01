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

.. code-block:: sh

    python s000installer.py

"""
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from getpass import getpass
from msnoise_table_def import *
from default import default
from database_tools import create_database_inifile

if __name__ == "__main__":
    print "Welcome to MSNoise"
    print
    
    print "What database technology do you want to use?"
    print " [1] sqlite"
    print " [2] mysql"
    tech = int(raw_input('Choice:'))
    
    if tech == 1:
        a = raw_input('Filename: [msnoise.sqlite]: ')
        hostname = a if len(a) != 0 else "msnoise.sqlite"
        engine = create_engine('sqlite:///%s'%hostname, echo=False)
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
        engine = create_engine('mysql://%s:%s@%s/%s'%(username, password, hostname, 
                                                      database), echo=False)
    
    create_database_inifile(tech, hostname, database, username, password)
    
    
    # create tables
    Base.metadata.create_all(engine)
    
    Session = sessionmaker(bind=engine)
    session = Session()
    
    configs = []
    for name in default.keys():
        session.add(Config(name=name,value=default[name][-1]))
    
    # session.add_all(configs)
    session.commit()
    session.close()
    print "Installation Done! - Go to Configuration Step!"