import cPickle
from database_tools import *
from getpass import getpass

print "Welcome to MsNoise"
print

try:
    f = open("db.ini","r")
    hostname,database,username,password = cPickle.load(f)
    a = raw_input("Parameters for %s@%s defined, use them? [Y]/N"%(username,hostname))
    use = True if a in ['y','Y',''] else False 
except:
    use = False
    
if not use:
    a = raw_input('mysql hostname [localhost]:')
    hostname = a if len(a) != 0 else "localhost"

    a = raw_input('mysql database [msnoise-test]:')
    database = a if len(a) != 0 else "msnoise-test"

    a = raw_input('mysql username [msnoise]:')
    username = a if len(a) != 0 else "msnoise"

    a = getpass('mysql password [msnoise]:')
    password = a if len(a) != 0 else "msnoise"
    
    create_database_inifile(hostname,database,username,password)

try:
    db = connect()
except:
    a = raw_input("This database doesn't exist, do you want to create it? [Y]/N")
    createdb = True if a in ['y','Y','']  else False
    if createdb:
        a = raw_input("Who has CREATE DATABASE privilege to this host? [root]:")
        root = a if len(a) != 0 else "root"
        
        a = getpass("His password?:")
        passwd = a if len(a) != 0 else ""
        
        create_database(root, passwd)
        db = connect()
        create_stations_table(db)
        create_filters_table(db)
        create_config_table(db)
        create_data_availability_table(db)
        create_jobs_table(db)
        create_dtt_table(db)

## THIS IS BRAND NEW, HELL YEAH
if get_config(db,'data_folder') == "":
    print "*"*70
    print "Now edit the configuration (data_structure and data_folder) using either \
configurator.py or phpMyAdmin and \
launch the intaller_populatestationtable.py script for populating the station table"
    print "*"*70