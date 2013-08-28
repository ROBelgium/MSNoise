# create DB

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

from msnoise_table_def import *
from default import default

engine = create_engine('sqlite:///msnoise.sqlite', echo=True)
# engine = create_engine('mysql://user:ass@server/rc1',echo=True)

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
    