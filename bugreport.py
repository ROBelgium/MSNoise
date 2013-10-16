# -*- coding: utf-8 -*-
#

def ispresent(module):
    try:
        mod = __import__(module)
        if hasattr(mod, '__version__'):
            print "ok",
        try:
            print "[X] %s: %s"%(module,mod.__version__)
        except:
            print "[X] %s: present (no version)"%module
    except:
        print "[ ] %s: not found"%module

import os

if __name__ == "__main__":

    print "************* Computer Report *************"
    
    print 
    print "----------------+PYTHON+-------------------"
    import sys
    print "Python",sys.version, "on",  sys.platform
    print 
    print "---------------+MODULES+-------------------"
    ispresent('numpy')
    ispresent('scipy')
    ispresent('pandas')
    ispresent('matplotlib')
    ispresent('statsmodels')
    ispresent('sqlalchemy')
    ispresent('traitsui')
    ispresent('traits')
    ispresent('enable')
    ispresent('scikits.samplerate')
    ispresent('obspy')
    
    ispresent('setuptools')
    ispresent('jinja2')
    ispresent('sphinx')
    ispresent('reportlab')
    ispresent('configobj')
    ispresent('MySQLdb')
    ispresent('PIL')
    ispresent('py2exe')
    ispresent('pkg_resources')
    ispresent('paramiko')
    ispresent('ctypes')
    ispresent('pyparsing')
    ispresent('distutils')
    ispresent('IPython')
    ispresent('vtk')
    ispresent('wx')
    ispresent('PyQt4')    
    
    print
    print "------------------+SYS+--------------------"
    
    
    for key in os.environ.keys():
        print key
        for value in os.environ[key].split(';'):
            if os.path.isdir(value) or os.path.isfile(value) :
                dir = "[X]"
            else:
                dir = "[ ]"
            print " ", dir, value
    
    
    
    
