# -*- coding: utf-8 -*-
#

import os
import platform
import sys
import argparse

def ispresent(module, how=None):
    try:
        mod = __import__(module)
        if hasattr(mod, '__version__'):
            print "[X] %s: %s"%(module,mod.__version__)
        else:
            print "[X] %s: present (no version)"%module
    except:
        print "[ ] %s: not found (install via %s)"% (module, how)


def main(system=False, modules=False, env=False, all=False):
    #~ parser = argparse.ArgumentParser(description='Helps determining what didn\'t work')
    #~ parser.add_argument('-s', '--sys', action="store_true",
                        #~ help='Outputs System info',
                        #~ default=False)
    #~ parser.add_argument('-m', '--modules', action="store_true",
                        #~ help='Outputs Python Modules Presence/Version',
                        #~ default=True)
    #~ parser.add_argument('-e', '--env', action="store_true",
                        #~ help='Outputs System Environment Variables',
                        #~ default=False)
    #~ parser.add_argument('-a', '--all', action="store_true",
                        #~ help='Outputs all of the above',
                        #~ default=False)
                                                
    #~ args = parser.parse_args()
    #~ if args.modules:
        #~ modules=True

    print "************* Computer Report *************"
    
    if system or all:
        print 
        print "----------------+SYSTEM+-------------------"
        print "\n".join(platform.uname())
        if platform.system() == "Linux":
            print " - ".join(platform.linux_distribution())
        print 

    print "----------------+PYTHON+-------------------"
    print "Python:",sys.version
    print
    if modules or all:
        print "---------------+MODULES+-------------------"
        print
        print "Required:"
        ispresent('setuptools')
        ispresent('click', 'easy_install click')
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
        
        ispresent('sphinx')
        ispresent('jinja2')
        
        ispresent('flask')
        ispresent('flask.ext.admin', 'easy_install flask-admin')
        ispresent('wtforms')
        ispresent('json')
        ispresent('psutil')
        
        print
        print "Backends: (at least one is required)"
        ispresent('wx')
        ispresent('PyQt4')
        ispresent('PySide')
        
        print
        print "Not required, just checking:"
        ispresent('reportlab')
        ispresent('configobj')
        ispresent('pkg_resources')
        ispresent('paramiko')
        ispresent('ctypes')
        ispresent('pyparsing')
        ispresent('distutils')
        ispresent('IPython')
        ispresent('vtk')
        
        print
    
    if env or all:    
        print "------------------+ENV+--------------------"
        
        
        for key in os.environ.keys():
            print key
            for value in os.environ[key].split(';'):
                if os.path.isdir(value) or os.path.isfile(value) :
                    dir = "[X]"
                else:
                    dir = "[ ]"
                print " ", dir, value

    
if __name__ == "__main__":
    main()
