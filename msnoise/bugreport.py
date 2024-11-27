# -*- coding: utf-8 -*-
#
import os
import platform
import sys
import importlib

def ispresent(module, how=None):
    try:
        mod = __import__(module)
        ver = importlib.metadata.version(module)
        if ver:
            return "[X] %s: %s"%(module, ver)
        else:
            return "[X] %s: present (no version)"%module
    except:
        return "[ ] %s: not found"% (module)


def main(system=False, modules=False, env=False, all=False, show=True):
    output = []
    output += "\n" + "************* Computer Report *************"
    
    if system or all:
        output += "\n"
        output += "\n" + "----------------+SYSTEM+-------------------"
        output += "\n" + "\n".join(platform.uname())
        output += "\n"

    output += "\n" + "----------------+PYTHON+-------------------"
    output += "\n" + "Python:" + sys.version
    output += "\n"
    output += "\n" + "This script is at " + os.path.abspath(__file__)
    output += "\n"
    if modules or all:
        output += "\n" + "---------------+MODULES+-------------------"
        output += "\n"
        output += "\n" + "Required:"
        output += "\n" + ispresent('setuptools')
        output += "\n" + ispresent('numpy')
        output += "\n" + ispresent('scipy')
        output += "\n" + ispresent('pandas')
        output += "\n" + ispresent('matplotlib')
        output += "\n" + ispresent('sqlalchemy')
        output += "\n" + ispresent('sqlalchemy_utils')
        output += "\n" + ispresent('obspy')
        output += "\n" + ispresent('click')
        output += "\n" + ispresent('pymysql')
        output += "\n" + ispresent('flask')
        output += "\n" + ispresent('flask_admin')
        output += "\n" + ispresent('markdown')
        output += "\n" + ispresent('wtforms')
        output += "\n" + ispresent('folium')
        output += "\n" + ispresent('jinja2')
        output += "\n" + ispresent('tables') + " (pytables)"
        output += "\n" + ispresent('xarray')
        output += "\n" + ispresent('logbook')
        output += "\n" + ispresent('pycwt')

        output += "\n"
        output += "\n" + "Only necessary if you plan to build the doc locally:"
        output += "\n" + ispresent('sphinx')
        output += "\n" + ispresent('sphinx_bootstrap_theme')
        output += "\n" + ispresent('sphinx_rtd_theme')
        output += "\n" + ispresent('sphinx_gallery')
        output += "\n" + ispresent('numpydoc')

        # output += "\n"
        # output += "\n" + "Graphical Backends:"
        # output += "\n" + ispresent('wx')
        # output += "\n" + ispresent('qt')
        # output += "\n" + ispresent('qt4')
        # output += "\n" + ispresent('qt5')
        # output += "\n" + ispresent('pyqt')
        # output += "\n" + ispresent('PyQt4')
        # output += "\n" + ispresent('PyQt5')
        # output += "\n" + ispresent('PySide')
        
        output += "\n"
        output += "\n" + "Not required, just checking:"
        output += "\n" + ispresent('json')
        output += "\n" + ispresent('psutil')
        output += "\n" + ispresent('reportlab')
        output += "\n" + ispresent('configobj')
        output += "\n" + ispresent('pkg_resources')
        output += "\n" + ispresent('paramiko')
        output += "\n" + ispresent('ctypes')
        output += "\n" + ispresent('pyparsing')
        output += "\n" + ispresent('distutils')
        output += "\n" + ispresent('IPython')
        output += "\n" + ispresent('notebook')
        # output += "\n" + ispresent('vtk')
        # output += "\n" + ispresent('enable')
        # output += "\n" + ispresent('traitsui')
        # output += "\n" + ispresent('traits')
        # output += "\n" + ispresent('scikits.samplerate')
        
        output += "\n"
    
    if env or all:    
        output += "\n" + "------------------+ENV+--------------------"

        for key in os.environ.keys():
            output += "\n" + key
            for value in os.environ[key].split(';'):
                if os.path.isdir(value) or os.path.isfile(value) :
                    dir = "[X]"
                else:
                    dir = "[ ]"
                output += "\n" + " ", dir, value
    output = "".join(output)
    if show:
        print(output)
    else:
        return output.split("\n")
    
if __name__ == "__main__":
    main()
