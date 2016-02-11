# -*- coding: utf-8 -*-
#

import os
import platform
import sys

def ispresent(module, how=None):
    try:
        mod = __import__(module)
        if hasattr(mod, '__version__'):
            return "[X] %s: %s"%(module,mod.__version__)
        else:
            return "[X] %s: present (no version)"%module
    except:
        return "[ ] %s: not found (install via %s)"% (module, how)


def main(system=False, modules=False, env=False, all=False, show=True):
    output = []
    output += "\n" + "************* Computer Report *************"
    
    if system or all:
        output += "\n"
        output += "\n" + "----------------+SYSTEM+-------------------"
        output += "\n" + "\n".join(platform.uname())
        if platform.system() == "Linux":
            output += "\n" + " - ".join(platform.linux_distribution())
        output += "\n"

    output += "\n" + "----------------+PYTHON+-------------------"
    output += "\n" + "Python:" + sys.version
    output += "\n"
    if modules or all:
        output += "\n" + "---------------+MODULES+-------------------"
        output += "\n"
        output += "\n" + "Required:"
        output += "\n" + ispresent('setuptools')
        output += "\n" + ispresent('click', 'easy_install click')
        output += "\n" + ispresent('numpy')
        output += "\n" + ispresent('scipy')
        output += "\n" + ispresent('pandas')
        output += "\n" + ispresent('matplotlib')
        output += "\n" + ispresent('statsmodels')
        output += "\n" + ispresent('sqlalchemy')
        output += "\n" + ispresent('traitsui')
        output += "\n" + ispresent('traits')
        output += "\n" + ispresent('scikits.samplerate')
        output += "\n" + ispresent('obspy')
        output += "\n" + ispresent('flask')
        output += "\n" + ispresent('flask.ext.admin', 'easy_install flask-admin')
        output += "\n" + ispresent('flask_admin')
        output += "\n" + ispresent('wtforms')
        output += "\n" + ispresent('bokeh')

        output += "\n"
        output += "\n" + "Only necessary if you plan to build the doc locally:"
        output += "\n" + ispresent('sphinx')
        output += "\n" + ispresent('jinja2')

        output += "\n"
        output += "\n" + "Graphical Backends: (at least one is required)"
        output += "\n" + ispresent('wx')
        output += "\n" + ispresent('PyQt4')
        output += "\n" + ispresent('PySide')
        
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
        output += "\n" + ispresent('vtk')
        output += "\n" + ispresent('enable')
        
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
