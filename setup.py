from setuptools import setup, find_packages
import os
import sys
import inspect


SETUP_DIRECTORY = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

UTIL_PATH = os.path.join(SETUP_DIRECTORY, "msnoise")
sys.path.insert(0, UTIL_PATH)

if UTIL_PATH:   # To avoid PEP8 E402
    from _version import get_git_version  # @UnresolvedImport

sys.path.pop(0)


setup(version=get_git_version(),
      name='msnoise',
      packages=find_packages(),
      include_package_data=True,
      install_requires=[
          'setuptools',
          'numpy>=1.0',
          'scipy',
          'pandas',
          'matplotlib',
          'sqlalchemy<2', # TEMP needed as it breaks flask-admin
          'obspy',
          'click',
          'pymysql',
          'flask',
          'flask-admin',
          'flask-wtf',
          'markdown',
          'folium',
          'wtforms',
          'jinja2',
          'scandir',  # useful for python < 3.5
          'logbook',
          'xarray'  # new in 2.0
      ],
      extras_require={
          'doc': [
              'sphinx',
              'sphinx_bootstrap_theme',
              'numpydoc',
              'sphinx_gallery'
              ],
          },
      entry_points='''
          [console_scripts]
          msnoise=msnoise.scripts.msnoise:run
      ''',
      author="Thomas Lecocq & MSNoise dev team",
      author_email="Thomas.Lecocq@seismology.be",
      description="A Python Package for Monitoring Seismic Velocity Changes"
                  " using Ambient Seismic Noise",
      license="EUPL-1.1",
      url="http://www.msnoise.org",
      keywords="noise monitoring seismic velocity change dvv dtt doublet"
               " stretching cross-correlation acoustics seismology",
      zip_safe=False,
      platforms='OS Independent',
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics'],
      )
