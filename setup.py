from setuptools import setup, find_packages
import versioneer

setup( version=versioneer.get_version(),
       cmdclass=versioneer.get_cmdclass(),
    name='msnoise',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'setuptools',
        'numpy>=1.0',
        'scipy',
        'pandas',
        'matplotlib',
        'statsmodels',
        'sqlalchemy',
        'traits',       
        'traitsui',
        'obspy',
        'click',
        'scikits.samplerate',
        'pymysql',
        'flask',
        'flask-admin',
        'bokeh',
        'multiprocessing_logging'
    ],
    entry_points='''
        [console_scripts]
        msnoise=msnoise.scripts.msnoise:run
    ''',
    author = "Thomas Lecocq & MSNoise dev team",
    author_email = "Thomas.Lecocq@seismology.be",
    description = "A Python Package for Monitoring Seismic Velocity Changes using Ambient Seismic Noise",
    license = "EUPL-1.1",
    url = "http://www.msnoise.org",
    keywords="noise monitoring seismic velocity change dvv dtt doublet stretching cross-correlation acoustics seismology"
)