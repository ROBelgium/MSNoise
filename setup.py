from setuptools import setup, find_packages

setup(version="1.5a",
      name='msnoise',
      packages=find_packages(),
      include_package_data=True,
      install_requires=[
          'setuptools',
          'numpy>=1.0',
          'scipy',
          'pandas',
          'matplotlib',
          'sqlalchemy',
          'obspy',
          'click',
          'pymysql',
          'flask',
          'flask-admin',
          'multiprocessing_logging',
          'markdown',
          'folium',
          'wtforms',
          'jinja2'
      ],
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
)