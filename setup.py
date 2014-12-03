from setuptools import setup, find_packages

setup(
    name='msnoise',
    version='1.3',
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
        'enable',
        'obspy',
        'Click',
    ],
    entry_points='''
        [console_scripts]
        msnoise=msnoise.scripts.msnoise:run
    ''',
)