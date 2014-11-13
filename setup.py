from setuptools import setup, find_packages

setup(
    name='msnoise',
    version='1.3',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        msnoise=msnoise.scripts.msnoise:run
    ''',
)