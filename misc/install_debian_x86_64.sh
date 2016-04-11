lsb_release -a # get info on the operating system
sudo apt-get update
sudo apt-get install libsamplerate0 libsamplerate0-dev
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
source .bashrc
# Update conda itself
conda update --yes conda
# Installed required packages
conda install --yes pip numpy scipy matplotlib statsmodels pandas
conda install --yes sqlalchemy sphinx jinja2 pymysql click
conda install --yes flask lxml
# Update remaining packages not available via conda
pip install flask-admin folium multiprocessing_logging markdown
pip install obspy

# Special install for scikits.samplerate
cd ..
sudo dpkg -L libsamplerate0
sudo dpkg -L libsamplerate0-dev
wget https://pypi.python.org/packages/source/s/scikits.samplerate/scikits.samplerate-0.3.3.tar.gz#md5=96c8d8ba3aa95a9db15994f78792efb4
tar -xvf scikits.samplerate-0.3.3.tar.gz
cd scikits.samplerate-0.3.3
echo "[samplerate]" >> site.cfg
echo "library_dirs=/usr/lib/x86_64-linux-gnu" >> site.cfg
echo "include_dirs=/usr/include" >> site.cfg
python setup.py build
python setup.py install
# Finished procedure, going to root folder
cd ..
ls -la