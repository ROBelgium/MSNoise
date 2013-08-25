PYTHON="/home/seisvolc/epd-7.3-2-rh5-x86_64/bin/python"
cd /data/thomas/MSNoise
$PYTHON 01.scan_archive_threaded.py
$PYTHON 02.new_jobs.py
$PYTHON 03.compute_cc.py
$PYTHON 04.stack.py
$PYTHON 05.compute_mwcs.py
$PYTHON 06.compute_dtt.py 
$PYTHON 07.plot_dtt.py
mutt -s "MSnoise Run" email@domain.com -a dtt_allmovstacksNEW_M0.png < MSNoise.log
