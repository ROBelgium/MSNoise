MSNoise
=======
A Python Package for Monitoring Seismic Velocity Changes using Ambient Seismic Noise.

[![Join the chat at https://gitter.im/ROBelgium/MSNoise](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/ROBelgium/MSNoise?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/ROBelgium/MSNoise.png)](https://travis-ci.org/ROBelgium/MSNoise)
[![Build status](https://ci.appveyor.com/api/projects/status/82k4dw58jhadfung?svg=true)](https://ci.appveyor.com/project/ThomasLecocq/msnoise)
[![codecov](https://codecov.io/gh/ROBelgium/MSNoise/branch/master/graph/badge.svg)](https://codecov.io/gh/ROBelgium/MSNoise)
[![Analytics](https://ga-beacon.appspot.com/UA-55331253-1/MSNoise/readme)](https://github.com/ROBelgium/MSNoise)


MSNoise is the first complete software package for computing and monitoring relative velocity variations using ambient seismic noise. 
MSNoise is a fully-integrated solution that automatically scans data archives and determines which jobs need to be done whenever the scheduled task is executed. 

MSNoise is developed by Thomas Lecocq (Royal Observatory of Belgium, ROB). Corentin Caudron used MSNoise during his PhD at ROB and still continuously provides invaluable debug information.
The group of active users (providing questions, feedback, snippets of code) is growing and the full list of Contributors is available here: http://msnoise.org/doc/contributors.html. 


History
-------

* 2010: MSNoise is based on Matlab, c++, csh and fortran codes developped at ISTerre/Univ. Grenoble and IPGP in the framework of the [ERC Whisper project](https://whisper.obs.ujf-grenoble.fr/).
* 2011/12: MSNoise is tested on Undervolc data, and used by Corentin for his PhD thesis.
* 2013: First release of MSNoise for the IAVCEI 2013 in Kagoshima ([Release Notes](http://msnoise.org/doc/releasenotes/msnoise-1.0.html)).
* 2014: Improvements and bugfixes, release 1.2.5. Publication of the [SRL article](http://srl.geoscienceworld.org/content/85/3/715.full) ([Release Notes](http://msnoise.org/doc/releasenotes/msnoise-1.2.5.html)).
* 2015: MSNoise 1.3: MSNoise is real python package, with a documented API and new plots ([Release Notes](http://msnoise.org/doc/releasenotes/msnoise-1.3.html)).
* 2016: MSNoise 1.4: new web admin interface, plugin support, instrument response removal and phase weighted stacking ([Release Notes](http://msnoise.org/doc/releasenotes/msnoise-1.4.html)).
* 2017: MSNoise 1.5: Autocorrelation / Single Station correlation support, rewritten preprocessing, new_jobs and scan_archive for more performance, better instrument response preloading ([Release Notes](http://msnoise.org/doc/releasenotes/msnoise-1.5.html)).
* 2019: MSNoise 1.6: Optimisation of the workflow (one job type per step), HPC support, faster *compute_cc* step, PSD-whitening, DB optimisations ([Release Notes](http://msnoise.org/doc/releasenotes/msnoise-1.6.html))

Documentation
-------------
The full documentation can be found on: http://www.msnoise.org.


Installation
------------

Please follow the instructions in the documentation: http://msnoise.org/doc/installation.html

Remember, always consider the current GitHub *master* as not stable!


Getting Help
------------
The best way to get help is to subscribe to the Mailing List and ask your question directly there. It is available on 
http://mailman-as.oma.be/mailman/listinfo/msnoise and the archive is http://mailman-as.oma.be/pipermail/msnoise/ or https://www.mail-archive.com/msnoise@mailman-as.oma.be/.

Citing MSNoise
--------------

If you use MSNoise, even a small part of it, for your research and publications, please consider citing it:

**Lecocq, T., C. Caudron, et F. Brenguier (2014)**, MSNoise, a Python Package
for Monitoring Seismic Velocity Changes Using Ambient Seismic Noise,
*Seismological Research Letters*, 85(3), 715‑726, doi:10.1785/0220130073.

Thanks to all [who already did so](http://www.msnoise.org/they-cite-msnoise)! 

Disclaimer
----------

Although we have cross-checked the whole code, we cannot warranty it is exempt of bugs. The package is provided as-is, we will not be held responsible for any use you make of it, nor for the results and conclusions you may find using MSNoise.



Licence
-------

MSNoise is released under EUPL v1.1
