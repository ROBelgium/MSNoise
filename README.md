MSNoise
=======
A Python Package for Monitoring Seismic Velocity Changes using Ambient Seismic Noise.

[![Join the chat at https://gitter.im/ROBelgium/MSNoise](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/ROBelgium/MSNoise?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/ROBelgium/MSNoise.png)](https://travis-ci.org/ROBelgium/MSNoise)
[![Build status](https://ci.appveyor.com/api/projects/status/82k4dw58jhadfung?svg=true)](https://ci.appveyor.com/project/ThomasLecocq/msnoise)
[![Analytics](https://ga-beacon.appspot.com/UA-55331253-1/MSNoise/readme)](https://github.com/ROBelgium/MSNoise)


MSNoise is the first complete software package for computing and monitoring relative velocity variations using ambient seismic noise. 
MSNoise is a fully-integrated solution that automatically scans data archives and determines which jobs need to be done whenever the scheduled task is executed. 

MSNoise is developed by Thomas Lecocq (Royal Observatory of Belgium, ROB) and Corentin Caudron (previously at ROB and EOSingapore, now at University of Cambridge).
The group of active users (providing questions, feedback, snippets of code) is growing and the full list of Contributors is available here: http://msnoise.org/doc/contributors.html. 

History
-------

* 2010: MSNoise is created based on Matlab, c++, csh and fortran codes from Florent Brenguier (IPGP, U. Grenoble-Alpes) and Daniel Clarke (IPGP).
* 2011/12: MSNoise is tested on Undervolc data, and used by Corentin for his PhD thesis.
* 2013: First release of MSNoise for the IAVCEI 2013 in Kagoshima ([Release Notes](http://msnoise.org/doc/releasenotes/msnoise-1.0.html)).
* 2014: Improvements and bugfixes, release 1.2.5. Publication of the SRL article: http://srl.geoscienceworld.org/content/85/3/715.full ([Release Notes](http://msnoise.org/doc/releasenotes/msnoise-1.2.5.html)).
* 2015: MSNoise 1.3: MSNoise is real python package ([Release Notes](http://msnoise.org/doc/releasenotes/msnoise-1.3.html)).
* 2016: MSNoise 1.4: new web-admin, plugin support, instrument response removal and phase weighted stacking ([Release Notes](http://msnoise.org/doc/releasenotes/msnoise-1.4.html)).


Documentation
-------------
The full documentation can be found on: http://www.msnoise.org.


Installation
------------

Please follow the instructions in the documentation: http://msnoise.org/doc/installation.html

Remember, always consider the the current *master* as not stable!


Getting Help
------------
The best way to get help is to subscribe to the Mailing List and ask your question directly there. It is available on 
http://mailman-as.oma.be/mailman/listinfo/msnoise and the archive is either http://mailman-as.oma.be/pipermail/msnoise/ or 
on http://news.gmane.org/gmane.science.geophysics.msnoise for a nice view.


Disclaimer
----------

Although we have cross-checked the whole code, we cannot warranty it is exempt of bugs. The package is provided as-is, we will not be held responsible for any use you make of it, nor for the results and conclusions you may find using MSNoise.


Licence
-------

Although MSNoise is provided freely on the internet and every individual may use it, we do not allow anyone to provide commercial service (read -paid in any way-), like training or teaching without contacting the original authors first. The authors and the contact email address are in the __init__.py file of the package.
