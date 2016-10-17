.. include:: ../configs.hrst

MSNoise 1.5
===========

Release date: XX XXXXX XXXX


Release type: major


Release notes:

.. contents::
    :local:

Introduction
------------
Just over XXX after the last major release (:doc:`msnoise-1.4`) we are proud
to announce the new :doc:`msnoise-1.5`. It is a **major** release, with a
massive amount of work since the last one: in `GitHub numbers
<https://github.com/ROBelgium/MSNoise/graphs/contributors?from=2016-04-20&to=2017-04-20&type=c>`_
, it's over XXX commits and about XXXX new lines of code and documentation
added!

MSNoise 1.5 introduces **XXXX major new features** :

This version has benefited from outputs/ideas/pull requests/questions from
several users/friends:

* X
* X
* all others (don't be mad :-) )


Thanks to all for using MSNoise, and please, let us know why/how you use it
(and please cite it!)!

To date, we found/are aware of XXX publications using MSNoise ! That's the best
validation of our project ever ! See the full list on the
`MSNoise website <http://www.msnoise.org/they-cite-msnoise/>`_.


*Thomas & Corentin*


~~~~

PS: if you use MSNoise for your research and prepare publications, please
consider citing it:

**Lecocq, T., C. Caudron, et F. Brenguier (2014)**, MSNoise, a Python Package
for Monitoring Seismic Velocity Changes Using Ambient Seismic Noise,
*Seismological Research Letters*, 85(3), 715‑726, doi:10.1785/0220130073.


Web-based Admin Interface Changes
---------------------------------


See :ref:`webadmin` for more details !


Plugin support
--------------




Command Line changes
--------------------


* msnoise config --sync

All commands are now documented: :doc:`../clickhelp/msnoise`.




New plots
---------

* X


Performance improvements
------------------------

Improvements in terms of performances have also been done for MSNoise 1.5:

* Added fftpack optimized nfft



Upgrading an existing project to MSNoise 1.5
--------------------------------------------

Some users will want to keep their current project without recomputing
everything. This requires adding a few configuration parameters to the database

Running the following command will take care of the upgrade from 1.4 to 1.5:

.. code-block:: sh

    msnoise upgrade_db





A final note about development pace and choices
-----------------------------------------------

* MSNoise team is

  * **1 developper** (Thomas)
  * 1 dedicated debugger (Corentin)
  * less than 5 really *active* users, providing feedback and/or lines of codes
    (Esteban, Raphaël, Aurélien, Carmelo, ...)

* All software engineering ideas are coming from too infrequent beerstormings
  between Thomas & others
* The web-interface and the plugin support were developed during Thomas'
  holidays

If you need help, please ask your question on the mailing list. Don't be afraid
to ask. If you have ideas, please share them. If you develop codes to supplement
MSNoise, please share them, even if very small, even if you don't master gitHub.
If you have complaints, post them too, but remember that the package you are
using has been coded by 1 person, and that it's not his full time job. So
MSNoise is provided "as-is", carefully written and tested, but there will be
bugs, issues, incompatibility with certain python installations, OS or module
versions. If you **want or need** developments made, contact Thomas via email
directly. If these developments are within the focus of the developers'
research, then a collaboration, i.e. resulting in a co-authored peer reviewed
publication, can be an option. Otherwise, you can contract us for
paid-developments.
