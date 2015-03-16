MSNoise 1.2.4
=============

Release date: 28 April 2014

Release notes:


Thanks to some great early MSNoise adopters who reported problems using the mailing-list, we have identified a few bugs and tricky situations where some steps failed. This is the case for

* scan_archive: there were major issues using Threads, this step uses Process for multiprocessing, which is much safer. The only problem remaining is that there is no more console-logging of the found/identified files. To be corrected soon.

* new_jobs: mostly rewritten, all jobs are properly identified now. There should no more lost jobs !

I've also reflected those changes in the documentation.
