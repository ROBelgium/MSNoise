MSNoise 1.2.2
=============

Release date: 6 November 2013

Release type: bugfix

Release notes:

BugFix in DatabaseTools
-----------------------

In s06compute_dtt.py, when updated_days_for_dates is called with pair as "%", it doesn't work in database_tools.py.
There is no query result back when calling

.. code-block:: python

    days = session.query(Job).filter(Job.pair == pair).filter(Job.day >= date1).filter(Job.day <= date2).filter(Job.type == type).filter(Job.lastmod >= lastmod).group_by(Job.day).order_by(Job.day).all().

It probably work under windows, but it doesn't work in Linux.

Thanks for reporting, Xiao, it's now corrected.