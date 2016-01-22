MSNoise 1.3.1
=============

Release date: 26 March 2015

Release type: bugfix

Release notes:

BugFixes
--------

When running the new_jobs procedure, the jobs are "inserted" in the database,
even if they already existed (they should be "updated"). This is because of the
complete rewrite of the code to optimize the operation.

As this optimization is mostly useful upon first run, I've added a parameter
`--init` to the command. If provided, the "massive insert" procedure is used,
if not, then the classic "insert or update if existing" is used.

So, upon first run : ```msnoise new_jobs --init```
And afterwards (in cron, e.g.): ```msnoise new_jobs```

Users who have already run MSNoise 1.3 on their archive need to clean the jobs
table in the database. The buggy jobs are those "CC" jobs which are still marked
"I"n progress after the compute_cc procedure, and with a "lastmod" = "NULL".
They can be easily identified and removed.

a classic SQL command would be:

```DELETE from jobs WHERE lastmod is NULL;```
