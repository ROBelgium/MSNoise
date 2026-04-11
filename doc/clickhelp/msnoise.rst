Help on the msnoise commands
============================

This page shows all the commands accessible from the command line.

Commands with an _old suffix are only visible in the console if a .old file is present in the folder. 

msnoise admin
-------------
.. code-block:: sh

    msnoise admin --help

    Usage:  [OPTIONS]
    
      Starts the Web Admin on http://localhost:5000 by default
    
    Options:
      -p, --port INTEGER  Port to open
      --help              Show this message and exit.


msnoise db init
---------------
.. code-block:: sh

    msnoise db init --help

    Usage:  [OPTIONS]
    
      This command initializes the current folder to be a MSNoise Project by
      creating a database and a db.ini file.
    
    Options:
      --tech TEXT      Database technology: 1=SQLite 2=MySQL/MariaDB 3=PostgreSQL
      --auto-workflow  Automatically create all default config sets, workflow
                       steps and links without prompting.
      --help           Show this message and exit.


msnoise db update_loc_chan
--------------------------
.. code-block:: sh

    msnoise db update_loc_chan --help

    Usage:  [OPTIONS]
    
      Populates the Location & Channel from the Data Availability table. Warning:
      rewrites automatically, no confirmation.
    
    Options:
      --help  Show this message and exit.


msnoise db execute
------------------
.. code-block:: sh

    msnoise db execute --help

    Usage:  [OPTIONS] SQL_COMMAND
    
      EXPERT MODE: Executes 'sql_command' on the database. Use this command at
      your own risk!!
    
    Options:
      -o, --outfile TEXT  Output filename (?="request.csv")
      -s, --show BOOLEAN  Show output (in case of SELECT statement)?
      --help              Show this message and exit.


msnoise db upgrade
------------------
.. code-block:: sh

    msnoise db upgrade --help

    Usage:  [OPTIONS]
    
      Upgrade the database from a previous version.
    
      Ensures every parameter defined in the config CSV files is present in the
      database with its default value.  Covers all categories (global, cc, mwcs,
      psd, …) for every config set already in the DB — not just global params.
    
      Safe to run on an already up-to-date project: existing values are never
      overwritten.
    
    Options:
      --help  Show this message and exit.


msnoise db clean_duplicates
---------------------------
.. code-block:: sh

    msnoise db clean_duplicates --help

    Usage:  [OPTIONS]
    
      Checks the Jobs table and deletes duplicate entries
    
    Options:
      --help  Show this message and exit.


msnoise db dump
---------------
.. code-block:: sh

    msnoise db dump --help

    Usage:  [OPTIONS]
    
      Dumps the complete database in formatted files, defaults to CSV.
    
    Options:
      --format TEXT
      --help         Show this message and exit.


msnoise db import
-----------------
.. code-block:: sh

    msnoise db import --help

    Usage:  [OPTIONS] TABLE
    
      Imports msnoise tables from formatted files (CSV).
    
    Options:
      --format TEXT
      --force
      --help         Show this message and exit.


msnoise info
------------
.. code-block:: sh

    msnoise info --help

    Usage:  [OPTIONS]
    
      Outputs general information about the current install and config, plus
      information about jobs and their status.
    
    Options:
      -j, --jobs  Jobs Info only
      --help      Show this message and exit.


msnoise config sync
-------------------
.. code-block:: sh

    msnoise config sync --help

    Usage:  [OPTIONS]
    
      Synchronise station metadata from inventory/dataless.
    
    Options:
      --help  Show this message and exit.


msnoise config set
------------------
.. code-block:: sh

    msnoise config set --help

    Usage:  [OPTIONS] KEY VALUE
    
      Set a configuration parameter using dot notation.
    
      KEY format:
        name                   →  global parameter (e.g. output_folder)
        category.name          →  category set 1  (e.g. cc.cc_sampling_rate)
        category.N.name        →  explicit set N  (e.g. mwcs.2.mwcs_wlen)
    
      Examples:
        msnoise config set output_folder /data/output
        msnoise config set cc.cc_sampling_rate 25
        msnoise config set cc.2.cc_sampling_rate 25
        msnoise config set mwcs.2.mwcs_wlen 10
        msnoise config set global.hpc Y
    
    Options:
      --help  Show this message and exit.


msnoise config get
------------------
.. code-block:: sh

    msnoise config get --help

    Usage:  [OPTIONS] KEY
    
      Get a configuration parameter value using dot notation.
    
      KEY format:
        name                   →  global parameter
        category.name          →  category set 1
        category.N.name        →  explicit set N
    
      Examples:
        msnoise config get output_folder
        msnoise config get cc.cc_sampling_rate
        msnoise config get mwcs.2.mwcs_wlen
    
    Options:
      --help  Show this message and exit.


msnoise config list
-------------------
.. code-block:: sh

    msnoise config list --help

    Usage:  [OPTIONS] [CATEGORY[.N]]
    
      List configuration parameters, optionally filtered.
    
      Examples:
        msnoise config list              # all categories, all sets
        msnoise config list cc           # all cc sets
        msnoise config list mwcs.2       # mwcs set 2 only
    
      Non-default values are marked with *.
    
    Options:
      --help  Show this message and exit.


msnoise config reset
--------------------
.. code-block:: sh

    msnoise config reset --help

    Usage:  [OPTIONS] KEY
    
      Reset a configuration parameter to its default value.
    
      KEY format: same dot notation as 'config set'.
    
      Examples:
        msnoise config reset cc.cc_sampling_rate
        msnoise config reset mwcs.2.mwcs_wlen
    
    Options:
      --help  Show this message and exit.


msnoise config create-set
-------------------------
.. code-block:: sh

    msnoise config create-set --help

    Usage:  [OPTIONS] SET_NAME
    
      Create a configuration set for a workflow step.
    
      SET_NAME: Name of the workflow step (e.g., mwcs, mwcs_dtt, etc.)
    
    Options:
      --help  Show this message and exit.


msnoise config delete-set
-------------------------
.. code-block:: sh

    msnoise config delete-set --help

    Usage:  [OPTIONS] SET_NAME SET_NUMBER
    
      Delete a configuration set.
    
      SET_NAME: The category name (e.g., mwcs, mwcs_dtt) SET_NUMBER: The set
      number to delete
    
    Options:
      --confirm  Skip confirmation prompt
      --help     Show this message and exit.


msnoise config list-sets
------------------------
.. code-block:: sh

    msnoise config list-sets --help

    Usage:  [OPTIONS]
    
      List all configuration sets.
    
    Options:
      -c, --category TEXT  Filter by category name
      --help               Show this message and exit.


msnoise config show-set
-----------------------
.. code-block:: sh

    msnoise config show-set --help

    Usage:  [OPTIONS] SET_NAME SET_NUMBER
    
      Show details of a configuration set.
    
      SET_NAME: The category name (e.g., mwcs, mwcs_dtt) SET_NUMBER: The set
      number to show
    
    Options:
      --help  Show this message and exit.


msnoise config copy-set
-----------------------
.. code-block:: sh

    msnoise config copy-set --help

    Usage:  [OPTIONS] OLD_SET_NAME OLD_SET_NUMBER NEW_SET_NAME NEW_SET_NUMBER
    
      Copy a configuration set to a new set.
    
      Useful for creating variations of existing configurations.
    
    Options:
      --help  Show this message and exit.


msnoise config create-all-sets
------------------------------
.. code-block:: sh

    msnoise config create-all-sets --help

    Usage:  [OPTIONS]
    
      Create one configuration set for each workflow category.
    
      This command creates a complete set of workflow configurations for each
      category.
    
    Options:
      --force    Force creation even if config sets already exist
      --dry-run  Show what would be created without actually creating it
      --help     Show this message and exit.


msnoise config create-workflow-step
-----------------------------------
.. code-block:: sh

    msnoise config create-workflow-step --help

    Usage:  [OPTIONS]
    
      Create a new workflow step interactively
    
    Options:
      --help  Show this message and exit.


msnoise config create-workflow-steps-from-configs
-------------------------------------------------
.. code-block:: sh

    msnoise config create-workflow-steps-from-configs --help

    Usage:  [OPTIONS]
    
      Create workflow steps automatically from all existing config sets.
    
      This command scans all configuration sets in the database and creates
      corresponding workflow steps, sorted by natural workflow order.
    
    Options:
      -v, --verbose  Show detailed output
      --help         Show this message and exit.


msnoise config list-workflow-steps
----------------------------------
.. code-block:: sh

    msnoise config list-workflow-steps --help

    Usage:  [OPTIONS]
    
      List all workflow steps
    
    Options:
      --help  Show this message and exit.


msnoise config show-workflow-graph
----------------------------------
.. code-block:: sh

    msnoise config show-workflow-graph --help

    Usage:  [OPTIONS]
    
      Show workflow graph
    
    Options:
      --help  Show this message and exit.


msnoise config create-workflow-links
------------------------------------
.. code-block:: sh

    msnoise config create-workflow-links --help

    Usage:  [OPTIONS]
    
      Create workflow links automatically between existing workflow steps.
    
      This command creates links between workflow steps following the natural
      workflow progression: preprocess -> cc -> filter -> stack ->
      mwcs/stretching/wavelet -> mwcs_dtt/wavelet_dtt.
    
      Links are created based on matching set numbers and workflow logic.
    
    Options:
      -v, --verbose  Show detailed output
      --help         Show this message and exit.


msnoise reset
-------------
.. code-block:: sh

    msnoise reset --help

    Usage:  [OPTIONS] JOBTYPE
    
      Resets the jobs to "T"odo. JOBTYPE is the step name (e.g. cc_1, stack_1). By
      default only resets jobs "I"n progress. --all resets all jobs, whatever the
      flag value.
    
    Options:
      -a, --all        Reset all jobs
      -r, --rule TEXT  Reset job that match this SQL rule
      --help           Show this message and exit.


msnoise populate
----------------
.. code-block:: sh

    msnoise populate --help

    Usage:  [OPTIONS]
    
      Rapidly scan the archive filenames and find Network/Stations, only works
      with known archive structures, or with a custom code provided by the user.
    
    Options:
      --fromDA  Populates the station table using network and station codes found
                in the data_availability table, overrides the default workflow
                step.
      --help    Show this message and exit.


msnoise scan_archive
--------------------
.. code-block:: sh

    msnoise scan_archive --help

    Usage:  [OPTIONS]
    
      Scan the archive and insert into the Data Availability table.
    
    Options:
      -i, --init         First run ?
      --path TEXT        Scan all files in specific folder, overrides the default
                         workflow step.
      -r, --recursively  When scanning a path, walk subfolders automatically ?
      --crondays TEXT    Number of past days to monitor, typically used in cron
                         jobs (overrides the 'crondays' configuration value). Must
                         be a float representing a number of days, or designate
                         weeks, days, and/or hours using the format 'Xw Xd Xh'.
      --help             Show this message and exit.


msnoise plot data_availability
------------------------------
.. code-block:: sh

    msnoise plot data_availability --help

    Usage:  [OPTIONS]
    
      Plots the Data Availability vs time
    
    Options:
      -c, --chan TEXT     Channel, you can use the ? wildcard, e.g. '?HZ'
                          (default) or 'HH?', etc.
      -s, --show BOOLEAN  Show figure interactively?
      -o, --outfile TEXT  Output filename (?=auto). Supports any matplotlib
                          format, e.g. ?.pdf for PDF with automatic naming.
      --help              Show this message and exit.


msnoise plot station_map
------------------------
.. code-block:: sh

    msnoise plot station_map --help

    Usage:  [OPTIONS]
    
      Plots the station map (very very basic)
    
    Options:
      -s, --show BOOLEAN  Show figure interactively?
      -o, --outfile TEXT  Output filename (?=auto). Supports any matplotlib
                          format, e.g. ?.pdf for PDF with automatic naming.
      --help              Show this message and exit.


msnoise new_jobs
----------------
.. code-block:: sh

    msnoise new_jobs --help

    Usage:  [OPTIONS]
    
      Determines if new jobs are to be defined
    
    Options:
      -i, --init    First run ? This disables the check for existing jobs.
      --nocc        Disable the creation of CC jobs.
      --after TEXT  Create the next runnable jobs in the workflow based on DONE
                    jobs of the given config-set type/category (e.g. 'cc',
                    'stack', 'mwcs'), skipping filter steps in between. Example:
                    'msnoise new_jobs --after cc' will create STACK jobs from CC
                    jobs marked D (via CC -> filter -> stack).
      --help        Show this message and exit.


msnoise cc preprocess
---------------------
.. code-block:: sh

    msnoise cc preprocess --help

    Usage:  [OPTIONS]
    
      Run preprocessing computations on workflow jobs
    
    Options:
      --threads INTEGER  Number of threads to use for processing
      -v, --verbose      Increase verbosity
      --help             Show this message and exit.


msnoise cc compute_cc
---------------------
.. code-block:: sh

    msnoise cc compute_cc --help

    Usage:  [OPTIONS]
    
      Computes the CC jobs (based on the "New Jobs" identified)
    
    Options:
      --chunk-size INTEGER  Max pairs to process per worker per day (0 = all,
                            default). Set >0 to share a day across multiple
                            parallel workers without write conflicts (recommended
                            for >50 stations).
      --help                Show this message and exit.


msnoise cc compute_cc_rot
-------------------------
.. code-block:: sh

    msnoise cc compute_cc_rot --help

    Usage:  [OPTIONS]
    
      Computes the CC jobs too (allows for R or T components)
    
    Options:
      --help  Show this message and exit.


msnoise cc stack
----------------
.. code-block:: sh

    msnoise cc stack --help

    Usage:  [OPTIONS]
    
      Stacks the [REF] or [MOV] windows. Computes the STACK jobs.
    
    Options:
      -r, --ref   Compute the REF Stack
      -m, --mov   Compute the MOV Stacks
      -s, --step  Compute the STEP Stacks
      --help      Show this message and exit.


msnoise cc stack_refstack
-------------------------
.. code-block:: sh

    msnoise cc stack_refstack --help

    Usage:  [OPTIONS]
    
      Compute REF stacks for all pending refstack configset jobs.
    
      Reads ref_begin/ref_end from the refstack configset:
    
      - Absolute date / 1970-01-01  -> Mode A: writes a mean REF file to disk.
    
      - Negative integer (e.g. -5)  -> Mode B: no file written; rolling reference
      computed on-the-fly inside the MWCS/stretching/WCT workers.
    
    Options:
      --help  Show this message and exit.


msnoise cc plot distance
------------------------
.. code-block:: sh

    msnoise cc plot distance --help

    Usage:  [OPTIONS] [EXTRA_ARGS]...
    
      Plots the REFs of all pairs vs distance
    
    Options:
      -f, --filterid INTEGER  Filter ID
      -c, --comp TEXT         Components (ZZ, ZE, NZ, 1E,...). Defaults to ZZ
      -a, --ampli FLOAT       Amplification of the individual lines on the
                              vertical axis (default=1)
      -s, --show BOOLEAN      Show figure interactively?
      -o, --outfile TEXT      Output filename (?=auto). Supports any matplotlib
                              format, e.g. ?.pdf for PDF with automatic naming.
      -r, --refilter TEXT     Refilter CCFs before plotting (e.g. 4:8 for
                              filtering CCFs between 4.0 and 8.0 Hz. This will
                              update the plot title.
      --virtual-source TEXT   Use only pairs including this station. Format must
                              be NET.STA
      --help                  Show this message and exit.


msnoise cc plot interferogram
-----------------------------
.. code-block:: sh

    msnoise cc plot interferogram --help

    Usage:  [OPTIONS] STA1 STA2 [EXTRA_ARGS]...
    
      Plots the interferogram between sta1 and sta2 (parses the CCFs) STA1 and
      STA2 must be provided with this format: NET.STA.LOC !
    
    Options:
      -p, --preprocessid INTEGER   Preprocessing step ID
      -cc, --ccid INTEGER          CC step ID
      -f, --filterid INTEGER       Filter ID
      -m, --stackid INTEGER        Stack step ID
      -mi, --stackid_item INTEGER  Mov Stack item within that Stack step ID
      -rs, --refstackid INTEGER    REF Stack step ID
      -c, --comp TEXT              Components (ZZ, ZE, NZ, 1E,...). Defaults to ZZ
      -s, --show BOOLEAN           Show figure interactively?
      -o, --outfile TEXT           Output filename (?=auto). Supports any
                                   matplotlib format, e.g. ?.pdf for PDF with
                                   automatic naming.
      -r, --refilter TEXT          Refilter CCFs before plotting (e.g. 4:8 for
                                   filtering CCFs between 4.0 and 8.0 Hz. This
                                   will update the plot title.
      --help                       Show this message and exit.


msnoise cc plot ccftime
-----------------------
.. code-block:: sh

    msnoise cc plot ccftime --help

    Usage:  [OPTIONS] STA1 STA2 [EXTRA_ARGS]...
    
      Plots the ccf vs time between sta1 and sta2 STA1 and STA2 must be provided
      with this format: NET.STA.LOC !
    
    Options:
      -p, --preprocessid INTEGER   Preprocessing step ID
      -cc, --ccid INTEGER          CC step ID
      -f, --filterid INTEGER       Filter ID
      -m, --stackid INTEGER        Stack step ID
      -mi, --stackid_item INTEGER  Mov Stack item within that Stack step ID
      -rs, --refstackid INTEGER    REF Stack step ID
      -c, --comp TEXT              Components (ZZ, ZE, NZ, 1E,...). Defaults to ZZ
      -a, --ampli FLOAT            Amplification of the individual lines on the
                                   vertical axis (default=1)
      -S, --seismic                Seismic style: fill the space between the zero
                                   and the positive wiggles
      -s, --show BOOLEAN           Show figure interactively?
      -o, --outfile TEXT           Output filename (?=auto). Supports any
                                   matplotlib format, e.g. ?.pdf for PDF with
                                   automatic naming.
      -e, --envelope               Plot envelope instead of time series
      -r, --refilter TEXT          Refilter CCFs before plotting (e.g. 4:8 for
                                   filtering CCFs between 4.0 and 8.0 Hz. This
                                   will update the plot title.
      --normalize TEXT
      --help                       Show this message and exit.


msnoise cc plot spectime
------------------------
.. code-block:: sh

    msnoise cc plot spectime --help

    Usage:  [OPTIONS] STA1 STA2 [EXTRA_ARGS]...
    
      Plots the ccf's spectrum vs time between sta1 and sta2 STA1 and STA2 must be
      provided with this format: NET.STA.LOC !
    
    Options:
      -p, --preprocessid INTEGER   Preprocessing step ID
      -cc, --ccid INTEGER          CC step ID
      -f, --filterid INTEGER       Filter ID
      -m, --stackid INTEGER        Stack step ID
      -mi, --stackid_item INTEGER  Mov Stack item within that Stack step ID
      -rs, --refstackid INTEGER    REF Stack step ID
      -c, --comp TEXT              Components (ZZ, ZE, NZ, 1E,...). Defaults to ZZ
      -a, --ampli FLOAT            Amplification of the individual lines on the
                                   vertical axis (default=1)
      -s, --show BOOLEAN           Show figure interactively?
      -o, --outfile TEXT           Output filename (?=auto). Supports any
                                   matplotlib format, e.g. ?.pdf for PDF with
                                   automatic naming.
      -r, --refilter TEXT          Refilter CCFs before plotting (e.g. 4:8 for
                                   filtering CCFs between 4.0 and 8.0 Hz. This
                                   will update the plot title.
      --help                       Show this message and exit.


msnoise cc dtt compute_mwcs
---------------------------
.. code-block:: sh

    msnoise cc dtt compute_mwcs --help

    Usage:  [OPTIONS]
    
      Computes the MWCS jobs
    
    Options:
      --help  Show this message and exit.


msnoise cc dtt compute_mwcs_dtt
-------------------------------
.. code-block:: sh

    msnoise cc dtt compute_mwcs_dtt --help

    Usage:  [OPTIONS]
    
      Computes the dt/t jobs based on the new MWCS data
    
    Options:
      --help  Show this message and exit.


msnoise cc dtt compute_stretching
---------------------------------
.. code-block:: sh

    msnoise cc dtt compute_stretching --help

    Usage:  [OPTIONS]
    
      Computes the stretching based on the new stacked data
    
    Options:
      --help  Show this message and exit.


msnoise cc dtt compute_wct
--------------------------
.. code-block:: sh

    msnoise cc dtt compute_wct --help

    Usage:  [OPTIONS]
    
      Computes the wavelet jobs based on the new STACK data
    
    Options:
      --help  Show this message and exit.


msnoise cc dtt compute_wct_dtt
------------------------------
.. code-block:: sh

    msnoise cc dtt compute_wct_dtt --help

    Usage:  [OPTIONS]
    
      Computes dv/v from WCT results (wavelet_dtt step, lineage-based)
    
    Options:
      --help  Show this message and exit.


msnoise cc dtt plot mwcs
------------------------
.. code-block:: sh

    msnoise cc dtt plot mwcs --help

    Usage:  [OPTIONS] STA1 STA2
    
      Plots the MWCS dt and coherence images for a station pair. STA1 and STA2
      must be provided as NET.STA.LOC. Lineage: -p preprocess, -cc cc, -f filter,
      -m stack, -mi stack_item, -w mwcs.
    
    Options:
      -p, --preprocessid INTEGER   Preprocessing step ID
      -cc, --ccid INTEGER          CC step ID
      -f, --filterid INTEGER       Filter ID
      -m, --stackid INTEGER        Stack step ID
      -mi, --stackid_item INTEGER  Mov Stack item within that Stack step ID
      -rs, --refstackid INTEGER    REF Stack step ID
      -w, --mwcsid INTEGER         MWCS step set number
      -d, --dttid INTEGER          MWCS-DTT step set number
      -c, --comp TEXT              Components (ZZ, ZE, NZ, 1E,...). Defaults to ZZ
      -s, --show BOOLEAN           Show figure interactively?
      -o, --outfile TEXT           Output filename (?=auto). Supports any
                                   matplotlib format, e.g. ?.pdf for PDF with
                                   automatic naming.
      --help                       Show this message and exit.


msnoise cc dtt plot mwcs_dtt_day
--------------------------------
.. code-block:: sh

    msnoise cc dtt plot mwcs_dtt_day --help

    Usage:  [OPTIONS] STA1 STA2 DAY
    
      Plots dt against t (scatter + regression) for a single day. STA1, STA2:
      NET.STA.LOC. DAY: YYYY-MM-DD. Lineage: -p preprocess, -cc cc, -f filter, -m
      stack, -mi stack_item, -w mwcs, -d dtt.
    
    Options:
      -p, --preprocessid INTEGER   Preprocessing step ID
      -cc, --ccid INTEGER          CC step ID
      -f, --filterid INTEGER       Filter ID
      -m, --stackid INTEGER        Stack step ID
      -mi, --stackid_item INTEGER  Mov Stack item within that Stack step ID
      -rs, --refstackid INTEGER    REF Stack step ID
      -w, --mwcsid INTEGER         MWCS step set number
      -d, --dttid INTEGER          MWCS-DTT step set number
      -c, --comp TEXT              Components (ZZ, ZE, NZ, 1E,...). Defaults to ZZ
      -s, --show BOOLEAN           Show figure interactively?
      -o, --outfile TEXT           Output filename (?=auto). Supports any
                                   matplotlib format, e.g. ?.pdf for PDF with
                                   automatic naming.
      --help                       Show this message and exit.


msnoise cc dtt plot mwcs_dtt_timing
-----------------------------------
.. code-block:: sh

    msnoise cc dtt plot mwcs_dtt_timing --help

    Usage:  [OPTIONS]
    
      Plots network-mean dt/t timeseries from MWCS-DTT results. Optionally
      highlight specific pairs with -p NET.STA.LOC:NET.STA.LOC.
    
    Options:
      -f, --filterid INTEGER   Filter ID
      -w, --mwcsid INTEGER     MWCS step set number
      -d, --dttid INTEGER      MWCS-DTT step set number
      -c, --comp TEXT          Components (ZZ, ZR,...). Defaults to ZZ
      -m, --mov_stack INTEGER  Plot specific mov stack (1-based index, 0=all)
      -p, --pair TEXT          Highlight a specific pair (NET.STA.LOC:NET.STA.LOC)
      -M, --dttname TEXT       DTT column: m (slope=dt/t) or m0 (zero-intercept)
      -s, --show BOOLEAN       Show figure interactively?
      -o, --outfile TEXT       Output filename (?=auto). Supports any matplotlib
                               format, e.g. ?.pdf for PDF with automatic naming.
      --help                   Show this message and exit.


msnoise cc dtt dvv compute_mwcs_dtt_dvv
---------------------------------------
.. code-block:: sh

    msnoise cc dtt dvv compute_mwcs_dtt_dvv --help

    Usage:  [OPTIONS]
    
      Aggregate MWCS dv/v across station pairs (mwcs_dtt_dvv step).
    
    Options:
      --help  Show this message and exit.


msnoise cc dtt dvv compute_stretching_dvv
-----------------------------------------
.. code-block:: sh

    msnoise cc dtt dvv compute_stretching_dvv --help

    Usage:  [OPTIONS]
    
      Aggregate Stretching dv/v across station pairs (stretching_dvv step).
    
    Options:
      --help  Show this message and exit.


msnoise cc dtt dvv compute_wavelet_dtt_dvv
------------------------------------------
.. code-block:: sh

    msnoise cc dtt dvv compute_wavelet_dtt_dvv --help

    Usage:  [OPTIONS]
    
      Aggregate WCT dv/v across station pairs (wavelet_dtt_dvv step). Supports
      multi-band extraction.
    
    Options:
      --help  Show this message and exit.


msnoise cc dtt dvv plot mwcs_dvv
--------------------------------
.. code-block:: sh

    msnoise cc dtt dvv plot mwcs_dvv --help

    Usage:  [OPTIONS]
    
      Plot dv/v from MWCS-DTT aggregate. Requires mwcs_dtt_dvv step.
    
    Options:
      -p, --preprocessid INTEGER   Preprocessing step ID
      -cc, --ccid INTEGER          CC step ID
      -f, --filterid INTEGER       Filter ID
      -m, --stackid INTEGER        Stack step ID
      -mi, --stackid_item INTEGER  Mov Stack item within that Stack step ID
      -rs, --refstackid INTEGER    REF Stack step ID
      -w, --mwcsid INTEGER         MWCS config set number
      -c, --comp TEXT              Components (ZZ, ZE, NZ, 1E,...). Defaults to ZZ
      -wi, --mwcsdttid INTEGER     MWCS-DTT config set number
      -D, --dvvid INTEGER          DVV aggregate config set number
      -M, --dttname TEXT           DTT column: m (slope) or m0 (zero-intercept
                                   slope)
      -p, --pair_type [CC|SC|AC]   Pair type to plot (CC/SC/AC). Default: CC
      -s, --show BOOLEAN           Show figure interactively?
      -o, --outfile TEXT           Output filename (?=auto). Supports any
                                   matplotlib format, e.g. ?.pdf for PDF with
                                   automatic naming.
      --help                       Show this message and exit.


msnoise cc dtt dvv plot stretching_dvv
--------------------------------------
.. code-block:: sh

    msnoise cc dtt dvv plot stretching_dvv --help

    Usage:  [OPTIONS]
    
      Plot dv/v from Stretching aggregate. Requires stretching_dvv step.
    
    Options:
      -f, --filterid INTEGER      Filter ID
      -S, --stretchingid INTEGER  Stretching config set number
      -D, --dvvid INTEGER         DVV aggregate config set number
      -c, --comp TEXT             Components (ZZ, ZE, NZ, 1E,...). Defaults to ZZ
      -m, --mov_stack INTEGER     Plot specific mov stack (1-based index, 0=all)
      -p, --pair_type [CC|SC|AC]  Pair type to plot (CC/SC/AC). Default: CC
      -s, --show BOOLEAN          Show figure interactively?
      -o, --outfile TEXT          Output filename (?=auto). Supports any
                                  matplotlib format, e.g. ?.pdf for PDF with
                                  automatic naming.
      --help                      Show this message and exit.


msnoise cc dtt dvv plot wavelet_dvv
-----------------------------------
.. code-block:: sh

    msnoise cc dtt dvv plot wavelet_dvv --help

    Usage:  [OPTIONS]
    
      Plot dv/v from WCT-DTT aggregate. Use -v heatmap for per-pair frequency
      view.
    
    Options:
      -f, --filterid INTEGER          Filter ID
      -c, --comp TEXT                 Components (ZZ, ZE, NZ, 1E,...). Defaults to
                                      ZZ
      -m, --mov_stack INTEGER         Plot specific mov stack (1-based index,
                                      0=all)
      -w, --wctid INTEGER             WCT config set number
      -d, --wctdttid INTEGER          WCT-DTT config set number
      -D, --dvvid INTEGER             DVV aggregate config set number
      -p, --pair_type [CC|SC|AC]      Pair type to plot (CC/SC/AC). Default: CC
      -v, --visualize [timeseries|heatmap]
                                      Plot style: timeseries (uses dvv aggregate)
                                      or heatmap (uses per-pair data)
      -r, --ranges TEXT               Frequency ranges for band averaging (first
                                      range used for timeseries)
      -s, --show BOOLEAN              Show figure interactively?
      -o, --outfile TEXT              Output filename (?=auto). Supports any
                                      matplotlib format, e.g. ?.pdf for PDF with
                                      automatic naming.
      --help                          Show this message and exit.


msnoise qc compute_psd
----------------------
.. code-block:: sh

    msnoise qc compute_psd --help

    Usage:  [OPTIONS]
    
      Computes the PSD jobs, saves results as NetCDF files. Based on New or
      Modified files identified by the new_jobs step.
    
    Options:
      --chunk-size INTEGER  Max stations to process per worker per day (0 = all,
                            default). Set >0 to share a day across multiple
                            parallel workers without write conflicts (recommended
                            for >50 stations).
      --help                Show this message and exit.


msnoise qc plot_psd
-------------------
.. code-block:: sh

    msnoise qc plot_psd --help

    Usage:  [OPTIONS] SEED_ID
    
      Plots the PSD and spectrogram based on NPZ files
    
    Options:
      --help  Show this message and exit.


msnoise qc compute_psd_rms
--------------------------
.. code-block:: sh

    msnoise qc compute_psd_rms --help

    Usage:  [OPTIONS]
    
      Computes the RMS from PSD NetCDF files.
    
    Options:
      --help  Show this message and exit.


msnoise qc plot_psd_rms
-----------------------
.. code-block:: sh

    msnoise qc plot_psd_rms --help

    Usage:  [OPTIONS] [SEED_ID...]
    
      Plot PSD-RMS results: time-series, clock plots, hour-maps, grid-maps.
    
      SEED_ID can be zero or more SEED identifiers (NET.STA.LOC.CHAN).  When
      omitted, all stations found on disk for the given psd/psd_rms lineage are
      plotted.  Uses MSNoiseResult for lineage-correct path resolution.
    
          msnoise qc plot_psd_rms BE.UCC..HHZ
          msnoise qc plot_psd_rms                           # all stations
          msnoise qc plot_psd_rms BE.UCC..HHZ BE.MEM..HHZ --type timeseries
          msnoise qc plot_psd_rms BE.UCC..HHZ --type clockplot --timezone Europe/Brussels
          msnoise qc plot_psd_rms BE.UCC..HHZ --type all --outfile noise.pdf
          msnoise qc plot_psd_rms --psd-id 2 --psd-rms-id 2 BE.UCC..HHZ
          msnoise qc plot_psd_rms BE.UCC..HHZ --day-start 8 --day-end 20
          msnoise qc plot_psd_rms BE.UCC..HHZ --day-start 0 --day-end 24  # disable
          msnoise qc plot_psd_rms BE.UCC..HHZ --resample-freq 1h --agg-func median
          msnoise qc plot_psd_rms BE.UCC..HHZ -R 15min -A max
    
    Options:
      --psd-id INTEGER                Config-set number for the psd step (matches
                                      psd_N in workflow)  [default: 1]
      --psd-rms-id INTEGER            Config-set number for the psd_rms step
                                      (matches psd_rms_N in workflow)  [default:
                                      1]
      -t, --type [timeseries|clockplot|hourmap|gridmap|dailyplot|all]
                                      Plot type to produce  [default: timeseries]
      -b, --band TEXT                 Frequency band label, e.g. "1.0-10.0"
                                      (default: first available)
      --scale FLOAT                   Amplitude scale factor applied before
                                      display  [default: 1000000000.0]
      -u, --unit TEXT                 Amplitude unit label for axis/colorbar
                                      [default: nm]
      -z, --timezone TEXT             Time zone for local-time plots (e.g.
                                      Europe/Brussels)  [default: UTC]
      --day-start FLOAT               Start of daytime in decimal hours (local
                                      time), e.g. 7.0 = 07:00  [default: 7.0]
      --day-end FLOAT                 End of daytime in decimal hours (local
                                      time), e.g. 19.0 = 19:00  [default: 19.0]
      --split-date TEXT               ISO date (YYYY-MM-DD) for before/after split
                                      (clockplot, dailyplot)
      -a, --annotate TEXT             Event annotations as DATE=LABEL,DATE=LABEL
      -R, --resample-freq TEXT        Pandas offset string for time-bin resampling
                                      (e.g. "15min", "1h", "2h").  Default: 30min
                                      [default: 30min]
      -A, --agg-func TEXT             Aggregation method for each resampled bin:
                                      mean, median, max, min, std, sum, or any
                                      pandas agg string.  Default: mean  [default:
                                      mean]
      -o, --outfile TEXT              Base filename for saved figure(s). Type
                                      suffix is appended.
      --no-show                       Do not call plt.show() (useful in non-
                                      interactive environments)
      --help                          Show this message and exit.


msnoise utils bugreport
-----------------------
.. code-block:: sh

    msnoise utils bugreport --help

    Usage:  [OPTIONS]
    
      This command launches the Bug Report script.
    
    Options:
      -s, --sys      System Info
      -m, --modules  Modules Info
      -e, --env      Environment Info
      -a, --all      All Info
      --help         Show this message and exit.


msnoise utils test
------------------
.. code-block:: sh

    msnoise utils test --help

    Usage:  [OPTIONS]
    
      Runs the test suite in a temporary folder.
    
      Use --fast for a quick smoke test that exercises the full workflow lifecycle
      using stub compute functions (no seismic data required, completes in < 30
      seconds).
    
    Options:
      -p, --prefix TEXT  Prefix for tables
      --tech INTEGER     Test using (1) SQLite or (2) MariaDB
      -c, --content      Run content tests instead of standard tests
      --fast             Run the fast smoke test suite (no real data needed)
      --help             Show this message and exit.


msnoise utils jupyter
---------------------
.. code-block:: sh

    msnoise utils jupyter --help

    Usage:  [OPTIONS]
    
      Launches an jupyter notebook in the current folder
    
    Options:
      --help  Show this message and exit.


msnoise utils export-params
---------------------------
.. code-block:: sh

    msnoise utils export-params --help

    Usage:  [OPTIONS]
    
      Export the full layered parameter chain for a lineage to YAML.
    
      Either supply --lineage as a slash-separated string, or build it from
      individual step IDs.  The resulting YAML contains one block per config
      category in lineage order — no key collisions, fully self-describing for
      reproducibility.
    
      Examples::
    
          msnoise utils export-params --lineage preprocess_1/cc_1/filter_1/stack_1
          msnoise utils export-params -p 1 -cc 1 -f 1 -m 1 -w 1 -wi 1     msnoise
          utils export-params -p 1 -cc 1 -f 1 -m 1 -S 1 -sd 1     msnoise utils
          export-params -p 1 -cc 1 -f 1 -m 1 -W 1 -wd 1 -wdv 1
    
    Options:
      -l, --lineage TEXT              Lineage string e.g. "preprocess_1/cc_1/filte
                                      r_1/stack_1/mwcs_1/mwcs_dtt_1"
      -p, --preprocessid INTEGER      Preprocessing step ID
      -cc, --ccid INTEGER             CC step ID
      -f, --filterid INTEGER          Filter ID
      -m, --stackid INTEGER           Stack step ID
      -rs, --refstackid INTEGER       REF Stack step ID (optional)
      -w, --mwcsid INTEGER            MWCS step ID (optional)
      -wi, --mwcsdttid INTEGER        MWCS-DTT step ID (optional)
      -S, --stretchingid INTEGER      Stretching step ID (optional)
      -sd, --stretchingdvvid INTEGER  Stretching-DVV step ID (optional)
      -W, --wctid INTEGER             WCT step ID (optional)
      -wd, --wctdttid INTEGER         WCT-DTT step ID (optional)
      -wdv, --waveletdvvid INTEGER    Wavelet-DVV step ID (optional)
      -o, --output TEXT               Output YAML path (default:
                                      params_<lineage>.yaml in current dir)
      --help                          Show this message and exit.


msnoise utils export-dvv
------------------------
.. code-block:: sh

    msnoise utils export-dvv --help

    Usage:  [OPTIONS]
    
      Export dv/v time series as NetCDF with embedded parameter provenance.
    
      Each output file contains dv/v statistics (mean, std, median, n_pairs,
      weighted/trimmed variants) on a ``times`` dimension, with the complete
      processing parameter chain embedded as a global YAML attribute for full
      reproducibility.
    
      Either supply --lineage as a slash-separated string ending at a DVV step, or
      build it from individual step IDs.
    
      Examples::
    
          # MWCS dv/v, all components and mov_stacks, output to current directory
          msnoise utils export-dvv -p 1 -cc 1 -f 1 -s 1 -r 1 -m 1 -md 1 -mdd 1
    
          # Stretching dv/v, ZZ only, custom output directory     msnoise utils
          export-dvv -p 1 -cc 1 -f 1 -s 1 -r 1 -S 1 -sd 1             --components
          ZZ --output /data/release/
    
          # Using a lineage string     msnoise utils export-dvv
          --lineage preprocess_1/cc_1/filter_1/stack_1/refstack_1/mwcs_1/mwcs_dtt_
          1/mwcs_dtt_dvv_1
    
          # Reload the result in Python     # >>> import xarray as xr, yaml     #
          >>> ds = xr.open_dataset("dvv_CC_ZZ__pre1-...nc")     # >>> params =
          yaml.safe_load(ds.attrs["msnoise_params"])     # >>>
          print(params["mwcs"]["mwcs_wlen"])
    
    Options:
      -l, --lineage TEXT              Slash-separated lineage string ending at a
                                      DVV step
      -p, --preprocessid INTEGER
      -cc, --ccid INTEGER
      -f, --filterid INTEGER
      -s, --stackid INTEGER
      -r, --refstackid INTEGER
      -m, --mwcsid INTEGER
      -md, --mwcsdttid INTEGER
      -mdd, --mwcsdvvid INTEGER
      -S, --stretchingid INTEGER
      -sd, --stretchingdvvid INTEGER
      -W, --wctid INTEGER
      -wd, --wctdttid INTEGER
      -wdv, --waveletdvvid INTEGER
      --pair-type TEXT                [default: CC]
      --components TEXT               Component to export (e.g. ZZ). Default: all.
      --mov-stack TEXT                Moving stack as window_step (e.g. 1D_1D).
                                      Default: all.
      -o, --output TEXT               Output directory (or full .nc path for a
                                      single export)  [default: .]
      --help                          Show this message and exit.


msnoise utils create_preprocess_jobs
------------------------------------
.. code-block:: sh

    msnoise utils create_preprocess_jobs --help

    Usage:  [OPTIONS]
    
      Create preprocess T jobs for FDSN/EIDA sources (bypasses scan_archive).
    
      For local/SDS sources, only creates jobs where DataAvailability records
      already exist for the requested date(s).  For FDSN/EIDA sources, creates
      jobs for all active stations regardless (the preprocess worker fetches on
      demand).
    
      Examples::
    
          msnoise utils create_preprocess_jobs --date 2026-03-28
    
          msnoise utils create_preprocess_jobs --date_range 2026-01-01 2026-03-28
    
    Options:
      --date TEXT             Single date to create jobs for (YYYY-MM-DD).
      --date_range START END  Date range: START END (inclusive, YYYY-MM-DD).
      --set-number INTEGER    Preprocess config set number (default 1).
      --help                  Show this message and exit.


msnoise utils create_psd_jobs
-----------------------------
.. code-block:: sh

    msnoise utils create_psd_jobs --help

    Usage:  [OPTIONS]
    
      Create psd T jobs for FDSN/EIDA sources (bypasses scan_archive).
    
      For FDSN/EIDA sources, creates jobs for all active stations for the
      requested date(s) regardless of DataAvailability (the psd worker fetches on
      demand).  For local/SDS sources, only creates jobs where DataAvailability
      records already exist for the requested date(s).
    
      Examples::
    
          msnoise utils create_psd_jobs --date 2026-03-28
    
          msnoise utils create_psd_jobs --date_range 2026-01-01 2026-03-28
    
    Options:
      --date TEXT             Single date to create jobs for (YYYY-MM-DD).
      --date_range START END  Date range: START END (inclusive, YYYY-MM-DD).
      --set-number INTEGER    PSD config set number (default 1).
      --help                  Show this message and exit.


msnoise utils import-stationxml
-------------------------------
.. code-block:: sh

    msnoise utils import-stationxml --help

    Usage:  [OPTIONS] SOURCE
    
      Import stations from a StationXML file or FDSN URL.
    
      SOURCE can be a local file path or a URL (e.g. an FDSN station web service
      query with level=channel).
    
      The parsed inventory is written to the project's ``response_path`` directory
      by default, making instrument responses immediately available to the
      preprocessing step.  Use --no-save to skip this.
    
      Examples:
        msnoise utils import-stationxml inventory.xml
        msnoise utils import-stationxml https://eida.ethz.ch/fdsnws/station/1/query?...
        msnoise utils import-stationxml inventory.xml --data-source-id 2 --no-save
    
    Options:
      -d, --data-source-id INTEGER  DataSource ID to assign to imported stations
                                    (default: project default).
      --no-save                     Do NOT save the inventory to response_path
                                    after import.
      --help                        Show this message and exit.


