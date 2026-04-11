.. include:: ../configs.hrst

.. _workflow_dvv:

**************************************
Workflow: Relative Velocity Variations
**************************************

This section presents the three methods available to compute dv/v,
assuming the :ref:`workflow_ccf` steps have been completed.

Each method is self-contained.  Run whichever chain(s) suit your use case.


Method 1 — MWCS (3 steps)
---------------------------

Moving-Window Cross-Spectrum → dt/t → dv/v aggregate.

.. toctree::
   :maxdepth: 1

   ../workflow_dvv/007_compute_mwcs
   ../workflow_dvv/007b_plot_mwcs
   ../workflow_dvv/008_compute_dtt
   ../workflow_dvv/009_compute_dvv
   ../workflow_dvv/009b_plot_dvv


Method 2 — Stretching (1 step)
--------------------------------

Direct dv/v by cross-correlating stretched reference traces.

.. toctree::
   :maxdepth: 1

   ../workflow_dvv/011_compute_stretching
   ../workflow_dvv/011b_plot_stretching


Method 3 — Wavelet Coherence Transform (2 steps, or fused in 1)
-----------------------------------------------------------------

WCT → dt/t → dv/v aggregate.  When ``|wavelet.wct_compute_dtt|`` is
``True`` (default), the dt/t step runs inline inside the WCT step and no
intermediate files are written.

.. toctree::
   :maxdepth: 1

   ../workflow_dvv/010_compute_wct
   ../workflow_dvv/010b_compute_wct_dtt
   ../workflow_dvv/010c_plot_wct
