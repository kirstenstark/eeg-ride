Quickstart
==========

The pipeline provides a single high-level function, ``group_pipeline()``, to carry out a full EEG analysis on a group of participants.

Here is a fairly minimal example for a (fictional) N400/P600 experiment with two experimental factors: ``semantics`` (e.g., related versus unrelated words) and emotional ``context`` (e.g., emotionally negative versus neutral).

Code block 1: Read in data

Code block 2: define cfg

Code block 3: Call ride 


.. code-block:: python

    from pipeline import group_pipeline

    cfg = RideCfg(comp_name=['s', 'r'],
                  comp_twd=[[0, 600], [-300, 300]],
                  comp_latency=[0, rt],
                  sfreq=500)

            

In this example we have specified:

- ``raw_files``, ``log_files``, ``output_dir``, ``besa_files``: The paths to the raw EEG data, to the behavioral log files, to the desired output directory, and to the BESA files for ocular correction

- ``triggers``: The four different numerical EEG trigger codes corresponding to each of the four cells in the 2 × 2 design

- ``skip_log_conditions``: Our log files may contain additional trials from a "filler" condition without corresponding EEG trials/triggers. These filler trials are marked with the condition label ``'filler'`` in the log file column ``semantics``

- ``components``: The *a priori* defined time windows and regions of interest for the relevant ERP components (N400 and P600)

- ``average_by``: The relevant groupings of trials for which by-participant averaged waveforms should be created. The keys (e.g., ``'related_negative'``) are custom labels of our choice; the values are the corresponding logical conditions that must be met for a trial to be included in the average.

For (way) more options, see :doc:`Pipeline inputs <inputs_py>`.
