Tips & tricks
=============

This page contains a collection of tips and tricks for using the RIDE toolbox.

Preprocessing and RIDE
----------------------

The EEG epochs that are fed into RIDE should be preprocessed and as "clean" as possible [#]_.
For instance, if RIDE is used to remove speech artifacts related to some verbal response, prior steps should be taken to remove any non-speech-related artifacts before using RIDE.
Typically, this involves steps like frequency-domain filtering, re-referencing, and removing ocular artifacts (e.g., using independent component analysis).
See [#]_, [#]_, and [#]_ for further information on how to preprocess EEG data.

One step where the preprocessing prior to RIDE may differ from your usual preprocessing pipeline is when choosing the length of your epochs during segmentation.
Oftentimes, researchers use epoch durations that end at 0.8, 1, or 2 s after stimulus onset, since this is usually long enough to capture the ERP of interest.
However, some speech artifacts will last significantly longer than this.
Since it is helpful for the performance of RIDE to capture the entire speech artifact, you may choose an epoch duration that is longer than what is typically used (e.g., ending at 3 s after stimulus onset).
You can use a time course or butterfly plot of your epochs to determine how long the speech artifacts are in your own data (that is, at what point the signal at most channels is back at baseline).

To save disk space and speed up subsequent analysis, you can shorten your epochs to the usual length (e.g., ending at 1 s after stimulus onset) *after* you have performed the RIDE correction (e.g., using MNE's ``Epochs.crop()`` method).

RIDE and artifact rejection
---------------------------

Many EEG preprocessing pipelines involve the rejection (i.e., deletion) of "bad" epochs, typically using some kind of amplitude-based threshold.
For instance, the ``reject`` argument in MNE's ``Epochs`` class uses a peak-to-peak threshold that, if exceeded at any one channel, marks the epoch as "bad" in order to exclude it from the downstream analysis.

You typically want to perform artifact rejection *before* feeding the epochs into RIDE, since RIDE works best with epochs that have already been cleaned from all non-speech-related artifacts (as discussed above).
However, note that the speech artifacts themselves can be rather large.
Therefore, you may want to use a rather lenient artifact rejection threshold (e.g., 200/250 µV) before RIDE correction, or otherwise you will lose many epochs and therefore reduce your statistical power---even though RIDE would have been able to "rescue" these epochs by removing the speech artifacts.

One way to deal with this is to perform amplitude-based rejection twice when using RIDE: once *before* the RIDE correction, using a rather lenient threshold to keep as many epochs in as possible, and once *after* RIDE correction, using a more stringent threshold to remove any remaining non-speech-related artifacts.

Another possibility is to run artifact rejection before RIDE correction and then use the results of RIDE to correct *all* the available trials, including the rejected ones.
Then you can repeat the rejection step and you will likely end up with more trials than before because RIDE has removed the speech artifacts from some of the rejected trials.
We have implemented this procedure in the ``hu-neuro-pipeline``` package (more details on that below).

Time windows of RIDE components
-------------------------------

You as the user need to specify the time window (onset and offset) in which RIDE will look for stimulus- and response-related components.
It is important that all relevant activity is contained within this time window.
Especially for the ``'r'`` component, it is important that the time window is long enough to capture the speech artifact entirely.

The default time windows from the RIDE toolbox in MATLAB are:

* ``[0, 600]`` for the ``'s'`` component, that is, 0 to 600 ms after to stimulus onset and

* ``[-300, 300]`` for the ``'r'`` component, that is, 300 ms before to 300 ms after response onset.

However, if the speech utterances in your own data are relatively long, you may need to specify a---sometimes significantly---longer time window for the ``'r'`` component.

Differences to the MATLAB version
---------------------------------

We tried to keep the current Python implementation as similar as possible to the original MATLAB toolbox [#]_.
However, there are currently a few differences that might be worth noting:

* The Python version currently only supports decomposition into ``s`` and ``r`` components (that is, stimulus- and response-related components), whereas the MATLAB version also supports one or more intermediate components (``c`` components).
  We plan to add support for ``c`` components in a future version of the package.

* The Python version also currently lacks support for microsaccade identification, which is included in the MATLAB version.

* The Python version contains a new function (``ride.correct_trials``) that can be used to correct spech artifacts on the single trial level.
  Specifically, the function takes as its inputs the fitted results from the RIDE correction and a set of EEG epochs (typically the same that were fed into RIDE).
  It then subtracts the ``'r'`` component (thought to captrue the speech artifact) from each epoch, shifted by the specific response onset for that epoch.

Integration in other packages
-----------------------------

RIDE is included in the ``hu-neuro-pipeline`` Python package [#]_ [#]_, which provides a user-friendly end-to-end pipeline for processing EEG data.
This package contains some reasonable default parameters for RIDE and implements all of the best practices discussed above.
It therefore might be a good starting point if you want to use RIDE for your own EEG data and do not yet have your own preprocessing pipeline in place.

Notes
-----

.. [#] Ouyang, G., Sommer, W., & Zhou, C. (2015). A toolbox for residue iteration decomposition (RIDE)—A method for the decomposition, reconstruction, and single trial analysis of event related potentials. *Journal of Neuroscience Methods*, 250, 7–21. https://doi.org/10.1016/j.jneumeth.2014.10.009
.. [#] Luck, S. J. (2014). *An Introduction to the Event-Related Potential Technique*. 2 ed. Cambridge, Massachusetts: MIT Press.
.. [#] Delorme, A. (2023). EEG is better left alone. *Scientific Reports*, 13(1), 2372.
.. [#] https://alexenge.github.io/intro-to-eeg/ipynb/2-preprocessing.html
.. [#] https://cns.hkbu.edu.hk/RIDE.htm
.. [#] https://hu-neuro-pipeline.readthedocs.io
.. [#] https://github.com/alexenge/hu-neuro-pipeline
