.. eeg-ride documentation master file, created by
   sphinx-quickstart on Mon Jul 22 17:49:34 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

eeg-ride
========

.. image:: https://img.shields.io/pypi/v/eeg-ride
   :target: https://pypi.org/project/eeg-ride
   :alt: Latest Version

.. image:: https://img.shields.io/pypi/pyversions/eeg-ride.svg
   :target: https://img.shields.io/pypi/pyversions/eeg-ride
   :alt: PyPI - Python Version

.. image:: https://img.shields.io/github/license/kirstenstark/eeg-ride
   :target: https://github.com/kirstenstark/eeg-ride/blob/main/LICENSE
   :alt: License

|

Separating EEG data into stimulus- and response-related components using Residue Iteration Decomposition (RIDE).

Based on the `MATLAB toolbox <https://cns.hkbu.edu.hk/RIDE.htm>`_ described in:

  Ouyang, G., Sommer W., & Zhou, C. (2015).
  A toolbox for residue iteration decomposition (RIDE)--A method for the decomposition, reconstruction, and single trial analysis of event related potentials.
  *Journal of Neuroscience Methods*, *250*, 7-21. 
  `https://doi.org/10.1016/j.jneumeth.2014.10.009 <https://doi.org/10.1016/j.jneumeth.2014.10.009>`_

One typical application is for correction of speech artifacts as described and validated in:

  Ouyang, G., Sommer, W., Zhou, C., Aristei, S., Pinkpank, T., & Abdel Rahman, R. (2016). 
  Articulation artifacts during overt language production in event-related brain potentials: Description and correction. 
  *Brain topography*, *29*, 791-813. 
  `https://doi.org/10.1007/s10548-016-0515-1 <https://doi.org/10.1007/s10548-016-0515-1>`_
  

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
