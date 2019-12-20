HIPSTER: HIgh-k Power SpecTrum EstimatoR
=========================================

Overview
----------

HIPSTER is a code to quickly compute small-scale power spectra for galaxy surveys and N-body simulations, based on the work of Philcox & Eisenstein (2019, accepted by MNRAS, `arXiv <https://arxiv.org/abs/1912.01010>`_) and Philcox (2020, in prep.). This computes the Legendre multipoles of the power spectrum, :math:`P_\ell(k)` in *configuration space*, by computing weighted pair counts over the survey or simulation box truncated at some maximum radius. This fully accounts for window function effects, does not include shot-noise, and is optimized for small-scale power spectrum computation.

The code can be run either in 'aperiodic' or 'periodic' mode, for galaxy surveys or N-body simulations respectively. The 'periodic' mode contains various optimizations relating to the periodic geometry, as detailed in the second paper.

The source code is publicly available on `Github <https://github.com/oliverphilcox/HIPSTER>`_, and contains many modules modified from the `RascalC <https://github.com/oliverphilcox/HIPSTER>`_ covariance matrix code.

To compute a 'periodic' matter power spectrum up to :math:`\ell=L` from a particles in a simulation box (``data.dat), with pair counts truncated at radius ``R0`` with :math:`k`-space binning file ``binning.csv`` we simply run::

  ./hipster_wrapper_periodic.sh --dat data.dat --l_max L -R0 R0 -k_bin binning.csv

For a galaxy power spectrum with a non-trivial survey geometry, we also need random particle files (``randoms.dat``) and the usage is simply::

    ./hipster_wrapper.sh --dat data.dat --ran_DR randoms.dat --ran_RR randoms.dat -l_max L -R0 R0 -k_bin binning.csv

This is described in detail in the accompanying pages.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage/installation
   usage/pre-processing
   usage/basic-usage
   usage/advanced-usage

For any queries regarding the code please contact `Oliver Philcox  <mailto:ohep2@alumni.cam.ac.uk>`_.


.. Indices and tables
.. ==================
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
