HIPSTER: HIgh-k Power SpecTrum and bispectrum EstimatoR
========================================================

Overview
----------

HIPSTER is a code to quickly compute small-scale power spectra and isotropic bispectra for galaxy surveys and cosmological simulations, based on the work of Philcox & Eisenstein (2019, MNRAS, `arXiv <https://arxiv.org/abs/1912.01010>`_) and Philcox (2020, submitted, `arXiv <https://arxiv.org/pdf/2005.01739.pdf>`_). This computes the Legendre multipoles of the power spectrum, :math:`P_\ell(k)` or bispectrum :math:`B_\ell(k_1,k_2)` in *configuration space*, by computing weighted pair counts over the survey or simulation box truncated at some maximum radius. This does not include shot-noise or aliasing, fully accounts for window function effects and is optimized for small-scale power spectrum and bispectrum computation. Note that the bispectrum code algorithm has the same complexity as the power spectrum algorithm, thanks to spherical harmonic decompositions. By combining HIPSTER with conventional FFT-based methods, spectra can be measured efficiency across a vast range of wavenumbers.

The code can be run either in *aperiodic* or *periodic* mode, for galaxy surveys or N-body simulations respectively. The *periodic* mode contains various optimizations relating to the periodic geometry, as detailed in the second paper. For the bispectrum, only *periodic* mode is currently supported, though the alternative is easily added. Generalization to anisotropic bispectra is straightforward (and requires no additional computing time) and can be added on request.

The source code is publicly available on `Github <https://github.com/oliverphilcox/HIPSTER>`_, and contains many modules modified from the `RascalC <https://github.com/oliverphilcox/HIPSTER>`_ covariance matrix code.

To compute a periodic matter power spectrum up to :math:`\ell=L` from a particles in a simulation box (``data.dat``), with pair counts truncated at radius ``R0`` with :math:`k`-space binning file ``binning.csv`` on 4 CPU-cores, we simply run::

  ./hipster_wrapper_periodic.sh --dat data.dat --l_max L -R0 R0 -k_bin binning.csv --nthreads 4

For a galaxy power spectrum with a non-trivial survey geometry, we also need random particle files (``randoms.dat``) and the usage is simply::

    ./hipster_wrapper.sh --dat data.dat --ran_DR randoms.dat --ran_RR randoms.dat -l_max L -R0 R0 -k_bin binning.csv --nthreads 4

To compute an isotropic bispectrum up to :math:`\ell=L` from a particles in a simulation box (``data.dat``), with pair counts truncated at radius ``R0`` with :math:`k`-space binning file ``binning.csv`` on 4 CPU-cores, we simply run::

      ./hipster_wrapper_bispectrum.sh --dat data.dat --l_max L -R0 R0 -k_bin binning.csv --nthreads 4

This is described in detail in the accompanying pages.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage/installation
   usage/pre-processing
   usage/basic-usage
   usage/advanced-usage

For any queries regarding the code please contact `Oliver Philcox  <mailto:ohep2@cantab.ac.uk>`_.


.. Indices and tables
.. ==================
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
