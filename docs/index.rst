HIPSTER: High-k Power SpecTrum EstimatoR
========================================

Overview
----------

HIPSTER is a code to quickly compute small-scale galaxy (or halo) power spectra in arbitary survey geometries, based on the work of Philcox & Eisenstein (2019b, submitted). This computes the Legendre multipoles of the power spectrum in *configuration space*, by computing weighted pair counts over the survey truncated at some maximum radius. This fully accounts for window function effects, does not include shot-noise, and is optimized for small-scale power spectrum computation.

The source code is publicly available on `Github <https://github.com/oliverphilcox/HIPSTER>`_, and contains many modules modified from the `RascalC` <https://github.com/oliverphilcox/HIPSTER>`_ covariance matrix code. Full documentation is provided in the accompanying pages.

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
