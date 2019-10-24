Package Installation
=====================

To install HIPSTER, simply clone the Github repository and compile the C++ code (see also the :ref:`dependencies` below). This is done as follows::

    git clone https://github.com/oliverphilcox/HIPSTER.git
    cd HIPSTER
    make

**NB**: HIPSTER can be run for either *periodic* or *aperiodic* data-sets. This is specified by adding compiler flags in the ``Makefile``. See :doc:`basic-usage` for more information.

.. todo:: add links here

Once HIPSTER is installed, it can be run via a simple bash wrapper, as described in the :doc:`basic-usage` section. This requires the inputs specified in :doc:`pre-processing`. For advanced usage, the routines in :doc:`advanced-usage` can be used.

.. _dependencies:

Dependencies
-------------

Hipster requires the following packages:

- [C compiler](https://gcc.gnu.org/): Tested with gcc 5.4.O
- [Gnu Scientific Library (GSL)](https://www.gnu.org/software/gsl/doc/html/index.html): Any recent version
- [Corrfunc](https://corrfunc.readthedocs.io): 2.0 or later (required for aperiodic surveys to compute geometry correction)
- [OpenMP](https://www.openmp.org/): Any recent version (optional, but required for parallelization)
- [Python](https://www.python.org/): 2.7 or later, 3.4 or later (required for pre- and post-processing)

Corrfunc can be installed using ``pip install corrfunc`` and is used for efficient pair counting.

Note that many of the code modules and convenience functions are based on those of [RascalC](https://RascalC.readthedocs.io), developed by Oliver Philcox, Daniel Eisenstein, Ross O'Connell and Alexander Wiegand.

Acknowledgements
-----------------

Authors:

- Oliver H. E. Philcox (Department of Astrophysical Sciences, Princeton)
- Daniel J. Eisenstein (Center for Astrophysics | Harvard & Smithsonian)

Citing:

.. todo:: add citing!!
