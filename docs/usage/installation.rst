Package Installation
=====================

To install HIPSTER on a Linux machine, simply clone the Github repository::

    git clone https://github.com/oliverphilcox/HIPSTER.git

Assuming the :ref:`dependencies` are satisfied, HIPSTER is now ready to run using the simple bash wrappers ``hipster_wrapper.sh``, ``hipster_wrapper_periodic.sh`` and ``hipster_wrapper_bispectrum.sh``, as described in the :doc:`basic-usage` section. This requires the inputs specified in :doc:`pre-processing`. For advanced usage, the routines in :doc:`advanced-usage` can be used.

For advanced users, the C++ code must also be compiled via::

    cd HIPSTER
    make [Periodic=-DPERIODIC] [Bispectrum=-DBISPECTRUM]

HIPSTER can be run for either *periodic* or *aperiodic* data-sets. This is specified by adding compiler flags in the Makefile, or can simply be activated by adding the line ``Periodic=-DPERIODIC`` to the ``make`` command. See :ref:`periodicity-note` for more information. Similarly the bispectrum mode is activated with ``Bispectrum=-DBISPECTRUM``. Note that compilation into the correct format is done automatically if the main wrappers are used.

**Note for Mac Users**: HIPSTER is primarily designed for Linux machines, though running on a Mac is also possible. To do so, you must have a recent version of GCC and `gnu-getopt <(http://macappstore.org/gnu-getopt/)>`_ (emulating the ``getopt`` Linux script). Furthermore the code must be compiled without OpenMP, by removing the ``-DOPENMP`` and ``-lgomp`` flags from the Makefile.

.. _dependencies:

Dependencies
-------------

HIPSTER requires the following (often pre-installed) packages:

- `C compiler <https://gcc.gnu.org/>`_: Tested with gcc 5.4.O
- `Gnu Scientific Library (GSL) <https://www.gnu.org/software/gsl/doc/html/index.html>`_: Any recent version above 5.1 (needed for C++11)
- `OpenMP <https://www.openmp.org/>`_: Any recent version (optional, but required for parallelization)
- `Python <(https://www.python.org/>`_: 2.7 or later, 3.4 or later (required for pre- and post-processing)

Additionally, for *periodic* surveys, we require `Corrfunc <https://corrfunc.readthedocs.io>`_ (v2.0 or later) to compute geometry correction functions via efficient pair counting. This can be installed using ``pip install corrfunc``.

Acknowledgements
-----------------

Main Authors:

- Oliver H. E. Philcox (Princeton / Harvard)
- Daniel J. Eisenstein (Harvard)

Additional Collaborators:

- David N. Spergel (CCA / Princeton)
- Francisco Villaescusa-Navarro (CCA / Princeton)
- Lehman Garrison (CCA / Harvard)
- Zachary Slepian (Florida)

Please cite the initial theory paper (Philcox & Eisenstein 2019, accepted by MNRAS, `arXiv <https://arxiv.org/abs/1912.01010>`_) and the periodic power spectrum and bispectrum paper (Philcox 2020, submitted to MNRAS, `arXiv <https://arxiv.org/pdf/2005.01739.pdf>`_) when using this code in your research.

Note that many of the code modules and convenience functions are based on those of `RascalC <https://RascalC.readthedocs.io>`_, developed by Oliver Philcox, Daniel Eisenstein, Ross O'Connell and Alexander Wiegand.
