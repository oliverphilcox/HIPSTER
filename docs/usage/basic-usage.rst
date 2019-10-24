Covariance Matrix Estimation
=============================

Overview
----------

This is the main section of RascalC, where 2PCF or 3PCF covariance matrix estimates are computed via Monte Carlo integration from a given set of input particles. For the 2PCF, depending on the number of input fields the code will compute either components for a single covariance matrix or all required components for 6 cross-covariance matrices (i.e. for multi-tracer analysis).

**Prerequisites**:
- In JACKKNIFE mode, the jackknife weights and binned pair counts must be computed via the :doc:`jackknife-weights` script before the C++ code is run.
- In DEFAULT mode, the RR pair counts must be computed via the :doc:`geometry-correction` script before the C++ code is run.
- In LEGENDRE and 3PCF modes, the survey correction functions :math:`\Phi` must be computed via the :doc:`geometry-correction` script before the C++ code is run.

.. _particle-grid:

Particle Grid and Cells
~~~~~~~~~~~~~~~~~~~~~~~~~



# Wrapper
Note that the code must be compiled in the correct manner before the wrapper is used.



# Note on choice of R0
Add time scaling and errors



The Makefile may need to be modified depending on the specific computational system. In addition, by specifying the ``-DPERIODIC`` flag in the Makefile, the code will assume periodic boundary conditions (as appropriate for periodic box simulations), measuring the angle $\mu$ from the z-axis. Furthermore, by removing the ``-DOPENMP`` flag, the code will run single threaded, without use of OPENMP.
