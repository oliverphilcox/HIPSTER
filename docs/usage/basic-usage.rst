Computing Configuration-Space Power Spectra
============================================

Overview
--------

For standard usage, we provide a simple bash wrapper for HIPSTER, which computes the power spectrum given a set of galaxy positions, a set of random particle positions and a :math:`k`-space binning file (with formats discussed in :doc:`pre-processing`). This is a basic wrapper around the C++ and Python scripts, which, for more advanced usage, can be run separately, as discussed in :doc:`advanced-usage`.

The basic structure of the wrapper is as follows:

  1) Compute the geometry correction function :math:`\Phi^{-1}` and fit it to a smooth model.
  2) Compute the weighted random-random (RR) pair counts.
  3) Compute the weighted data-random (DR) pair counts.
  4) Compute the weighted data-data (DD) pair counts.
  5) Combine the pair counts and output the power spectrum estimates.

Note that steps (1) and (2) depend only on the survey geometry and random particle files, thus, if multiple mocks are being analyzed, they need only be run once. If the relevant option is specified, the wrapper will look for pre-computed survey-correction functions and RR pair counts and only re-create them if they do not exist. For aperiodic surveys with large truncation radii or many random particles, these steps are slow, thus this provides a significant speed boost.

.. _periodicity-note:

Note on Periodicity
~~~~~~~~~~~~~~~~~~~~

HIPSTER can be run in either *periodic* or *aperiodic* mode. In the former, we assume the simulation takes the form of a cubic box and measure the angle :math:`\mu` from the Z-axis, as appropriate for many simulations. In the aperiodic case, we measure :math:`\mu` relative to the local line of sight, as appropriate for (non-uniform and non-cubic) surveys.

To specify periodicity, we can pass the ``--periodic`` argument to the bash script. In addition, the C++ code must be compiled with the ``-DPERIODIC`` flag which should be manually added to the Makefile. The C++ code will crash if this is not specified.

.. _truncation-radius-note

Note on choice of Truncation Radius and Bin Widths
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A key hyperparameter of the code is the power spectrum estimation is the *truncation radius* :math:`R_0`. This is the maximum radius up to which particle counts are computed and sets the computation time of the algorithm (which scales as :math:`R_0^3`). As discussed in the introductory paper, the effect of :math:`R_0` is to convolve the true power spectrum with a window function of characteristic scale :math:`3/R_0`, giving a small bias which is important at low-:math:`k`, but negligible on small-scales. Considering moments up to :math:`\ell=4`, we find :math:`R_0=50h^{-1}\mathrm{Mpc}` to be sufficient for measuring :math:`k\gtrsim 0.5h\,\mathrm{Mpc}^{-1}` and :math:`R_0=100h^{-1}\mathrm{Mpc}` to be sufficient for :math:`k\gtrsim 0.25h\,\mathrm{Mpc}^{-1}`. For fixed truncation error, :math:`R_0` scales inversely with the minimum :math:`k`-bin of interest.

The choice of :math:`R_0` also sets the :math:`k`-binning scale via :math:`Delta k\gtrsim 3/R_0` (assuming linear binning). Using narrow :math:`k`-bins will not give additional information, but lead to the :math:`k`-bins becoming more correlated.

Using the HIPSTER Wrapper
--------------------------

The HIPSTER wrapper can be run simply via ``./hipster_wrapper.sh``. The following arguments are required:

    - ``--dat``: Data file in (x,y,z,weight) co-ordinates.
    - ``--ran_DR``: Random file for DR pair counting.
    - ``--ran_RR``: Random file for RR pair counting (and survey correction function estimation).
    - ``--l_max``: Maximum Legendre multipole.
    - ``--R0``: Pair count truncation radius (see note above).
    - ``--k_bin``: :math:`k`-space binning file (created in :doc:`pre-processing` or user-defined).

A number of additional arguments are possible:

    - ``--string``: (Optional): Identification string for output file names.
    - ``--nthreads``: (Optional): Number of CPU threads on which to run. Default: 10.
    - ``--periodic``: If set, assume periodic boundary conditions. C++ code must also be compiled with the ``-DPERIODIC`` flag, as discussed above.
    - ``--load_RR``: If set, load previously computed RR pair counts and survey correction functions for a large speed boost. If these are not found, they will be recomputed.
    - ``-h``: Display the command line options.

Note that two different random catalogs can be provided; one to compute the DR counts and one to compute the RR pair counts. It is usually preferable to use a larger random catalog for the DR pair counts to reduce noise. We recommend :math:`\sim 50`x randoms for DR counts and :math:`\sim 10`x randoms for the RR counts.

As an example, consider computing the isotropic (:math:`\ell=0`) power spectrum cut at :math:`R_0=50h^{-1}\mathrm{Mpc}` from a single set of galaxies (``galaxies.dat``) and randoms (``randoms.dat``), given some :math:`k`-binning file  ``binning.csv``::

    ./hipster_wrapper.sh --dat galaxies.dat --ran_DR randoms.dat --ran_RR randoms.dat -l_max 0 -R0 50 -k_bin binning.csv

The output of the wrapper is saved as ``output/{STRING}_power_spectrum_n{K_BINS}_l{MAX_L}.txt`` where {STRING} is the identification string described above, {MAX_L} is the maximum Legendre multipole and {K_BINS} is the number of :math:`k` bins in the binning file. The output file contains power spectrum estimates for each :math:`k`-bin on a separate line, with the column indicating the (even) Legendre multipole.
