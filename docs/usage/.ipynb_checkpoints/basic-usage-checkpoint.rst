Computing Configuration-Space Spectra
=======================================

Overview
--------

For standard usage, we provide simple bash wrappers for HIPSTER, which compute the power spectrum or bispectrum given (a) a set of galaxy or simulation particle positions, and (b) a :math:`k`-space binning file (with formats discussed in :doc:`pre-processing`). For non-periodic data-sets (e.g. galaxy surveys), a set of random particle positions is also required for the power spectrum. (Though the periodic bispectrum also requires random particles, these are created internally.) These are basic wrappers around the C++ and Python scripts, which, for more advanced usage, can be run separately, as discussed in :doc:`advanced-usage`.

The basic structure of the wrappers is as follows.

**Power Spectrum**:
  1) *(Optional)* Subsample the data (and random) files, such that the analysis uses only a fraction of the full dataset. This can be useful to speed up slow pair counts.
  2) *(Non-Periodic Only)* Compute the geometry correction function :math:`\Phi^{-1}` and fit it to a smooth model.
  3) *(Non-Periodic Only)* Compute the weighted random-random (RR) pair counts.
  4) *(Non-Periodic Only)* Compute the weighted data-random (DR) pair counts.
  5) Compute the weighted data-data (DD) pair counts.
  6) Combine the pair counts and output the power spectrum estimates.

**Bispectrum**:
  1) *(Optional)* Subsample the data file, as above.
  2) Compute the weighted pair counts and construct the bispectrum estimate.

Note that, for the power spectrum, steps (2) and (3) depend only on the survey geometry and random particle files, thus, if multiple mocks are being analyzed, they need only be run once. If the relevant option is specified, the wrapper will look for pre-computed survey-correction functions and RR pair counts and only re-create them if they do not exist. For aperiodic surveys with large truncation radii or many random particles, these steps are slow, thus this provides a significant speed boost. For aperiodic surveys, step (6) is done via a Python script, whilst for periodic simulations and bispectra, it takes place in the main C++ code.

Using the HIPSTER Wrapper
--------------------------

The HIPSTER wrappers can be run simply via ``./hipster_wrapper.sh``, ``./hipster_wrapper_periodic.sh`` or ``./hipster_wrapper_bispectrum.sh``. The following arguments are required:

    - ``--dat``: Data file in (x,y,z,weight) co-ordinates.
    - ``--l_max``: Maximum Legendre multipole. For the power spectrum, HIPSTER currently only uses even multipoles, but all multipoles are used for the bispectrum.
    - ``--R0``: Pair count truncation radius (see note below).
    - ``--f_rand``: *(Bispectrum Only)* Ratio of random particles to data points (see note below).
    - ``--k_bin``: :math:`k`-space binning file (created in :doc:`pre-processing` or user-defined).
    - ``--ran_DR``: *(Non-Periodic Power Only)* Random file for DR pair counting.
    - ``--ran_RR``: *(Non-Periodic Power Only)* Random file for RR pair counting (and survey correction function estimation).

A number of additional arguments are possible:

    - ``--string`` (Optional): Identification string for output file names. Default: hipster.
    - ``--nthreads`` (Optional): Number of CPU threads on which to run. Default: 10.
    - ``--subsample`` (Optional):  Factor by which to sub-sample the data. Default: 1 (no subsampling)
    - ``--load_RR``: *(Non-Periodic Power Only)* If set, load previously computed RR pair counts and survey correction functions for a large speed boost. If these are not found, they will be recomputed.
    - ``-h``: Display the command line options.

Note that, for the power spectrum of non-periodic surveys, two different random catalogs can be provided; one to compute the DR counts and one to compute the RR pair counts. It is usually preferable to use a larger random catalog for the DR pair counts to reduce noise. We recommend around 50x randoms for DR counts and :math:`\sim 10`x randoms for the RR counts.

As an example, consider computing the isotropic (:math:`\ell=0`) power spectrum cut at :math:`R_0=50h^{-1}\mathrm{Mpc}` from a single set of galaxies (``galaxies.dat``) and randoms (``randoms.dat``), given some :math:`k`-binning file  ``binning.csv``::

    ./hipster_wrapper.sh --dat galaxies.dat --ran_DR randoms.dat --ran_RR randoms.dat -l_max 0 -R0 50 -k_bin binning.csv --nthread 4

We've specified that the code should run on 4 cores here.

Similarly, for a simulation with periodic boundary conditions with data-file ``data.dat`` containing particle positions::

    ./hipster_wrapper_periodic.sh --dat data.dat -l_max 0 -R0 50 -k_bin binning.csv --subsample 2

Here, we've set the subsampling to 2, meaning that we'll use (a randomly selected) half of the available data, to get faster computation.

The output of the wrapper is saved in the HIPSTER directory as ``output/{STRING}_power_spectrum_n{K_BINS}_l{MAX_L}_R0{R0}.txt`` where {STRING} is the identification string described above, {MAX_L} is the maximum Legendre multipole, {K_BINS} is the number of :math:`k` bins in the binning file and {R0} is the truncation radius. The output file contains power spectrum estimates for each :math:`k`-bin on a separate line, with the column indicating the (even) Legendre multipole.

For the bispectrum, we can run the following to get the spectrum up to :math:`\ell=4` (assuming that the input file ``data.dat`` has periodic boundary conditions)::

    ./hipster_wrapper_bispectrum.sh --dat data.dat -l_max 4 -R0 50 -k_bin binning.csv --subsample 2 --f_rand 3

We've chosen to use three times as many randoms as data points, which is usually sufficient. This outputs a bispectrum in the HIPSTER directory as ``output/{STRING}_bispectrum_n{K_BINS}_l{MAX_L}_R0{R0}.txt`` where {STRING} is the identification string described above, {MAX_L} is the maximum Legendre multipole, {K_BINS} is the number of :math:`k` bins in the binning file and {R0} is the truncation radius. The output file contains bispectrum estimates for each Legendre multipole in a separate column with the row indicating the :math:`k_1,k_2` bins (in :math:`h\,\mathrm{Mpc}^{-1}` units), using the indexing :math:`n_\mathrm{collapsed} = in_\mathrm{bin}+j` for the :math:`i`-th :math:`k_1` and :math:`j`-th :math:`k_2` bin, with :math:`n_\mathrm{bin}` total bins.

.. _periodicity-note:

Note on Periodicity
~~~~~~~~~~~~~~~~~~~~

For the power spectrum, HIPSTER can be run in either *periodic* or *aperiodic* mode. In the former, we assume the simulation takes the form of a cubic box and measure the angle :math:`\mu` from the Z-axis, as appropriate for most simulations. In the aperiodic case, we measure :math:`\mu` relative to the local line of sight, as appropriate for (non-uniform and non-cubic) surveys. The periodic wrapper runs many times faster than the non-periodic one; this is as a result of many simplifications in the underlying equations. Currently the bispectrum is only supported in periodic mode.

To specify periodicity when using the C++ code alone (without the bash wrapper), we can pass the ``-perbox`` argument to the C++ code, which must be compiled with the ``-DPERIODIC`` flag (that can be manually added to the Makefile). The C++ code will crash if this is not specified.

.. _truncation-radius-note

Note on choice of Truncation Radius and Bin Widths
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A key hyperparameter of the code is the power spectrum estimation is the *truncation radius* :math:`R_0`. This is the maximum radius up to which particle counts are computed and sets the computation time of the algorithm (which scales as :math:`R_0^3`). As discussed in the introductory paper, the effect of :math:`R_0` is to convolve the true power spectrum with a window function of characteristic scale :math:`3/R_0`, giving a small bias which is important at low-:math:`k`, but negligible on small-scales. Considering moments up to :math:`\ell=4`, we find :math:`R_0=50h^{-1}\mathrm{Mpc}` to be sufficient for measuring :math:`k\gtrsim 0.5h\,\mathrm{Mpc}^{-1}` and :math:`R_0=100h^{-1}\mathrm{Mpc}` to be sufficient for :math:`k\gtrsim 0.25h\,\mathrm{Mpc}^{-1}`. For fixed truncation error, :math:`R_0` scales inversely with the minimum :math:`k`-bin of interest.

The choice of :math:`R_0` also sets the :math:`k`-binning scale via :math:`\Delta k\gtrsim 3/R_0` (assuming linear binning). Using narrow :math:`k`-bins will not give additional information, but lead to the :math:`k`-bins becoming more correlated.

.. _bispectrum-randoms-note

Note on Random Particles in the Bispectrum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whilst it is possible to compute the periodic-box bispectrum without any use of random particles (by performing all random particle integrals analytically), this turns out to be very computationally intensive for the bispectrum. As detailed in the second HIPSTER paper, there is one particular term (labelled :math:`\widetilde{DDR}^{II}`) that is difficult to compute, thus we elect to compute it via pair counts with a random catalog which is created internally via HIPSTER. The HIPSTER parameter :math:`f_\mathrm{rand}` is the ratio of random particles to galaxies (after subsampling, if applied), and controls this effect. Generally a ratio of order a few gives little sampling noise, but this can be easily experimented with. The runtime of the code scales in proportion to :math:`(1+f_\mathrm{rand})`.

.. _particle-weights-note

Note on Particle Weights
~~~~~~~~~~~~~~~~~~~~~~~~~

In the input data file to HIPSTER, particle weights can optionally be specified. For aperiodic computations, these can take any form, for example FKP weights for the power spectrum. For periodic surveys, weights are also supported, which are conventionally used to compute polyspectra from multiple-tracer simulations, for example with dark matter and neutrino particles. In this case, each tracer particle carries a weight proportional to its mass. These weights have the additional constraint that, if the particles are unclustered, the weights must not be varying across the survey, though can vary between sets of particles. In other words, HIPSTER can use a collection of different particle sets at once, each set of which has a (different) uniform weight. If unspecified, weights are set to unity.
