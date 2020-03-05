Advanced Usage
===============

Here we describe the individual code modules used to compute the power spectrum multipoles, :math:`P_\ell(k)` or bispectrum multipoles :math:`B_\ell(k_1,k_2)`. For basic usage, the HIPSTER wrappers (:doc:`basic-usage`) can be used, which simply run all the relevant packages in sequence connecting the relevant inputs and outputs. If the user requires more flexibility (e.g. to use multiple sets of random particle files or to run particular sections of the code on an HPC cluster) the modules can be run separately, as described below.

.. _survey-correction-function:

Survey Correction Function
---------------------------

***Only Relevant for Non-Periodic Power Spectra***

An important ingredient in the weighted power spectrum pair counts is the 'survey correction function' :math:`\Phi`, defined as the ratio of ideal and true (unweighted) RR pair counts. For periodic data, this is simply unity, and does not need to be computed. In this package, we compute :math:`\Phi` using the `Corrfunc <https://Corrfunc.readthedocs.io>`_ pair counting routines (for aperiodic data) and store the results as quadratic fits to the first few multipoles of :math:`\Phi^{-1}`. Note that a normalization is also performed for later convenience.

This can be computed using the ``compute_correction_function.py`` script::

    python python/compute_correction_function.py {RANDOM_PARTICLE_FILE} {OUTFILE} {R_MAX} {N_R_BINS} {N_MU_BINS} {NTHREADS}

with inputs:

    - {RANDOM_PARTICLE_FILE}: File containing random particles in the survey geometry. Since the correction function is being fit to a smooth function, we can use a relatively small random catalog here.
    - {OUTFILE}: Location of output ASCII file. This is automatically read in by the C++ code
    - {PERIODIC}: Whether to assume periodic boundary conditions
    - {R_MAX}: Radius (in :math:`h^{-1}\mathrm{Mpc}`) up to which to count pairs and fit the correction function. This should be at least as large as the truncation radius (:math:`R_0`) used for the power spectrum computation. Note that the computation time scales as :math:`R_\mathrm{max}^3`.
    - {N_R_BINS}: Number of radial bins in the pair counting. We recommend using :math:`\sim 1\,h^{-1}\mathrm{Mpc}` radial binning here.
    - {N_MU_BINS}: Number of angular bins used in the pair counting. We recommend 100 :math:`\mu` bins.
    - {N_THREADS}: Number of CPU threads on which to perform pair counts.

This may take some time to compute, but only needs to be computed once for a given survey. Whilst the user can specify the number of radial and angular bins used by the pair counts, this does not have a significant affect on the output function, assuming a moderately fine binning is used.

.. _main-c-code:

Computing Weighted Pair Counts
-------------------------------

The main functionality of HIPSTER is to compute weighted pair counts across a survey, and combine these into a power spectrum or bispectrum. This is done via the ``grid_power.cpp`` C++ code, which must be compiled before use.

Compilation
~~~~~~~~~~~~

To compile the C++ code, simply run the following in the installation directory::

    bash clean
    make [Periodic=-DPERIODIC] [Bispectrum=-DBISPECTRUM]

This creates the ``./power`` executable in the working directory. If the ``Periodic=-DPERIODIC`` statement is included, the code is compiled assuming periodic boundary conditions (measuring the angle :math:`\mu` from the z-axis), as appropriate for N-body simulations. If the ``Bispectrum=-DBISPECTRUM`` statement is added, the code will compute the bispectrum rather than power spectrum.

The Makefile may need to be modified depending on the specific computational system. In particular, Mac users may need to remove the OpenMP references (``-DOPENMP`` and ``-lgomp``) to ensure compilation occurs. Note that this will force the code to run single threaded.

Running the Code
~~~~~~~~~~~~~~~~~

Configuration space power spectra are computed by estimating DD (and, for non-periodic surveys RR and DR) pair counts (analogous to 2PCF computation) weighted by a :math:`k`-bin-dependent kernel. For the bispectra, we estimate DDD and DDR counts, via spherical harmonic decomposition and pair counting, with the remaining counts performed semi-analytically.

To compute power spectrum pair counts for two input fields we simply run the ``./power`` executable, specifying the input parameters on the command line. (Parameters can also be altered in the ``power_mod/power_parameters.h`` file, but this requires the code to be re-compiled after each modification.). As an example let's compute a set of data-random (DR) counts from input files ``galaxy_positions.txt`` and ``random_positions.txt``::

    ./power -in galaxy_positions.txt -in2 random_positions.txt -binfile binfile.csv -output output/ -out_string DR -max_l 2 -R0 100 -inv_phi_file inv_phi_coefficients.txt -nthread 10

This runs in a few seconds to hours (depending on the catalog size and computational resources available).

Analogously, for the bispectrum we can run the ``./power`` executable, which, in this case, will compute all required pair counts and save the bispectrum::

    ./power -in galaxy_positions.txt -binfile binfile.csv -output output/ -out_string DR -max_l 2 -R0 100 -inv_phi_file inv_phi_coefficients.txt -nthread 10 -f_rand 3

The code uses the following main arguments:

    - ``-in``: First input ASCII file containing space separated (x,y,z,weight) positions of particles in comoving :math:`h^{-1}\mathrm{Mpc}` units. Note that the weight column is optional, and will be set to unity if not included (see :ref:`particle-weights-note`).
    - ``-in2``: (*Aperiodic power spectrum only*) Second input ASCII file, with format as above.
    - ``-binfile``: :math:`k`-space ASCII binning file, as described above.
    - ``-output``: Directory in which to house output products. This will be created if not already in existence.
    - ``-out_string``: String to include in output filename for identification (e.g. RR, DR, DD or a simulation name)
    - ``-max_l``: Maximum Legendre multipole required (must be even for the power spectrum). Currently, only multipoles up to the tetrahexacontapole (:math:`\ell = 6`) have been implemented for the power spectrum or :math:`\ell = 10` for the bispectrum, but more can be added if required. For power spectra (but not bispectra), all multipoles are even; if there is a need for odd multipoles, please contact the author and these can be easily added in.
    - ``-R0``: Truncation radius in Mpc/h units (default: :math:`100\,h^{-1}\mathrm{Mpc}`). See :ref:`truncation-radius-note`.
    - ``-inv_phi_file``: (*Non-Periodic Power spectrum only*) Location of survey geometry correction file, as produced above.
    - ``-nthread``: Number of CPU threads to use for the computation.
    - ``-perbox``: This flag must be set if we require the code to be run with *periodic* boundary conditions, measuring :math:`\mu` from the z-axis. The code must also be compiled with the -DPERIODIC flag.
    - ``f_rand``: (*Bispectrum only*) Ratio of random particles to galaxies used for DDR counts. See :ref:`bispectrum-randoms-note`.

Note that a full list of command line options to the executable can be shown by running ``./power`` without any arguments.

For the power spectrum, the code creates the output file ``{OUT_STRING}_power_counts_n{N_BINS}_l{MAX_L}_R0{R0}.txt``, specifying the ``out_string`` parameter, the number of radial bins and the maximum Legendre multipole. Each line of the output file has the (weighted) pair count with the column specifying the Legendre multipole. If the code has been run in periodic mode, it additionally outputs ``{OUT_STRING}_analyt_RR_power_counts_n{N_BINS}_l{MAX_L}_R0{R0}.txt`` containing the RR counts (computed from a 1-dimensional Hankel transform) and ``{OUT_STRING}_power_spectrum_n{N_BINS}_l{MAX_L}_R0{R0}.txt`` containing the full power spectrum estimate. This is the main output of the code.

To compute the full power spectra for non-periodic surveys, the data-data (DD), data-random (DR) and random-random (RR) pair counts must be computed. (For periodic surveys, we require only the data-data counts). We do *not* have to use the same sized random catalogs for the DR and RR counts. It is usually preferable to use a larger random catalog for the DR pair counts to reduce noise. We recommend :math:`\sim 50\times` randoms for DR counts and :math:`\sim 10\times` for the RR counts. Note that the RR counts are the most computationally intensive procedure, but they only need be computed for each survey once (i.e. when analyzing mock data, the RR pair counts are the same for each mock).

For the bispectrum, the code instead outputs the files ``{OUT_STRING}_{TYPE}_n{N_BINS}_l{MAX_L}_{R0}R0.txt`` where {TYPE} is ``bispectrum``, ``DDD_counts``, ``DDR_I_counts``, ``DDR_II_counts`` and ``analyt_RRR_counts``, giving the full bispectrum and various components. Each line of the output file has the (weighted) pair count in a combination of :math:`k_1,k_2` bins with the column specifying the Legendre multipole. The :math:`i`-th :math:`k_1` and :math:`j`-th :math:`k_2` bin is indexed as :math:`in_\mathrm{bins}+j`.

.. _power-spectrum-reconstruction

Reconstructing the Power Spectrum
----------------------------------

***Only Relevant for Non-Periodic Surveys***

Once the pair counts have been computed, it is straightforward to reconstruct the power spectrum. This can be done via a simple Python script::

    python python/reconstruct_power.py {DD_FILE} {DR_FILE} {RR_FILE} {GAL_FILE} {N_RAND_RR} {N_RAND_DR} {OUTFILE}

where {DD_FILE}, {DR_FILE} and {RR_FILE} give the locations of the DD, DR and RR weighted pair counts, {GAL_FILE} gives the input galaxy file (needed for normalization), {N_RAND_RR} and {N_RAND_DR} give the number of random particles used for RR and DR counts. {PERIODIC} is unity if the code is computed with periodic boundary conditions and zero else. The output power spectrum is given in ASCII format in the specified {OUTFILE}, with the power spectrum estimates for each :math:`k`-bin on a separate line, with the column indicating the (even) Legendre multipole.

For periodic simulations, the full power spectrum is created inside the C++ code, as described above.
