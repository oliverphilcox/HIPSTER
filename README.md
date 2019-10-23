# ConfigPowerSpectra

Code to compute configuration-space power spectra for arbitrary survey geometries, based on the work of Philcox & Eisenstein (2019). This computes the Legendre multipoles of the power spectrum, $P_\ell(k)$ by computing weighted pair counts over the survey, truncated at some maximum radius $R_0$. This fully accounts for window function effects, does not include shot-noise, and is optimized for small-scale power spectrum computation in real- or redshift-space.

# Code Requirements

- C++ Compiler (g++ used by default)
- OpenMP Installation (optional, but required for parallelization)
- Python (optional, but useful for pre- and post-processing.)

# Computing Binning Functions

We provide two simply Python routines to produce $k$-space binning files. These are either in linear or logarithmic (base e) binning and are run as::

    python python/compute_binning_file_linear.py {N_LOG_BINS} {MIN_K} {MAX_K} {OUTPUT_FILE}
    python python/compute_binning_file_log.py {N_LINEAR_BINS} {MIN_K} {MAX_K} {OUTPUT_FILE}

with the output file saved to the specified destination. The binning file can also be manually specified as an ASCII file with each line specifying the upper and lower coordinates of each $k$-bin (in $h\,\mathrm{Mpc}^{-1}$ units). Note that the bins are required to be contiguous (i.e. the upper limit of the $n$th bin should equal the lower limit of the $(n+1)$th bin.

# Inputs

# Phi Files

# Note on choice of R0
Add time scaling and errors

# Code compilation

To compile the C++ code, simply run the following in the installation directory::

    bash clean
    make

This creates the ``./power`` executable in the working directory.

The Makefile may need to be modified depending on the specific computational system. In addition, by specifying the ``-DPERIODIC`` flag in the Makefile, the code will assume periodic boundary conditions (as appropriate for periodic box simulations), measuring the angle $\mu$ from the z-axis. Furthermore, by removing the ``-DOPENMP`` flag, the code will run single threaded, without use of OPENMP.


# Compute Weighted Pair Counts

Configuration space power spectra are computed by estimating RR, DR and DD pair counts (analogous to 2PCF computation) weighted by a $k$-bin-dependent kernel. To compute pair counts for two input fields we simply run the ``./power`` executable, specifying the input parameters on the command line. (Parameters can also be altered in the ``power_mod/ppower_parameters.h`` file, but this requires the code to be re-compiled after each modification.). As an example let's compute a set of data-random (DR) counts from input files ``galaxy_positions.txt`` and ``random_positions.txt``.

    ./power -in galaxy_positions.txt -in2 random_positions.txt -binfile binfile.csv -output output/ -out_string DR -max_l 2 -R0 100 -inv_phi_file inv_phi_coefficients.txt -nthread 10

This runs in a few minutes to hours (depending on the catalog size and computational resources available) and uses the following main arguments:

    -``-in``: First input ASCII file containing space separated (x,y,z,weight) positions of particles in comoving Mpc/h units.
    -``-in2``: Second input ASCII file, with format as above.
    -``-binfile``: $k$-space ASCII binning file, as described above.
    -``-output``: Directory in which to house output products. This will be created if not already in existence.
    -``-out_string``: String to include in output filename for identification (e.g. RR, DR or DD)
    -``-max_l``: Maximum Legendre multipole required (must be even). Maximum 10.
    -``-R0``: Truncation radius in Mpc/h units (default: $100\,h^{-1}\mathrm{Mpc}$). See note above.
    -``-inv_phi_file``: Location of survey geometry correction file, as produced above.
    -``-nthread``: Number of CPU threads to use for the computation.
    -``-periodic``: This flag must be set if we require the code to be run with *periodic* boundary conditions, measuring $\mu$ from the z-axis.

Note that a full list of command line options to the executable can be shown by running ``./power`` without any arguments. The code creates the output file ``{OUT_STRING}_power_counts_n{N_BINS}_m{MAX_L}_full.txt``, specifying the ``out_string`` parameter, the number of radial bins and the maximum Legendre multipole. Each line of the output file has the (weighted) pair count with the column specifying the Legendre multipole.

To compute the full power spectra, the data-data (DD), data-random (DR) and random-random (RR) pair counts must be computed. We do *not* have to use the same sized random catalogs for the DR and RR counts. It is usually preferable to use a larger random catalog for the DR pair counts to reduce noise. We recommend $\sim$50x randoms for DR counts and $\sim$10x randoms for the DD counts. Note that the RR counts are the most computationally intensive procedure, but they only need be computed for each survey once (i.e. when analyzing mock data, the RR pair counts are the same for each mock).

# Reconstruct Power Spectrum

Once the pair counts have been computed, it is straightforward to reconstruct the power spectrum. This can be done via a simple Python script;

    python python/reconstruct_power.py {DD_FILE} {DR_FILE} {RR_FILE} {GAL_FILE} {N_RAND_RR} {N_RAND_DR} {PERIODIC} {OUTFILE}

where {DD_FILE}, {DR_FILE} and {RR_FILE} give the locations of the DD, DR and RR weighted pair counts, {GAL_FILE} gives the input galaxy file (needed for normalization), {N_RAND_RR} and {N_RAND_DR} give the number of random particles used for RR and DR counts. {PERIODIC} is unity if the code is computed with periodic boundary conditions and zero else. The output power spectrum is given in ASCII format in the specified {OUTFILE}, with the power spectrum estimates for each $k$-bin on a separate line, with the column indicating the (even) Legendre multipole.
