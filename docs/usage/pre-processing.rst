File Inputs
============

The HIPSTER code requires two main inputs; (a) files containing galaxy, simulation particle or random particle positions and (b) a file specifying the desired :math:`k`-space binning strategy. The file types of these are described below. We additionally provide a number of utility functions to assist with the creation of these files. These are located in the ``python/`` directory.

.. _particle_file_inputs

Data and Random Particle Files
---------------------------------

The main inputs to the HIPSTER code are files containing the locations and weights of galaxies or simulation particles (i.e. 'data') and random particles (i.e. 'randoms', only for non-periodic surveys). The random particles are of the same form as those used in real-space correlation function analyses, and we expect their distribution to match that of unclustered galaxies in the survey. These files are usually provided by clustering teams or can be simply created in the case of a periodic box geometry. The file format is a list of particle positions in space-separated (x,y,z,weight) format, with the co-ordinates given in comoving :math:`h^{-1}\mathrm{Mpc}`` units. For periodic simulation box data, weights are assumed to be unity everywhere and do not need to be specified.

We provide a convenience function to convert galaxy/random files in comoving (RA,Dec,z,weight) co-ordinates to the required format (using a simple WCDM co-ordinate converter by Daniel Eisenstein)::

    python python/convert_to_xyz.py {INFILE} {OUTFILE} {OMEGA_M} {OMEGA_K} {W_DARK_ENERGY}

where {INFILE} and {OUTFILE} are the filenames for the (RA,Dec,redshift,weight) and the remaining parameters specify the (present-day) cosmology. If these are not specified a cosmology of :math:`\{\Omega_m = 0.31,\Omega_k = 0,w_\Lambda = -1\}` is assumed by default.

For pair-count analysis, we usually require random particle files larger than the data; we suggest using :math:`N_\mathrm{rand}\sim 50N_\mathrm{gal}` for DR counts and :math:`N_\mathrm{rand}\sim 10N_\mathrm{gal}` for RR counts to minimize random noise. For this reason we provide a further convenience function to draw a random subset of a given particle file (in (x,y,z,weight) co-ordinates). This is run as follows (where {N_PARTICLES} specifies the size of the output file)::

    python python/take_subset_of_particles/py {INPUT_FILE} {OUTPUT_FILE} {N_PARTICLES}

For periodic boxes, we do not need to generate random particle files, but the above script can be useful for subsampling the data, for faster runtimes.

.. _binning_function_input

Binning Functions
------------------

In addition to the sets of data/random positions, we require a file to specify the desired :math:`k`-space binning of the output power spectra. Two Python routines are provided to produce the relevant files in linear or logarithmic (base :math:`e`) binning and are run as::

        python python/compute_binning_file_linear.py {N_LOG_BINS} {MIN_K} {MAX_K} {OUTPUT_FILE}
        python python/compute_binning_file_log.py {N_LINEAR_BINS} {MIN_K} {MAX_K} {OUTPUT_FILE}

with the output file saved to the specified destination.

The binning file can also be manually specified as an ASCII file with each line specifying the upper and lower coordinates of each :math:`k`-bin (in comoving :math:`h\,\mathrm{Mpc}^{-1}`` units). Note that the bins are required to be contiguous (i.e. the upper limit of the :math:`n`-th bin should equal the lower limit of the :math:`(n+1)`-th bin.
