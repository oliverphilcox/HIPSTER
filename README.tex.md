# HIPSTER

## HIgh-k Power Spectrum EstimatoR

Code to compute small-scale power spectra and isotropic bispectra for cosmological simulations and galaxy surveys of arbitrary shape, based on the configuration space estimators of Philcox & Eisenstein (2019, MNRAS, [arXiv](https://arxiv.org/abs/1912.01010)) and Philcox (2020, submitted, [arXiv](https://arxiv.org/pdf/2005.01739.pdf)). This computes the Legendre multipoles of the power spectrum, $P_\ell(k)$ or bispectrum, $B_\ell(k_1,k_2)$ by computing weighted pair counts over the simulation box or survey, truncated at some maximum radius $R_0$. This fully accounts for window function effects, does not include shot-noise or aliasing, and is optimized for small-scale spectrum computation in real- or redshift-space. Both the power spectrum and bispectrum estimators have efficiency $\mathcal{O}\left(N^2\right)$ for $N$ particles, and become faster at large $k$. We additionally include a *Lyman-alpha* power spectrum estimator, which performs an efficient real-space summation over different observational skewers. In this case, the estimator does *not* remove the window function effects (which can be large, since the geometry is highly non-spherical).

The code can be run either in 'aperiodic' or 'periodic' mode, for galaxy surveys or cosmological simulations respectively. The 'periodic' mode contains various optimizations relating to the periodic geometry, as detailed in the second paper. HIPSTER also supports weighted spectra, for example when tracer particles are weighted by their mass in a multi-species simulation. Generalization to anisotropic bispectra is straightforward (and requires no additional computing time) and can be added on request. 

Full documentation of HIPSTER is available on [ReadTheDocs](https://HIPSTER.readthedocs.io).

### Basic Usage

To compute a power spectrum from particles in a *periodic* simulation box (``data.dat``) up to $\ell=L$ with pair-counts truncated at radius $R_0$ using $k$-space binning file ``binning.csv`` and 4 CPU cores, run:

    ./hipster_wrapper_periodic.sh --dat data.dat --l_max L -R0 R0 -k_bin binning.csv --nthreads 4

To compute a power spectrum from galaxies in a *non-periodic* survey (``data.dat``), defined by a set of randoms (``randoms.dat``), up to $\ell=L$, truncating pair-counts at $R_0$ and using $k$-space binning file ``binning.csv``, with 4 CPU-cores, run:

    ./hipster_wrapper.sh --dat galaxies.dat --ran_DR randoms.dat --ran_RR randoms.dat -l_max L -R0 R0 -k_bin binning.csv --nthreads 4

To compute an isotropic bispectrum from particles in a *periodic* simulation box (``data.dat``) up to $\ell=L$ with pair-counts truncated at radius $R_0$ using $k$-space binning file ``binning.csv`` and 4 CPU cores, using 3 times as many random points as data points, run:

    ./hipster_wrapper_bispectrum.sh --dat data.dat --l_max L -R0 R0 -k_bin binning.csv --nthreads 4 --f_rand 3
        
To compute a 3D Lyman-alpha power spectrum from Lya pixels in a *non-periodic* survey (``lya-data.dat``, including the x, y, z positions, (weighted) fractional overdensity fields and skewer IDs), up to $\ell=L$, truncating pair-counts at $R_0$ and using $k$-space binning file ``binning.csv``, with 4 CPU-cores, run:

    ./hipster_wrapper_lya.sh --dat lya-data.dat --ran_RR lya-pixels.dat -l_max L -R0 R0 -k_bin binning.csv --nthreads 4

For Lyman-alpha, one must convolve the desired theory models with an appropriate window function. This can be computed with the ```python/compute_window_lya.py``` script.

For any queries regarding the code please contact [Oliver Philcox](mailto:ohep2@cantab.ac.uk).

**New for version 2**: Optimizations for periodic N-body simulations

**New for version 3**: A new pair-count estimator for the periodic bispectrum

**New for version 4**: A new pair-count estimator for Lyman-alpha power spectra
