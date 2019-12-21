# HIPSTER

## HIgh-k Power Spectrum EstimatoR

Code to compute small-scale power spectra for N-body simulations and galaxy surveys of arbitrary shape, based on the configuration space estimators of Philcox & Eisenstein (2019, accepted by MNRAS, [arXiv](https://arxiv.org/abs/1912.01010)) and Philcox 2020 (in prep.). This computes the Legendre multipoles of the power spectrum, <img src="/tex/a939b8abbd34a6a7097130a860c9ebc2.svg?invert_in_darkmode&sanitize=true" align=middle width=38.738704949999985pt height=24.65753399999998pt/> by computing weighted pair counts over the simulation box or survey, truncated at some maximum radius <img src="/tex/12d208b4b5de7762e00b1b8fb5c66641.svg?invert_in_darkmode&sanitize=true" align=middle width=19.034022149999988pt height=22.465723500000017pt/>. This fully accounts for window function effects, does not include shot-noise, and is optimized for small-scale power spectrum computation in real- or redshift-space.

The code can be run either in 'aperiodic' or 'periodic' mode, for galaxy surveys or N-body simulations respectively. The 'periodic' mode contains various optimizations relating to the periodic geometry, as detailed in the second paper.

Full documentation of HIPSTER is available on [ReadTheDocs](https://HIPSTER.readthedocs.io).

### Basic Usage

To compute a power spectrum from particles in a *periodic* simulation box (``data.dat``) up to <img src="/tex/720b52da688c892f252bc47ce206b36d.svg?invert_in_darkmode&sanitize=true" align=middle width=39.95424014999999pt height=22.831056599999986pt/> with pair-counts truncated at radius <img src="/tex/12d208b4b5de7762e00b1b8fb5c66641.svg?invert_in_darkmode&sanitize=true" align=middle width=19.034022149999988pt height=22.465723500000017pt/> using <img src="/tex/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode&sanitize=true" align=middle width=9.075367949999992pt height=22.831056599999986pt/>-space binning file ``binning.csv`` and 4 CPU cores, run:

    ./hipster_wrapper_periodic.sh --dat data.dat --l_max L -R0 R0 -k_bin binning.csv --nthreads 4

To compute a power spectrum from galaxies in a *non-periodic* survey (``data.dat``), defined by a set of randoms (``randoms.dat``), up to <img src="/tex/720b52da688c892f252bc47ce206b36d.svg?invert_in_darkmode&sanitize=true" align=middle width=39.95424014999999pt height=22.831056599999986pt/>, truncating pair-counts at <img src="/tex/12d208b4b5de7762e00b1b8fb5c66641.svg?invert_in_darkmode&sanitize=true" align=middle width=19.034022149999988pt height=22.465723500000017pt/> and using <img src="/tex/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode&sanitize=true" align=middle width=9.075367949999992pt height=22.831056599999986pt/>-space binning file ``binning.csv``, with 4 CPU-cores, run:

    ./hipster_wrapper.sh --dat galaxies.dat --ran_DR randoms.dat --ran_RR randoms.dat -l_max L -R0 R0 -k_bin binning.csv --nthreads 4

For any queries regarding the code please contact [Oliver Philcox](mailto:ohep2@alumni.cam.ac.uk).

**New for version 2**: Optimizations for periodic N-body simulations
