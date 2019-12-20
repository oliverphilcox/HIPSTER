# HIPSTER

## HIgh-k Power Spectrum EstimatoR

Code to compute configuration-space power spectra for N-body simulations and galaxy surveys of arbitrary shape, based on the work of Philcox & Eisenstein (2019, accepted by MNRAS). This computes the Legendre multipoles of the power spectrum, $P_\ell(k)$ by computing weighted pair counts over the survey, truncated at some maximum radius $R_0$. This fully accounts for window function effects, does not include shot-noise, and is optimized for small-scale power spectrum computation in real- or redshift-space.

Full documentation of HIPSTER is available on [ReadTheDocs](https://HIPSTER.readthedocs.io). To compute a galaxy power spectrum up to $\ell=L$ with weighted pair counts from galaxies (``galaxies.dat``) and randoms (``randoms.dat``) truncated at radius $R_0$ with $k$-binning file ``binning.csv``, the basic usage is simply:

    ./hipster_wrapper.sh --dat galaxies.dat --ran_DR randoms.dat --ran_RR randoms.dat -l_max L -R0 R0 -k_bin binning.csv

For any queries regarding the code please contact [Oliver Philcox](mailto:ohep2@alumni.cam.ac.uk).

New for version 2: Optimizations for Periodic N-body Simulations
