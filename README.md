# HIPSTER

## HIgh-k Power Spectrum EstimatoR

Code to compute configuration-space power spectra for galaxy surveys of arbitrary shape, based on the work of Philcox & Eisenstein (2019, accepted by MNRAS). This computes the Legendre multipoles of the power spectrum, <img src="/tex/a939b8abbd34a6a7097130a860c9ebc2.svg?invert_in_darkmode&sanitize=true" align=middle width=38.738704949999985pt height=24.65753399999998pt/> by computing weighted pair counts over the survey, truncated at some maximum radius <img src="/tex/12d208b4b5de7762e00b1b8fb5c66641.svg?invert_in_darkmode&sanitize=true" align=middle width=19.034022149999988pt height=22.465723500000017pt/>. This fully accounts for window function effects, does not include shot-noise, and is optimized for small-scale power spectrum computation in real- or redshift-space.

Full documentation of HIPSTER is available on [ReadTheDocs](HIPSTER.readthedocs.io). To compute a galaxy power spectrum up to :math:`\ell=L` with weighted pair counts from galaxies (``galaxies.dat``) and randoms (``randoms.dat``) truncated at radius ``R0`` with :math:`k`-binning file ``binning.csv``, the basic usage is simply:

    ./hipster_wrapper.sh --dat galaxies.dat --ran_DR randoms.dat --ran_RR randoms.dat -l_max L -R0 R0 -k_bin binning.csv
