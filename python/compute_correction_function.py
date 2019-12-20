## Function to compute the survey correction function for a given survey geometry (specified an input random particle file).
## This is based on a simple wrapper for Corrfunc, to compute RR counts which are compared to an idealized model
## NB: This requires an aperiodic survey.
## Results are stored as fits to the multipoles of 1/Phi = RR_true / RR_model, which is read by the C++ code

import sys
import numpy as np

# PARAMETERS
if len(sys.argv)!=4:
    if len(sys.argv)!=8:
        print("Usage: python RR_counts.py {RANDOM_PARTICLE_FILE} {OUTFILE} {PERIODIC} [{R_MAX} {N_R_BINS} {N_MU_BINS} {NTHREADS}]")
        sys.exit()
fname = str(sys.argv[1])
outfile = str(sys.argv[2])
periodic = int(sys.argv[3])

if not periodic:
    if len(sys.argv)!=8:
        print("Usage: python RR_counts.py {RANDOM_PARTICLE_FILE} {OUTFILE} {PERIODIC} [{R_MAX} {N_R_BINS} {N_MU_BINS} {NTHREADS}]")
        sys.exit()
    r_max = float(sys.argv[4])
    nrbins = int(sys.argv[5])
    nmu_bins = int(sys.argv[6])
    nthreads = int(sys.argv[7])
    mu_max = 1.;
else:
    print("We don't need to compute a correction function for periodic data! Exiting.")
    sys.exit()
    
## First read in weights and positions:
dtype = np.double

print("Counting lines in file")
total_lines=0
for n, line in enumerate(open(fname, 'r')):
    total_lines+=1

X,Y,Z,W=[np.zeros(total_lines) for _ in range(4)]

print("Reading in data");
for n, line in enumerate(open(fname, 'r')):
    if n%1000000==0:
        print("Reading line %d of %d" %(n,total_lines))
    split_line=np.array(line.split(" "), dtype=float)
    X[n]=split_line[0];
    Y[n]=split_line[1];
    Z[n]=split_line[2];
    W[n]=split_line[3];

N = len(X) # number of particles

print("Number of random particles %.1e"%N)

if not periodic:
    print("Using non-periodic input data. This requires computation of RR counts");
    print('Computing pair counts up to a maximum radius of %.2f'%r_max)

    binfile = np.linspace(0,r_max,nrbins+1)
    binfile[0]=1e-4 # to avoid zero errors
    r_lo = binfile[1:]
    r_hi = binfile[:-1]
    r_cen = 0.5*(r_lo+r_hi)

    # Compute RR counts for the non-periodic case (measuring mu from the radial direction)
    def coord_transform(x,y,z):
        # Convert the X,Y,Z coordinates into Ra,Dec,comoving_distance (for use in corrfunc)
        # Shamelessly stolen from astropy
        xsq = x ** 2.
        ysq = y ** 2.
        zsq = z ** 2.

        com_dist = (xsq + ysq + zsq) ** 0.5
        s = (xsq + ysq) ** 0.5

        if np.isscalar(x) and np.isscalar(y) and np.isscalar(z):
            Ra = math.atan2(y, x)*180./np.pi
            Dec = math.atan2(z, s)*180./np.pi
        else:
            Ra = np.arctan2(y, x)*180./np.pi+180.
            Dec = np.arctan2(z, s)*180./np.pi

        return com_dist, Ra, Dec

    # Convert coordinates to spherical coordinates
    com_dist,Ra,Dec = coord_transform(X,Y,Z);

    # Now compute RR counts
    from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks

    print('Computing unweighted RR pair counts')
    RR=DDsmu_mocks(1,2,nthreads,mu_max,nmu_bins,binfile,Ra,Dec,com_dist,weights1=W,weight_type='pair_product',
                   verbose=False,is_comoving_dist=True)
    # Weight by average particle weighting
    RR_counts=(RR[:]['npairs']*RR[:]['weightavg']).reshape((nrbins,nmu_bins))

    # Now compute ideal model for RR Counts
    print("Compute correction function model")
    mu_cen = np.arange(1/(2*nmu_bins),1.+1/(2*nmu_bins),1/nmu_bins)
    delta_mu = (mu_cen[-1]-mu_cen[-2])
    delta_mu_all = delta_mu*np.ones_like(mu_cen).reshape(1,-1)
    norm = np.sum(W**2.) # partial normalization - divided by np.sum(W_gal**2) in reconstruction script
    RR_model = 4.*np.pi*(r_hi**3.-r_lo**3.).reshape(-1,1)*delta_mu_all*norm/3.

    # Compute inverse Phi function and multipoles
    inv_phi = RR_counts/RR_model
    l_max = 4
    from scipy.special import legendre
    inv_Phi_multipoles = np.zeros([l_max//2+1,len(inv_phi)])
    for i in range(len(RR_counts)):
        for l_i,ell in enumerate(np.arange(0,l_max+2,2)):
            inv_Phi_multipoles[l_i,i]=(2.*ell+1.)*delta_mu*np.sum(legendre(ell)(mu_cen)*inv_phi[i,:])

    # Now fit to a smooth model
    def inv_phi_ell_model(r,*par):
        return par[0]+par[1]*r+par[2]*r**2.
    from scipy.optimize import curve_fit
    all_ell = np.arange(0,l_max+2,2)
    coeff = np.asarray([curve_fit(inv_phi_ell_model,r_cen[1:],inv_Phi_multipoles[ell//2][1:],
                                  p0=[0 for _ in range(3)])[0] for ell in all_ell])

    np.savetxt(outfile,coeff,delimiter="\t")
    print("Saved correction function to %s"%outfile)

else:
    print('Assuming periodic data, so no RR counts need to be computed')

    # This is much simpler for the periodic case
    # The partially-normalized inverse correction function is simply the number density of random particles
    x_range = max(X)-min(X)
    y_range = max(Y)-min(Y)
    z_range = max(Z)-min(Z)

    if np.abs(x_range/y_range-1.)>0.01 or np.abs(x_range/z_range-1.)>0.01:
        print('Is this data periodic? Different x, y and z axis dimensions found. Exiting')
        sys.exit();

    V = x_range*y_range*z_range
    n_rand = len(X)/V # number density of randomds
    l_max = 4
    coeff = np.zeros([3,3])
    coeff[0,0] = n_rand # only need coefficient of monopole without r dependence
    np.savetxt(outfile,coeff,delimiter="\t")
    print("Saved correction function to %s"%outfile)
