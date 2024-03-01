## Function to compute the window function multipoles for a given survey geometry (specified an input random particle file).
## This is based on a simple wrapper for Corrfunc, to compute RR counts which are compared to an idealized model
## NB: This assumes an aperiodic survey.
## We store the multipoles Phi_ell = RR_ell / RR_ideal which can be used to forward model the spectra, e.g.
## P_0(k) = FT[xi_0(r)W(r;R_0)Phi_0(r) + 1/5 xi_2(r)W(r;R_0)Phi_2(r) + ...](k)

import sys
import numpy as np
from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks

# PARAMETERS
if len(sys.argv)!=8:
    print("Usage: python compute_window_function_lya.py {SKEWER_PARTICLE_FILE} {OUTFILE} {R_MAX} {L_MAX} {N_R_BINS} {N_MU_BINS} {NTHREADS}")
    sys.exit()
    
fname = str(sys.argv[1])
outfile = str(sys.argv[2])
r_max = float(sys.argv[3])
l_max = int(sys.argv[4])
nrbins = int(sys.argv[5])
nmubins = int(sys.argv[6]) 
nthreads = int(sys.argv[7])

## First read in weights and positions:
dtype = np.double

print("Counting lines in file")
total_lines=0
for n, line in enumerate(open(fname, 'r')):
    total_lines+=1

X,Y,Z,W,TID=[np.zeros(total_lines) for _ in range(5)]

print("Reading in data");
for n, line in enumerate(open(fname, 'r')):
    if n%1000000==0:
        print("Reading line %d of %d" %(n,total_lines))
    split_line=np.array(line.split(" "), dtype=float)
    X[n]=split_line[0];
    Y[n]=split_line[1];
    Z[n]=split_line[2];
    W[n]=split_line[3];
    TID[n]=split_line[4];
    
N = len(X) # number of particles

# Estimate volume
Vol = (np.max(Z)-np.min(Z))*(np.max(Y)-np.min(Y))*(np.max(X)-np.min(X))

print("Number of skewer points %.1e with density %.2e"%(N,N/Vol))

print('Computing pair counts up to a maximum radius of %.2f with %d radial and %d angular bins'%(r_max,nrbins,nmubins))

# Now compute RR counts (with full LoS angles)
#from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks
from Corrfunc.theory.DDsmu import DDsmu

binfile = np.linspace(0.,r_max,nrbins+1) # define binning 
binfile[0] += 1e-6 # to eliminate self-counts
r_hi = binfile[1:]
r_lo = binfile[:-1]
r_cen = 0.5*(r_lo+r_hi)
mu_max = 1. # to avoid self-skewers

assert nmubins>=8*l_max, "Need sufficient angular resolution to measure higher-multipoles"
assert l_max%2==0, "l_max should be even"

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

# Compute Corrfunc counts
RR = DDsmu_mocks(1,2,nthreads,mu_max,nmubins,binfile,Ra,Dec,com_dist,weights1=W,weight_type='pair_product',
            verbose=True,is_comoving_dist=True,output_savg=True)

# Assemble counts
RR_counts = (RR[:]['npairs']*RR[:]['weightavg']).reshape((nrbins,nmubins))
r_av = 0.5*(binfile[1:]+binfile[:-1])

# Subtract off self counts
all_tids = np.unique(TID)
skewerRR_counts = 0.
for tid in all_tids:
    if tid%50==0: print("Removing self-counts from skewer %d"%tid)
    filt = (TID==tid)
    if np.sum(filt)==1: continue
    skewerRR = DDsmu_mocks(1,2,nthreads,mu_max,nmubins,binfile,Ra[filt],Dec[filt],com_dist[filt],weights1=W[filt],weight_type='pair_product',is_comoving_dist=True,verbose=False)

    skewerRR_counts += (skewerRR[:]['npairs']*skewerRR[:]['weightavg']).reshape((nrbins,nmubins))
    
RR_counts -= skewerRR_counts
    
# Compute ideal RR counts
RR_th = 4.*np.pi/3.*(binfile[1:]**3-binfile[:-1]**3)*1./nmubins

# Define Legendre polynomials
from scipy.special import legendre
mu_bins = np.linspace(0.,mu_max-1./nmubins,nmubins)+0.5/nmubins
legs = [legendre(l)(mu_bins) for l in np.arange(0,l_max+1,2)]
    
# Compute window function multipoles
Phis = np.asarray([(2.*l+1.)*(legs[i]*RR_counts).mean(axis=1)/RR_th for i,l in enumerate(np.arange(0,l_max+1,2))])
# Add radial axis
rPhis = np.vstack([r_av,Phis])
names = 'r\t'
for l in np.arange(0,l_max+1,2):
    names += 'Phi_%d\t'%l

np.savetxt(outfile,rPhis,delimiter="\t",header=names)
print("Saved window function multipoles to %s"%outfile)
