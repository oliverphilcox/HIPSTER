## Function to compute the survey correction function for a given survey geometry (specified an input random particle file).
## This is based on a simple wrapper for Corrfunc, to compute RR counts which are compared to an idealized model
## NB: This requires an aperiodic survey.
## Results are stored as fits to the multipoles of 1/Phi = RR_true / RR_model, which is read by the C++ code

import sys
import numpy as np
from Corrfunc.mocks import DDrppi_mocks

# PARAMETERS
if len(sys.argv)!=8:
    print("Usage: python compute_correction_function_lya.py {RANDOM_PARTICLE_FILE} {OUTFILE} {R_MAX} {N_R_BINS} {PI_MAX} {DELTA_PI} {NTHREADS}")
    sys.exit()
    
fname = str(sys.argv[1])
outfile = str(sys.argv[2])
r_max = float(sys.argv[3])
nrbins = int(sys.argv[4])
pimax = int(sys.argv[5]) # this should be an integer (for reasons known only to Corrfunc) 
Delta_pi = int(sys.argv[6]) # this should be an integer (for reasons known only to Corrfunc) 
nthreads = int(sys.argv[7])

## Check divisibility of Delta_pi and adjust as appropriate
if pimax%Delta_pi!=0: 
    pimax = (pimax//Delta_pi+1)*Delta_pi
print("Resetting pi_max to %d"%pimax)

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

print('Computing pair counts up to a maximum radius of %.2f'%r_max)

binfile = np.linspace(0.01,r_max,nrbins+1) # define binning (fixing r_min = 0.01)
r_hi = binfile[1:]
r_lo = binfile[:-1]
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

# Compute Corrfunc counts
RRp = DDrppi_mocks(1,2,nthreads,pimax,binfile,Ra,Dec,com_dist,weights1=W,weight_type='pair_product',
            verbose=True,is_comoving_dist=True,output_rpavg=True)

# Assemble counts
RRp_counts_raw=(RRp[:]['npairs']*RRp[:]['weightavg']).reshape((nrbins,-1))
rp_av = 0.5*(binfile[1:]+binfile[:-1])
pi_max_raw = RRp[:]['pimax'].reshape((nrbins,-1)).mean(axis=0)

# Average together pi bins, as desired
RRp_counts = np.sum([RRp_counts_raw[:,i::Delta_pi] for i in range(Delta_pi)],axis=0)
pi_max = pi_max_raw[Delta_pi-1::Delta_pi]

# Now compute ideal model for RR Counts
print("Computing correction function model")
pi_cen = np.arange(1/(2*len(pi_max)),max(pi_max),max(pi_max)/len(pi_max))

delta_pi = (pi_cen[-1]-pi_cen[-2])
delta_pi_all = delta_pi*np.ones_like(pi_cen).reshape(1,-1)
norm = np.sum(W)**2. # partial normalization

RRp_model = 2*np.pi*(r_hi**2.-r_lo**2.).reshape(-1,1)*delta_pi_all*norm

# Compute inverse Phi function
inv_phi = RRp_counts/RRp_model

# Now fit to a smooth model
all_rp = rp_av[:,None]*np.ones((1,len(pi_cen)))
all_pi = pi_cen[None,:]

scale = np.mean(inv_phi)

def inv_phi_model(input_par,scale=scale):
    par = np.asarray(input_par)*scale
    output = 0.
    output += (par[0]/(1+all_rp)+par[1]+par[2]*all_rp)
    output += (par[3]/(1+all_rp)+par[4]+par[5]*all_rp)*all_pi
    output += (par[6]/(1+all_rp)+par[7]+par[8]*all_rp)/(1+all_pi)
    output = output*(all_rp>1)+(par[9]+par[10]*all_pi)*(all_rp<1)
    return output
def to_min(par):
    return np.sum((inv_phi_model(par, scale)-inv_phi)**2./scale**2)

from scipy.optimize import minimize
coeff = minimize(to_min,np.asarray([1 for _ in range(11)])).x*scale
np.savetxt(outfile,coeff,delimiter="\t")

print("Saved correction function to %s"%outfile)
