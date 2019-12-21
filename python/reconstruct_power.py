## Utility function to compute the power spectrum from measured (weighted) RR, DR and DD counts.
## We assume data is NOT periodic (for periodic data, the power spectrum is created inside the C++ code).

import numpy as np
import sys

if len(sys.argv)<8:
    print("Please specify input parameters in the form {DD_FILE} {DR_FILE} {RR_FILE} {GAL_FILE} {N_RAND_RR} {N_RAND_DR} {OUTFILE}")
    sys.exit()

DD_input = str(sys.argv[1])
DR_input = str(sys.argv[2])
RR_input = str(sys.argv[3])
gal_file = str(sys.argv[4])
N_rand_RR = int(sys.argv[5])
N_rand_DR = int(sys.argv[6])
outfile = str(sys.argv[7])

# Load pair counts
DD = np.loadtxt(DD_input)
DR = np.loadtxt(DR_input)
RR = np.loadtxt(RR_input)

# Compute normalization factor
print('Computing normalization factor')
norm = 0
N_gal = 0
with open(gal_file) as infile:
    for l,line in enumerate(infile):
        this_w = float(line.split()[3])
        norm+=this_w**2.
        N_gal+=1
print("Norm = %.2f"%norm)

# Compute number of pairs (for normalization)
DR_pairs = N_gal*N_rand_DR
DD_pairs = N_gal*(N_gal-1)
RR_pairs = N_rand_RR*(N_rand_RR-1.)

# Compute power
power = (DD-2.*DR/(DR_pairs/DD_pairs)+RR/(RR_pairs/DD_pairs))/norm

# Save to file
np.savetxt(outfile,power)
print('Saved power spectrum estimate to %s'%outfile)
