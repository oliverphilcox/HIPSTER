## Utility function to compute the power spectrum from measured (weighted) RR, DR and DD counts.
## We assume data is NOT periodic (for periodic data, the power spectrum is created inside the C++ code).

import numpy as np
import sys

if len(sys.argv)<8:
#    print("Please specify input parameters in the form {DD_FILE} {DR_FILE} {RR_FILE} {GAL_FILE} {N_RAND_RR} {N_RAND_DR} {OUTFILE}")
    print("Please specify input parameters in the form {DD_FILE} {DR_FILE} {RR_FILE} {GAL_FILE} {RAN_FILE_DR} {RAN_FILE_RR} {OUTFILE}")
    sys.exit()

DD_input = str(sys.argv[1])
DR_input = str(sys.argv[2])
RR_input = str(sys.argv[3])
gal_file = str(sys.argv[4])
ran_file_DR = str(sys.argv[5])
ran_file_RR = str(sys.argv[6])
#N_rand_RR = int(sys.argv[5])
#N_rand_DR = int(sys.argv[6])
outfile = str(sys.argv[7])

# Load pair counts
DD = np.loadtxt(DD_input)
DR = np.loadtxt(DR_input)
RR = np.loadtxt(RR_input)

# Compute normalization factor
print('Computing normalization factors')
N_gal = 0
with open(gal_file) as infile:
    for l,line in enumerate(infile):
        this_w = float(line.split()[3])
        N_gal += this_w
N_ran_DR = 0
with open(ran_file_DR) as infile:
    for l,line in enumerate(infile):
        this_w = float(line.split()[3])
        N_ran_DR += this_w
N_ran_RR = 0
with open(ran_file_RR) as infile:
    for l,line in enumerate(infile):
        this_w = float(line.split()[3])
        N_ran_RR += this_w

norm = N_gal**2. # to correct for Phi normalization 
print("Norm = %.2f"%norm)

# Compute number of pairs (for normalization)
DD_pairs = N_gal**2.
DR_pairs = N_gal*N_ran_DR
RR_pairs = N_ran_RR**2.

# Compute power
power = (DD-2.*DR/(DR_pairs/DD_pairs)+RR/(RR_pairs/DD_pairs))/norm

# Save to file
np.savetxt(outfile,power)
print('Saved power spectrum estimate to %s'%outfile)
