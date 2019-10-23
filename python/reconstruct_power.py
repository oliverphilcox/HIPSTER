## Utility function to compute the power spectrum from measured (weighted) RR, DR and DD counts.

import numpy as np
import sys

if len(sys.argv)<9:
    print("Please specify input parameters in the form {DD_FILE} {DR_FILE} {RR_FILE} {GAL_FILE} {N_RAND_RR} {N_RAND_DR} {PERIODIC} {OUTFILE}")
    sys.exit()

DD_input = str(sys.argv[1])
DR_input = str(sys.argv[2])
RR_input = str(sys.argv[3])
gal_file = str(sys.argv[4])
N_rand_RR = int(sys.argv[5])
N_rand_DR = int(sys.argv[6])
periodic = int(sys.argv[7])
outfile = str(sys.argv[8])

if periodic:
    print('Assuming periodic boundary conditions')
else:
    print('Assuming non-periodic boundary conditions (and FKP weights)')

# Load pair counts
DD = np.loadtxt(DD_input)
DR = np.loadtxt(DR_input)
RR = np.loadtxt(RR_input)

# Compute normalization factor
print('Computing normalization factor')
if periodic:
    norm = 0.
    gal_x,gal_y,gal_z,gal_w = np.loadtxt(gal_file).T
    N_gal = len(gal_x)
    x_range = max(gal_x)-min(gal_x)
    y_range = max(gal_y)-min(gal_y)
    z_range = max(gal_z)-min(gal_z)

    if np.abs(x_range/y_range-1.)>0.01 or np.abs(x_range/z_range-1.)>0.01:
        print('Is this data periodic? Different x, y and z axis dimensions found. Exiting')
        sys.exit();

    V = x_range*y_range*z_range
    norm = np.sum(gal_w**2.)
    print('Norm = %.2f'%norm)
else:
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
