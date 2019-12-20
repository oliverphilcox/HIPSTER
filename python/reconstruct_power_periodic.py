## Utility function to compute the power spectrum from measured (weighted) DD counts, assuming periodicity.

import numpy as np
import sys

if len(sys.argv)<5:
    print("Please specify input parameters in the form {DD_FILE} {RR_FILE} {GAL_FILE} {OUTFILE}")
    sys.exit()

DD_input = str(sys.argv[1])
RR_input = str(sys.argv[2])
gal_file = str(sys.argv[3])
outfile = str(sys.argv[4])

# Load pair counts
DD = np.loadtxt(DD_input)
RR_analyt = np.loadtxt(RR_input)

# Compute normalization factor
print('Computing normalization factor')
norm = 0.
gal_x,gal_y,gal_z = np.loadtxt(gal_file).T[:3]
print('Assigning unit weight to all particles.')
N_gal = len(gal_x)
x_range = max(gal_x)-min(gal_x)
y_range = max(gal_y)-min(gal_y)
z_range = max(gal_z)-min(gal_z)

if np.abs(x_range/y_range-1.)>0.01 or np.abs(x_range/z_range-1.)>0.01:
    print('Is this data periodic? Different x, y and z axis dimensions found. Exiting')
    sys.exit();

V = x_range*y_range*z_range
norm = N_gal*N_gal/V
print('Norm = %.2f'%norm)

# Compute power
power = DD/norm-RR_analyt

# Save to file
np.savetxt(outfile,power)
print('Saved power spectrum estimate to %s'%outfile)
