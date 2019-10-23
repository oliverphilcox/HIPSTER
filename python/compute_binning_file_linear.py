## Utility function to create a k-binning scheme for power spectrum computation
## (Taken from RascalC code)

import sys
import numpy as np

if len(sys.argv)<5:
    print("Please specify input parameters in the form {N_LINEAR_BINS} {MIN_K} {MAX_K} {OUTPUT_FILE}.")
    sys.exit()
nkbins = int(sys.argv[1])
k_min = float(sys.argv[2])
k_max = float(sys.argv[3])
out_file = str(sys.argv[4])

print("Using LINEAR binning");

# Define radial bins
kbins = np.linspace(k_min,k_max,nkbins+1)

if kbins[0]<=1e-4:
    kbins[0]=1e-4 # for stability

# PRINT binning:
with open(out_file,'w+') as writefile:
    for i in range(nkbins-1):
        writefile.write("%.8f\t%.8f\n" %(kbins[i],kbins[i+1]))
    writefile.write("%.8f\t%.8f" %(kbins[nkbins-1],kbins[nkbins]))
print("Binning file '%s' written successfully."%out_file)
