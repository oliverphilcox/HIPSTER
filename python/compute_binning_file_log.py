## Utility function to create a k-binning scheme for power spectrum computation (for logarithmic binning with base e)
## (Taken from RascalC code)

import sys
import numpy as np

if len(sys.argv)<5:
    print("Please specify input parameters in the form {N_LOG_BINS} {MIN_K} {MAX_K} {OUTPUT_FILE}.")
    sys.exit()
nkbins = int(sys.argv[1])
k_min = float(sys.argv[2])
assert k_min>0,'Minimum k must be greater than zero to take logarithm'
k_max = float(sys.argv[3])
out_file = str(sys.argv[4])

print("Using LOG binning");

# Define radial bins
kbins = np.logspace(np.log(k_min),np.log(k_max),nkbins+1,base=np.e)

# PRINT binning:
with open(out_file,'w+') as writefile:
    for i in range(nkbins-1):
        writefile.write("%.8f\t%.8f\n" %(kbins[i],kbins[i+1]))
    writefile.write("%.8f\t%.8f" %(kbins[nkbins-1],kbins[nkbins]))
print("Binning file '%s' written successfully."%out_file)
