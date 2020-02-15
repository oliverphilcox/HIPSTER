import numpy as np

k_bins = np.loadtxt('/tigress/ophilcox/hipster_data/k_bin_quijote.csv')#quijote.csv')#k_bin_small3.csv')
k_mean = np.mean(k_bins,axis=1)

import readfof,MAS_library as MASL
import Pk_library as PKL

import sys
N = int(sys.argv[1])

n1 = N%len(k_mean)
n2 = N//(len(k_mean))
print(n1,n2)

# input files
snapdir = '/projects/QUIJOTE/Halos/fiducial/1/' #folder hosting the catalogue
snapnum = 4                                            #redshift 0

# determine the redshift of the catalogue
z_dict = {4:0.0, 3:0.5, 2:1.0, 1:2.0, 0:3.0}
redshift = z_dict[snapnum]

# read the halo catalogue
FoF = readfof.FoF_catalog(snapdir, snapnum, long_ids=False,
                          swap=False, SFR=False, read_IDs=False)

# get the properties of the halos
pos_h = FoF.GroupPos/1e3            #Halo positions in Mpc/h                                                                                                                                                                       
mass  = FoF.GroupMass*1e10          #Halo masses in Msun/h

good_pos = pos_h[mass>3.1e13]

# input parameters
grid    = 512  
BoxSize = 1000. #Mpc/h
MAS     = 'CIC'

# define the array hosting the density field
delta = np.zeros((grid,grid,grid), dtype=np.float32)

good_pos = good_pos.astype(np.float32)   #pos should be a numpy float array

# compute density field
MASL.MA(good_pos,delta,BoxSize,MAS)

# compute overdensity field
delta /= np.mean(delta, dtype=np.float64);  delta -= 1.0 

mu_grid = np.arange(-1.,1.,0.02)
#theta_grid = np.arange(0.,2.*np.pi,np.pi/10.)
theta_grid = np.arccos(mu_grid)

k1 = k_mean[n1]
k2 = k_mean[n2]

print("Computing bispectrum for k1,k2 = (%.2f, %.2f)"%(k1,k2))
Bk = PKL.Bk(delta, BoxSize, k1, k2, theta_grid, MAS='CIC', threads=20)
#all_Bk.append(Bk)
#all_k1k2.append([k1,k2])

np.savez('/tigress/ophilcox/pylians_bispectra/pylians_bispectrum_%d.npz'%N,k1=k1,k2=k2,Bk=Bk,theta=theta_grid,mu=mu_grid)
#2=all_k1k2,all_Bk=all_Bk,all_theta=theta_grid)
