// grid_power.cpp -- Oliver Philcox - based on Oliver Philcox's RascalC code

#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex>
#include <assert.h>
#include "threevector.hh"
#include "STimer.cc"
#include <sys/stat.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

// For multi-threading:
#ifdef OPENMP
#include <omp.h>
#include <sched.h>
#endif

// In order to not print the output matrices

#define PAGE 4096     // To force some memory alignment.

typedef unsigned long long int uint64;
typedef unsigned int uint;
#ifdef BISPECTRUM
typedef std::complex<double> Complex;
#endif

// Could swap between single and double precision here.
typedef double Float;
typedef double3 Float3;

// BISPECTRUM SPECIFIC CODE
#ifdef BISPECTRUM
	#define MAXORDER 10
	#define NLM_MAX ((MAXORDER+1)*(MAXORDER+2)/2)

	// Some global constants for the a_lm normalizations.
	// From: http://en.wikipedia.org/wiki/Table_of_spherical_harmonics
	// Including an extra factor of 2 in all m!=0 cases.

	// All factors are of the form a*sqrt(b/pi), so let's use that:
	#define YNORM(a,b) (2.0*(1.0*a)*(1.0*a)*(1.0*b)/M_PI)
	// This includes the factor of 2, so divide the m=0 by 2.

	static Float almnorm[NLM_MAX] = {
	    YNORM(1/2,1)/2.0,

	    YNORM(1/2,3)/2.0,
	    YNORM(1/2,3/2),

	    YNORM(1/4,5)/2.0,
	    YNORM(1/2,15/2),
	    YNORM(1/4,15/2),

	    YNORM(1/4,7)/2.0,
	    YNORM(1/8,21),
	    YNORM(1/4,105/2),
	    YNORM(1/8,35),

	    YNORM(3/16,1)/2.0,
	    YNORM(3/8,5),
	    YNORM(3/8,5/2),
	    YNORM(3/8, 35),
	    YNORM(3/16, 35/2),

	    YNORM(1/16, 11)/2.0,
	    YNORM(1/16, 165/2),
	    YNORM(1/8, 1155/2),
	    YNORM(1/32, 385),
	    YNORM(3/16, 385/2),
	    YNORM(3/32, 77),

	    YNORM(1/32, 13)/2.0,
	    YNORM(1/16, 273/2),
	    YNORM(1/64, 1365),
	    YNORM(1/32, 1365),
	    YNORM(3/32, 91/2),
	    YNORM(3/32, 1001),
	    YNORM(1/64, 3003),

	    YNORM(1/32, 15)/2.0,
	    YNORM(1/64, 105/2),
	    YNORM(3/64, 35),
	    YNORM(3/64, 35/2),
	    YNORM(3/32, 385/2),
	    YNORM(3/64, 385/2),
	    YNORM(3/64, 5005),
	    YNORM(3/64, 715/2),

	    YNORM(1/256, 17)/2.0,
	    YNORM(3/64, 17/2),
	    YNORM(3/128, 595),
	    YNORM(1/64, 19635/2),
	    YNORM(3/128, 1309/2),
	    YNORM(3/64, 17017/2),
	    YNORM(1/128, 7293),
	    YNORM(3/64, 12155/2),
	    YNORM(3/256, 12155/2),

	    YNORM(1/256, 19)/2.0,
	    YNORM(3/256, 95/2),
	    YNORM(3/128, 1045),
	    YNORM(1/256, 21945),
	    YNORM(3/128, 95095/2),
	    YNORM(3/256, 2717),
	    YNORM(1/128, 40755),
	    YNORM(3/512, 13585),
	    YNORM(3/256, 230945/2),
	    YNORM(1/512, 230945),

	    YNORM(1/512, 21)/2.0,
	    YNORM(1/256, 1155/2),
	    YNORM(3/512, 385/2),
	    YNORM(3/256, 5005),
	    YNORM(3/256, 5005/2),
	    YNORM(3/256, 1001),
	    YNORM(3/1024, 5005),
	    YNORM(3/512, 85085),
	    YNORM(1/512, 255255/2),
	    YNORM(1/512, 4849845),
	    YNORM(1/1024, 969969)
	};
#endif

// Define module files
#include "power_mod/power_parameters.h"
#include "power_mod/cell_utilities.h"
#include "power_mod/grid.h"
#include "power_mod/driver.h"
#include "power_mod/survey_correction_legendre.h"
#ifdef BISPECTRUM
		#include "power_mod/kernel_interp_bispectrum.h"
		#include "power_mod/triple_counter.h"
		#include "power_mod/bispectrum_counts.h"
#else
		#include "power_mod/kernel_interp.h"
		#include "power_mod/pair_counter.h"
		#include "power_mod/power_counts.h"
#endif

STimer TotalTime;



// =========================== main() ==================================

int main(int argc, char *argv[]) {

	Parameters par=Parameters(argc,argv);

    // Now read in particles to grid:
#ifdef BISPECTRUM
		int n_grid = 3;
#else
		int n_grid = 2;
#endif
		Grid all_grid[n_grid]; // create empty grids
    for(int index=0;index<n_grid;index++){
        Float3 shift;
        Particle *orig_p;
        if (!par.make_random){
            char *filename;
            if(index==0) filename=par.fname;
						else if(index==1) filename=par.fname2;
            else filename=par.fname3;
            orig_p = read_particles(par.rescale, &par.np, filename, par.rstart, par.nmax);
            assert(par.np>0);
            par.perbox = compute_bounding_box(orig_p, par.np, par.rect_boxsize, par.cellsize, par.rmax, shift, par.nside);
        } else {
        // If you want to just make random particles instead:
        assert(par.np>0);
        orig_p = make_particles(par.rect_boxsize, par.np);
        par.cellsize = par.rect_boxsize.x/float(par.nside);
        // set as periodic if we make the random particles
        par.perbox = true;
        }

        if (par.qinvert) invert_weights(orig_p, par.np);
        if (par.qbalance) balance_weights(orig_p, par.np);

        // Now ready to compute!
        // Sort the particles into the grid.
        Grid tmp_grid(orig_p, par.np, par.rect_boxsize, par.cellsize, par.nside, shift, 1.);

        Float grid_density = (double)par.np/tmp_grid.nf;
        printf("\n RANDOM CATALOG %d DIAGNOSTICS:\n",index+1);
        printf("Grid cell-size = %.2fMpc/h\n", tmp_grid.cellsize);
        printf("Average number of particles per grid cell = %6.2f\n", grid_density);
        Float max_density = 200.0;
        if (grid_density>max_density){
            fprintf(stderr,"Average particle density exceeds maximum advised particle density (%.0f particles per cell) - exiting.\n",max_density);
            exit(1);
        }
        if (grid_density<0.01){
            printf("#\n# WARNING: grid appears inefficiently fine; exiting.\n#\n");
            exit(1);
        }

        printf("# Done gridding the particles\n");
        printf("# %d particles in use, %d with positive weight\n", tmp_grid.np, tmp_grid.np_pos);
        printf("# Weights: Positive particles sum to %f\n", tmp_grid.sumw_pos);
        printf("#          Negative particles sum to %f\n", tmp_grid.sumw_neg);

        // Now save grid to global memory:
        all_grid[index].copy(&tmp_grid);

        free(orig_p); // Particles are now only stored in grid

        fflush(NULL);
    }

    // Read in survey correction function
    SurveyCorrection sc(&par,1,1);

    // Count number of second/third field cells enclosed by the maximum truncation radius
    Float cellsize = all_grid[1].cellsize;
    Float filled_vol = 4./3.*M_PI*pow(par.R0+2.*cellsize,3);
    int n_close = ceil(filled_vol/pow(cellsize,3)); // number of close cells

    // Define cell separations (dimensionless) within truncation radius (same for power spectrum or bispectrum)
    Float3 cell_sep_close_tmp[n_close];
    int len = ceil((par.R0+cellsize/2.)/cellsize);
    int len_cell_sep_close=0.; // counter
    integer3 this_pos;

    for(int i=-1*len;i<=len;i++){
        for(int j=-1*len;j<=len;j++){
            for(int k=-1*len;k<=len;k++){
                if((sqrt(pow(i,2)+pow(j,2)+pow(k,2))*cellsize)<(par.R0+cellsize)){
                    this_pos = {i,j,k};
                    cell_sep_close_tmp[len_cell_sep_close] = this_pos;
                    len_cell_sep_close++;
                }
            }
        }
    }
    assert(len_cell_sep_close<=n_close);

    integer3 cell_sep_close[len_cell_sep_close]; // proper array to house cell separations of correct length
    for(int i=0;i<len_cell_sep_close;i++) cell_sep_close[i] = cell_sep_close_tmp[i];

    // Compute kernel interpolation functions
    printf("Creating kernel interpolator function\n");
    KernelInterp interp_func(par.R0,par.rmin,par.rmax);

    // RUN Pair/Triple Counter
#ifdef BISPECTRUM
		group_counter(&all_grid[0],&all_grid[1],&all_grid[2],&par,&sc,&interp_func,cell_sep_close,len_cell_sep_close);
#else
    group_counter(&all_grid[0],&all_grid[1],&par,&sc,&interp_func,cell_sep_close,len_cell_sep_close);
#endif

    rusage ru;
    getrusage(RUSAGE_SELF, &ru);

    fprintf(stderr,"# user CPU time used: %ld \n# system CPU time used: %ld \n# maximum resident set size: %ld \n# integral shared memory size: %ld \n# integral unshared data size: %ld \n# integral unshared stack size: %ld \n# page reclaims (soft page faults): %ld \n# page faults (hard page faults): %ld \n# swaps: %ld \n# block input operations: %ld \n# block output operations: %ld \n# IPC messages sent: %ld \n# IPC messages received: %ld \n# signals received: %ld \n# voluntary context switches: %ld \n# involuntary context switches: %ld \n",ru.ru_utime.tv_sec,ru.ru_stime.tv_sec,ru.ru_maxrss,ru.ru_ixrss,ru.ru_idrss,ru.ru_isrss,ru.ru_minflt,ru.ru_majflt,ru.ru_nswap,ru.ru_inblock,ru.ru_oublock,ru.ru_msgsnd,ru.ru_msgrcv,ru.ru_nsignals,ru.ru_nvcsw,ru.ru_nivcsw);

    return 0;
}
