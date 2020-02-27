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

// Define module files
#include "power_mod/power_parameters.h"
#include "power_mod/cell_utilities.h"
#include "power_mod/grid.h"
#include "power_mod/driver.h"
#ifdef BISPECTRUM
		#include "power_mod/kernel_interp_bispectrum.h"
		#include "power_mod/triple_counter.h"
		#include "power_mod/bispectrum_counts.h"
#else
		#include "power_mod/survey_correction_legendre.h"
		#include "power_mod/kernel_interp.h"
		#include "power_mod/pair_counter.h"
		#include "power_mod/power_counts.h"
#endif

STimer TotalTime;



// =========================== main() ==================================

int main(int argc, char *argv[]) {

	Parameters par=Parameters(argc,argv);

    // Now read in particles to grid:
		int n_grid = 2;

		Grid all_grid[n_grid]; // create empty grids
    for(int index=0;index<n_grid;index++){
        Float3 shift;
        Particle *orig_p;
				Float this_np;
				// Load in particles.
				// NB: For bispectrum the second grid is a random particle file created internally
#ifdef BISPECTRUM
				if (index==0){
#else
				if(true){
#endif
            char *filename;
#ifndef BISPECTRUM
            if(index==0) filename=par.fname;
						else filename=par.fname2;
#else
						filename=par.fname;
#endif
						if(index>0) if(strcmp(filename,par.fname)==0){
								all_grid[index].copy(&all_grid[0]);
								continue;
						}
	          orig_p = read_particles(par.rescale, &par.np, filename, par.rstart, par.nmax);
            assert(par.np>0);
						this_np = par.np;
            par.perbox = compute_bounding_box(orig_p, par.np, par.rect_boxsize, par.cellsize, par.rmax, shift, par.nside);
        }
#ifdef BISPECTRUM
				else {
        // If you want to just make random particles instead (as appropriate for bispectra)
        assert(par.np>0);
				// Generate some number of randoms (scaled to number of data points)
				this_np = par.np*par.f_rand;
				orig_p = make_particles(par.rect_boxsize, this_np);
        }
#endif

        if (par.qinvert) invert_weights(orig_p, this_np);
        if (par.qbalance) balance_weights(orig_p, this_np);

        // Now ready to compute!
        // Sort the particles into the grid.
        Grid tmp_grid(orig_p, this_np, par.rect_boxsize, par.cellsize, par.nside, shift, 1.);

        Float grid_density = (double)this_np/tmp_grid.nf;
        printf("\nCATALOG %d DIAGNOSTICS:\n",index+1);
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
        printf("#          Negative particles sum to %f\n\n", tmp_grid.sumw_neg);

        // Now save grid to global memory:
        all_grid[index].copy(&tmp_grid);


        free(orig_p); // Particles are now only stored in grid

        fflush(NULL);
		}

		Float max_sep = par.R0;

    // Count number of second/third field cells enclosed by the maximum truncation radius
		// This only depends on the cellsize, so works for our random grid too
    Float cellsize = all_grid[1].cellsize;
    Float filled_vol = 4./3.*M_PI*pow(max_sep+2.*cellsize,3);
    int n_close = ceil(filled_vol/pow(cellsize,3)); // number of close cells

    // Define cell separations (dimensionless) within truncation radius
    Float3 cell_sep_close_tmp[n_close];
    int len = ceil((max_sep+cellsize/2.)/cellsize);
    int len_cell_sep_close=0.; // counter
    integer3 this_pos;

    for(int i=-1*len;i<=len;i++){
        for(int j=-1*len;j<=len;j++){
            for(int k=-1*len;k<=len;k++){
                if((sqrt(pow(i,2)+pow(j,2)+pow(k,2))*cellsize)<(max_sep+cellsize)){
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
    KernelInterp interp_func(&par);

    // RUN Pair/Triple Counter
#ifdef BISPECTRUM
		Float DDRII_counts[par.nbin*par.nbin*par.mbin]; // these are the random counts we will generate below + pass to data-data counter

		// First compute the counts with randoms
		group_counter(&all_grid[1],&all_grid[0],&par,&interp_func,cell_sep_close,len_cell_sep_close,true,DDRII_counts);

		// Now compute all the counts involving only data points
		group_counter(&all_grid[0],&all_grid[0],&par,&interp_func,cell_sep_close,len_cell_sep_close,false,DDRII_counts);

#else
		// Read in survey correction function
		SurveyCorrection sc(&par,1,1);

    group_counter(&all_grid[0],&all_grid[1],&par,&sc,&interp_func,cell_sep_close,len_cell_sep_close);
#endif

    rusage ru;
    getrusage(RUSAGE_SELF, &ru);

    fprintf(stderr,"# user CPU time used: %ld \n# system CPU time used: %ld \n# maximum resident set size: %ld \n# integral shared memory size: %ld \n# integral unshared data size: %ld \n# integral unshared stack size: %ld \n# page reclaims (soft page faults): %ld \n# page faults (hard page faults): %ld \n# swaps: %ld \n# block input operations: %ld \n# block output operations: %ld \n# IPC messages sent: %ld \n# IPC messages received: %ld \n# signals received: %ld \n# voluntary context switches: %ld \n# involuntary context switches: %ld \n",ru.ru_utime.tv_sec,ru.ru_stime.tv_sec,ru.ru_maxrss,ru.ru_ixrss,ru.ru_idrss,ru.ru_isrss,ru.ru_minflt,ru.ru_majflt,ru.ru_nswap,ru.ru_inblock,ru.ru_oublock,ru.ru_msgsnd,ru.ru_msgrcv,ru.ru_nsignals,ru.ru_nvcsw,ru.ru_nivcsw);

    return 0;
}
