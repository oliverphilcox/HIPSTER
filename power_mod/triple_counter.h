// triple_counter.h - Module to run the triple counting for grid_power.cpp in BISPECTRUM mode - Oliver Philcox 2019

#ifndef GROUP_COUNTER_H
#define GROUP_COUNTER_H

#include "bispectrum_counts.h"

class group_counter{

private:
    int nbin,mbin;

public:
    void check_threads(Parameters *par,int print){
        // Set up OPENMP and define which threads to use
#ifdef OPENMP
        cpu_set_t mask[par->nthread+1];
        int tnum=0;
        sched_getaffinity(0, sizeof(cpu_set_t), &mask[par->nthread]);
        if(print==1) fprintf(stderr, "CPUs used are: ");
        for(int ii=0;ii<64;ii++){
            if(CPU_ISSET(ii, &mask[par->nthread])){
                if(print==1) fprintf(stderr,"%d ", ii);
                CPU_ZERO(&mask[tnum]);
                CPU_SET(ii,&mask[tnum]);
                tnum++;
            }
        }
        fprintf(stderr,"\n");
#endif
    }

public:

    group_counter(Grid *grid1, Grid *grid2, Parameters *par, KernelInterp *kernel_interp, integer3 *cell_sep, int len_cell_sep, bool random_counts, Float *DDRII_counts){
        // if random_counts = true, we just need to compute a DDR term with R at the center of the triangle
        // this is the same as for the DDD counts

        // Define parameters
        nbin = par->nbin; // number of radial bins
        mbin = par->mbin; // number of Legendre bins
        Float percent_counter=0.,used_cells=0,cell_attempts=0;
        uint64 used_particles=0;

        STimer initial, TotalTime; // time initialization
        initial.Start();

        //-----------INITIALIZE OPENMP + CLASSES----------
        BispectrumCounts global_counts(par,kernel_interp,random_counts);

        check_threads(par,1); // define threads

        fprintf(stderr, "Init time: %g s\n",initial.Elapsed());
        printf("# 1st grid filled cells: %d\n",grid1->nf);
        printf("# All 1st grid points in use: %d\n",grid1->np);
        printf("# Max points in one cell in grid 1%d\n",grid1->maxnp);
        fflush(NULL);

        TotalTime.Start(); // Start timer

        if(random_counts) printf("\n### COMPUTING THE DATA+RANDOM PARTICLE COUNTS\n\n");
        else printf("\n### COMPUTING THE DATA+DATA PARTICLE COUNTS\n\n");

#ifdef OPENMP
        #pragma omp parallel firstprivate(par,grid1,grid2,kernel_interp) shared(global_counts,TotalTime) reduction(+:percent_counter,used_cells,used_particles)
        { // start parallel loop
        // Decide which thread we are in
        int thread = omp_get_thread_num();
        assert(omp_get_num_threads()<=par->nthread);
        if (thread==0) printf("# Starting power-spectrum particle-count computation on %d threads.\n\n", omp_get_num_threads());
#else
        int thread = 0;
        printf("# Starting power-spectrum particle-count computation single threaded.\n\n");
        { // start loop
#endif

            // DEFINE LOCAL THREAD VARIABLES
            Particle *prim_list; // list of particles in first cell
            Float3 separation; // separation vector
            int* prim_ids; // list of particle IDs in first cell
            int prim_id_1D,sec_id_1D; // 1D cell ids
            integer3 prim_id,delta,sec_id; // particle positions in grid
            int mnp = grid1->maxnp; // max number of particles in grid1 cell
            Cell prim_cell,sec_cell; // cell objects
            Particle particle_i; // primary particle
            Float3* sep_register; // list of separation vectors between primary and secondary
            int register_index; // index of register holding particle separations

            BispectrumCounts loc_counts(par,kernel_interp,random_counts);
            // Assign memory for intermediate steps
            int ec=0;
            ec+=posix_memalign((void **) &prim_list, PAGE, sizeof(Particle)*mnp);
            ec+=posix_memalign((void **) &prim_ids, PAGE, sizeof(int)*mnp);
            ec+=posix_memalign((void **) &sep_register, PAGE, sizeof(Float3)*mnp*len_cell_sep); // hold maximum number of close particles
            assert(ec==0);

#ifdef OPENMP
#pragma omp for schedule(dynamic)
#endif
            // Loop over all filled n1 cells
            for(int n1=0;n1<grid1->nf;n1++){

                // Print time left
                if((float(n1)/float(grid1->nf)*100)>=percent_counter){
                    if(percent_counter>0){
                        if(random_counts) printf("Data+Random Counts: Counting cell %d of %d on thread %d: %.0f percent complete\n",n1+1,grid1->nf,thread,percent_counter);
                        else printf("Data+Data Counts: Counting cell %d of %d on thread %d: %.0f percent complete\n",n1+1,grid1->nf,thread,percent_counter);
                    }
                    percent_counter+=5.;
                }

                // Pick first cell
                prim_id_1D = grid1->filled[n1]; // 1d ID for cell i
                prim_id = grid1->cell_id_from_1d(prim_id_1D); // define first cell
                prim_cell = grid1->c[prim_id_1D];
                if(prim_cell.np==0) continue; // skip if empty

                cell_attempts+=len_cell_sep; // total number of cells used

                // Select one particle from this cell at a time
                for(int i=prim_cell.start;i<(prim_cell.start+prim_cell.np);i++){
                    particle_i = grid1->p[i];
                    register_index=0;

                    // Now find all secondary particles in allowed region around this primary particle
                    // Iterate over all nearby second cells
                    for(int n2=0;n2<len_cell_sep;n2++){ // second cell index
                        delta = cell_sep[n2]; // difference in cell positions
                        sec_id = prim_id + delta;
                        sec_id_1D = grid2->test_cell(sec_id);
                    #ifdef PERIODIC
                        separation = grid2->cell_sep(delta); // separation vector between grid cells
                    #else
                        separation = {0,0,0}; // zero vector (unused)
                    #endif
                        if(sec_id_1D<0) continue; // if cell not in grid
                        //if((one_grid==1)&&(sec_id_1D<prim_id_1D)) continue; // already counted this pair of cells!
                        sec_cell = grid2->c[sec_id_1D];
                        if(sec_cell.np==0) continue; // if empty cell
                        if(i==prim_cell.start) used_cells++; // update number of cells used on first time round i-loop

                        // Now we have a primary particle and a cell full of secondary particles
                        // Iterate over secondary particles in the cell

                        for(int j=sec_cell.start;j<(sec_cell.start+sec_cell.np);j++){
                            if(j==i) continue; // skip if identical particles
                            //if((one_grid==1)&&(j<=i)) continue; // skip if already counted or identical particles (for same grids only)
                            used_particles++;
                            // Now save the distances of the cells to the register
                            sep_register[register_index++]=grid2->p[j].pos+separation-particle_i.pos;
                        }
                    }

                    // Now time to fill up the bispectrum counts
                    loc_counts.fill_up_counts(sep_register, register_index);
                }
            }
#ifdef OPENMP
#pragma omp critical // only one processor at once
#endif
        {
            // Sum up power-sums
            global_counts.sum_counts(&loc_counts);
            loc_counts.reset();
        }

    } // end OPENMP loop


    // ----- REPORT AND SAVE OUTPUT ------------
    TotalTime.Stop();

    int runtime = TotalTime.Elapsed();
    fflush(NULL);
    if(random_counts) printf("\nDATA+RANDOM TRIPLE COUNTS COMPLETE\n");
    else printf("\nDATA+DATA TRIPLE COUNTS COMPLETE\n");
    printf("\nTotal process time for %.2e particle pairs: %d s, i.e. %2.2d:%2.2d:%2.2d hms\n", double(global_counts.used_pairs),runtime, runtime/3600,runtime/60%60,runtime%60);
    printf("We tried %.2e pairs of cells and accepted %.2e pairs of cells.\n", double(cell_attempts),double(used_cells));
    printf("Cell acceptance ratio is %.3f.\n",(double)used_cells/cell_attempts);
    printf("We tried %.2e pairs of particles and accepted %.2e pairs of particles.\n", double(used_particles),double(global_counts.used_pairs));
    printf("Particle acceptance ratio is %.3f.\n",(double)global_counts.used_pairs/used_particles);
    printf("Average of %.2f pairs accepted per primary particle.\n\n",(Float)global_counts.used_pairs/grid1->np);
    printf("Trial speed: %.2e cell pairs per core per second\n",double(used_cells)/(runtime*double(par->nthread)));
    printf("Acceptance speed: %.2e particle pairs per core per second\n",double(global_counts.used_pairs)/(runtime*double(par->nthread)));

  if(random_counts){
    global_counts.save_counts2();
    printf("DDRII counts saved to %s/%s_DRR_II_counts_n%d_l%d_R0%d.txt\n", par->out_file,par->out_string,nbin, (mbin-1),int(par->R0));
    // Now copy in DDR counts to local array
    for(int i=0;i<nbin*nbin*mbin;i++) DDRII_counts[i] = global_counts.bispectrum_counts[i];
  }
  // Compute RRR analytic counts (only need to do this once)
  else{
    Float* Wka; // array to hold W_a(R0) functions (with RRR_ab = W_a(R0) * W_b(R0))
    int ec=0;
    ec+=posix_memalign((void **) &Wka, PAGE, sizeof(Float)*nbin);
    assert(ec==0);
    global_counts.randoms_analytic(Wka);

    global_counts.save_counts(Wka);
    printf("DDD counts saved to %s/%s_DDD_counts_n%d_l%d_R0%d.txt\n", par->out_file,par->out_string,nbin, (mbin-1),int(par->R0));
    global_counts.save_spectrum(Wka,DDRII_counts);
    printf("Printed full bispectrum to file as %s/%s_bispectrum_n%d_l%d_R0%d.txt\n\n", par->out_file,par->out_string,nbin, (mbin-1),int(par->R0));
  }
}
};

#endif
