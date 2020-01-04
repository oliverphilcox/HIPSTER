// bispectrum_counts.h for the grid_power.cpp code - Oliver Philcox 2019

#ifndef BISPECTRUM_COUNTS_H
#define BISPECTRUM_COUNTS_H

class BispectrumCounts{

private:
    int nbin,mbin,n_lm, n_mult;
    Float rmin,rmax,R0; //Ranges in r and truncation radius
    Float *r_high, *r_low; // Max and min of each radial bin
    Float *bispectrum_counts; // Power counts
    bool box; // to decide if we have a periodic box
    char *out_file, *out_string;
    int max_legendre; // maximum order of Legendre polynomials needed
    SurveyCorrection *sc; // survey correction function
    KernelInterp *kernel_interp;
#ifdef PERIODIC
    Float bispectrum_norm;
#define MAXORDER 10
    int map[MAXORDER+1][MAXORDER+1][MAXORDER+1];   // The multipole index of x^a y^b z^c

#endif

public:
    uint64 used_triples; // total number of pairs used

public:
    void sum_counts(BispectrumCounts *bc){
        // Add counts accumulated in different threads
        for(int i=0;i<nbin*nbin*mbin;i++) bispectrum_counts[i]+=bc->bispectrum_counts[i];
        used_triples+=bc->used_triples;
    }

public:
    BispectrumCounts(Parameters *par, SurveyCorrection *_sc, KernelInterp *_kernel_interp){
        nbin = par->nbin; // number of radial bins
        mbin = par->mbin; // number of Legendre bins
        out_file = par->out_file; // output directory
        out_string = par->out_string; // type specifier string
        R0 = par-> R0; // truncation radius

        kernel_interp = new KernelInterp(_kernel_interp);

        make_map(); // Generate spherical harmonic coefficients

        sc = _sc;
        max_legendre = par->max_l;

        n_lm = (max_legendre+1)*(max_legendre+2)/2; // number of Y_lm bins up to maximum ell
        n_mult = ((max_legendre+1)*(max_legendre+2)*(max_legendre+3)/6) // total number of Cartesian multipoles, satisfying a+b+c<=max-legendre

        // define power spectrum normalization if periodic = n^3 V = N^3 / V^2
        bispectrum_norm = pow(par->np,3.)/pow(par->rect_boxsize[0]*par->rect_boxsize[1]*par->rect_boxsize[2],2.);

        int ec=0;
        ec+=posix_memalign((void **) &bispectrum_counts, PAGE, sizeof(double)*nbin*nbin*mbin);
        assert(ec==0);

        reset();

        box=par->perbox;

        rmax=par->rmax;
        rmin=par->rmin;

        r_high = par->radial_bins_high;
        r_low = par->radial_bins_low;
    }

    void reset(){
        for(int j=0;j<nbin*nbin*mbin;j++){
            bispectrum_counts[j]=0.;
        }
        used_triples=0.;
    }

    ~BispectrumCounts(){
        free(bispectrum_counts);
    }

#ifdef PERIODIC
    void randoms_analytic(Float* RRR_analyt){
      // Compute the analytic randoms term from numerical integration of the window function in each bin.
      // NB: This is the outer product of two W(k) functions in each bin

      printf("\nComputing analytic RRR term");
      fflush(NULL);
      // Define r
      Float integrand,tmp_r,old_kr,old_kr3,old_kernel,new_kernel,new_kr,new_kr3,diff_kr;
      int npoint=100000;
      Float delta_r = R0/double(npoint);
      Float RRR_analyt_tmp[nbin];

      for(int i=0;i<nbin;i++) RRR_analyt_tmp[i]=0; // initialize array which will hold W(k) in correct bin.

      for(int j=0;j<npoint;j++){ // iterate over all r in [0,R0]
          tmp_r = double(j+1)/double(npoint)*R0;

          integrand = 4*M_PI*pow(tmp_r,2)*pair_weight(tmp_r)*delta_r;  // 4 pi r^2 W(r) dr

          for(int i=0;i<nbin;i++){ // iterate over all k-bins
              // Now compute multipole contributions assuming k-bins are contiguous
              if(i==0){
                  old_kr = tmp_r*r_low[0];
                  old_kr3 = pow(old_kr,3);
                  old_kernel = kernel_interp->kernel(0,old_kr);
              }
              new_kr = tmp_r*r_high[i];
              new_kr3 = pow(new_kr,3);
              diff_kr = new_kr3 - old_kr3;

              // Now compute kernel at ell=0
              new_kernel = kernel_interp->kernel(0,new_kr);

              // Add contribution to integral
              RRR_analyt_tmp[i]+=integrand*(new_kernel-old_kernel)/diff_kr; // ie 4 pi r^2 W(r) j_0^a(r) dr

              // Update values for next time
              old_kr = new_kr;
              old_kr3 = new_kr3;
              old_kernel = new_kernel;
          }
      }

    // Now save output
    char RRR_periodic_name[1000];
    Float tmp_out;

    snprintf(RRR_periodic_name, sizeof RRR_periodic_name, "%s/%s_analyt_RRR_counts_n%d_l%d.txt", out_file, out_string, nbin, 2*(mbin-1));
    FILE * RRR_periodic_file = fopen(RRR_periodic_name,"w");

    for (int i1=0;i1<nbin;i1++){
        for (int i2=0;i2<nbin;i2++){
            RRR_analyt[i1*nbin+i2]=RRR_analyt_tmp[i1]*RRR_analyt_tmp[i2];
            for(int j=0;j<mbin;j++){
                if(j==0) tmp_out = RRR_analyt_tmp[i1]*RRR_analyt_tmp[i2];
                else tmp_out = 0.;
                fprintf(RRR_periodic_file,"%le\t",tmp_out);
            }
            fprintf(RRR_periodic_file,"\n");
        }
    }
    fflush(NULL);

    // Close open files
    fclose(RRR_periodic_file);

    printf("\nComputation complete, and saved to %s.",RRR_periodic_name);

    }
#endif

    inline void fill_up_counts(Float3 *separation_register, int register_size){
        // Count triples of particles with some weighting function
        // Input is the list of all separations of a primary particle with secondaries and length of this list.
        // Full bispectrum kernel is made from this.

        // TODO: Define this earlier + allocate memory?
        Float partial_kernel[register_size*nbin]; // NB: partial kernel will be overwritten for each ell bins
        Float new_kr, old_kr, new_kr3, old_kr3, old_kernel, new_kernel,tmp_counts,particle_sep,k1;
        Float fi, fij, fijk,L_ell;
        Float m[n_mult]; // powers of x,y,z used for Y_lms
        Complex alm[n_lm*register_size];   // Y_lm (for m>0) for each particle

        // First iterate over all particles in register
        for(int n=0;n<register_size;n++){
            if(separation_register[n]==0){
                  fprintf(stderr,"Zero found in register; this suggests a self-count problem.\n");
                  exit(1);
            }

            // Compute separation length and check if <R0;
            particle_sep = separation_register.norm();
            if(particle_sep<R0){
                partial_kernel[n*nbin]=-1;
                continue;
            }
            sep_weight = pair_weight(particle_sep); // compute W(r;R_0)

            //// STORE Y_LM COEFFICIENTS

            // First consider all possible powers of separation vector x,y,z coordinates
            norm_sep = separation_register[n]/particle_sep;

            Float *mm = m; // grab a pointer to m
            fi = 1.;
            for(int i=0;i<=max_legendre;i++) {
                fij = fi;
                for(int j=0;j<=max_legendre-i;j++) {
                    fijk = fij;
                    for(int k=0;k<=max_legendre-i-j;k++) {
                        *mm += fijk;
                		    fijk *= norm_sep.z;
                		    mm++;
                    }
                    fij *= norm_sep.y;
                }
                fi *= norm_sep.x;
            }

#define CM(a,b,c) m[map[a][b][c]]
            Complex *almbin = &(alm[n][0]);

            // Next fill up Y_lms (using tedious spherical harmonic code in another file)
            #include "spherical_harmonics.cpp"

            // We now have Y_lm values (unnormalized) for each particle pair.

            // Now let's compute the k-binning kernels in each k-bin and ell
            for(int ell=0;ell<max_legendre;ell++){

                // Now iterate over k-bin (assuming they are contiguous)
                for(int i=0;i<nbin;i++){

                    if(i==0){
                        old_kr = particle_sep*r_low[0];
                        old_kr3 = pow(old_kr,3);
                        old_kernel = kernel_interp->kernel(ell,old_kr);
                    }

                    new_kr = particle_sep*r_high[i];
                    new_kr3 = pow(new_kr,3);
                    new_kernel = kernel_interp->kernel(ell,new_kr);

                    // Fill up array of (binned) j_ell(kr)W(kr) values
                    partial_kernel[n*nbin+i] = sep_weight*(new_kernel-old_kernel)/(new_kr3-old_kr3);

                    // Save the kernel functions for the next bin
                    old_kr = new_kr;
                    old_kr3 = new_kr3;
                    old_kernel = new_kernel;
                }
            }

        // We now have Y_lm and the kernels, i.e. the full kernel for each pair of particles.
        // Now take outer product of these and add to relevant bins

#define RealProduct(a,b) (a.real()*b.real()+a.imag()*b.imag())

        // Now loop over both sets of pairs of particles (excluding self-counts)
        for(int n=0;n<register_size-1;n++){
            if(partial_kernel[n*nbin]==-1) continue; // outside R0
            k1 = partial_kernel[n*nbin+i];
            for(int n2=n+1;n2<register_size;n2++){
                if(partial_kernel[n2*nbin]==-1) continue; // outside R0

                // Now loop over Legendre multipoles
                for (int ell=0, nn=0; ell<=max_legendre; ell++) {
                    L_ell = 0;
                    // Now compute the L_ell part from the pair;
                    for (int mm=0; mm<=ell; mm++, nn++) {
                    // nn counts the (ell,m)
                    // Our definition is that the power is 4*pi/(2*l+1) times alm alm*.
                    L_ell += RealProduct(alm[n][nn],alm[n2][nn])*almnorm[nn]* 4.0 * M_PI / (2*ell+1.0);

                    // Now loop over k-bins
                    for(int i=0;i<nbin;i++){
                        k1 = partial_kernel[n1*nbin+i];
                        for(int i2=0;i2<nbin;i2++){
                            bispectrum_counts[(i*nbin+i2)*mbin+ell] += L_ell*partial_kernel[n2*nbin+i2];
                        }
                    }
                }
            }
        }
    }

        // Now add in the Legendre moments (can be done separately)

                        used_triples++; // update number of pairs used



                // NB: we implicitly assume survey correction function is unity here.
                //for(int l_i=0;l_i<sc->l_bins;l_i++) tmp_phi_inv+=legendre[l_i];

                w_ijk = pi.w*pj.w*pk.w*pair_weight(r_ij)*pair_weight(r_ik);///tmp_phi_inv;

                // NB: This assumes bins are contiguous
                for(int i=0;i<nbin;i++){

                    // Now compute multipole contributions
                    if(i==0){
                        old_kr = r_ij*r_low[0];
                        old_kr3 = pow(old_kr,3);
                    }
                    new_kr = r_ij*r_high[i];
                    new_kr3 = pow(new_kr,3);
                    diff_kr = new_kr3 - old_kr3;

                    for(int i2=0;i2<nbin;i2++){

                        // Now compute multipole contributions
                        if(i2==0){
                            old_kr2 = r_ik*r_low[0];
                            old_kr32 = pow(old_kr2,3);
                        }
                        new_kr2 = r_ik*r_high[i2];
                        new_kr32 = pow(new_kr2,3);
                        diff_kr2 = new_kr32 - old_kr32;

                        for(int j=0;j<mbin;j++){
                            //TODO: Should put a lot of these functions before the second loop
                            if(i==0) old_kernel[j] = kernel_interp->kernel(j*2,old_kr);
                            if(i2==0) old_kernel2[j] = kernel_interp->kernel(j*2,old_kr2);
                            new_kernel = kernel_interp->kernel(j*2,new_kr);
                            new_kernel2 = kernel_interp->kernel(j*2,new_kr2);
                            bispectrum_counts[(i*nbin+i2)*mbin+j]+=w_ijk*legendre[j]*(new_kernel-old_kernel[j])/diff_kr*(new_kernel2-old_kernel2[j])/diff_kr2;
                            old_kernel[j] = new_kernel;
                            old_kernel2[j] = new_kernel2;
                        }
                        // Update values for next time

                        old_kr2 = new_kr2;
                        old_kr32 = new_kr32;
                      }
                  old_kr = new_kr;
                  old_kr3 = new_kr3;
                }
            }


    inline Float pair_weight(Float sep){
        // Compute weight function W(r;R_0)
        if(sep<R0/2) return 1.;
        else{
            Float x = sep/R0;
            if(sep<3*R0/4) return 1.-8.*pow(2*x-1,3)+8.*pow(2*x-1,4);
            else return -64.*pow(x-1,3)-128.*pow(x-1,4);
        }
    }

    void save_counts(int one_grid) {
        // Print bispectrum-count output to file.
        // Create output files

        char count_name[1000];
  #ifdef PERIODIC
        snprintf(count_name, sizeof count_name, "%s/%s_DDD_counts_n%d_l%d.txt", out_file,out_string,nbin, (mbin-1));
  #else
        snprintf(count_name, sizeof count_name, "%s/%s_bispectrum_counts_n%d_l%d.txt", out_file,out_string,nbin, (mbin-1));
  #endif
        FILE * CountFile = fopen(count_name,"w");

        for (int i=0;i<nbin;i++){
            for(int i2=0;i2<nbin;i2++){
                for(int j=0;j<mbin;j++){
                    //if(one_grid==1) power_counts[(i*nbin+i2)*mbin+j]*=2.; // since we ignore i-j switches
                    fprintf(CountFile,"%le\t",bispectrum_counts[(i*nbin+i2)*mbin+j]);
                }
                fprintf(CountFile,"\n");
            }
        }

        fflush(NULL);

        // Close open files
        fclose(CountFile);
      }

#ifdef PERIODIC
  void save_spectrum(Float* RRR_analytic){
        // If periodic, we can output the whole bispectrum estimate here
        char bkk_name[1000];
        Float output_bkk;
        printf("Norm = %.2e\n",bispectrum_norm);
        snprintf(bkk_name, sizeof bkk_name, "%s/%s_bispectrum_n%d_l%d.txt", out_file,out_string,nbin, (mbin-1));
        FILE * BkkFile = fopen(bkk_name,"w");

        for (int i=0;i<nbin;i++){
            for (int i2=0;i2<nbin;i2++){
                for(int j=0;j<mbin;j++){
                    output_bkk = bispectrum_counts[(i*nbin+i2)*mbin+j]/bispectrum_norm;
                    if(j==0) output_bkk+=2*RRR_analytic[i];
                fprintf(BkkFile,"%le\t",output_bkk);
                }
            fprintf(BkkFile,"\n");
            }
        }

        fflush(NULL);

        // Close open files
        fclose(BkkFile);
}
#endif

// SPHERICAL HARMONIC INFORMATION


void make_map() {
	// Construct the index number in our multipoles for x^a y^b z^c
        for(int i=0;i<=MAXORDER;i++)
            for(int j=0;j<=MAXORDER-i;j++)
                for(int k=0;k<=MAXORDER-i-j;k++) map[i][j][k] = 0;
	int n=0;
        for(int i=0;i<=ORDER;i++)
            for(int j=0;j<=ORDER-i;j++)
                for(int k=0;k<=ORDER-i-j;k++) {
		    map[i][j][k] = n; n++;
		}
	return;
}

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

};

#endif
