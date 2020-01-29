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

    // Triple count arrays
    int map[MAXORDER+1][MAXORDER+1][MAXORDER+1]; // The multipole index of x^a y^b z^c
    Complex *Aa_lm; // A^a_lm vector
    Float *Cab_l; // C^{ab}_lm vector
    Float *j0a_W; // j_0^a(|x_i-x_j|)W(|x_i-x_j|)
    Float *DDR_II_sum; // E_l^{II,ab}(|x_i-x_j|) sum
    Float *kernel_arr; // j^a_ell(kr)W(r) for specific ell

    Float *m; // powers of x,y,z used for Y_lms
    Complex *alm;   // Y_lm (for m>0) for a particle

#ifdef PERIODIC
    Float bispectrum_norm,bispectrum_norm2;
#endif

public:
    uint64 used_pairs; // total number of pairs used

public:
    void sum_counts(BispectrumCounts *bc){
        // Add counts accumulated in different threads
        for(int i=0;i<nbin*nbin*mbin;i++) bispectrum_counts[i] += bc->bispectrum_counts[i];
        for(int i=0;i<nbin;i++) j0a_W[i] += bc->j0a_W[i];
        for(int i=0;i<nbin*nbin*mbin;i++) DDR_II_sum[i] += bc->DDR_II_sum[i];
        used_pairs+=bc->used_pairs;
    }

public:
    BispectrumCounts(Parameters *par, SurveyCorrection *_sc, KernelInterp *_kernel_interp){
        nbin = par->nbin; // number of radial bins
        mbin = par->mbin; // number of Legendre bins
        out_file = par->out_file; // output directory
        out_string = par->out_string; // type specifier string
        R0 = par-> R0; // truncation radius

        kernel_interp = new KernelInterp(_kernel_interp);

        sc = _sc;
        max_legendre = par->max_l;

        n_lm = (max_legendre+1)*(max_legendre+2)/2; // number of Y_lm bins up to maximum ell
        n_mult = ((max_legendre+1)*(max_legendre+2)*(max_legendre+3)/6); // total number of Cartesian multipoles, satisfying a+b+c<=max-legendre

        make_map(); // Generate spherical harmonic coefficients

        // define power spectrum normalization if periodic = n^3 V = N^3 / V^2, and n^2V = N^2/V for DDR term
        bispectrum_norm = pow(par->np,3.)/pow(par->rect_boxsize[0]*par->rect_boxsize[1]*par->rect_boxsize[2],2.);
        bispectrum_norm2 = pow(par->np,2.)/pow(par->rect_boxsize[0]*par->rect_boxsize[1]*par->rect_boxsize[2],1.);

        // Generate necessary arrays
        int ec=0;
        ec+=posix_memalign((void **) &bispectrum_counts, PAGE, sizeof(double)*nbin*nbin*mbin);
        ec+=posix_memalign((void **) &Aa_lm, PAGE, sizeof(Complex)*n_lm*nbin);
        ec+=posix_memalign((void **) &Cab_l, PAGE, sizeof(Float)*nbin*nbin*mbin);
        ec+=posix_memalign((void **) &DDR_II_sum, PAGE, sizeof(Float)*nbin*nbin*mbin);
        ec+=posix_memalign((void **) &j0a_W, PAGE, sizeof(Float)*nbin);
        ec+=posix_memalign((void **) &kernel_arr, PAGE, sizeof(Float)*nbin);

        ec+=posix_memalign((void **) &m, PAGE, sizeof(Float)*n_mult); // powers of x,y,z used for Y_lms
        ec+=posix_memalign((void **) &alm, PAGE, sizeof(Complex)*n_lm);   // Y_lm (for m>0) for a particle
        assert(ec==0);

        // Zero relevant arrays;
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
            DDR_II_sum[j]=0.;
        }
        for(int i=0;i<nbin;i++) j0a_W[i]=0;
        used_pairs=0.;
    }

    ~BispectrumCounts(){
        free(bispectrum_counts);
        free(Aa_lm);
        free(Cab_l);
        free(j0a_W);
        free(DDR_II_sum);
        free(kernel_arr);
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

    snprintf(RRR_periodic_name, sizeof RRR_periodic_name, "%s/%s_analyt_RRR_counts_n%d_l%d.txt", out_file, out_string, nbin, (mbin-1));
    FILE * RRR_periodic_file = fopen(RRR_periodic_name,"w");

    for (int i1=0;i1<nbin;i1++){
        for (int i2=0;i2<nbin;i2++){
            RRR_analyt[i1*nbin+i2]=RRR_analyt_tmp[i1]*RRR_analyt_tmp[i2];
            for(int j=0;j<mbin;j++){
                if(j==0) tmp_out = RRR_analyt[i1*nbin+i2];
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

        if(register_size==0) return;
        Float sep_weight,particle_sep, tmp_kernel,tmp_E,tmp_count,this_kernel,old_kr,old_kr3,new_kr,new_kr3,old_kernel,new_kernel;
        Float3 norm_sep;
        int start_n;

        // Must first zero useful arrays - ARE ALL THESE NECESSARY?
        for(int mm=0;mm<n_lm*nbin;mm++) Aa_lm[mm]=0;
        for(int mm=0;mm<nbin*nbin*mbin;mm++) Cab_l[mm]=0;

        // First iterate over all particles in register and compute W(r;R_0)j_ell^a and Y_lm;
        for(int n=0;n<register_size;n++){
            used_pairs++; // update number of pairs used

            // Compute separation length
            particle_sep = separation_register[n].norm();

            // Compute the E_II sum for all particles up to 2*R0;
            if((particle_sep>=2*R0)) continue;
            for(int n1=0;n1<nbin;n1++){
                for(int n2=n1;n2<nbin;n2++){
                    for(int ell=0;ell<=max_legendre;ell++){
                        tmp_E = kernel_interp->kernel_E_II(n1,n2,ell,particle_sep);
                        DDR_II_sum[(n1*nbin+n2)*mbin+ell] += tmp_E;
                        // Use symmetry to fill other bin
                        if(n1!=n2) DDR_II_sum[(n2*nbin+n1)*mbin+ell] += tmp_E;
                    }
                }
            }

            // Only need separation up to R0 for all other terms
            if((particle_sep>=R0)) continue;
            if((particle_sep==0)){
                  fprintf(stderr,"Zero found in register; this suggests a self-count problem.\n");
                  exit(1);
            }

            // Zero arrays - Needed?
            for(int mm=0;mm<n_mult;mm++) m[mm]=0; // powers of x,y,z used for Y_lms
            for(int mm=0;mm<n_lm;mm++) alm[mm]=0;

            sep_weight = pair_weight(particle_sep); // compute W(r;R_0)

            // Compute Y_lm for particle pair in the alm variable
            norm_sep = separation_register[n]/particle_sep;
            compute_Ylm(norm_sep, alm, m);
            //
            // // Compute L_ell for checking
            // #define RealProduct(a,b) (a.real()*b.real()+a.imag()*b.imag())
            // int nnn=0;
            // for(int elll=0;elll<=max_legendre;elll++){
            //     tmp_count=0.;
            //     for(int mmm=0;mmm<=elll;mmm++,nnn++){
            //         printf("%.2e ",RealProduct(alm[nnn],alm[nnn]));
            //         tmp_count+=RealProduct(alm[nnn],alm[nnn])*almnorm[nnn];
            //     }
            //     Float Lell = tmp_count*4*M_PI/(2*elll+1);
            //
            //     Float mu_ang = norm_sep.z/norm_sep.norm(); // z-axis angle
            //     printf("l = %d, L_l = %.2e vs %.2e\n",elll,Lell,Pn(mu_ang,elll));
            //  }

            //for(int nn=0;nn<n_lm;nn++) printf("%.2e %.2e\n",alm[nn].real(),alm[nn].imag());

            // Now let's compute the k-binning kernels in each k-bin and ell
            start_n=0; // counts starting location in Y_lm array for certain ell
            for(int ell=0;ell<=max_legendre;ell++){

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

                    // Compute W(r)j^a_ell(r) values
                    tmp_kernel = sep_weight*(new_kernel-old_kernel)/(new_kr3-old_kr3);

                    // Fill up array of (binned) j^a_ell(r)W(r)
                    kernel_arr[i]=tmp_kernel;

                    // Assign j0a_W to array
                    if(ell==0) j0a_W[i]+=tmp_kernel;

                    // Fill up array of (binned) j^a_ell(r)W(r)Y_lm(r) values
                    for(int mm=0;mm<=ell;mm++){
                      Aa_lm[i*n_lm+start_n+mm]+= tmp_kernel*alm[start_n+mm];
                    }
                    // Save the kernel functions for the next bin
                    old_kr = new_kr;
                    old_kr3 = new_kr3;
                    old_kernel = new_kernel;

                }
                start_n+=ell+1;

                // Now these are filled we can compute Cab_l;
                for(int i=0;i<nbin;i++){
                  this_kernel = kernel_arr[i];
                    for(int i2=0;i2<nbin;i2++){
                        Cab_l[(i*nbin+i2)*mbin+ell]+=this_kernel*kernel_arr[i2];
                    }
                }
            }
        }

        // We can now perform the sum over m to compute A^a_lm A^b*_lm term;
#define RealProduct(a,b) (a.real()*b.real()+a.imag()*b.imag())

        // Loop over bispectrum bins
        for(int i=0;i<nbin;i++){
            for(int i2=0;i2<nbin;i2++){
                for(int ell=0, nn=0;ell<=max_legendre;ell++){
                    tmp_count=0;
                    for (int mm=0; mm<=ell; mm++, nn++){ // Loop over m
                        tmp_count+=RealProduct(Aa_lm[i*n_lm+nn],Aa_lm[i2*n_lm+nn])*almnorm[nn];
                    }
                    // Added to binned estimate and remove self-count
                    bispectrum_counts[(i*nbin+i2)*mbin+ell]+= tmp_count * 4.0 * M_PI / (2*ell+1.0) - Cab_l[(i*nbin+i2)*mbin+ell];
                }
            }
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
        // DDD contains unnormalized bispectrum counts
        // DDR-I contains sum over j_0^a(r)W(r;R_0)  - must be multiplied by -2 delta^K_{ell,0}/n^2V * tilde{W}^a for spectral contribution
        // DDR-II contains sum over E^{II,ab}_ell(r) - unnormalized

        char count_name[1000], count_name2[1000], count_name3[1000];
        snprintf(count_name, sizeof count_name, "%s/%s_DDD_counts_n%d_l%d.txt", out_file,out_string,nbin, (mbin-1));
        snprintf(count_name2, sizeof count_name2, "%s/%s_DDR_I_counts_n%d_l%d.txt", out_file,out_string,nbin, (mbin-1));
        snprintf(count_name3, sizeof count_name3, "%s/%s_DDR_II_counts_n%d_l%d.txt", out_file,out_string,nbin, (mbin-1));

        FILE * CountFile = fopen(count_name,"w");
        FILE * CountFile2 = fopen(count_name2,"w");
        FILE * CountFile3 = fopen(count_name3,"w");

        for (int i=0;i<nbin;i++){
            fprintf(CountFile2,"%le\n",j0a_W[i]); // print the DDR-I term stochastic component
            for(int i2=0;i2<nbin;i2++){
                for(int j=0;j<mbin;j++){
                    //if(one_grid==1) power_counts[(i*nbin+i2)*mbin+j]*=2.; // since we ignore i-j switches
                    fprintf(CountFile,"%le\t",bispectrum_counts[(i*nbin+i2)*mbin+j]);
                    fprintf(CountFile3,"%le\t",DDR_II_sum[(i*nbin+i2)*mbin+j]);
                }
                fprintf(CountFile,"\n");
                fprintf(CountFile3,"\n");
            }
        }

        fflush(NULL);

        // Close open files
        fclose(CountFile);
        fclose(CountFile2);
        fclose(CountFile3);
      }

#ifdef PERIODIC
  void save_spectrum(Float* RRR_analytic){
        // If periodic, we can output the whole bispectrum estimate here
        char bkk_name[1000];
        Float output_bkk, tmp_Wka;
        printf("Norm = %.2e\n",bispectrum_norm);
        snprintf(bkk_name, sizeof bkk_name, "%s/%s_bispectrum_n%d_l%d.txt", out_file,out_string,nbin, (mbin-1));
        FILE * BkkFile = fopen(bkk_name,"w");

        for (int i=0;i<nbin;i++){
            // Compute W(k) in bin i (this is square-root of diagonal of RRR_analytic)
            tmp_Wka = pow(RRR_analytic[i*nbin+i],0.5);
            for (int i2=0;i2<nbin;i2++){
                for(int j=0;j<mbin;j++){
                    // Compute DDD term
                    output_bkk = bispectrum_counts[(i*nbin+i2)*mbin+j]/bispectrum_norm;

                    // Add DDR-I term
                    if(j==0) output_bkk-=2.*j0a_W[i2]*tmp_Wka/bispectrum_norm2;

                    // Add DDR-II term
                    output_bkk -= DDR_II_sum[(i*nbin+i2)*mbin+j]/bispectrum_norm2;

                    if(j==0) output_bkk+=2*RRR_analytic[i*nbin+i2];
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

inline void compute_Ylm(Float3 norm_sep, Complex* alm, Float* m){
      // Compute Y_lm coefficients given input normalized vector and holder of coefficients + empty array for powers.
      Float *mm = m;
      Float fi, fij, fijk;
      // First consider all possible powers of separation vector x,y,z coordinates
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

while (true){
#define CM(a,b,c) m[map[a][b][c]]
      Complex *almbin = &(alm[0]);

      // Next fill up Y_lms (using tedious spherical harmonic code in another file)
      #include "spherical_harmonics.h"
}
      //for(int mm=0;mm<n_lm;mm++) printf("%d %.2e %.2e\n",mm,alm[mm].real(),alm[mm].imag());
}

inline void make_map() {
  // Construct the index number in our multipoles for x^a y^b z^c
        for(int i=0;i<=MAXORDER;i++)
            for(int j=0;j<=MAXORDER-i;j++)
                for(int k=0;k<=MAXORDER-i-j;k++) map[i][j][k] = 0;
  int n=0;
        for(int i=0;i<=max_legendre;i++)
            for(int j=0;j<=max_legendre-i;j++)
                for(int k=0;k<=max_legendre-i-j;k++) {
        map[i][j][k] = n; n++;
		}
	return;
}

};


#endif
