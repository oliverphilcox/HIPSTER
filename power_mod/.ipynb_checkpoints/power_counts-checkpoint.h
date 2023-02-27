// power_counts.h for the grid_power.cpp code - Oliver Philcox 2019

#ifndef POWER_COUNTS_H
#define POWER_COUNTS_H

class PowerCounts{

private:
    int nbin,mbin;
    Float rmin,rmax,R0; //Ranges in r and truncation radius
    Float *r_high, *r_low; // Max and min of each radial bin
    Float *power_counts; // Power counts
    bool box; // to decide if we have a periodic box
    char *out_file, *out_string;
    int max_legendre; // maximum order of Legendre polynomials needed
    SurveyCorrection *sc; // survey correction function
    KernelInterp *kernel_interp;
#ifdef PERIODIC
    Float power_norm;
#endif

public:
    uint64 used_pairs; // total number of pairs used

    void cleanup_l(Float3 p1,Float3 p2,Float& norm,Float& mu){
        Float3 pos=p1-p2;
        norm = pos.norm();
#ifndef PERIODIC
        Float3 los=p1+p2; // No 1/2 as normalized anyway below
        mu = fabs(pos.dot(los)/norm/los.norm());
#else
        // In the periodic case use z-direction for mu
        mu = fabs(pos.z/norm);
#endif
        }

public:
    void sum_counts(PowerCounts *pc){
        // Add counts accumulated in different threads
        for(int i=0;i<nbin*mbin;i++) power_counts[i]+=pc->power_counts[i];
        used_pairs+=pc->used_pairs;
    }

public:
    PowerCounts(Parameters *par, SurveyCorrection *_sc, KernelInterp *_kernel_interp){
        nbin = par->nbin; // number of radial bins
        mbin = par->mbin; // number of Legendre bins
        out_file = par->out_file; // output directory
        out_string = par->out_string; // type specifier string
        R0 = par-> R0; // truncation radius

        kernel_interp = new KernelInterp(_kernel_interp);

        sc = _sc;
        max_legendre = fmax((sc->l_bins-1)*2,par->max_l);

#ifdef PERIODIC
        // define power spectrum normalization if periodic = (sum_weights_1)*(sum_weights_2)/V
        power_norm = (par->sum_w1*par->sum_w2)/(par->rect_boxsize[0]*par->rect_boxsize[1]*par->rect_boxsize[2]);
#endif

        int ec=0;
        ec+=posix_memalign((void **) &power_counts, PAGE, sizeof(double)*nbin*mbin);
        assert(ec==0);

        reset();

        box=par->perbox;

        rmax=par->rmax;
        rmin=par->rmin;

        r_high = par->radial_bins_high;
        r_low = par->radial_bins_low;
    }

    void reset(){
        for(int j=0;j<nbin*mbin;j++){
            power_counts[j]=0.;
        }
        used_pairs=0.;
    }

    ~PowerCounts(){
        free(power_counts);
    }

#ifdef PERIODIC
    void randoms_analytic(Float* RR_analyt){
      // Compute the analytic randoms term from numerical integration of the window function in each bin.

      printf("\nComputing analytic RR term");
      fflush(NULL);
      // Define r
      Float integrand,tmp_r,old_kr,old_kr3,old_kernel,new_kernel,new_kr,new_kr3,diff_kr;
      int npoint=200000;
      Float delta_r = R0/double(npoint);

      for(int i=0;i<nbin;i++) RR_analyt[i]=0; // initialize array

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
              RR_analyt[i]+=integrand*(new_kernel-old_kernel)/diff_kr; // ie 4 pi r^2 W(r) j_0^a(r) dr
              
              // Update values for next time
              old_kr = new_kr;
              old_kr3 = new_kr3;
              old_kernel = new_kernel;
          }
      }

    // Now save output
    char RR_periodic_name[1000];
    Float tmp_out;

    snprintf(RR_periodic_name, sizeof RR_periodic_name, "%s/%s_analyt_RR_power_counts_n%d_l%d_R0%d.txt", out_file, out_string, nbin, 2*(mbin-1),int(R0));
    FILE * RR_periodic_file = fopen(RR_periodic_name,"w");

    for (int i=0;i<nbin;i++){
        for(int j=0;j<mbin;j++){
            if(j==0) tmp_out = RR_analyt[i];
            else tmp_out = 0.;
            fprintf(RR_periodic_file,"%le\t",tmp_out);
        }
        fprintf(RR_periodic_file,"\n");
    }
    fflush(NULL);

    // Close open files
    fclose(RR_periodic_file);

    printf("\nComputation complete, and saved to %s.",RR_periodic_name);

    }
#endif

    inline void count_pairs(Particle pi, Particle pj,Float3 separation){
        // Count pairs of particles with some weighting function
        // separation is the separation of the cell centers (needed for periodic datasets)
        Float r_ij, mu_ij, w_ij, tmp_phi_inv=0, new_kr, old_kr,new_kernel,diff_kr,new_kr3, old_kr3;
        Float legendre[max_legendre/2+1]; // Even-ell Legendre coefficients
        Float old_kernel[mbin]; // saved kernel functions

#ifdef PERIODIC
        pj.pos+=separation; // add on separation of cell centers (since periodic particles are defined relative to cell center)
#endif
        cleanup_l(pi.pos,pj.pos,r_ij,mu_ij); // define distance and angle

        if(r_ij>=R0) return; // outside correct size
        if(r_ij==0){
            printf("Zero separation! This indicates a bug.\n");
            exit(1);
        }
        if(r_ij==0) return;

        used_pairs++; // update number of pairs used

        // First define all required legendre moments used for kernel and Phi function
        legendre[0]=1.;
        if(max_legendre>1){
            Float mu2 = mu_ij*mu_ij;
            legendre[1]=0.5*(3.*mu2-1.);
            if(max_legendre>3){
                Float mu4 = mu2*mu2;
                legendre[2]=1./8.*(35.*mu4-30.*mu2+3.);
                if(max_legendre>5){
                    Float mu6 = mu4*mu2;
                    legendre[3] = 1./16.*(231.*mu6-315.*mu4+105.*mu2-5.);
                    if(max_legendre>7){
                        Float mu8 = mu6*mu2;
                        legendre[4] = 1./128.*(6435.*mu8-12012.*mu6+6930.*mu4-1260.*mu2+35.);
                        if(max_legendre>9){
                            legendre[5] = 1./256.*(46189.*mu8*mu2-109395.*mu8+90090.*mu6-30030.*mu4+3465.*mu2-63.);
                        }
                    }
                }
            }
        }
        
    #ifdef LYA
        // Define rp and pi distances
        Float3 pos = pi.pos-pj.pos;
        Float3 los=pi.pos+pj.pos; // No 1/2 as normalized anyway below
        Float pi_ij = abs(pos.dot(los)/los.norm()); // note we only need abs(pi) here by parity
        Float rp_ij = sqrt(r_ij*r_ij-pi_ij*pi_ij);
        
        // Compute phi inverse 
        tmp_phi_inv = sc->inv_correction_function(rp_ij, pi_ij);
    #else
        for(int l_i=0;l_i<sc->l_bins;l_i++) tmp_phi_inv+=legendre[l_i]*sc->inv_correction_function(l_i*2,r_ij);
    #endif
        
        w_ij = pi.w*pj.w*pair_weight(r_ij)/tmp_phi_inv;
        
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
            for(int j=0;j<mbin;j++){
                if(i==0) old_kernel[j] = kernel_interp->kernel(j*2,old_kr);
                new_kernel = kernel_interp->kernel(j*2,new_kr);
                power_counts[i*mbin+j]+=w_ij*legendre[j]*(new_kernel-old_kernel[j])/diff_kr;
                old_kernel[j] = new_kernel;
            }
            // Update values for next time
            old_kr = new_kr;
            old_kr3 = new_kr3;
        }

    }

    inline Float pair_weight(Float sep){
        // Compute weight function W(r;R_0)
        // NB: input r should never be greater than R0!
        if(sep<R0/2.) return 1.;
        if(sep>R0){
            printf("\nParticle separation should not be larger than R0, this indicates a bug!");
            exit(1);
        }
        else{
            Float x = sep/R0;
            if(sep<3.*R0/4.) return 1.-8.*pow(2.*x-1.,3.)+8.*pow(2.*x-1.,4.);
            else return -64.*pow(x-1.,3.)-128.*pow(x-1.,4.);
        }
    }

    void save_counts(int one_grid) {
        // Print power-count output to file.
        // Create output files

        char pow_name[1000];
  #ifdef PERIODIC
        snprintf(pow_name, sizeof pow_name, "%s/%s_DD_power_counts_n%d_l%d_R0%d.txt", out_file,out_string,nbin, 2*(mbin-1),int(R0));
  #else
        snprintf(pow_name, sizeof pow_name, "%s/%s_power_counts_n%d_l%d_R0%d.txt", out_file,out_string,nbin, 2*(mbin-1),int(R0));
  #endif
        FILE * PowFile = fopen(pow_name,"w");

        for (int i=0;i<nbin;i++){
            for(int j=0;j<mbin;j++){
                if(one_grid==1) power_counts[i*mbin+j]*=2.; // since we ignore i-j switches
                fprintf(PowFile,"%le\t",power_counts[i*mbin+j]);
            }
            fprintf(PowFile,"\n");
        }

        fflush(NULL);

        // Close open files
        fclose(PowFile);
      }

#ifdef PERIODIC
  void save_spectrum(Float* RR_analytic){
        // If periodic, we can output the whole power spectrum estimate here
        char pk_name[1000];
        Float output_pk;
        printf("Norm = %.2e\n",power_norm);
        snprintf(pk_name, sizeof pk_name, "%s/%s_power_spectrum_n%d_l%d_R0%d.txt", out_file,out_string,nbin, 2*(mbin-1),int(R0));
        FILE * PkFile = fopen(pk_name,"w");

        for (int i=0;i<nbin;i++){
            for(int j=0;j<mbin;j++){
                output_pk = power_counts[i*mbin+j]/power_norm;
                if(j==0) output_pk-=RR_analytic[i];
                fprintf(PkFile,"%le\t",output_pk);
            }
            fprintf(PkFile,"\n");
        }

        fflush(NULL);

        // Close open files
        fclose(PkFile);
}
#endif

};

#endif
