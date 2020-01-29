// kernel_interp_bispectrum.h - Interpolator for kernel functions for bispectrum computation
// This stores 3(-1)^ell D_ell(kw)

#ifndef KERNEL_INTERP_H
#define KERNEL_INTERP_H


class KernelInterpSingle{
    // Store an interpolator for a single E_II^{ab}_l function

  public:
      int npoint;
      Float R0;
      gsl_interp_accel *sep_a;
      gsl_spline *interp_func;
      double *r_vals, *kernel_vals;

  public:

      // empty initializer
      KernelInterpSingle(){}

      KernelInterpSingle(KernelInterpSingle *kern){
          // Copy constructor
          kern->copy(npoint,R0,&r_vals,&kernel_vals);

          interpolate();
      }


      void initialize(Float _R0, Float *_r_vals, Float *_kernel_vals, int _npoint){
          // Create the initializer

          //TODO - all use same accel vector?
        R0 = _R0;
        npoint = _npoint;

        r_vals = (double *)malloc(sizeof(double)*npoint);
        kernel_vals = (double *)malloc(sizeof(double)*npoint);

        for(int i=0;i<npoint;i++){
            r_vals[i] = _r_vals[i];
            kernel_vals[i] = _kernel_vals[i];
        }

        // TODO: do we actually need to carry all these arrays around?

        interpolate();

      }

      void interpolate(){
        // Construct interpolation object
        interp_func = gsl_spline_alloc(gsl_interp_cspline, npoint);
        gsl_spline_init(interp_func,r_vals, kernel_vals, npoint);
        sep_a = gsl_interp_accel_alloc();
      }

      ~KernelInterpSingle(){
          // Destructor
          free(kernel_vals);
          free(r_vals);
          gsl_spline_free(interp_func);
          gsl_interp_accel_free(sep_a);
      }

      void copy(int& _npoint, Float& _R0, double **_r_vals, double** _kernel_vals){
          _npoint = npoint;
          _R0 = R0;

          *_r_vals = (double *)malloc(sizeof(double)*npoint);
          *_kernel_vals = (double *)malloc(sizeof(double)*npoint);

          for(int i=0;i<npoint;i++){
            (*_r_vals)[i]=r_vals[i];
            (*_kernel_vals)[i] = kernel_vals[i];
          }

      }

      inline double kernel(double this_sep){
          // Compute kernel function accelerated by interpolation
          if(this_sep>=2*R0) return 0.; // kernel is identically zero beyond this
          return gsl_spline_eval(interp_func, this_sep, sep_a);

      }


};


class KernelInterp{
    // Compute a interpolation kernel function

private:
    Float R0,k_max,k_min;
    Float *k_high, *k_low; // max and min of each radial bin
    int max_legendre; // maximum order of Legendre polynomials needed
    char *out_file, *out_string;

public:
    int nbin,mbin;
    int npoint=1000;
    Float min_val,max_val;
    double *sep, *kernel_vals_0,*kernel_vals_1,*kernel_vals_2,*kernel_vals_3,*kernel_vals_4,*kernel_vals_5;
    double *kernel_vals_6,*kernel_vals_7,*kernel_vals_8,*kernel_vals_9,*kernel_vals_10;
    gsl_interp_accel *sep_a;
    gsl_spline *interp_kernel_0, *interp_kernel_1, *interp_kernel_2,*interp_kernel_3;
    gsl_spline *interp_kernel_4, *interp_kernel_5, *interp_kernel_6,*interp_kernel_7;
    gsl_spline *interp_kernel_8, *interp_kernel_9, *interp_kernel_10;

    // E_II functions
    Float *kernel_II_vals, *r_vals;
    // Create interpolators
    KernelInterpSingle *E_II_interpolators;

public:
    inline double kernel(int ell, double this_sep){
        // Kernel function accelerated by interpolation
        if(this_sep<=min_val) return 0.;
        else{
            if (ell==0) return gsl_spline_eval(interp_kernel_0,this_sep,sep_a);
            else if (ell==1) return gsl_spline_eval(interp_kernel_1,this_sep,sep_a);
            else if (ell==2) return gsl_spline_eval(interp_kernel_2,this_sep,sep_a);
            else if (ell==3) return gsl_spline_eval(interp_kernel_3,this_sep,sep_a);
            else if (ell==4) return gsl_spline_eval(interp_kernel_4,this_sep,sep_a);
            else if (ell==5) return gsl_spline_eval(interp_kernel_5,this_sep,sep_a);
            else if (ell==6) return gsl_spline_eval(interp_kernel_6,this_sep,sep_a);
            else if (ell==7) return gsl_spline_eval(interp_kernel_7,this_sep,sep_a);
            else if (ell==8) return gsl_spline_eval(interp_kernel_8,this_sep,sep_a);
            else if (ell==9) return gsl_spline_eval(interp_kernel_9,this_sep,sep_a);
            else if (ell==10) return gsl_spline_eval(interp_kernel_10,this_sep,sep_a);
            else{
                printf("Multipoles greater than ell = 10 not yet implemented");
                exit(1);
            }
        }
    }

public:
    void copy_function(KernelInterp *kern){
        // Copy pre-existing kernel function into object
        // NB: This isn't currently used, so may have memory leaks.
        npoint = kern->npoint;
        min_val = kern->min_val;
        nbin = kern->nbin;
        mbin = kern->mbin;

        // Allocate memory
        sep = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_0 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_1 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_2 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_3 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_4 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_5 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_6 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_7 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_8 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_9 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_10 = (double *)malloc(sizeof(double)*npoint);

        kernel_II_vals = (double *)malloc(sizeof(double)*npoint*nbin*nbin*mbin);
        r_vals = (double *)malloc(sizeof(double)*npoint);

        // Read in grids
        for(int i=0;i<npoint;i++){
            sep[i]=kern->sep[i];
            kernel_vals_0[i]=kern->kernel_vals_0[i];
            kernel_vals_1[i]=kern->kernel_vals_1[i];
            kernel_vals_2[i]=kern->kernel_vals_2[i];
            kernel_vals_3[i]=kern->kernel_vals_3[i];
            kernel_vals_4[i]=kern->kernel_vals_4[i];
            kernel_vals_5[i]=kern->kernel_vals_5[i];
            kernel_vals_6[i]=kern->kernel_vals_6[i];
            kernel_vals_7[i]=kern->kernel_vals_7[i];
            kernel_vals_8[i]=kern->kernel_vals_8[i];
            kernel_vals_9[i]=kern->kernel_vals_9[i];
            kernel_vals_10[i]=kern->kernel_vals_10[i];

            r_vals[i]=kern->r_vals[i];
        }

        for(int i=0;i<npoint*nbin*nbin*mbin;i++) kernel_II_vals[i] = kern->kernel_II_vals[i];

        // Remake E_II interpolators
        for(int i=0;i<nbin*nbin*mbin;i++){
            E_II_interpolators[i] = KernelInterpSingle(kern->E_II_interpolators[i]);
        }

        // activate interpolator
        interpolate();
    }

private:
    void copy(int& _npoint, int& _nbin, int&  _mbin, Float& _min_val, double **_sep, double **_kernel_vals_0,double **_kernel_vals_1,double **_kernel_vals_2,double **_kernel_vals_3,double **_kernel_vals_4,double **_kernel_vals_5,double **_kernel_vals_6,double **_kernel_vals_7,double **_kernel_vals_8,double **_kernel_vals_9,double **_kernel_vals_10,
    KernelInterpSingle **_E_II_interpolators){
        _npoint = npoint;
        _min_val = min_val;
        _nbin = nbin;
        _mbin = mbin;

        *_sep = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_0 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_1 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_2 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_3 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_4 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_5 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_6 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_7 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_8 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_9 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_10 = (double *)malloc(sizeof(double)*npoint);

        *_E_II_interpolators = (KernelInterpSingle *)malloc(sizeof(KernelInterpSingle)*nbin*nbin*mbin);

        for(int i=0;i<npoint;i++){
            (*_sep)[i]=sep[i];
            (*_kernel_vals_0)[i]=kernel_vals_0[i];
            (*_kernel_vals_1)[i]=kernel_vals_1[i];
            (*_kernel_vals_2)[i]=kernel_vals_2[i];
            (*_kernel_vals_3)[i]=kernel_vals_3[i];
            (*_kernel_vals_4)[i]=kernel_vals_4[i];
            (*_kernel_vals_5)[i]=kernel_vals_5[i];
            (*_kernel_vals_6)[i]=kernel_vals_6[i];
            (*_kernel_vals_7)[i]=kernel_vals_7[i];
            (*_kernel_vals_8)[i]=kernel_vals_8[i];
            (*_kernel_vals_9)[i]=kernel_vals_9[i];
            (*_kernel_vals_10)[i]=kernel_vals_10[i];
        }

        for(int i=0;i<nbin*nbin*mbin;i++){
            (*_E_II_interpolators)[i] = KernelInterpSingle(E_II_interpolators[i]);
        }

;    }


    void interpolate(){

        interp_kernel_0 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_1 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_2 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_3 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_4 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_5 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_6 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_7 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_8 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_9 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_10 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        gsl_spline_init(interp_kernel_0,sep,kernel_vals_0,npoint);
        gsl_spline_init(interp_kernel_1,sep,kernel_vals_1,npoint);
        gsl_spline_init(interp_kernel_2,sep,kernel_vals_2,npoint);
        gsl_spline_init(interp_kernel_3,sep,kernel_vals_3,npoint);
        gsl_spline_init(interp_kernel_4,sep,kernel_vals_4,npoint);
        gsl_spline_init(interp_kernel_5,sep,kernel_vals_5,npoint);
        gsl_spline_init(interp_kernel_6,sep,kernel_vals_6,npoint);
        gsl_spline_init(interp_kernel_7,sep,kernel_vals_7,npoint);
        gsl_spline_init(interp_kernel_8,sep,kernel_vals_8,npoint);
        gsl_spline_init(interp_kernel_9,sep,kernel_vals_9,npoint);
        gsl_spline_init(interp_kernel_10,sep,kernel_vals_10,npoint);
        sep_a = gsl_interp_accel_alloc();
    }

public:
    KernelInterp(){
        // Empty constructor
    }

    KernelInterp(KernelInterp *kern){
        // Copy constructor

        kern->copy(npoint,nbin, mbin, min_val,&sep,&kernel_vals_0,&kernel_vals_1,&kernel_vals_2,&kernel_vals_3,&kernel_vals_4,&kernel_vals_5,&kernel_vals_6,&kernel_vals_7,&kernel_vals_8,&kernel_vals_9,&kernel_vals_10,&E_II_interpolators);
        interpolate();
    }

    KernelInterp(Parameters *par){
        // Construct interpolation object
        Float tmp_kw,Si_int,tmp_bessel, tmp_sin, tmp_cos, tmp_kernel;

        R0 = par->R0; // truncation radius
        k_max = par->rmax;
        k_min = par->rmin;
        nbin = par->nbin; // number of radial bins
        mbin = par->mbin; // number of Legendre bins
        max_legendre = par->max_l;

        k_high = par->radial_bins_high;
        k_low = par->radial_bins_low;

        out_file = par->out_file;
        out_string = par->out_string;

        // Allocate memory
        sep = (double *)malloc(sizeof(double)*npoint); // k*|r_i-r_j|
        kernel_vals_0 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_1 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_2 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_3 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_4 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_5 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_6 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_7 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_8 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_9 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_10 = (double *)malloc(sizeof(double)*npoint);

        r_vals = (double *)malloc(sizeof(double)*npoint);
        kernel_II_vals = (double *)malloc(sizeof(double)*npoint*nbin*nbin*mbin);


        min_val = 0.1*(R0*k_min)/double(npoint); // minimum interpolation value
        max_val = 2.1*(R0*k_max); // maximum interpolation value

        for(int i=0;i<npoint;i++){
            tmp_kw = min_val + double(i)/double(npoint-1)*(max_val-min_val);
            sep[i] = tmp_kw;
            tmp_sin = sin(tmp_kw);
            tmp_cos = cos(tmp_kw);

            Si_int = gsl_sf_Si(tmp_kw);

            // Compute multipole contributions
            for(int ell=0;ell<=max_legendre;ell+=1){
                tmp_kernel = 3.;
                if(ell==0){
                    tmp_bessel = tmp_sin - tmp_kw*tmp_cos;
                    kernel_vals_0[i] = tmp_bessel*tmp_kernel;
                }
                else if(ell==1){
                    tmp_bessel = -2.*tmp_cos-tmp_kw*tmp_sin;
                    kernel_vals_1[i] = tmp_bessel*tmp_kernel;
                }
                else if(ell==2){
                    tmp_bessel = tmp_kw*tmp_cos - 4*tmp_sin+3*Si_int;
                    kernel_vals_2[i] = tmp_bessel*tmp_kernel;
                }
                else if(ell==3){
                    tmp_bessel = 7.*tmp_cos - 15.*tmp_sin/tmp_kw+tmp_kw*tmp_sin;
                    kernel_vals_3[i] = tmp_bessel*tmp_kernel;
                }
                else if(ell==4){
                    tmp_bessel = 0.5*((105./tmp_kw-2.*tmp_kw)*tmp_cos+(22.-105./pow(tmp_kw,2))*tmp_sin+15*Si_int);
                    kernel_vals_4[i] = tmp_bessel*tmp_kernel;
                }
                else if(ell==5){
                    tmp_bessel = (tmp_kw*(315 - 16*pow(tmp_kw,2))*tmp_cos - (315 - 105 *pow(tmp_kw,2) + pow(tmp_kw,4))*tmp_sin)/pow(tmp_kw,3.);
                    kernel_vals_5[i] = tmp_bessel*tmp_kernel;
                }
                else if(ell==6){
                    tmp_bessel = 1./(8.*pow(tmp_kw,4))*(tmp_kw*(20790. - 1575.*pow(tmp_kw,2) + 8.*pow(tmp_kw,4.))*tmp_cos + (-20790. + 8505.*pow(tmp_kw,2) - 176.*pow(tmp_kw,4.))*tmp_sin + 105.*pow(tmp_kw,4.)*Si_int);
                    kernel_vals_6[i] = tmp_bessel*tmp_kernel;
                }
                else if(ell==7){
                    tmp_bessel = 1./pow(tmp_kw,5.)*(tmp_kw*(27027 - 2772*pow(tmp_kw,2) + 29*pow(tmp_kw,4))*tmp_cos + (-27027 + 11781*pow(tmp_kw,2.) - 378*pow(tmp_kw,4.) + pow(tmp_kw,6.))*tmp_sin);
                    kernel_vals_7[i] = tmp_bessel*tmp_kernel;
                }
                else if(ell==8){
                    tmp_bessel = 1./(16*pow(tmp_kw,6.))*(tmp_kw*(5405400 - 630630*pow(tmp_kw,2) + 10395*pow(tmp_kw,4) -  16.*pow(tmp_kw,6))*tmp_cos + (-5405400 + 2432430*pow(tmp_kw,2) - 100485*pow(tmp_kw,4) + 592*pow(tmp_kw,6))*tmp_sin + 315*pow(tmp_kw,6)*Si_int);
                    kernel_vals_8[i] = tmp_bessel*tmp_kernel;
                }
                else if(ell==9){
                    tmp_bessel = 1./pow(tmp_kw,7.)*(tmp_kw*(4922775 - 617760*pow(tmp_kw,2.) + 12870*pow(tmp_kw,4.) - 46*pow(tmp_kw,6.))*tmp_cos - (4922775 - 2258685*pow(tmp_kw,2.) + 109395*pow(tmp_kw,4.) - 990*pow(tmp_kw,6.) + pow(tmp_kw,8.))*tmp_sin);
                    kernel_vals_9[i] = tmp_bessel*tmp_kernel;
                }
                else if(ell==10){
                    tmp_bessel = 1./(128.*pow(tmp_kw,8.))*(tmp_kw*(10475665200. - 1378377000*pow(tmp_kw,2.) + 34144110*pow(tmp_kw,4.) - 186615*pow(tmp_kw,6.) + 128*pow(tmp_kw,8.))*tmp_cos - 7.*(1496523600 - 695752200*pow(tmp_kw,2.) + 37258650*pow(tmp_kw,4.) - 444015*pow(tmp_kw,6.) + 1024*pow(tmp_kw,8.))*tmp_sin +
  3465.*pow(tmp_kw,8.)*Si_int);
                    kernel_vals_10[i] = tmp_bessel*tmp_kernel;
                }
                else{
                    printf("\nOnly ell = 0 to 10 implemented\n");
                    exit(1);
                }
            }
        }

        // activate interpolator function
        interpolate();

        // Now compute the E_II kernel
        load_E_II_interpolators();
    }

    inline Float pair_weight(Float sep){
        // Compute weight function W(r;R_0)
        if(sep<R0/2.) return 1.;
        else{
            Float x = sep/R0;
            if(sep<3.*R0/4.) return 1.-8.*pow(2.*x-1.,3.)+8.*pow(2.*x-1.,4.);
            else return -64.*pow(x-1.,3.)-128.*pow(x-1.,4.);
        }
    }

    void load_E_II_interpolators(){
        // Load the E_ell^{II,ab}(r;R_0) kernel from file or compute if not yet written.

        // Load output kernel grid if it exists.
        char output_name[1000];
        snprintf(output_name, sizeof output_name, "%s/E_II_kernel_values_n%d_l%d_R0%d.txt",out_file, nbin, (mbin-1),int(R0));
        int ec = read_interpolation_grid(output_name);
        if (ec==1){
            fprintf(stderr,"\nInterpolation grid has not been computed for this choice of k-bins. This will be recomputed.\n");
            compute_E_II_kernel();
        }

        // Now we have the kernel values for each a,b pair of k values and ell
        // Next, we create the interpolators
        printf("Loading interpolators\n");

        // Allocate memory to interpolators
        E_II_interpolators = (KernelInterpSingle *)malloc(sizeof(KernelInterpSingle)*nbin*nbin*mbin);

        Float tmp_kernel[npoint];
        for(int ell=0;ell<=max_legendre;ell++){
            for(int n1=0;n1<nbin;n1++){
                for(int n2=0;n2<nbin;n2++){
                    // copy kernel into local grid
                    for(int pi=0;pi<npoint;pi++) tmp_kernel[pi] = kernel_II_vals[((n1*nbin+n2)*mbin+ell)*npoint+pi];
                    E_II_interpolators[(n1*nbin+n2)*mbin+ell].initialize(R0, r_vals, tmp_kernel, npoint);
                  }
            }
      }
        // now done with this array
        free(kernel_II_vals);
        free(r_vals);
    }

    void compute_E_II_kernel(){
      // Compute the E_ell^{II,ab}(r;R_0) kernel for each radial bin pair and ell
      // First, compute the omega_ell^a functions for each radial bin and ell
      // This writes the interpolation grid to file which can be reread easily.

      printf("Computing E_II_kernel from scratch\n");

      Float dr = R0/(Float(npoint)-1); // step size in r array
      Float dx = 2.*R0/(Float(npoint)-1); // step-size in x array
      Float max_p = 2.01*k_max;
      Float dp = max_p/(Float(npoint)-1); // step-size in p array
      Float old_kr, old_kr3, old_kernel,new_kr,new_kr3,new_kernel,r_value, tmp_r2_W, p_value, tmp_int;
      Float p2, j0_tmp, om_1;
      Float bessel_array[mbin]; // array of j_ell(pr) values

      //TODO: memalign this?
      Float omega_ell_a[nbin*mbin*npoint]; // array to hold omega_ell values
      Float tmp_j_ell_a[nbin*mbin]; // array to hold temporate j_ell^a for single r

      // First initialize omega vector
      for(int i=0;i<nbin*mbin*npoint;i++) omega_ell_a[i]=0;

      int percent=0;
      for(int j=0;j<npoint;j++){

          if (Float(j)/Float(npoint)*100.>=Float(percent)){
              printf("Computing numerical integral 1 of 2: %d%% complete.\n",percent);
              percent+=10;
          }

          // Create real-space array and do ell, k independent pieces first
          r_value = Float(j)*dr+dr/100.; // avoid zero errors
          tmp_r2_W = 4.*M_PI*pow(r_value,2)*pair_weight(r_value)*dr;

          // Compute j_l^a(r) functions
          for(int ell=0;ell<=max_legendre;ell++){
              // Iterate over k-bins (assuming contiguous)
              for(int i=0;i<nbin;i++){
                  if(i==0){
                      old_kr = r_value*k_low[0];
                      old_kr3 = pow(old_kr,3);
                      old_kernel = kernel(ell,old_kr);
                  }
                  new_kr = r_value*k_high[i];
                  new_kr3 = pow(new_kr,3);
                  new_kernel = kernel(ell,new_kr);

                  // compute j_ell^a
                  tmp_j_ell_a[i*mbin+ell] = (new_kernel-old_kernel)/(new_kr3-old_kr3);

                  //Save kernel functions for next time
                  old_kr = new_kr;
                  old_kr3 = new_kr3;
                  old_kernel = new_kernel;
              }
          }

          // Iterate over p array
          for(int n=0;n<npoint;n++){
              p_value = dp*Float(n) + dp/100.; // avoid zero errors

              // Compute j_l(pr) array
              gsl_sf_bessel_jl_array(max_legendre, p_value*r_value, bessel_array);

              // Now fill up array with integral contribution
              for(int ell=0;ell<=max_legendre;ell++){
                  tmp_int = tmp_r2_W*bessel_array[ell];
                  for(int i=0;i<nbin;i++){
                      //TODO: fix horrible indexing here
                      omega_ell_a[(i*mbin+ell)*npoint+n] += tmp_int * tmp_j_ell_a[i*mbin+ell];
                  }
              }
          }
      }

      // NB: omega functions tested against python and are correct
      printf("\n");

      // We now have omega_ell^a(p) functions for all multipoles and k-bins
      // NB: probably don't need to store this for each ell separately - just compute each ell at once

      // TODO: Only store kernel for a <= b since symmetric

      // Zero array
      for(int i=0;i<npoint*(nbin+1)/2*nbin*mbin;i++) kernel_II_vals[i]=0;

      percent = 0;
      // First iterate over p array
      for(int pi=0;pi<npoint;pi++){
          if (Float(pi)/Float(npoint)*100.>=Float(percent)){
              printf("Computing numerical integral 2 of 2: %d%% complete.\n",percent);
              percent+=10;
          }
          p_value = dp*pi + dp/100.;
          p2 = pow(p_value,2)/(2.*M_PI)*dp;

          // Now iterate over real-space array
          for(int xi=0;xi<npoint;xi++){
              if(pi==0) r_vals[xi] = xi*dx + dx/100.;
              j0_tmp = gsl_sf_bessel_j0(p_value*r_vals[xi]); // j0(pr)

              for(int ell=0;ell<=max_legendre;ell++){
                  for(int n1=0;n1<nbin;n1++){
                      om_1 = omega_ell_a[(n1*mbin+ell)*npoint+pi];
                      for(int n2=0;n2<nbin;n2++){
                          // TODO: indexing is horrible again - p should be on outside
                          kernel_II_vals[((n1*nbin+n2)*mbin+ell)*npoint+xi] +=  p2 * j0_tmp * om_1 * omega_ell_a[(n2*mbin+ell)*npoint+pi];
                      }
                  }
              }
          }
      }
      
      // Save output kernel grid
      char output_name[1000];
      snprintf(output_name, sizeof output_name, "%s/E_II_kernel_values_n%d_l%d_R0%d.txt",out_file, nbin, (mbin-1),int(R0));
      FILE *output_file = fopen(output_name,"w");

      printf("Saving interpolation grid to %s\n",output_name);
      for(int i=0;i<npoint;i++){
          for(int n1=0;n1<nbin;n1++){
              for(int n2=0;n2<nbin;n2++){
                  fprintf(output_file,"%le\t%le\t%le\t%le\t%le\t",k_low[n1],k_high[n1],k_low[n2],k_high[n2],r_vals[i]);
                  for(int ell=0;ell<=max_legendre;ell++) fprintf(output_file,"%le\t",kernel_II_vals[((n1*nbin+n2)*mbin+ell)*npoint+i]);
                  fprintf(output_file,"\n");
              }
          }
      }
  }

    Float kernel_E_II(int n1, int n2, int ell, Float r){
        // Compute the E_l^{II,ab}(r) kernel from prebuilt interpolators
        // n1,n2 specify the radial bins
        return E_II_interpolators[(n1*nbin+n2)*mbin+ell].kernel(r);
    }

    ~KernelInterp() {
        //Destructor
        gsl_spline_free(interp_kernel_0);
        gsl_spline_free(interp_kernel_1);
        gsl_spline_free(interp_kernel_2);
        gsl_spline_free(interp_kernel_3);
        gsl_spline_free(interp_kernel_4);
        gsl_spline_free(interp_kernel_5);
        gsl_spline_free(interp_kernel_6);
        gsl_spline_free(interp_kernel_7);
        gsl_spline_free(interp_kernel_8);
        gsl_spline_free(interp_kernel_9);
        gsl_spline_free(interp_kernel_10);
        gsl_interp_accel_free(sep_a);

        free(E_II_interpolators);

        free(sep);
        free(kernel_vals_0);
        free(kernel_vals_1);
        free(kernel_vals_2);
        free(kernel_vals_3);
        free(kernel_vals_4);
        free(kernel_vals_5);
        free(kernel_vals_6);
        free(kernel_vals_7);
        free(kernel_vals_8);
        free(kernel_vals_9);
        free(kernel_vals_10);
    }

    int read_interpolation_grid(char* interp_name){
        // Read the saved interpolation grid
        // Returns 1 if fails, in which case grid will be remade
        char line[100000];

        FILE *fp;
        fp = fopen(interp_name,"r");
        if (fp==NULL){
            fprintf(stderr,"Interpolation grid file %s not found\n",interp_name);
            return 1;
        }
        fprintf(stderr,"\nReading interpolation grid file %s\n",interp_name);

        int line_count=0; // line counter
        int counter=0; // counts which element in line
        int i=0; //counts index of interpolation vector
        // Read in values to file
        while (fgets(line,100000,fp)!=NULL) {
            // Select required lines in file
            if (line[0]=='#') continue;
            if (line[0]=='\n') continue;

            // Split into variables
            char * split_string;
            split_string = strtok(line, "\t");
            counter=0;

            // Define which n1,n2,i,ell bins we are on
            int n2 = line_count%nbin;
            int n1 = (line_count/nbin)%nbin;
            i = (line_count/nbin)/nbin;
            int ell;
            // Iterate over line
            while (split_string!=NULL){
                // Check we have the same radial bins
                if(counter==0) if(atof(split_string)!=k_low[n1]) return 1;
                if(counter==1) if(atof(split_string)!=k_high[n1]) return 1;
                if(counter==2) if(atof(split_string)!=k_low[n2]) return 1;
                if(counter==3) if(atof(split_string)!=k_high[n2]) return 1;
                // Read in radial bin on first iteration only
                if((counter==4)&&(n1==0)&&(n2==0)) r_vals[i] = atof(split_string);
                // Read in angular bins
                if((counter>4)&&(counter<=max_legendre+5)){
                    ell=counter-5;
                    kernel_II_vals[((n1*nbin+n2)*mbin+ell)*npoint+i] = atof(split_string);
                }
                split_string = strtok(NULL,"\t");
                counter++;
            }
            if((i==0)&&(n1==0)&&(n2==0)){
                if(counter<(max_legendre+6)) return 1;
            }
            line_count++;

        }
        if(npoint!=(i+1)) return 1;
        printf("Read in interpolation grid with %d points\n",npoint);
        return 0; // if all went well
    }

};


#endif
