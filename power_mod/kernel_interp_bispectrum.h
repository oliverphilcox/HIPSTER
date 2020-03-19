// kernel_interp_bispectrum.h - Interpolator for kernel functions for bispectrum computation
// This stores 3(-1)^ell D_ell(kw)

#ifndef KERNEL_INTERP_H
#define KERNEL_INTERP_H

// BISPECTRUM SPECIFIC CODE
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

class KernelInterp{
    // Compute a interpolation kernel function

private:
    Float R0,k_max,k_min;
    Float *k_high, *k_low; // max and min of each radial bin
    int max_legendre; // maximum order of Legendre polynomials needed
    int nthread;
    char *out_file, *out_string;

public:
    int nbin,mbin;
    int npoint=100000;
    Float min_val,max_val;
    double *sep, *kernel_vals_0,*kernel_vals_1,*kernel_vals_2,*kernel_vals_3,*kernel_vals_4,*kernel_vals_5;
    double *kernel_vals_6,*kernel_vals_7,*kernel_vals_8,*kernel_vals_9,*kernel_vals_10;
    gsl_interp_accel *sep_a;
    gsl_spline *interp_kernel_0, *interp_kernel_1, *interp_kernel_2,*interp_kernel_3;
    gsl_spline *interp_kernel_4, *interp_kernel_5, *interp_kernel_6,*interp_kernel_7;
    gsl_spline *interp_kernel_8, *interp_kernel_9, *interp_kernel_10;

public:
    inline double kernel(int ell, double this_sep){
        // Kernel function accelerated by interpolation
        if(this_sep<=min_val) this_sep=min_val;
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

public:
    void copy_function(KernelInterp *kern){
        // Copy pre-existing kernel function into object
        // NB: This isn't currently used, so may have memory leaks.
        npoint = kern->npoint;
        min_val = kern->min_val;
        max_val = kern->max_val;
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
        }

        // activate interpolator
        interpolate();
    }

private:
    void copy(int& _npoint, int& _nbin, int&  _mbin, Float& _min_val, Float& _max_val, double **_sep, double **_kernel_vals_0,double **_kernel_vals_1,double **_kernel_vals_2,double **_kernel_vals_3,double **_kernel_vals_4,double **_kernel_vals_5,double **_kernel_vals_6,double **_kernel_vals_7,double **_kernel_vals_8,double **_kernel_vals_9,double **_kernel_vals_10){
        _npoint = npoint;
        _min_val = min_val;
        _max_val = max_val;
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
    }

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

        kern->copy(npoint,nbin, mbin, min_val,max_val,&sep,&kernel_vals_0,&kernel_vals_1,&kernel_vals_2,&kernel_vals_3,&kernel_vals_4,&kernel_vals_5,&kernel_vals_6,&kernel_vals_7,&kernel_vals_8,&kernel_vals_9,&kernel_vals_10);
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
        nthread = par->nthread;

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

        // This minimum value should be sufficient for multipoles up to ell = 10
        min_val = 0.1*(R0*k_min)/double(npoint);
        max_val = 2.1*(R0*k_max); // maximum interpolation value

        for(int i=0;i<npoint;i++){
            tmp_kw = min_val + double(i)/double(npoint-1)*(max_val-min_val);
            sep[i] = tmp_kw;
            tmp_sin = sin(tmp_kw);
            tmp_cos = cos(tmp_kw);

            Si_int = gsl_sf_Si(tmp_kw);

            // Insert kw = 0 limit for low k to avoid blow-up
            if(tmp_kw<0.1){
                for(int ell=0;ell<=max_legendre;ell+=1){
                    if(ell==0) kernel_vals_0[i] = 0;
                    else if(ell==1) kernel_vals_1[i] = -6;
                    else if(ell==2) kernel_vals_2[i] = 0;
                    else if(ell==3) kernel_vals_3[i] = -24;
                    else if(ell==4) kernel_vals_4[i] = 0;
                    else if(ell==5) kernel_vals_5[i] = -48;
                    else if(ell==6) kernel_vals_6[i] = 0;
                    else if(ell==7) kernel_vals_7[i] = -384./5.;
                    else if(ell==8) kernel_vals_8[i] = 0;
                    else if(ell==9) kernel_vals_9[i] = -768./7.;
                    else if(ell==10) kernel_vals_10[i] = 0;
                    else{
                        printf("\nOnly ell = 0 to 10 implemented\n");
                        exit(1);
                    }
                }
            }
            else{
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
                        tmp_bessel = 1./(pow(tmp_kw,5.))*(tmp_kw*(27027. - 2772.*pow(tmp_kw,2) + 29.*pow(tmp_kw,4.))*tmp_cos + (-27027. + 11781.*pow(tmp_kw,2.) - 378.*pow(tmp_kw,4.) + pow(tmp_kw,6.))*tmp_sin);
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
        }

        for(int ii=0;ii<100;ii++) printf("2: %.2e %.4e\n",sep[ii],kernel_vals_2[ii]);
        for(int ii=0;ii<100;ii++) printf("9: %.2e %.4e\n",sep[ii],kernel_vals_9[ii]);
        for(int ii=0;ii<100;ii++) printf("10: %.2e %.4e\n",sep[ii],kernel_vals_10[ii]);
        exit(1);

        // activate interpolator function
        interpolate();

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
};


#endif
