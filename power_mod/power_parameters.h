
// parameter function file for grid_power.cpp

#ifndef PARAMETERS_H
#define PARAMETERS_H

class Parameters{

public:
	// Important variables to set!  Here are the defaults:

    //---------- ESSENTIAL PARAMETERS -----------------

    // Name of the first particle field
    char *fname = NULL;
    const char default_fname[500] = "";

    // Name of the radial binning .csv file in k-space
    char *radial_bin_file = NULL;
    const char default_radial_bin_file[500] = "";

    // Output directory
    char *out_file = NULL;
    const char default_out_file[500] = "";

    // Optional File Prefix for output
    char *out_string = NULL;
    const char default_out_string[500] = "";

    // The number of threads to run on
	  int nthread = 10;

    // Whether or not we are using a periodic box
  	bool perbox = false;

    // Maximum Legendre moment (must be even for power spectra)
    int max_l = 4; // max Legendre moment

    // Kernel truncation radius (in Mpc/h)
    Float R0 = 100;

#ifndef BISPECTRUM
    //---------- POWER SPECTRUM PARAMETERS -----------------

    // Name of second particle field (power spectrum only)
    char *fname2 = NULL;
    const char default_fname2[500] = "";

    // Survey correction function coefficient file
    char *inv_phi_file = NULL;
    const char default_inv_phi_file[500] = "";

#else
  //---------- BISPECTRUM PARAMETERS -----------------

    // Ratio of randoms to data particles (needed for DDRII count)
    Float f_rand = 3;

#endif

    //-------- OTHER PARAMETERS ----------------------------------------------

    // The grid size, which is used to allocate particles internally.
    int nside = 50;

  	// The maximum number of points to read
  	uint64 nmax = 1000000000000;

  	// Whether to balance the weights or multiply them by -1
  	int qinvert = 0, qbalance = 0;

	  //---------------- INTERNAL PARAMETERS -----------------------------------
    // (no more user defined parameters below this line)

    // Number of particles
    // (will be overwritten if reading from a file)
    int np = -1;

    // Summed weights (only used if periodic)
#ifdef PERIODIC
    Float sum_w1,sum_w2;
#endif
  	// The periodicity of the position-space cube.
  	Float boxsize = 1000; // this is only used if the input particles are made randomly

  	// The particles will be read from the unit cube, but then scaled by boxsize.
  	Float rescale = 1.;   // If left zero or negative, set rescale=boxsize

#ifndef BISPECTRUM
    // For consistency with other modules
    char *inv_phi_file2 = NULL; // Survey correction function coefficient file
    char *inv_phi_file12 = NULL; // Survey correction function coefficient file
#endif

	   // The periodicity of the position-space cuboid in 3D.
    Float3 rect_boxsize = {boxsize,boxsize,boxsize}; // this is overwritten on particle read-in

		Float cellsize;

    // Radial binning parameters (will be set from file)
    int nbin=0,mbin;
    Float rmin, rmax;
    Float * radial_bins_low;
    Float * radial_bins_high;

    // Variable to decide if we are using multiple tracers:
    bool multi_tracers;

    // Constructor
	  Parameters(int argc, char *argv[]){

	    if (argc==1) usage();
	    int i=1;
	    while (i<argc) {
      if (!strcmp(argv[i],"-fname")) fname = argv[++i];
      else if (!strcmp(argv[i],"-rescale")) rescale = atof(argv[++i]);
	    else if (!strcmp(argv[i],"-R0")) R0 = atof(argv[++i]);
      else if (!strcmp(argv[i],"-nside")) nside = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-in")) fname = argv[++i];
#ifndef BISPECTRUM
      else if (!strcmp(argv[i],"-in2")) fname2 = argv[++i];
      else if (!strcmp(argv[i],"-inv_phi_file")) inv_phi_file=argv[++i];
#else
      else if (!strcmp(argv[i],"-f_rand")) f_rand=atof(argv[++i]);
#endif
		else if (!strcmp(argv[i],"-nmax")) nmax = atoll(argv[++i]);
		else if (!strcmp(argv[i],"-nthread")) nthread = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-balance")) qbalance = 1;
		else if (!strcmp(argv[i],"-invert")) qinvert = 1;
    else if (!strcmp(argv[i],"-output")) out_file = argv[++i];
    else if (!strcmp(argv[i],"-out_string")) out_string = argv[++i];
    else if (!strcmp(argv[i],"-binfile")) radial_bin_file=argv[++i];
    else if (!strcmp(argv[i],"-max_l")) max_l=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-perbox")) perbox = 1;
		else if (!strcmp(argv[i],"-def")) { fname = NULL; }
		else {
		    fprintf(stderr, "Don't recognize %s\n", argv[i]);
		    usage();
		}
		i++;
	    }
#ifdef BISPECTRUM
#ifndef PERIODIC
      fprintf(stderr,"\nCode was compiled in BISPECTRUM mode without PERIODIC mode enabled. Support for this is not yet available. Exiting.\n\n");
      exit(1);
#endif
#endif

#ifdef PERIODIC
        if (perbox!=true){
            printf("\nC++ code compiled with periodic flag, but periodic box parameter is not set! Exiting.\n\n");
            exit(1);
        }
#else
        if (perbox==true){
            printf("\nC++ code not compiled with periodic flag, but periodic box parameter is set! Exiting.\n\n");
            exit(1);
        }
#endif

        if(R0<1){
            printf("\nTruncation radius (%.0f Mpc/h) is too small for accurate spectral computation. Exiting.\n\n",R0);
            exit(1);
        }
        if(R0>400){
            printf("\nTruncation radius (%.0f Mpc/h) is too large for efficient spectral computation. Exiting.\n\n",R0);
            exit(1);
        }
#ifdef BISPECTRUM
        if(f_rand>30){
            printf("\nYou have chosen a huge random catalog of %.2f x more particles than the data. This is too big for fast computation! Exiting.\n\n",f_rand);
            exit(1);
        }
        if(f_rand<1){
            printf("\nThe random catalog should be at least as large as the data. You've selected f_rand = %.2f. Exiting.\n\n",f_rand);
            exit(1);
        }
#endif
	    // compute smallest and largest boxsizes
	    Float box_min = fmin(fmin(rect_boxsize.x,rect_boxsize.y),rect_boxsize.z);
	    Float box_max = fmax(fmax(rect_boxsize.x,rect_boxsize.y),rect_boxsize.z);

	    assert(i==argc);  // For example, we might have omitted the last argument, causing disaster.

#ifndef BISPECTRUM
        assert(max_l%2==0); // check maximum ell is even
        assert(max_l<=6); // ell>6 not yet implemented!
        mbin = max_l/2+1; // number of angular bins is set to number of Legendre bins
        if (fname2==NULL) fname2 = (char *) default_fname2;   // No name was given
        if (inv_phi_file==NULL) {inv_phi_file = (char *) default_inv_phi_file;} // no phi file specified
#ifdef PERIODIC
        fprintf(stderr,"Survey correction function %s specified, but in PERIODIC mode. Survey correction function will be ignored.",inv_phi_file);
#endif
#else
        assert(max_l<=10); // ell>10 not yet implemented!
        mbin = max_l+1;
#endif
        if (rescale<=0.0) rescale = box_max;   // This would allow a unit cube to fill the periodic volume
	    if (out_file==NULL) out_file = (char *) default_out_file; // no output savefile
	    if (out_string==NULL) out_string = (char *) default_out_string; // no output string
	    if (radial_bin_file==NULL) {radial_bin_file = (char *) default_radial_bin_file;} // No radial binning

	    if (fname==NULL) fname = (char *) default_fname;   // No name was given
#ifndef BISPECTRUM
    if (!strcmp(fname2,"")) fname2 = fname; // read from first file if second not given
#endif
	    create_directory();

	    // Read in the radial binning
	    read_radial_binning(radial_bin_file);
        printf("Read in %d radial k-space bins in range (%.0f, %.0f) successfully.\n\n",nbin,rmin,rmax);

	    assert(box_min>0.0);
	    assert(rmax>0.0);
	    assert(nside>0);

#ifdef OPENMP
		omp_set_num_threads(nthread);
#else
		nthread=1;
#endif

		// Output for posterity
		printf("Box Size = {%6.5e,%6.5e,%6.5e}\n", rect_boxsize.x,rect_boxsize.y,rect_boxsize.z);
		printf("Grid = %d\n", nside);
		printf("Radial Bins = %d\n", nbin);
		printf("Radial k-space binning = {%6.5f, %6.5f} over %d bins (user-defined bin widths) \n",rmin,rmax,nbin);
		printf("Output directory: '%s'\n\n",out_file);

	}
private:
	void usage() {
	      fprintf(stderr, "\nUSAGE FOR HIPSTER\n\n");
        fprintf(stderr, "Key Parameters:\n");
        fprintf(stderr, "   -def: This allows one to accept the defaults without giving other entries.\n");
	      fprintf(stderr, "   -in <file>: The input file for particle-set 1 (space-separated x,y,z,[w]).\n");
        fprintf(stderr, "   -binfile <filename>: File containing the desired k-space radial bins\n");
        fprintf(stderr, "   -output: Directory to save output spectra into\n");
        fprintf(stderr, "   -out_string: (Optional) String to add to file name for identification\n");
        fprintf(stderr, "   -nthread <nthread>: The number of CPU threads ot use for parallelization.\n");
        fprintf(stderr, "   -perbox <perbox>: Boolean, whether the box is periodic is not\n");
        fprintf(stderr, "   -max_l <max_l>: Maximum legendre multipole (must be even for power spectra)\n");
        fprintf(stderr, "   -R0 <R0>: Truncation radius for pair-separation window (in Mpc/h)\n");
#ifndef BISPECTRUM
        fprintf(stderr, "   -in2 <file>: The input file for particle-set 2 (space-separated x,y,z,[w]).\n");
        fprintf(stderr, "   -inv_phi_file <filename>: Survey inverse correction function multipole coefficient file\n");
#else
        fprintf(stderr, "   -f_rand <f_rand>: Ratio of random particles to galaxies. Typically this should be order a few.\n");
#endif
        fprintf(stderr, "\nOther Parameters:\n");
        fprintf(stderr, "   -nside <nside>: The grid size for accelerating the pair count.  Default 250.\n");
	      fprintf(stderr, "          There are {nside} cells along the longest dimension of the periodic box.\n");
	      fprintf(stderr, "   -nmax <nmax>: The maximum number of particles to read in from the particle files. Default 1000000000000\n");
  	    fprintf(stderr, "   -invert: Multiply all the weights by -1.\n");
  	    fprintf(stderr, "   -balance: Rescale the negative weights so that the total weight is zero.\n");
        fprintf(stderr, "   -rs <rstart>:  If inverting particle weights, this sets the index from which to start weight inversion. Default 0\n");
  	    fprintf(stderr, "\n");
  	    fprintf(stderr, "\n");
	    exit(1);
	}

	void create_directory(){
        // Initialize output directory:
	    if (mkdir(out_file,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)==0){
            printf("\nCreating output directory\n");
        }
    }

    void read_radial_binning(char* binfile_name){
        // Read the radial binning file and determine the number of bins
        char line[100000];

        FILE *fp;
        fp = fopen(binfile_name,"r");
        if (fp==NULL){
            fprintf(stderr,"Radial binning file %s not found\n",binfile_name);
            abort();
        }
        fprintf(stderr,"\nReading radial binning file '%s'\n",binfile_name);

        // Count lines to construct the correct size
        while (fgets(line,10000,fp)!=NULL){
            if (line[0]=='#') continue; // comment line
            if (line[0]=='\n') continue;
                nbin++;
            }
            printf("# Found %d radial bins in the file\n",nbin);
            rewind(fp); // restart file

            // Now allocate memory to the weights array
            int ec=0;
            ec+=posix_memalign((void **) &radial_bins_low, PAGE, sizeof(Float)*nbin);
            ec+=posix_memalign((void **) &radial_bins_high, PAGE, sizeof(Float)*nbin);
            assert(ec==0);

            int line_count=0; // line counter
            int counter=0; // counts which element in line

            // Read in values to file
            while (fgets(line,100000,fp)!=NULL) {
                // Select required lines in file
                if (line[0]=='#') continue;
                if (line[0]=='\n') continue;

                // Split into variables
                char * split_string;
                split_string = strtok(line, "\t");
                counter=0;

                // Iterate over line
                while (split_string!=NULL){
                    if(counter==0){
                        radial_bins_low[line_count]=atof(split_string);
                        }
                    if(counter==1){
                        radial_bins_high[line_count]=atof(split_string);
                        }
                    if(counter>1){
                        fprintf(stderr,"Incorrect file format");
                        abort();
                    }
                    split_string = strtok(NULL,"\t");
                    counter++;
                }
                line_count++;
            }

            rmin = radial_bins_low[0];
            rmax = radial_bins_high[line_count-1];
            assert(line_count==nbin);
    }
};
#endif
