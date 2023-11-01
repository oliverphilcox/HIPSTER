// driver.h - this contains various c++ functions to create particles in random positions / read them in from file. Based on code by Alex Wiegand.
#include "cell_utilities.h"

#ifndef DRIVER_H
#define DRIVER_H

// ====================  The Driver ===========================

Particle *make_particles(Float3 rect_boxsize, int np) {
    // Make np random particles
    srand48(time(NULL));      // set to a fixed value for reproducibility
    Particle *p = (Particle *)malloc(sizeof(Particle)*np);
    for (int j=0; j<np; j++) {
        p[j].pos.x = drand48()*rect_boxsize.x;
        p[j].pos.y = drand48()*rect_boxsize.y;
        p[j].pos.z = drand48()*rect_boxsize.z;
        p[j].w = 1.0;
    }
    printf("# Done making %d random particles, periodically distributed.\n", np);
    return p;
}

Particle *read_particles(Float rescale, int *np, const char *filename, uint64 nmax) {
    // This will read particles from a file, space-separated x,y,z,w,JK for weight w, (jackknife region JK)
    // Particle positions will be rescaled by the variable 'rescale'.
    // For example, if rescale==boxsize, then inputting the unit cube will cover the periodic volume
    char line[1000];
    int j=0,n=0;
#ifdef LYA
    Float sum_weight;
#endif
    FILE *fp;
    int stat;
    double tmp[5];

    fp = fopen(filename, "r");
    if (fp==NULL) {
        fprintf(stderr,"File %s not found\n", filename); abort();
    }

    // Count lines to construct the correct size
    while (fgets(line,1000,fp)!=NULL&&(uint)n<nmax) {
        if (line[0]=='#') continue;
        if (line[0]=='\n') continue;
        n++;
    }
    rewind(fp);
    *np = n;

    Particle *p = (Particle *)malloc(sizeof(Particle)*n);
    printf("# Found %d particles from %s\n", n, filename);
    printf("# Rescaling input positions by factor %f\n", rescale);

    while (fgets(line,1000,fp)!=NULL&&j<n) {
        if (line[0]=='#') continue;
        if (line[0]=='\n') continue;
        stat=sscanf(line, "%lf %lf %lf %lf %lf", tmp, tmp+1, tmp+2, tmp+3, tmp+4);

        if (stat<3) {
        	fprintf(stderr,"Particle %d has bad format\n", j); // Not enough coordinates
        	abort();
        }

        p[j].pos.x = tmp[0]*rescale;
        p[j].pos.y = tmp[1]*rescale;
        p[j].pos.z = tmp[2]*rescale;
        // Get the weights from line 4 if present, else fill with +1
        if((stat!=4)&&(stat!=5)) p[j].w = 1.;
        else p[j].w = tmp[3];
#ifdef LYA
        p[j].tid = tmp[4];
        sum_weight += p[j].w;
#endif
        j++;
    }
    fclose(fp);
    printf("# Done reading the particles\n");
    
#ifdef LYA
    // Check if the mean weight is close to unity
    if ((abs(sum_weight)/j)>0.01) {
        fprintf(stderr, "Mean Lya weight should be approximately zero; exiting\n");
        abort();
    }
#endif

    return p;
}




bool compute_bounding_box(Particle *p, int np, Float3 &rect_boxsize, Float &cellsize, Float rmax, Float3& pmin, int nside) {
    // Compute the boxsize of the bounding cuboid box and determine whether we are periodic
    Float3 pmax;
    bool box=false;
    pmin.x = pmin.y = pmin.z = 1e30;
    pmax.x = pmax.y = pmax.z = -1e30;
    for (int j=0; j<np; j++) {
        pmin.x = fmin(pmin.x, p[j].pos.x);
        pmin.y = fmin(pmin.y, p[j].pos.y);
        pmin.z = fmin(pmin.z, p[j].pos.z);
        pmax.x = fmax(pmax.x, p[j].pos.x);
        pmax.y = fmax(pmax.y, p[j].pos.y);
        pmax.z = fmax(pmax.z, p[j].pos.z);
    }
    printf("# Range of x positions are %6.2f to %6.2f\n", pmin.x, pmax.x);
    printf("# Range of y positions are %6.2f to %6.2f\n", pmin.y, pmax.y);
    printf("# Range of z positions are %6.2f to %6.2f\n", pmin.z, pmax.z);
    Float3 prange = pmax-pmin;
    Float         biggest = prange.x;
    biggest = fmax(biggest, prange.y);
    biggest = fmax(biggest, prange.z);
    if (prange.x>0.95*biggest && prange.y>0.95*biggest && prange.z>0.95*biggest) {
        // Probably using a cube of inputs, intended for a periodic box
    	box=true;
#ifndef PERIODIC
    	// fprintf(stderr,"#\n# WARNING: cubic input detected but you have not compiled with PERIODIC flag!\n#\n");
    	// printf("#\n# WARNING: cubic input detected but you have not compiled with PERIODIC flag!\n#\n");
    	// if(biggest+2*rmax>rect_boxsize.x){
    	// printf("#\n# WARNING: Box periodicity is too small; particles will overlap on periodic wrapping!");
    	// fprintf(stderr,"#\n# WARNING: Box periodicity is too small; particles will overlap on periodic wrapping!");
    	// }
    	// biggest = rect_boxsize.x; // just use the input boxsize
        // set max_boxsize to just enclose the biggest dimension plus r_max
        // NB: We natively wrap the grid (to allow for any position of the center of the grid)
        // Must add rmax to biggest to ensure there is no periodic overlap in this case.
        box=false;
        Float max_boxsize = 1.05*(biggest+2*rmax);
        cellsize = max_boxsize/nside; // compute the width of each cell
        // Update ranges to add rmax
        prange.x=1.05*(prange.x+2*rmax);
        prange.y=1.05*(prange.y+2*rmax);
        prange.z=1.05*(prange.z+2*rmax);
        printf("max: %.3f, rmax: %.3f prangex: %.3f\n",max_boxsize,rmax,prange.x);
        // Now compute the size of the box in every dimension
        rect_boxsize = ceil3(prange/cellsize)*cellsize; // to ensure we fit an integer number of cells in each direction
        printf("# Setting non-periodic box-size to {%6.2f,%6.2f,%6.2f}\n", rect_boxsize.x,rect_boxsize.y,rect_boxsize.z);
#else
        // Set boxsize to be the biggest dimension which allows for periodic overlap
        rect_boxsize= {biggest,biggest,biggest};
        cellsize = biggest/nside;
        printf("# Setting periodic box-size to %6.2f\n", biggest);
#endif
    } else {
        // Probably a non-periodic input (e.g. a real dataset)
        box=false;
#ifdef PERIODIC
    	fprintf(stderr,"#\n# WARNING: non-cubic input detected but you have compiled with PERIODIC flag!\n#\n");
    	printf("#\n# WARNING: non-cubic input detected but you have compiled with PERIODIC flag!\n#\n");
#endif
        // set max_boxsize to just enclose the biggest dimension plus r_max
        // NB: We natively wrap the grid (to allow for any position of the center of the grid)
        // Must add rmax to biggest to ensure there is no periodic overlap in this case.
        Float max_boxsize = 1.05*(biggest+2*rmax);
        cellsize = max_boxsize/nside; // compute the width of each cell
        // Update ranges to add rmax
        prange.x=1.05*(prange.x+2*rmax);
        prange.y=1.05*(prange.y+2*rmax);
        prange.z=1.05*(prange.z+2*rmax);
        // Now compute the size of the box in every dimension
        rect_boxsize = ceil3(prange/cellsize)*cellsize; // to ensure we fit an integer number of cells in each direction
        printf("# Setting non-periodic box-size to {%6.2f,%6.2f,%6.2f}\n", rect_boxsize.x,rect_boxsize.y,rect_boxsize.z);
    }

    return box;
}


void invert_weights(Particle *p, int np) {
    for (int j=0; j<np; j++) p[j].w *= -1.0;
    printf("# Multiplying all weights by -1\n");
}

void balance_weights(Particle *p, int np) {
    Float sumpos = 0.0, sumneg = 0.0;
    for (int j=0; j<np; j++)
	if (p[j].w>=0.0) sumpos += p[j].w;
	    else sumneg += p[j].w;
    if (sumneg==0.0 || sumpos==0.0) {
	fprintf(stderr,"Asked to rebalance weights, but there are not both positive and negative weights\n");
	abort();
    }
    Float rescale = sumpos/(-sumneg);
    printf("# Rescaling negative weights by %f\n", rescale);
    for (int j=0; j<np; j++)
	if (p[j].w<0.0) p[j].w *= rescale;
    return;
}

#endif
