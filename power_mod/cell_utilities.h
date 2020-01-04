// cell_utilities.h - This contains many cell-related functions including the cell and particle classes. Based on code by Alex Wiegand.

#ifndef CELL_UTILITIES_H
#define CELL_UTILITIES_H

// We need a vector floor3 function
Float3 floor3(float3 p) {
    return Float3(floor(p.x), floor(p.y), floor(p.z));
}

// we need a vector ceil3 function
Float3 ceil3(float3 p) {
    return Float3(ceil(p.x), ceil(p.y), ceil(p.z));
}

// =================== Particles ====================
// This is the info about each particle that we load in and store in the Grid.

class Particle {
  public:
    Float3 pos;
    Float w;  // The weight for each particle
    Float JK; // The Jackknife region ID for each particle (stored as float for ease)
    int rand_class; // Integer 0 or 1 which defines which set of random particles this is in.
    // This is set at read-in at random and used for the EE computation to avoid diagonal non-cancellation.
};


// ====================  The Cell and Grid classes ==================

/* The Grid class holds a new copy of the particles.
These are sorted into cells, and all positions are referenced to the
cell center.  That way, we can handle periodic wrapping transparently,
simply by the cell indexing.

For simplicity, we opt to flatten the index of the cells into a 1-d number.
For example, this makes multi-threading over cells simpler.
*/

class Cell {
  public:
    int start;	// The starting index of the particle list
    int np;
    int np1; // number of particles in cell in random-partition 1
    int np2; // number of particles in cell in random-partition 2
};


#endif

// LEGENDRE MULTIPOLES
// legendre[0]=1.;
// if(max_legendre>0){
//   legendre[1]=mu_ijk;
//   if(max_legendre>1){
//       Float mu2 = mu_ijk*mu_ijk;
//       legendre[2]=0.5*(3.*mu2-1.);
//       if(max_legendre>2){
//         Float mu3 = mu2*mu_ijk;
//         legendre[3] = 0.5*(-3.*mu_ijk+5.*mu3);
//         if(max_legendre>3){
//           Float mu4 = mu3*mu_ijk;
//           legendre[4]=1./8.*(35.*mu4-30.*mu2+3.);
//           if(max_legendre>4){
//             Float mu5 = mu4*mu_ijk;
//             legendre[5]=1./8.*(15.*mu_ijk-70.*mu3+63*mu5);
//             if(max_legendre>5){
//               Float mu6 = mu5*mu_ijk;
//               legendre[6] = 1./16.*(231.*mu6-315.*mu4+105.*mu2-5.);
//               if(max_legendre>6){
//                 Float mu7 = mu6*mu_ijk;
//                 legendre[7] = 1./16.*(-35.*mu_ijk+315*mu3-693.*mu5+429*mu7);
//                 if(max_legendre>7){
//                   Float mu8 = mu7*mu_ijk;
//                   legendre[8] = 1./128.*(6435.*mu8-12012.*mu6+6930.*mu4-1260.*mu2+35.);
//                   if(max_legendre>8){
//                     Float mu9 = mu8*mu_ijk;
//                     legendre[9] = 1./128.*(315.*mu_ijk-4620*mu3+18018.*mu5-25740.*mu7+12155.*mu9);
//                     if(max_legendre>9){
//                       legendre[5] = 1./256.*(46189.*mu8*mu2-109395.*mu8+90090.*mu6-30030.*mu4+3465.*mu2-63.);
//                     }
//                   }
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }
//
