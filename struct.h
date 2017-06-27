
struct elem {
  real q0, qx, qy, qz, q, q2;
  real d[3][3], e[3][4];
  real cs1, cs2, cs3, cb1, cb2, cb3;       // Modulus 
  real s01, s02, s03, b01, b02, b03;       // Reference strains
  real s1, s2, s3, b1, b2, b3;             // Strains
  real f1, f2, f3, t1, t2, t3;             // Force and couple
  real pe;                                 // Pot Energy
  real peT1, peT2, peT3, peR1, peR2, peR3; // pot energies of the elements

}; 					     
					     
struct node {				     
  real d[3][3], e[3][4];		     
  real rx, ry, rz, q0, qx, qy, qz;         // Coordinates
  real px, py, pz, lx, ly, lz;             // Momenta in Space frame 
  real p1, p2, p3, l1, l2, l3;             // Momenta in Body frame
  real f1, f2, f3, t1, t2, t3;             // Force and couple
  real fx, fy, fz, t0, tx, ty, tz;         // Forces in SFF
  real ke; 				     
  real keT1, keT2, keT3, keR1, keR2, keR3; // energies of the DoF */

};

struct rod {
  struct node *nodes;                        // Nodal data structure
  struct elem *elems;                        // Element data structure
  int  n_nodes, n_elems;                     // Number of nodes, elements
  real YI;                                   // Parameter inherited from the main code
  real dia, rho, Y, G, sigma;                // diameter, volume density, Y, G, sigma 
  real i1, i2, i3, xm, xa, le;               // Inertia, mass, area, length
  real cs1, cs2, cs3, cb1, cb2, cb3;         // Modulus 

  real At1, Bt1, Ct1, At2, Bt2, Ct2, At3, Bt3, Ct3; //OU process coefficients - translational
  real Ar1, Br1, Cr1, Ar2, Br2, Cr2, Ar3, Br3, Cr3; //OU process coefficients - rotational
  real zt1, zt2, zt3, zr1, zr2, zr3;                //Friction per unit mass coeff

  real keT1, keT2, keT3, keR1, keR2, keR3;     
  real peT1, peT2, peT3, peR1, peR2, peR3;     
					       
  real xcm, ycm, zcm;			     // Center of mass coordinates
  real pxcm, pycm, pzcm, lxcm, lycm, lzcm;   // Center of mass momenta

};

