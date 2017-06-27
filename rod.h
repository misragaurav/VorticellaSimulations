/*  Element properties for a homogenous rod  */

//#define real_float
#define real double
#define Pi   acos(-1.0)
#define q_normalization
#define exclude_ends

/* Properties of undeformed rod */
#define GAMMA0_1  0.0
#define GAMMA0_2  0.0
#define GAMMA0_3  1.0

#define OMEGA0_1  0.0
#define OMEGA0_2  0.0
#define OMEGA0_3  0.0

#define DIA      4.0    // unit - 1 micron
#define DENSITY  0.001  // volume density
#define AREA     Pi*(DIA/2.0)*(DIA/2.0)
#define MASS     DENSITY*AREA  //mass  per unit length
#define INERTIA  (1.0/4.0)*AREA*(DIA/2.0)*(DIA/2.0) //moment of area

/* Initial configuration */

#define GAMMA_1   0.0
#define GAMMA_2   0.0
#define GAMMA_3   1.0

#define OMEGA_1   0.0
#define OMEGA_2   0.0 //(2.0*Pi)/(128*1) 
#define OMEGA_3   0.0
