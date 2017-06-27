#include "struct.h"

// Functions needed across files
void init (struct rod *rod, int n_rods, char *datadir);        // function in init.c file for writing the config file
void pop(struct rod *rod, int n_rods, char *origdir);          // function in init.c file for populating the structure from config or checkpoint
void propagator(struct rod *rod, real dt, int n_rods);         // function in propagatorSympl.c for doing symplectic update of q
void force(struct rod *rod);                        // function in forces.c file
void normal_modes(struct rod *rod, real dt, int t, int tprint, int n_rods, char *datadir); // function in misc.c file
void file_name (char *name, char *work_dir, int task_number);                              // function in misc.c file
void com(struct rod *rod);                                                                 // function in misc.c file
void bc(struct rod *rod, char var);                                                        // function in misc.c file
void energy(struct rod *rod);
void friction(struct rod *rod, real dt, real factor);
void frictionSphere(struct rod *rod, real dt, real factor);
void rotmatrix_node(struct rod *rod);
void rotmatrix_elem(struct rod *rod);
void tip(struct rod *rod);
void Qinv(real Q[][8]);
void normalize(struct rod *rod);
void popspas(struct rod *rod);

// Variables needed across files

char datadir[64] ;
int N ;
int M0 ;
real LENGTH ;
real AR ;
real  dt ;
int t, nt, tprint ;
