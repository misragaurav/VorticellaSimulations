// The initialization file - it sets up the rod's coordinates and momenta according to the given strains and initial momenta
#include "proto.h"

#define X  	0.0
#define Y  	0.0
#define Z  	0.0
#define TH  	0.0
#define PH 	0.0
#define PS 	0.0

#define Y     0.001
#define G     Y/(2+2*1.0/3.0)

void init (struct rod *rod, int n_rods, char *datadir)  /* Initialize rods */
{
  static real i1, i2, i3, xm, xa, le;
  static real s01, s02, s03, b01, b02, b03;
  static real s1, s2, s3, b1, b2, b3;
  static real f0x, f0y, f0z, m0x, m0y, m0z;
  static real fLx, fLy, fLz, mLx, mLy, mLz;
  static int    v01, v02, v03, w01, w02, w03;
  static int    vL1, vL2, vL3, wL1, wL2, wL3;
  struct elem   *pe;
  struct node   *pn;

  real R[3][3], Q[4][8], Qplus[4][4];
  real rx, ry, rz, q0, qx, qy, qz, q;
  real dq0, dqx, dqy, dqz;
  int  n_nodes, n_elems, n, m, k, l, m0;
  FILE *filep=0;
  char file1[64];

   xm = MASS;  xa = AREA;  le = LENGTH;
   i1=INERTIA*DENSITY; i2=INERTIA*DENSITY; i3=(2*INERTIA)*DENSITY; // Moments of inertia per unit length

   //modulii
   rod->cs1 = xa*G;        rod->cs2 = xa*G;        rod->cs3 = xa*Y;
   rod->cb1 = INERTIA*Y;   rod->cb2 = INERTIA*Y;   rod->cb3 = (2*INERTIA)*G;

   //strains
  s01=GAMMA0_1;  s02=GAMMA0_2;  s03=GAMMA0_3; b01=OMEGA0_1;  b02=OMEGA0_2;  b03=OMEGA0_3;
  s1=GAMMA_1;  s2=GAMMA_2;  s3=GAMMA_3;  b1=OMEGA_1;  b2=OMEGA_2;  b3=OMEGA_3;

  f0x=F0x;  f0y=F0y;  f0z=F0_z;  m0x=M0x;  m0y=M0y;  m0z=M0z;
  fLx=FLx;  fLy=FLy;  fLz=FLz;   mLx=MLx;  mLy=MLy;  mLz=MLz;

  n_nodes = N+2;      n_elems = N+1;       m0 = M0;

  /* Material properties */

  for (n = 0; n < n_rods; n++)
    {
      rod[n].n_nodes = n_nodes;                                  /* Number of nodes/elements */
      rod[n].n_elems = n_elems;

      rod[n].v01 = v01;  rod[n].v02 = v02;  rod[n].v03 = v03;  /* Set boundary conditions */
      rod[n].w01 = w01;  rod[n].w02 = w02;  rod[n].w03 = w03;
      rod[n].vL1 = vL1;  rod[n].vL2 = vL2;  rod[n].vL3 = vL3;
      rod[n].wL1 = wL1;  rod[n].wL2 = wL2;  rod[n].wL3 = wL3;

      rod[n].f0x = f0x;  rod[n].f0y = f0y;  rod[n].f0z = f0z;
      rod[n].m0x = m0x;  rod[n].m0y = m0y;  rod[n].m0z = m0z;
      rod[n].fLx = fLx;  rod[n].fLy = fLy;  rod[n].fLz = fLz;
      rod[n].mLx = mLx;  rod[n].mLy = mLy;  rod[n].mLz = mLz;

      rod[n].xm  = xm;   rod[n].xa  = xa;   rod[n].le  = le;
      rod[n].i1  = i1;   rod[n].i2  = i2;   rod[n].i3  = i3;

      rod[n].nodes = (struct node *) calloc(rod[n].n_nodes, sizeof(struct node));
      rod[n].elems = (struct elem *) calloc(rod[n].n_elems, sizeof(struct elem));

      for (m = 0; m < rod[n].n_elems; m++)                       /* Set element variables */
	{
	  pe = &rod[n].elems[m];

	  pe->cs1 = rod->cs1;  pe->cs2 = rod->cs2;  pe->cs3 = rod->cs3;
	  pe->cb1 = rod->cb1;  pe->cb2 = rod->cb2;  pe->cb3 = rod->cb3;
	  
	  pe->s01 = s01;  pe->s02 = s02;  pe->s03 = s03;
	  pe->b01 = b01;  pe->b02 = b02;  pe->b03 = b03;
	  pe->s1  = s1;   pe->s2  = s2;   pe->s3  = s3;
	  pe->b1  = b1;   pe->b2  = b2;   pe->b3  = b3;
	}
      
      /* Nodal coordinates */

      rx = X;  ry = Y;  rz = Z;                                   /* 1st node (node m0) */
      q0 = cos(TH/2.0)*cos((PS+PH)/2.0); qx = sin(TH/2.0)*cos((PH-PS)/2.0);
      qy = sin(TH/2.0)*sin((PH-PS)/2.0); qz = cos(TH/2.0)*sin((PS+PH)/2.0);

      if(m0==0)
	m0=m0+1;  // To avoid the ghost node for m0=0
      
      if(m0==N+1)
	m0=m0-1;  // To avoid the ghost node for m0=N+1
      
      pn = &rod[n].nodes[m0];
      pn->rx = rx;  pn->ry = ry;  pn->rz = rz;
      pn->q0 = q0;  pn->qx = qx;  pn->qy = qy; pn->qz = qz;
      
      pn->px = 0.0; pn->py = 0.0; pn->pz = 0.0;
      pn->l1 = 0.0; pn->l2 = 0.0; pn->l3 = 0.0;

      for (m = m0-1; m >= 1; m--)                                 /* Calculate quaternions */
	{
	  pe = &rod[n].elems[m];
	  pn = &rod[n].nodes[m];
	  pn->px = 0.0; pn->py = 0.0; pn->pz = 0.0;
	  pn->l1 = 0.0; pn->l2 = 0.0; pn->l3 = 0.0;

	  Q[0][0] =  0.0;     Q[0][1] = -pe->b1;  Q[0][2] = -pe->b2;  Q[0][3] = -pe->b3;
	  Q[1][0] =  pe->b1;  Q[1][1] =  0.0;     Q[1][2] =  pe->b3;  Q[1][3] =  -pe->b2;
	  Q[2][0] =  pe->b2;  Q[2][1] = -pe->b3;  Q[2][2] =  0.0;     Q[2][3] = -pe->b1;
	  Q[3][0] =  pe->b3;  Q[3][1] = pe->b2;  Q[3][2] =  -pe->b1;  Q[3][3] =  0.0;

	  for (k = 0; k < 4; k++)
	    for (l = 0; l < 4; l++)
	      {
		Q[k][l] *= -0.25*rod[n].le;
		Qplus[k][l] = Q[k][l];
	      }

	  for (k = 0; k < 4; k++)
	    for (l = 0; l < 4; l++)
	      Q[k][l] *= -1.0;
	  for (k = 0; k < 4; k++)
	    {
	      Q[k][k] += 1.0;
	      Qplus[k][k] += 1.0;
	    }
	  Qinv (Q);

	  dq0 = Qplus[0][0]*q0 + Qplus[0][1]*qx + Qplus[0][2]*qy + Qplus[0][3]*qz;
	  dqx = Qplus[1][0]*q0 + Qplus[1][1]*qx + Qplus[1][2]*qy + Qplus[1][3]*qz;
	  dqy = Qplus[2][0]*q0 + Qplus[2][1]*qx + Qplus[2][2]*qy + Qplus[2][3]*qz;
	  dqz = Qplus[3][0]*q0 + Qplus[3][1]*qx + Qplus[3][2]*qy + Qplus[3][3]*qz;	  

	  q0 =  Q[0][4]*dq0 + Q[0][5]*dqx + Q[0][6]*dqy + Q[0][7]*dqz;
	  qx =  Q[1][4]*dq0 + Q[1][5]*dqx + Q[1][6]*dqy + Q[1][7]*dqz;
	  qy =  Q[2][4]*dq0 + Q[2][5]*dqx + Q[2][6]*dqy + Q[2][7]*dqz;
	  qz =  Q[3][4]*dq0 + Q[3][5]*dqx + Q[3][6]*dqy + Q[3][7]*dqz;

	  q  = sqrt( q0*q0 + qx*qx + qy*qy + qz*qz);
	  q0 /= q; qx /= q; qy /= q;  qz /= q; 
	  pn->q0 = q0; pn->qx = qx; pn->qy = qy; pn->qz = qz; 
	}

      q0 = cos(TH/2.0)*cos((PS+PH)/2.0); qx = sin(TH/2.0)*cos((PH-PS)/2.0);
      qy = sin(TH/2.0)*sin((PH-PS)/2.0); qz = cos(TH/2.0)*sin((PS+PH)/2.0);

      for (m = m0+1; m < rod[n].n_nodes-1; m++)                      /* Calculate quaternions */
	{
	  pe = &rod[n].elems[m];
	  pn = &rod[n].nodes[m];
	  pn->px = 0.0; pn->py = 0.0; pn->pz = 0.0;
	  pn->l1 = 0.0; pn->l2 = 0.0; pn->l3 = 0.0;

	  Q[0][0] =  0.0;     Q[0][1] = -pe->b1;  Q[0][2] = -pe->b2;  Q[0][3] = -pe->b3;
	  Q[1][0] =  pe->b1;  Q[1][1] =  0.0;     Q[1][2] =  pe->b3;  Q[1][3] =  -pe->b2;
	  Q[2][0] =  pe->b2;  Q[2][1] = -pe->b3;  Q[2][2] =  0.0;     Q[2][3] = -pe->b1;
	  Q[3][0] =  pe->b3;  Q[3][1] = pe->b2;   Q[3][2] = -pe->b1;  Q[3][3] =  0.0;

	  for (k = 0; k < 4; k++)
	    for (l = 0; l < 4; l++)
	      {
		Q[k][l] *= 0.25*rod[n].le;
		Qplus[k][l] = Q[k][l];
	      }

	  for (k = 0; k < 4; k++)
	    for (l = 0; l < 4; l++)
	      Q[k][l] *= -1.0;
	  for (k = 0; k < 4; k++)
	    {
	      Q[k][k] += 1.0;
	      Qplus[k][k] += 1.0;
	    }
	  Qinv (Q);

	  dq0 = Qplus[0][0]*q0 + Qplus[0][1]*qx + Qplus[0][2]*qy + Qplus[0][3]*qz;
	  dqx = Qplus[1][0]*q0 + Qplus[1][1]*qx + Qplus[1][2]*qy + Qplus[1][3]*qz;
	  dqy = Qplus[2][0]*q0 + Qplus[2][1]*qx + Qplus[2][2]*qy + Qplus[2][3]*qz;
	  dqz = Qplus[3][0]*q0 + Qplus[3][1]*qx + Qplus[3][2]*qy + Qplus[3][3]*qz;	  

	  q0 =  Q[0][4]*dq0 + Q[0][5]*dqx + Q[0][6]*dqy + Q[0][7]*dqz;
	  qx =  Q[1][4]*dq0 + Q[1][5]*dqx + Q[1][6]*dqy + Q[1][7]*dqz;
	  qy =  Q[2][4]*dq0 + Q[2][5]*dqx + Q[2][6]*dqy + Q[2][7]*dqz;
	  qz =  Q[3][4]*dq0 + Q[3][5]*dqx + Q[3][6]*dqy + Q[3][7]*dqz;

	  q  = sqrt( q0*q0 + qx*qx + qy*qy + qz*qz);
	  q0 /= q; qx /= q; qy /= q;  qz /= q; 
	  pn->q0 = q0; pn->qx = qx; pn->qy = qy; pn->qz = qz; 
	}
      
      rx = X;  ry = Y;  rz = Z;

      for (m = m0-1; m >= 1; m--)                                 /* Calculate coordinates */
	{
	  pe = &rod[n].elems[m];
	  pn = &rod[n].nodes[m];

	  q0 = 0.5*(rod[n].nodes[m+1].q0 + pn->q0);
	  qx = 0.5*(rod[n].nodes[m+1].qx + pn->qx);	 
	  qy = 0.5*(rod[n].nodes[m+1].qy + pn->qy);
	  qz = 0.5*(rod[n].nodes[m+1].qz + pn->qz);

	  q  = sqrt( q0*q0 + qx*qx + qy*qy + qz*qz);
	  q0 /= q; qx /= q; qy /= q;  qz /= q; 
	  
	  R[0][0] = -qy*qy+qx*qx-qz*qz+q0*q0; R[0][1] = 2*qy*qx+2*qz*q0; R[0][2] = -2*qy*q0+2*qx*qz;
	  R[1][0] = 2*qy*qx-2*qz*q0; R[1][1] =  qy*qy-qx*qx-qz*qz+q0*q0; R[1][2] =  2*qy*qz+2*qx*q0;
	  R[2][0] = 2*qy*q0+2*qx*qz; R[2][1] =  2*qy*qz-2*qx*q0; R[2][2] = -qy*qy-qx*qx+qz*qz+q0*q0;

	  rx -= le*(R[0][0]*pe->s1 + R[1][0]*pe->s2 + R[2][0]*pe->s3);
	  ry -= le*(R[0][1]*pe->s1 + R[1][1]*pe->s2 + R[2][1]*pe->s3);
	  rz -= le*(R[0][2]*pe->s1 + R[1][2]*pe->s2 + R[2][2]*pe->s3);
	  pn->rx = rx;  pn->ry = ry;  pn->rz = rz;
	}

      rx = X;  ry = Y;  rz = Z;

      for (m = m0+1; m < rod[n].n_nodes-1; m++)                      /* Calculate coordinates */
	{
	  pe = &rod[n].elems[m];
	  pn = &rod[n].nodes[m];

	  q0 = 0.5*(rod[n].nodes[m-1].q0 + pn->q0);
	  qx = 0.5*(rod[n].nodes[m-1].qx + pn->qx);	 
	  qy = 0.5*(rod[n].nodes[m-1].qy + pn->qy);
	  qz = 0.5*(rod[n].nodes[m-1].qz + pn->qz);

	  q  = sqrt( q0*q0 + qx*qx + qy*qy + qz*qz);
	  q0 /= q; qx /= q; qy /= q;  qz /= q; 

	  R[0][0] = -qy*qy+qx*qx-qz*qz+q0*q0; R[0][1] = 2*qy*qx+2*qz*q0; R[0][2] = -2*qy*q0+2*qx*qz;
	  R[1][0] = 2*qy*qx-2*qz*q0; R[1][1] =  qy*qy-qx*qx-qz*qz+q0*q0; R[1][2] =  2*qy*qz+2*qx*q0;
	  R[2][0] = 2*qy*q0+2*qx*qz; R[2][1] =  2*qy*qz-2*qx*q0; R[2][2] = -qy*qy-qx*qx+qz*qz+q0*q0;

	  rx += le*(R[0][0]*pe->s1 + R[1][0]*pe->s2 + R[2][0]*pe->s3);
	  ry += le*(R[0][1]*pe->s1 + R[1][1]*pe->s2 + R[2][1]*pe->s3);
	  rz += le*(R[0][2]*pe->s1 + R[1][2]*pe->s2 + R[2][2]*pe->s3);

	  pn->rx = rx;  pn->ry = ry;  pn->rz = rz;
	}

      tip(rod);
      //printf("%.6f \t  %.6f \t  %.6f \t  %.6f \t  %.6f \t  %.6f \t %.6f \n", rod->x0, rod->y0, rod->z0, rod->q00, rod->qx0, rod->qy0, rod->qz0);
      //printf("%.6f \t  %.6f \t  %.6f \t  %.6f \t  %.6f \t  %.6f \t %.6f \n", rod->xL, rod->yL, rod->zL, rod->q0L, rod->qxL, rod->qyL, rod->qzL);

      strcpy(file1,datadir);
      strcat(file1,"/config.ods");
      filep=fopen(file1,"w");
      for(k=1; k<n_nodes-1; k++)
	{
	  pn = &rod[0].nodes[k];
	  	  fprintf (filep, "%i \t %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \n", k, pn->rx, pn->ry, pn->rz, pn->q0, pn->qx, pn->qy, pn->qz, pn->px, pn->py, pn->pz, pn->l1, pn->l2, pn->l3);

	  //	  printf ("%i \t %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \t  %.16f \n", k, pn->rx, pn->ry, pn->rz, pn->q0, pn->qx, pn->qy, pn->qz, pn->px, pn->py, pn->pz, pn->l1, pn->l2, pn->l3);

	}      
      fclose(filep);
      
      for(k=0; k<n_nodes; k++)   // setting the structure elements (positions and momenta) back to zero
	{
	  pn = &rod[0].nodes[k];
	  pn->rx = 0;  pn->ry = 0;  pn->rz = 0; pn->qx = 0;  pn->qy = 0;  pn->qz = 0; pn->q0 = 0;
	  pn->px = 0;  pn->py = 0;  pn->pz = 0; pn->l1 = 0;  pn->l2 = 0;  pn->l3 = 0;
	}      
    }
}

void pop(struct rod *rod, int n_rods, char *origdir)
{
  struct node *pn;
  int k=0, node_num=123;
  FILE *filep=0;
  char file1[64];

#ifdef checkpoint_resume
  extern int t;
  extern unsigned short Xi[3];

  strcpy(file1,origdir);
  strcat(file1,"/checkpoint.ods");
  filep=fopen(file1,"r");
  fscanf(filep,"%i %i %i %i",&t,&Xi[0],&Xi[1],&Xi[2]);
  for(k=1; k<rod[0].n_nodes-1; k++)
    {
      pn = &rod[0].nodes[k];

      fscanf(filep, "%i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &node_num, &(pn->rx), &(pn->ry), &(pn->rz), &(pn->q0), &(pn->qx), &(pn->qy), &(pn->qz), &(pn->px), &(pn->py), &(pn->pz), &(pn->l1), &(pn->l2), &(pn->l3) );

      if (node_num != k)
	{
	  printf("node number not equal to k, something wrong, aborting \n");
	  exit(0);
	}
    }      
  fclose(filep);

  pn = &rod[0].nodes[0];
  pn->rx=0; pn->ry=0; pn->rz=0; pn->q0=1; pn->qx=0; pn->qy=0; pn->qz=0; pn->px=0; pn->py=0; pn->pz=0; pn->l1=0; pn->l2=0; pn->l3=0;
  pn = &rod[0].nodes[rod->n_nodes-1];
  pn->rx=0; pn->ry=0; pn->rz=0; pn->q0=1; pn->qx=0; pn->qy=0; pn->qz=0; pn->px=0; pn->py=0; pn->pz=0; pn->l1=0; pn->l2=0; pn->l3=0;

  rotmatrix_node(rod);
  bc(rod,'q');
  rotmatrix_elem(rod);
  bc(rod,'r');

#else

  strcpy(file1,origdir);
  strcat(file1,"/config.ods");
  filep=fopen(file1,"r");

  for(k=1; k<rod[0].n_nodes-1; k++)
    {
      pn = &rod[0].nodes[k];

#ifndef real_float
      fscanf(filep, "%i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &node_num, &(pn->rx), &(pn->ry), &(pn->rz), &(pn->q0), &(pn->qx), &(pn->qy), &(pn->qz), &(pn->px), &(pn->py), &(pn->pz), &(pn->l1), &(pn->l2), &(pn->l3) );
#else
      fscanf(filep, "%i %f %f %f %f %f %f %f %f %f %f %f %f %f", &node_num, &(pn->rx), &(pn->ry), &(pn->rz), &(pn->q0), &(pn->qx), &(pn->qy), &(pn->qz), &(pn->px), &(pn->py), &(pn->pz), &(pn->l1), &(pn->l2), &(pn->l3) );
#endif
      if (node_num != k)
	{
	  printf("node number not equal to k, something wrong, aborting \n");
	  exit(0);
	}
    }      
  fclose(filep);

  pn = &rod[0].nodes[0]; //// Setting virtual nodes - virtual nodes are not used.
  pn->rx=0; pn->ry=0; pn->rz=0; pn->q0=1; pn->qx=0; pn->qy=0; pn->qz=0; pn->py=0; pn->pz=0; pn->l1=0; pn->l2=0; pn->l3=0;
  pn = &rod[0].nodes[rod->n_nodes-1];				                                                         
  pn->rx=0; pn->ry=0; pn->rz=0; pn->q0=1; pn->qx=0; pn->qy=0; pn->qz=0; pn->py=0; pn->pz=0; pn->l1=0; pn->l2=0; pn->l3=0;

#ifdef real_float
  normalize(rod);
#endif

  rotmatrix_node(rod);
  bc(rod,'q');
  rotmatrix_elem(rod);
  bc(rod,'r');
#endif
}
