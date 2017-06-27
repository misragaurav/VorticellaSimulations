#include "proto.h"

void load (struct rod *rod)
{
  struct node *pn;
  real f1, f2, f3, t1, t2, t3;
  int m;

  for (m=0; m < rod->n_nodes; m++)
    {
      pn = &rod->nodes[m];
      
      f1 = pn->d[0][0]*Fx + pn->d[0][1]*Fy + pn->d[0][2]*Fz + F1; 
      f2 = pn->d[1][0]*Fx + pn->d[1][1]*Fy + pn->d[1][2]*Fz + F2; 
      f3 = pn->d[2][0]*Fx + pn->d[2][1]*Fy + pn->d[2][2]*Fz + F3; 
      t1 = pn->d[0][0]*Mx + pn->d[0][1]*My + pn->d[0][2]*Mz + M1; 
      t2 = pn->d[1][0]*Mx + pn->d[1][1]*My + pn->d[1][2]*Mz + M2; 
      t3 = pn->d[2][0]*Mx + pn->d[2][1]*My + pn->d[2][2]*Mz + M3; 
      
      pn->f1 += f1;	//force per unit length
      pn->f2 += f2;
      pn->f3 += f3;
      pn->t1 += t1;
      pn->t2 += t2;
      pn->t3 += t3;
    }

  //Boundary load at L
  pn = &rod->nodes[rod->n_nodes-1];

  f1 = pn->d[0][0]*FLx + pn->d[0][1]*FLy + pn->d[0][2]*FLz; 
  f2 = pn->d[1][0]*FLx + pn->d[1][1]*FLy + pn->d[1][2]*FLz; 
  f3 = pn->d[2][0]*FLx + pn->d[2][1]*FLy + pn->d[2][2]*FLz; 
  t1 = pn->d[0][0]*MLx + pn->d[0][1]*MLy + pn->d[0][2]*MLz; 
  t2 = pn->d[1][0]*MLx + pn->d[1][1]*MLy + pn->d[1][2]*MLz; 
  t3 = pn->d[2][0]*MLx + pn->d[2][1]*MLy + pn->d[2][2]*MLz; 

  pn->f1 += f1;
  pn->f2 += f2;
  pn->f3 += f3;
  pn->t1 += t1;
  pn->t2 += t2;
  pn->t3 += t3;
}

void force(struct rod *rod)  /* Accelerations in BFF */
{
  struct elem   *pe;
  struct node   *pn, *pn1, *pn2;
  real ep[3][4];
  real drx, dry, drz, dq0, dqx, dqy, dqz;
  real s1, s2, s3, b1, b2, b3;
  real f1, f2, f3, t1, t2, t3;
  real fx, fy, fz, t0, tx, ty, tz;
  real tx0, tx1, tx2, tx3, tdq;
  real le;
  int  m;

  le = rod->le;

  // Setting nodal forces to zero at the beginning of each time step
  for (m = 0; m < rod->n_nodes; m++)
  {
      pn = &rod->nodes[m];
      pn->fx = pn->fy = pn->fz = 0.0;
      pn->t0 = pn->tx = pn->ty = pn->tz = 0.0;
  }

  /*  Calculate nodal forces and torques  */

  for (m = 0; m < rod->n_elems; m++)
    {
      pe  = &rod->elems[m];                /* Element m touches nodes m and m+1 */
      pn1 = &rod->nodes[m];
      pn2 = &rod->nodes[m+1];

      drx = (pn2->rx - pn1->rx)/le;        /* Derivatives */
      dry = (pn2->ry - pn1->ry)/le;
      drz = (pn2->rz - pn1->rz)/le;
      dq0 = (pn2->q0 - pn1->q0)/le;
      dqx = (pn2->qx - pn1->qx)/le;
      dqy = (pn2->qy - pn1->qy)/le;
      dqz = (pn2->qz - pn1->qz)/le;

      ep[0][0] = -dqx;  ep[0][1] =  dq0; ep[0][2] =  dqz;  ep[0][3] = -dqy;
      ep[1][0] = -dqy;  ep[1][1] = -dqz; ep[1][2] =  dq0;  ep[1][3] =  dqx;
      ep[2][0] = -dqz;  ep[2][1] =  dqy; ep[2][2] = -dqx;  ep[2][3] =  dq0;

#ifdef q_normalization
      ep[0][0] /= pe->q;  ep[0][1] /= pe->q; ep[0][2] /= pe->q;  ep[0][3] /= pe->q;
      ep[1][0] /= pe->q;  ep[1][1] /= pe->q; ep[1][2] /= pe->q;  ep[1][3] /= pe->q;
      ep[2][0] /= pe->q;  ep[2][1] /= pe->q; ep[2][2] /= pe->q;  ep[2][3] /= pe->q;
#endif

      s1  = pe->d[0][0]*drx + pe->d[0][1]*dry + pe->d[0][2]*drz;  /* Strains */
      s2  = pe->d[1][0]*drx + pe->d[1][1]*dry + pe->d[1][2]*drz;
      s3  = pe->d[2][0]*drx + pe->d[2][1]*dry + pe->d[2][2]*drz;
      b1  = 2*(pe->e[0][0]*dq0 + pe->e[0][1]*dqx + pe->e[0][2]*dqy + pe->e[0][3]*dqz);
      b2  = 2*(pe->e[1][0]*dq0 + pe->e[1][1]*dqx + pe->e[1][2]*dqy + pe->e[1][3]*dqz);
      b3  = 2*(pe->e[2][0]*dq0 + pe->e[2][1]*dqx + pe->e[2][2]*dqy + pe->e[2][3]*dqz);

/* Direct shear and bending forces */
      f1  = pe->cs1*(s1-pe->s01);
      f2  = pe->cs2*(s2-pe->s02);
      f3  = pe->cs3*(s3-pe->s03);
      t1  = pe->cb1*(b1-pe->b01);
      t2  = pe->cb2*(b2-pe->b02);
      t3  = pe->cb3*(b3-pe->b03);

      //printf("m: %i \t f1: %f \t f2: %f \t f3: %f \t t1: %f \t t2: %f \t t3 %f \n", m, f1, f2, f3, t1, t2, t3);

/* Gamma X F */
      tx0 = s1*f1 + s2*f2 + s3*f3;                /* Include dot product for projection */
      tx1 = s2*f3 - s3*f2;
      tx2 = s3*f1 - s1*f3;
      tx3 = s1*f2 - s2*f1;

/* Rotate forces to natural coordinate systems */
      fx  = pe->d[0][0]*f1 + pe->d[1][0]*f2 + pe->d[2][0]*f3;                  /* dF/ds */
      fy  = pe->d[0][1]*f1 + pe->d[1][1]*f2 + pe->d[2][1]*f3;
      fz  = pe->d[0][2]*f1 + pe->d[1][2]*f2 + pe->d[2][2]*f3;

      pn1->fx += fx/le;
      pn1->fy += fy/le;
      pn1->fz += fz/le;
      pn2->fx -= fx/le;
      pn2->fy -= fy/le;
      pn2->fz -= fz/le;

      t0  = 2*(pe->e[0][0]*t1 + pe->e[1][0]*t2 + pe->e[2][0]*t3)/le;           /* dT/ds */
      tx  = 2*(pe->e[0][1]*t1 + pe->e[1][1]*t2 + pe->e[2][1]*t3)/le;
      ty  = 2*(pe->e[0][2]*t1 + pe->e[1][2]*t2 + pe->e[2][2]*t3)/le;
      tz  = 2*(pe->e[0][3]*t1 + pe->e[1][3]*t2 + pe->e[2][3]*t3)/le;

      pn1->t0 += t0;
      pn1->tx += tx;
      pn1->ty += ty;
      pn1->tz += tz;
      pn2->t0 -= t0;
      pn2->tx -= tx;
      pn2->ty -= ty;
      pn2->tz -= tz;

#ifdef q_normalization
      t0  = (pe->e[0][0]*tx1 + pe->e[1][0]*tx2 + pe->e[2][0]*tx3)/pe->q; // Gamma X F
      tx  = (pe->e[0][1]*tx1 + pe->e[1][1]*tx2 + pe->e[2][1]*tx3)/pe->q;
      ty  = (pe->e[0][2]*tx1 + pe->e[1][2]*tx2 + pe->e[2][2]*tx3)/pe->q;
      tz  = (pe->e[0][3]*tx1 + pe->e[1][3]*tx2 + pe->e[2][3]*tx3)/pe->q;
#else
      t0  = (pe->e[0][0]*tx1 + pe->e[1][0]*tx2 + pe->e[2][0]*tx3 - pe->q0*tx0)/pe->q2; // Gamma X F
      tx  = (pe->e[0][1]*tx1 + pe->e[1][1]*tx2 + pe->e[2][1]*tx3 - pe->qx*tx0)/pe->q2;
      ty  = (pe->e[0][2]*tx1 + pe->e[1][2]*tx2 + pe->e[2][2]*tx3 - pe->qy*tx0)/pe->q2;
      tz  = (pe->e[0][3]*tx1 + pe->e[1][3]*tx2 + pe->e[2][3]*tx3 - pe->qz*tx0)/pe->q2;
#endif

      pn1->t0 += t0;
      pn1->tx += tx;
      pn1->ty += ty;
      pn1->tz += tz;
      pn2->t0 += t0;
      pn2->tx += tx;
      pn2->ty += ty;
      pn2->tz += tz;

      t0  = ep[0][0]*t1 + ep[1][0]*t2 + ep[2][0]*t3;               /* Omega X T */
      tx  = ep[0][1]*t1 + ep[1][1]*t2 + ep[2][1]*t3;
      ty  = ep[0][2]*t1 + ep[1][2]*t2 + ep[2][2]*t3;
      tz  = ep[0][3]*t1 + ep[1][3]*t2 + ep[2][3]*t3;

#ifdef q_normalization
      tdq = t0*pe->q0 + tx*pe->qx + ty*pe->qy + tz*pe->qz;
      t0  -= pe->q0*tdq;
      tx  -= pe->qx*tdq;
      ty  -= pe->qy*tdq;
      tz  -= pe->qz*tdq;
#endif
      pn1->t0 += t0;
      pn1->tx += tx;
      pn1->ty += ty;
      pn1->tz += tz;
      pn2->t0 += t0;
      pn2->tx += tx;
      pn2->ty += ty;
      pn2->tz += tz;
    }
  
/* Rotate forces and torques to body frame */

  for (m = 0; m < rod->n_nodes; m++)
    {
      pn  = &rod->nodes[m];

      fx = pn->fx; 
      fy = pn->fy; 
      fz = pn->fz;
      t0 = pn->t0; 
      tx = pn->tx; 
      ty = pn->ty; 
      tz = pn->tz;

      pn->f1 =  pn->d[0][0]*fx + pn->d[0][1]*fy + pn->d[0][2]*fz;
      pn->f2 =  pn->d[1][0]*fx + pn->d[1][1]*fy + pn->d[1][2]*fz;
      pn->f3 =  pn->d[2][0]*fx + pn->d[2][1]*fy + pn->d[2][2]*fz;
      pn->t1 = (pn->e[0][0]*t0 + pn->e[0][1]*tx + pn->e[0][2]*ty + pn->e[0][3]*tz)*0.5;
      pn->t2 = (pn->e[1][0]*t0 + pn->e[1][1]*tx + pn->e[1][2]*ty + pn->e[1][3]*tz)*0.5;
      pn->t3 = (pn->e[2][0]*t0 + pn->e[2][1]*tx + pn->e[2][2]*ty + pn->e[2][3]*tz)*0.5;

      /*
  FILE *file_ptr=0;
      if (t%(tprint) == 0)    
	{
	  char forces[64];
	  strcpy(forces,datadir);
	  strcat(forces, "/forces.ods");
	  file_ptr = fopen (forces, "a");
	  fprintf (file_ptr,"%.3f \t %i \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \n", t*dt, m, pn->f1, pn->f2, pn->f3, pn->t1, pn->t2, pn->t3);
	  fclose (file_ptr);
	}
      */
      //      printf("m: %i \t f1: %f \t f2: %f \t f3: %f \t t1: %f \t t2: %f \t t3 %f \n", m, pn->f1, pn->f2, pn->f3, pn->t1, pn->t2, pn->t3);
    }

  //    load(rod);


  // BC at 1
  /*
      pn  = &rod->nodes[0];
      pn->f1 = 0;
      pn->f2 = 0;
      pn->f3 = 0;        
      pn->t1 = 0;
      pn->t2 = 0;
      pn->t3 = 0;

  // BC at  L
     
      pn  = &rod->nodes[rod->n_nodes-1];
      pn->f1 = 0;
      pn->f2 = 0;
      pn->f3 = 0;        
      pn->t1 = 0;
      pn->t2 = 0;
      pn->t3 = 0;
      */
  
  //  sumtandf(rod);

}
