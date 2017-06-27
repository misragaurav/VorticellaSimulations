#include "proto.h"
#define  rng   ran48
#define  ran32 ( (real)int_ran[(int)(12.0*rand()/((real)(RAND_MAX) + 1.0))] )
#define  ran48 ( (real)int_ran[(int)(12.0*erand48(Xi))] )
real int_ran[12] = {-2.0,-1.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,2.0};
extern unsigned short Xi[3];

int M = 1;

void q(struct rod *rod, real dt, int dir)
{
  struct node  *pn;
  real l[3], q[4];
  real i[3];
  real eq[4], c, s, c2, s2, halfangle;
  real l1_temp, l2_temp;           
  int  m, n, dir1, dir2;

  i[0] = rod->i1;  i[1] = rod->i2;  i[2] = rod->i3;

  for (m = 0; m < rod->n_nodes; m++)
    {
      pn = &rod->nodes[m];      
      q[0] = pn->q0;  q[1] = pn->qx;  q[2] = pn->qy; q[3] = pn->qz;
      l[0] = pn->l1;  l[1] = pn->l2;  l[2] = pn->l3;

      switch (dir)
	{
	case 0: eq[0] =-q[1]; eq[1] = q[0]; eq[2] = q[3]; eq[3] = -q[2]; break;
	case 1: eq[0] =-q[2]; eq[1] =-q[3]; eq[2] = q[0]; eq[3] =  q[1]; break;
	case 2: eq[0] =-q[3]; eq[1] = q[2]; eq[2] =-q[1]; eq[3] =  q[0]; break;
	}

      halfangle = 0.5*( (l[dir])/(i[dir]) )*dt;

      c = cos(halfangle);  s = sin(halfangle);

      for (n = 0; n < 4; n++)
	{ q[n] = c*q[n] + s*eq[n]; }
      pn->q0 = q[0];
      pn->qx = q[1];
      pn->qy = q[2];
      pn->qz = q[3];

      c2 = cos(2*halfangle); s2 = sin(2*halfangle);

      dir1 = (dir+1)%3;
      dir2 = (dir+2)%3;
      
      l1_temp = l[dir1];                   // stored old values of momenta
      l2_temp = l[dir2];
      l[dir2] = l2_temp*c2 - l1_temp*s2;   // Rotating the frames
      l[dir1] = l2_temp*s2 + l1_temp*c2;
      
      pn->l1 = l[0];
      pn->l2 = l[1];
      pn->l3 = l[2];
    }
}

void r(struct rod *rod, real dt)
{
  struct node *pn;
  real xm;
  int m;

  xm = rod->xm;
  for (m = 0; m < rod->n_nodes; m++)
    {      
      pn = &rod->nodes[m];
      pn->rx += pn->px*dt/xm;
      pn->ry += pn->py*dt/xm;
      pn->rz += pn->pz*dt/xm;
    }
}

void momenta(struct rod *rod, real dt)
{
  struct node *pn;
  int  m;

  for (m = 0; m < rod->n_nodes; m++)
    {
      pn = &rod->nodes[m];
      pn->p1 = pn->d[0][0]*pn->px + pn->d[0][1]*pn->py + pn->d[0][2]*pn->pz;
      pn->p2 = pn->d[1][0]*pn->px + pn->d[1][1]*pn->py + pn->d[1][2]*pn->pz;
      pn->p3 = pn->d[2][0]*pn->px + pn->d[2][1]*pn->py + pn->d[2][2]*pn->pz;
    }

  if (VISCOCITY==0)
    {
      for (m = 0; m < rod->n_nodes; m++)
	{
	  pn = &rod->nodes[m];
	  pn->p1 += pn->f1*dt;
	  pn->p2 += pn->f2*dt;
	  pn->p3 += pn->f3*dt;
	  pn->l1 += pn->t1*dt;
	  pn->l2 += pn->t2*dt;
	  pn->l3 += pn->t3*dt;
	}
    }
  else
    {
      friction(rod, dt, 1.0); // If friction changes with time then this call is needed
      
      for (m = 0; m < rod->n_nodes; m++)
	{
	  pn = &rod->nodes[m];
	  pn->p1 = pn->p1*rod->At1 + pn->f1*rod->Bt1 + rod->Ct1*rng;       
	  pn->p2 = pn->p2*rod->At2 + pn->f2*rod->Bt2 + rod->Ct2*rng;     
	  pn->p3 = pn->p3*rod->At3 + pn->f3*rod->Bt3 + rod->Ct3*rng;
	  pn->l1 = pn->l1*rod->Ar1 + pn->t1*rod->Br1 + rod->Cr1*rng;     
	  pn->l2 = pn->l2*rod->Ar2 + pn->t2*rod->Br2 + rod->Cr2*rng;     
	  pn->l3 = pn->l3*rod->Ar3 + pn->t3*rod->Br3 + rod->Cr3*rng;
	}      
    }

  for (m = 0; m < rod->n_nodes; m++)
    {
      pn = &rod->nodes[m];
      pn->px = pn->d[0][0]*pn->p1 + pn->d[1][0]*pn->p2 + pn->d[2][0]*pn->p3;
      pn->py = pn->d[0][1]*pn->p1 + pn->d[1][1]*pn->p2 + pn->d[2][1]*pn->p3;
      pn->pz = pn->d[0][2]*pn->p1 + pn->d[1][2]*pn->p2 + pn->d[2][2]*pn->p3;
    }
}

void propagator(struct rod *rod, real dt, int n_rods)
{
  int n, m;

  for (n = 0; n < n_rods; n++)
    {
      for (m = 0; m < M; m++) //Qs must be updated before Rs so that rotation of gamma0 is accurate
	{
	  q(&rod[n], 0.5*dt/(2*M), 2);  // Free rotation (L3)
	  q(&rod[n], 0.5*dt/(2*M), 1);  // Free rotation (L2)
	  q(&rod[n], 0.5*dt/(2*M), 0);  // Free rotation (L1)
	  q(&rod[n], 0.5*dt/(2*M), 0);  // Free rotation (L1)
	  q(&rod[n], 0.5*dt/(2*M), 1);  // Free rotation (L2)
	  q(&rod[n], 0.5*dt/(2*M), 2);  // Free rotation (L3)
	}

      rotmatrix_node(rod);
      rotmatrix_elem(rod);

      r(&rod[n], dt/2);
	   
      force(&rod[n]);			    
      momenta(&rod[n], dt);			    

      r(&rod[n], dt/2);

      for (m = 0; m < M; m++) 
	{
	  q(&rod[n], 0.5*dt/(2*M), 2);  // Free rotation (L3) //
    	  q(&rod[n], 0.5*dt/(2*M), 1);  // Free rotation (L2) // 
    	  q(&rod[n], 0.5*dt/(2*M), 0);  // Free rotation (L1) //
    	  q(&rod[n], 0.5*dt/(2*M), 0);  // Free rotation (L1) // 
    	  q(&rod[n], 0.5*dt/(2*M), 1);  // Free rotation (L2) // 
    	  q(&rod[n], 0.5*dt/(2*M), 2);  // Free rotation (L3) // 
	}

      rotmatrix_node(rod);
      rotmatrix_elem(rod);
    }   
}
