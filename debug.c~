#include "proto.h"

void normalize(struct rod *rod)
{
  struct node *pn;
  int m;
  real q0, qx, qy, qz, q2, q;

  for (m = 1; m < rod->n_nodes-1; m++)
    {	 
      pn = &rod->nodes[m];
      q0 = pn->q0;
      qx = pn->qx;
      qy = pn->qy;
      qz = pn->qz;
      q2 = q0*q0 + qx*qx + qy*qy + qz*qz;
      q   = sqrt(q2);
      pn->q0 /= q; qx /= q; qy /= q; qz /= q;
    }
}

void sumtandf(struct rod *rod) // prints total F and T and any given time
{
  struct node *pn;
  real tx, ty, tz, fx, fy, fz;
  real sumtx=0, sumty=0, sumtz=0;
  real sumfx=0, sumfy=0, sumfz=0;
  int m;

  com(rod); //populates CoM position at the current time point

  for (m = 1; m < rod->n_nodes-1; m++)
    {	 
      pn = &rod->nodes[m];
            
      fx = pn->d[0][0]*pn->f1 + pn->d[1][0]*pn->f2 + pn->d[2][0]*pn->f3;
      fy = pn->d[0][1]*pn->f1 + pn->d[1][1]*pn->f2 + pn->d[2][1]*pn->f3;
      fz = pn->d[0][2]*pn->f1 + pn->d[1][2]*pn->f2 + pn->d[2][2]*pn->f3;
      tx = pn->d[0][0]*pn->t1 + pn->d[1][0]*pn->t2 + pn->d[2][0]*pn->t3;
      ty = pn->d[0][1]*pn->t1 + pn->d[1][1]*pn->t2 + pn->d[2][1]*pn->t3;
      tz = pn->d[0][2]*pn->t1 + pn->d[1][2]*pn->t2 + pn->d[2][2]*pn->t3;
      
      sumfx += fx;
      sumfy += fy;
      sumfz += fz;
      sumtx += tx;
      sumty += ty;
      sumtz += tz;
      sumtx += fz*(pn->ry - rod->ycm) - fy*(pn->rz - rod->zcm);
      sumty += fx*(pn->rz - rod->zcm) - fz*(pn->rx - rod->xcm);
      sumtz += fy*(pn->rx - rod->xcm) - fx*(pn->ry - rod->ycm);

      //  printf("sumT+rXf: \t node %i \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \n", m, sumtx, sumty, sumtz, sumfx, sumfy, sumfz );

    }
    printf("sumT+rXf: \t node %i \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \n", m, sumtx, sumty, sumtz, sumfx, sumfy, sumfz );
}  

void suml(struct rod *rod) // prints total L and any given time
{
  struct node *pn;
  real lx, ly, lz ;
  real sumx=0, sumy=0, sumz=0;
  int m;

  com(rod); //populates com position

  for (m = 1; m < rod->n_nodes-1; m++)
    {	 
      pn = &rod->nodes[m];
            
      lx = pn->d[0][0]*pn->l1 + pn->d[1][0]*pn->l2 + pn->d[2][0]*pn->l3;
      ly = pn->d[0][1]*pn->l1 + pn->d[1][1]*pn->l2 + pn->d[2][1]*pn->l3;
      lz = pn->d[0][2]*pn->l1 + pn->d[1][2]*pn->l2 + pn->d[2][2]*pn->l3;

      sumx += lx;
      sumy += ly;
      sumz += lz;

      sumx += pn->pz*(pn->ry - rod->ycm) - pn->py*(pn->rz - rod->zcm);
      sumy += pn->px*(pn->rz - rod->zcm) - pn->pz*(pn->rx - rod->xcm);
      sumz += pn->py*(pn->rx - rod->xcm) - pn->px*(pn->ry - rod->ycm);
    }
  printf("%.16e \t %.16e \t %.16e \n", sumx, sumy, sumz);
}

void lalpha(struct rod *rod) //Rotates L from body to space frame
{
  struct node *pn;
  int m;

  for (m = 1; m < rod->n_nodes-1; m++)
    {	 
      pn = &rod->nodes[m];
      pn->lx = pn->d[0][0]*pn->l1 + pn->d[1][0]*pn->l2 + pn->d[2][0]*pn->l3;
      pn->ly = pn->d[0][1]*pn->l1 + pn->d[1][1]*pn->l2 + pn->d[2][1]*pn->l3;
      pn->lz = pn->d[0][2]*pn->l1 + pn->d[1][2]*pn->l2 + pn->d[2][2]*pn->l3;
    }
}

void li(struct rod *rod) //Rotates L from space to body frame
{
  struct node *pn;
  int m;

  for (m = 1; m < rod->n_nodes-1; m++)
    {	 
      pn = &rod->nodes[m];
      pn->l1 = pn->d[0][0]*pn->lx + pn->d[0][1]*pn->ly + pn->d[0][2]*pn->lz;
      pn->l2 = pn->d[1][0]*pn->lx + pn->d[1][1]*pn->ly + pn->d[1][2]*pn->lz;
      pn->l3 = pn->d[2][0]*pn->lx + pn->d[2][1]*pn->ly + pn->d[2][2]*pn->lz;
    }
}

void mm(int r1, int c1, real fir[][c1], int r2, int c2, real sec[][c2], real prod[][c2] ) // Matrix multiplication function
{
  int i,j,k;
  real sum=0;

  if(c1 != r2)
    {
      printf("\n Incompatible matrices, aborting now \n");
      exit(0);
    }

  for (i=0;i<r1;i++)
    {
      for (j=0;j<c2;j++)
	{
	  for (k=0;k<c1;k++)
	    {
	      sum = sum + fir[i][k]*sec[k][j];
	    }
	  prod[i][j]=sum;
	  sum=0;
	}
    }
  printf("\n");
  for (i=0;i<r1;i++)
    {
       printf("\n");
       for (j=0;j<c2;j++)
	{
	  printf("%.1f \t",prod[i][j]);
	}
    }
}
