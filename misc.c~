#include "proto.h"

void file_name (char *name, char *work_dir, int task_number)
{
  static char   name0[] = {'0','1','2','3','4','5','6','7','8','9'};
  int    index;

  if (work_dir != '\0')
    {
      for (index = strlen(name); index > 0; index--)
	name[index+strlen(work_dir)]=name[index-1];
      name[strlen(work_dir)] = '/';
      for (index = 0; index < strlen(work_dir); index++)
	name[index] = work_dir[index];
    }

  index = strlen(name);
  name[index]   = name0[task_number/100];
  task_number %= 100;
  name[index+1] = name0[task_number/10];
  task_number %= 10;
  name[index+2] = name0[task_number];
  name[index+3] = '.';
  name[index+4] = 'o';
  name[index+5] = 'd';
  name[index+6] = 's';
  name[index+7] = '\0';
}

void Qinv (real Q[4][8])  /* Invert Q matrix */
{
  real ratio;
  int    i, j, k;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      Q[i][j+4] = 0.0;
  for (i = 0; i < 4; i++)
    Q[i][i+4] = 1.0;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      {
	if (i == j)  continue;
	ratio = Q[j][i]/Q[i][i];
	for (k = 0; k < 8; k++)
	  Q[j][k] -= ratio*Q[i][k];
      }
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      Q[i][j+4] /= Q[i][i];
}

void rotmatrix_node(struct rod *rod)
{
  struct node *pn;
  real q0, qx, qy, qz;
  int m;

  for (m = 1; m < rod->n_nodes-1; m++)
    {
      pn  = &rod->nodes[m];
      q0 = pn->q0; qx = pn->qx; qy = pn->qy; qz = pn->qz;
      
      pn->d[0][0] = -qy*qy+qx*qx-qz*qz+q0*q0; pn->d[0][1] = 2*qy*qx+2*qz*q0;          pn->d[0][2] = -2*qy*q0+2*qx*qz;
      pn->d[1][0] = 2*qy*qx-2*qz*q0;          pn->d[1][1] = qy*qy-qx*qx-qz*qz+q0*q0;  pn->d[1][2] = 2*qy*qz+2*qx*q0;
      pn->d[2][0] = 2*qy*q0+2*qx*qz;          pn->d[2][1] = 2*qy*qz-2*qx*q0;          pn->d[2][2] = -qy*qy-qx*qx+qz*qz+q0*q0;
      pn->e[0][0] = -qx;  pn->e[0][1] =  q0; pn->e[0][2] =  qz;  pn->e[0][3] = -qy;
      pn->e[1][0] = -qy;  pn->e[1][1] = -qz; pn->e[1][2] =  q0;  pn->e[1][3] =  qx;
      pn->e[2][0] = -qz;  pn->e[2][1] =  qy; pn->e[2][2] = -qx;  pn->e[2][3] =  q0;
    }
}

void rotmatrix_elem(struct rod *rod)
{
  struct node *pn1, *pn2;
  struct elem *pe;
  real q0, qx, qy, qz, q2, q;
  int m;

  for (m = 0; m < rod->n_nodes-1; m++)
    {
      pe  = &rod->elems[m];               
      pn1 = &rod->nodes[m];
      pn2 = &rod->nodes[m+1];

      q0 = 0.5*(pn1->q0 + pn2->q0);
      qx = 0.5*(pn1->qx + pn2->qx);
      qy = 0.5*(pn1->qy + pn2->qy);
      qz = 0.5*(pn1->qz + pn2->qz);
      q2 = q0*q0 + qx*qx + qy*qy + qz*qz;

#ifdef q_normalization
      q   = sqrt(q2);
      q0 /= q; qx /= q; qy /= q; qz /= q;
#endif

      pe->q0=q0; pe->qx=qx; pe->qy=qy; pe->qz=qz;
      pe->q =q ; pe->q2=q2;

      pe->d[0][0] = -qy*qy+qx*qx-qz*qz+q0*q0; pe->d[0][1] = 2*qy*qx+2*qz*q0;          pe->d[0][2] = -2*qy*q0+2*qx*qz;
      pe->d[1][0] = 2*qy*qx-2*qz*q0;          pe->d[1][1] = qy*qy-qx*qx-qz*qz+q0*q0;  pe->d[1][2] = 2*qy*qz+2*qx*q0;
      pe->d[2][0] = 2*qy*q0+2*qx*qz;          pe->d[2][1] = 2*qy*qz-2*qx*q0;          pe->d[2][2] = -qy*qy-qx*qx+qz*qz+q0*q0;
      pe->e[0][0] = -qx;  pe->e[0][1] =  q0; pe->e[0][2] =  qz;  pe->e[0][3] = -qy;
      pe->e[1][0] = -qy;  pe->e[1][1] = -qz; pe->e[1][2] =  q0;  pe->e[1][3] =  qx;
      pe->e[2][0] = -qz;  pe->e[2][1] =  qy; pe->e[2][2] = -qx;  pe->e[2][3] =  q0;
    }
}

void frictionSphere(struct rod *rod, real dt, real factor)
{
  real taut1, taut2, taut3, taur1, taur2, taur3;
  real a = HeadRad;

  rod->zt1 = factor*(4.5*VISCOCITY/(DENSITY*a*a)); 
  rod->zt2 = factor*(4.5*VISCOCITY/(DENSITY*a*a)); 
  rod->zt3 = factor*(4.5*VISCOCITY/(DENSITY*a*a)); 
  rod->zr1 = factor*( 15*VISCOCITY/(DENSITY*a*a)); 
  rod->zr2 = factor*( 15*VISCOCITY/(DENSITY*a*a)); 
  rod->zr3 = factor*( 15*VISCOCITY/(DENSITY*a*a)); 

  taut1 =  1/rod->zt1;   rod->At1 = exp(-dt/taut1); 
  taut2 =  1/rod->zt2;   rod->At2 = exp(-dt/taut2); 
  taut3 =  1/rod->zt3;   rod->At3 = exp(-dt/taut3); 
  taur1 =  1/rod->zr1;   rod->Ar1 = exp(-dt/taur1); 
  taur2 =  1/rod->zr2;   rod->Ar2 = exp(-dt/taur2); 
  taur3 =  1/rod->zr3;   rod->Ar3 = exp(-dt/taur3); 

  rod->Bt1 = taut1*(1-rod->At1);    rod->Ct1 = 0; //no fluctuations
  rod->Bt2 = taut2*(1-rod->At2);    rod->Ct2 = 0;
  rod->Bt3 = taut3*(1-rod->At3);    rod->Ct3 = 0;
  rod->Br1 = taur1*(1-rod->Ar1);    rod->Cr1 = 0;
  rod->Br2 = taur2*(1-rod->Ar2);    rod->Cr2 = 0;
  rod->Br3 = taur3*(1-rod->Ar3);    rod->Cr3 = 0;

}

void friction(struct rod *rod, real dt, real factor)
{
  real taut1, taut2, taut3, taur1, taur2, taur3;
  real xm, i1, i2, i3;

  xm=rod->xm; i1=rod->i1; i2=rod->i2; i3=rod->i3;

  rod->zt1 = factor*(   16*VISCOCITY/(DENSITY*DIA*DIA*log(2*N*AR))); // reduced approx
  rod->zt2 = factor*(   16*VISCOCITY/(DENSITY*DIA*DIA*log(2*N*AR))); 
  rod->zt3 = factor*(    8*VISCOCITY/(DENSITY*DIA*DIA*log(2*N*AR))); 
  rod->zr1 = factor*(13.84*VISCOCITY/(DENSITY*DIA*DIA)  ); 
  rod->zr2 = factor*(13.84*VISCOCITY/(DENSITY*DIA*DIA)  ); 
  rod->zr3 = factor*(   32*VISCOCITY/(DENSITY*DIA*DIA)  ); 

  taut1 =  1/rod->zt1;    rod->At1 = exp(-dt/taut1); 
  taut2 =  1/rod->zt2;    rod->At2 = exp(-dt/taut2); 
  taut3 =  1/rod->zt3;    rod->At3 = exp(-dt/taut3); 
  taur1 =  1/rod->zr1;    rod->Ar1 = exp(-dt/taur1); 
  taur2 =  1/rod->zr2;    rod->Ar2 = exp(-dt/taur2); 
  taur3 =  1/rod->zr3;    rod->Ar3 = exp(-dt/taur3); 

  rod->Bt1 = taut1*(1-rod->At1);    rod->Ct1 = sqrt((1-rod->At1*rod->At1)*xm*T); //********* check if xm and i correct for fric per unit mass
  rod->Bt2 = taut2*(1-rod->At2);    rod->Ct2 = sqrt((1-rod->At2*rod->At2)*xm*T);
  rod->Bt3 = taut3*(1-rod->At3);    rod->Ct3 = sqrt((1-rod->At3*rod->At3)*xm*T);
  rod->Br1 = taur1*(1-rod->Ar1);    rod->Cr1 = sqrt((1-rod->Ar1*rod->Ar1)*i1*T);
  rod->Br2 = taur2*(1-rod->Ar2);    rod->Cr2 = sqrt((1-rod->Ar2*rod->Ar2)*i2*T);
  rod->Br3 = taur3*(1-rod->Ar3);    rod->Cr3 = sqrt((1-rod->Ar3*rod->Ar3)*i3*T);

}

void com(struct rod *rod)
{
  struct node *pn;
  real Rx=0, Ry=0, Rz=0;
  int k;

  rod->xcm  = 0.0;
  rod->ycm  = 0.0;
  rod->zcm  = 0.0;
  rod->pxcm = 0.0;
  rod->pycm = 0.0;
  rod->pzcm = 0.0;
  rod->lxcm = 0.0;
  rod->lycm = 0.0;
  rod->lzcm = 0.0;

  for (k = 1; k < rod->n_nodes-1; k++)
    {
      pn = &rod[0].nodes[k];
      rod->xcm += pn->rx/(rod->n_nodes-2);
      rod->ycm += pn->ry/(rod->n_nodes-2);
      rod->zcm += pn->rz/(rod->n_nodes-2);
    }
    
  for (k = 1; k < rod->n_nodes-1; k++)
    {
      pn = &rod->nodes[k];

      rod->pxcm += rod->le*pn->px;
      rod->pycm += rod->le*pn->py;
      rod->pzcm += rod->le*pn->pz;
    }

  for (k = 1; k < rod->n_nodes-1; k++)
    {
      pn = &rod->nodes[k];
      pn->lx = pn->d[0][0]*pn->l1 + pn->d[1][0]*pn->l2 + pn->d[2][0]*pn->l3;
      pn->ly = pn->d[0][1]*pn->l1 + pn->d[1][1]*pn->l2 + pn->d[2][1]*pn->l3;
      pn->lz = pn->d[0][2]*pn->l1 + pn->d[1][2]*pn->l2 + pn->d[2][2]*pn->l3;
    }

  for (k = 1; k < rod->n_nodes-1; k++)
    {
      pn = &rod->nodes[k];
	  
      Rx = pn->rx - rod->xcm;
      Ry = pn->ry - rod->ycm;
      Rz = pn->rz - rod->zcm;

      rod->lxcm += rod->le*(pn->lx + Ry*pn->pz - Rz*pn->py);
      rod->lycm += rod->le*(pn->ly + Rz*pn->px - Rx*pn->pz);
      rod->lzcm += rod->le*(pn->lz + Rx*pn->py - Ry*pn->px);
    }
}

void energy(struct rod *rod)
{
  struct node *pn1, *pn2;
  struct elem *pe;
  real drx, dry, drz, dq0, dqx, dqy, dqz;
  real s1, s2, s3, b1, b2, b3;
  real p1, p2, p3, l1, l2, l3;
  real xm, i1, i2, i3;
  int m;
  
  xm = rod->xm;  i1 = rod->i1;  i2 = rod->i2;  i3 = rod->i3;

  for (m = 1; m < rod->n_nodes-2; m++)
    {
      pn1 = &rod->nodes[m];

      p1 = pn1->d[0][0]*pn1->px + pn1->d[0][1]*pn1->py + pn1->d[0][2]*pn1->pz;
      p2 = pn1->d[1][0]*pn1->px + pn1->d[1][1]*pn1->py + pn1->d[1][2]*pn1->pz;
      p3 = pn1->d[2][0]*pn1->px + pn1->d[2][1]*pn1->py + pn1->d[2][2]*pn1->pz;
      l1 = pn1->l1;  
      l2 = pn1->l2;  
      l3 = pn1->l3;

      pn1->keT1 = p1*p1/(2*xm);
      pn1->keT2 = p2*p2/(2*xm);
      pn1->keT3 = p3*p3/(2*xm);	       
      pn1->keR1 = l1*l1/(2*i1);
      pn1->keR2 = l2*l2/(2*i2);
      pn1->keR3 = l3*l3/(2*i3);

      pn1->ke = pn1->keT1 + pn1->keT2 + pn1->keT3 + pn1->keR1 + pn1->keR2 + pn1->keR3 ;
    }

      pn1 = &rod->nodes[rod->n_nodes-2];

      p1 = pn1->d[0][0]*pn1->px + pn1->d[0][1]*pn1->py + pn1->d[0][2]*pn1->pz;
      p2 = pn1->d[1][0]*pn1->px + pn1->d[1][1]*pn1->py + pn1->d[1][2]*pn1->pz;
      p3 = pn1->d[2][0]*pn1->px + pn1->d[2][1]*pn1->py + pn1->d[2][2]*pn1->pz;
      l1 = pn1->l1;  
      l2 = pn1->l2;  
      l3 = pn1->l3;

      pn1->keT1 = p1*p1/(2*headmass*xm);
      pn1->keT2 = p2*p2/(2*headmass*xm);
      pn1->keT3 = p3*p3/(2*headmass*xm);	       
      pn1->keR1 = l1*l1/(2*headinertia*i1);
      pn1->keR2 = l2*l2/(2*headinertia*i2);
      pn1->keR3 = l3*l3/(2*headinertia*i3);

      pn1->ke = pn1->keT1 + pn1->keT2 + pn1->keT3 + pn1->keR1 + pn1->keR2 + pn1->keR3 ;

#ifdef exclude_ends
  for (m = 1; m < rod->n_elems-1; m++)
#else
  for (m = 0; m < rod->n_elems; m++)
#endif
    {
      pe  = &rod->elems[m];
      pn1 = &rod->nodes[m];
      pn2 = &rod->nodes[m+1];

      drx = (pn2->rx - pn1->rx)/rod->le;
      dry = (pn2->ry - pn1->ry)/rod->le;
      drz = (pn2->rz - pn1->rz)/rod->le;
      dq0 = (pn2->q0 - pn1->q0)/rod->le;
      dqx = (pn2->qx - pn1->qx)/rod->le;
      dqy = (pn2->qy - pn1->qy)/rod->le;
      dqz = (pn2->qz - pn1->qz)/rod->le;

      s1 = pe->d[0][0]*drx + pe->d[0][1]*dry + pe->d[0][2]*drz;
      s2 = pe->d[1][0]*drx + pe->d[1][1]*dry + pe->d[1][2]*drz;
      s3 = pe->d[2][0]*drx + pe->d[2][1]*dry + pe->d[2][2]*drz;
      b1 = 2*(pe->e[0][0]*dq0 + pe->e[0][1]*dqx + pe->e[0][2]*dqy + pe->e[0][3]*dqz);
      b2 = 2*(pe->e[1][0]*dq0 + pe->e[1][1]*dqx + pe->e[1][2]*dqy + pe->e[1][3]*dqz);
      b3 = 2*(pe->e[2][0]*dq0 + pe->e[2][1]*dqx + pe->e[2][2]*dqy + pe->e[2][3]*dqz);

      s1 -= pe->s01;
      s2 -= pe->s02;
      s3 -= pe->s03;
      b1 -= pe->b01;
      b2 -= pe->b02;
      b3 -= pe->b03;

      pe->peT1 = 0.5*s1*s1*pe->cs1;
      pe->peT2 = 0.5*s2*s2*pe->cs2;
      pe->peT3 = 0.5*s3*s3*pe->cs3;
      pe->peR1 = 0.5*b1*b1*pe->cb1;
      pe->peR2 = 0.5*b2*b2*pe->cb2;
      pe->peR3 = 0.5*b3*b3*pe->cb3;
      
      pe->pe = pe->peT1+pe->peT2+pe->peT3+pe->peR1+pe->peR2+pe->peR3;
    }

  for (m=0, rod->pe=0, rod->peT1=0, rod->peT2=0, rod->peT3=0, rod->peR1=0, rod->peR2=0, rod->peR3=0; m < rod->n_elems; m++)
    {
      pe = &rod->elems[m];
      rod->pe   += pe->pe;
      rod->peT1 += pe->peT1 ;
      rod->peT2 += pe->peT2 ;
      rod->peT3 += pe->peT3 ;
      rod->peR1 += pe->peR1 ;
      rod->peR2 += pe->peR2 ;
      rod->peR3 += pe->peR3 ;
    }

  for (m=1, rod->ke=0, rod->keT1=0, rod->keT2=0, rod->keT3=0, rod->keR1=0, rod->keR2=0, rod->keR3=0; m < rod->n_nodes-1; m++)
    {
      pn1 = &rod->nodes[m];
      rod->ke   += pn1->ke;
      rod->keT1 += pn1->keT1 ;
      rod->keT2 += pn1->keT2 ;
      rod->keT3 += pn1->keT3 ;
      rod->keR1 += pn1->keR1 ;
      rod->keR2 += pn1->keR2 ;
      rod->keR3 += pn1->keR3 ;
    }

  // Spasmoneme energy

  real strain=0;

  for (m=1; m < rod->n_nodes-1; m++)
    {
      pn1 = &rod->nodes[m];      
      pn1->xx = pn1->rx + pn1->d[0][0]*pn1->x1 + pn1->d[1][0]*pn1->x2 + pn1->d[2][0]*0;
      pn1->xy = pn1->ry + pn1->d[0][1]*pn1->x1 + pn1->d[1][1]*pn1->x2 + pn1->d[2][1]*0;
      pn1->xz = pn1->rz + pn1->d[0][2]*pn1->x1 + pn1->d[1][2]*pn1->x2 + pn1->d[2][2]*0;
    }

  for (m=1; m < rod->n_elems-1; m++)
    {
      pe  = &rod->elems[m];            
      pn1 = &rod->nodes[m];
      pn2 = &rod->nodes[m+1];
      strain = sqrt( (pn2->xx - pn1->xx)*(pn2->xx - pn1->xx)+(pn2->xy - pn1->xy)*(pn2->xy - pn1->xy)+(pn2->xz - pn1->xz)*(pn2->xz - pn1->xz) )/rod->delS - pe->X0;
      pe->peS = 0.5*pe->c*strain*strain;;
    }

  for (m=1, rod->peS=0; m < rod->n_elems-1; m++)
    {
      pe = &rod->elems[m];
      rod->peS += pe->peS;
    }

  //Free BC
  pe  = &rod->elems[0]; 
  pe->L = 0;
  pe  = &rod->elems[rod->n_elems-1]; 
  pe->L = 0;           

}

void bc(struct rod *rod, char var)
{
  struct node *pn1, *pn2;
  real q0, qx, qy, qz;
  real le;

  le = rod->le;

  switch (var)
    {
    case 'q':
      // Update of ghost node 0
      pn1 = &rod->nodes[0];
      pn2 = &rod->nodes[1];

      q0 = 2*rod->q00 - pn2->q0 ;
      qx = 2*rod->qx0 - pn2->qx ;
      qy = 2*rod->qy0 - pn2->qy ;
      qz = 2*rod->qz0 - pn2->qz ;
      /*      
#ifdef q_normalization
real q;
	q  = sqrt( q0*q0 + qx*qx + qy*qy + qz*qz);
	q0 /= q; qx /= q; qy /= q;  qz /= q; 
#endif
      */
      pn1->q0 = q0; pn1->qx = qx; pn1->qy = qy; pn1->qz = qz; 
      
      // Update of ghost node N+1
      pn1 = &rod->nodes[rod->n_nodes-1];
      pn2 = &rod->nodes[rod->n_nodes-2];

      q0 = 2*rod->q0L - pn2->q0 ;
      qx = 2*rod->qxL - pn2->qx ;
      qy = 2*rod->qyL - pn2->qy ;
      qz = 2*rod->qzL - pn2->qz ;
      /*      
#ifdef q_normalization
      q  = sqrt( q0*q0 + qx*qx + qy*qy + qz*qz);
      q0 /= q; qx /= q; qy /= q;  qz /= q; 
#endif
      */
      pn1->q0 = q0; pn1->qx = qx; pn1->qy = qy; pn1->qz = qz; 
      
      break;

    case 'r':
      // Update of ghost node 0
      pn1 = &rod->nodes[0];
      pn2 = &rod->nodes[1];

      pn1->rx  = 2*rod->x0 - pn2->rx;    
      pn1->ry  = 2*rod->y0 - pn2->ry;
      pn1->rz  = 2*rod->z0 - pn2->rz;

      // Update of ghost node N+1
      pn1 = &rod->nodes[rod->n_nodes-1];
      pn2 = &rod->nodes[rod->n_nodes-2];

      pn1->rx  = 2*rod->xL - pn2->rx;    
      pn1->ry  = 2*rod->yL - pn2->ry;
      pn1->rz  = 2*rod->zL - pn2->rz;
      break;

    case 'l':                  // Might be needed for clamped BC
      // Update of ghost node 0
      pn1 = &rod->nodes[0];
      pn2 = &rod->nodes[1];
      pn1->l1 = pn2->l1; 
      pn1->l2 = pn2->l2;
      pn1->l3 = pn2->l3;

      // Update of ghost node N+1
      pn1 = &rod->nodes[rod->n_nodes-1];
      pn2 = &rod->nodes[rod->n_nodes-2];
      pn1->l1 = pn2->l1; 
      pn1->l2 = pn2->l2;
      pn1->l3 = pn2->l3;
      break;

    case 'p':                  // Might be needed for clamped BC
      // Update of ghost node 0
      pn1 = &rod->nodes[0];
      pn2 = &rod->nodes[1];
      pn1->px = pn2->px ;
      pn1->py = pn2->py ;
      pn1->pz = pn2->pz ;
  
      // Update of ghost node N+1
      pn1 = &rod->nodes[rod->n_nodes-1];
      pn2 = &rod->nodes[rod->n_nodes-2];
      pn1->px = pn2->px ;
      pn1->py = pn2->py ;
      pn1->pz = pn2->pz ;
      break;
    }
}

void tip(struct rod *rod)
{
  struct elem *pe;
  struct node *pn;
  real d[3][3], Q[4][8], Qplus[4][4];
  real rx, ry, rz, q0, qx, qy, qz, q;
  real dq0, dqx, dqy, dqz;
  int k,l;

  pn = &rod[0].nodes[1];
  rx = pn->rx;  ry = pn->ry;  rz = pn->rz;
  q0 = pn->q0;  qx = pn->qx;  qy = pn->qy; qz = pn->qz;

  // Calculate quaternions at tip 0
  pe = &rod[0].elems[0]; // really I need the bend at the point between the tip and first real node
  pn = &rod[0].nodes[0];

  Q[0][0] =  0.0;     Q[0][1] = -pe->b1;  Q[0][2] = -pe->b2;  Q[0][3] = -pe->b3;
  Q[1][0] =  pe->b1;  Q[1][1] =  0.0;     Q[1][2] = -pe->b3;  Q[1][3] =  pe->b2;
  Q[2][0] =  pe->b2;  Q[2][1] =  pe->b3;  Q[2][2] =  0.0;     Q[2][3] = -pe->b1;
  Q[3][0] =  pe->b3;  Q[3][1] = -pe->b2;  Q[3][2] =  pe->b1;  Q[3][3] =  0.0;

  for (k = 0; k < 4; k++)
    for (l = 0; l < 4; l++)
      {
	Q[k][l] *= -0.25*rod[0].le/2; //half le
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
  rod->q00 = q0;  rod->qx0 = qx;  rod->qy0 = qy;  rod->qz0 = qz;
 
  // Calculate quaternions at tip L
  pn = &rod[0].nodes[rod->n_nodes-2];
  rx = pn->rx;  ry = pn->ry;  rz = pn->rz;
  q0 = pn->q0;  qx = pn->qx;  qy = pn->qy; qz = pn->qz;

  pe = &rod[0].elems[rod->n_elems-1];
  pn = &rod[0].nodes[rod->n_nodes-1];

  Q[0][0] =  0.0;     Q[0][1] = -pe->b1;  Q[0][2] = -pe->b2;  Q[0][3] = -pe->b3;
  Q[1][0] =  pe->b1;  Q[1][1] =  0.0;     Q[1][2] = -pe->b3;  Q[1][3] =  pe->b2;
  Q[2][0] =  pe->b2;  Q[2][1] =  pe->b3;  Q[2][2] =  0.0;     Q[2][3] = -pe->b1;
  Q[3][0] =  pe->b3;  Q[3][1] = -pe->b2;  Q[3][2] =  pe->b1;  Q[3][3] =  0.0;

  for (k = 0; k < 4; k++)
    for (l = 0; l < 4; l++)
      {
	Q[k][l] *= 0.25*rod[0].le/2; //half le
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
  rod->q0L = q0;  rod->qxL = qx;  rod->qyL = qy;  rod->qzL = qz;

  // Calculate coordinates at 0 tip
  pn = &rod[0].nodes[1];
  rx = pn->rx;  ry = pn->ry;  rz = pn->rz;  

  q0 = 0.5*(rod->q00 + pn->q0);
  qx = 0.5*(rod->qx0 + pn->qx);	 
  qy = 0.5*(rod->qy0 + pn->qy);
  qz = 0.5*(rod->qz0 + pn->qz);

  q  = sqrt( q0*q0 + qx*qx + qy*qy + qz*qz);
  q0 /= q; qx /= q; qy /= q;  qz /= q; 
	  
  d[0][0] = -qy*qy+qx*qx-qz*qz+q0*q0; d[0][1] = 2*qy*qx+2*qz*q0; d[0][2] = -2*qy*q0+2*qx*qz;
  d[1][0] = 2*qy*qx-2*qz*q0; d[1][1] =  qy*qy-qx*qx-qz*qz+q0*q0; d[1][2] =  2*qy*qz+2*qx*q0;
  d[2][0] = 2*qy*q0+2*qx*qz; d[2][1] =  2*qy*qz-2*qx*q0; d[2][2] = -qy*qy-qx*qx+qz*qz+q0*q0;

  pe = &rod[0].elems[0];

  rx -= rod[0].le*0.5*(d[0][0]*pe->s1 + d[1][0]*pe->s2 + d[2][0]*pe->s3); //half le
  ry -= rod[0].le*0.5*(d[0][1]*pe->s1 + d[1][1]*pe->s2 + d[2][1]*pe->s3);
  rz -= rod[0].le*0.5*(d[0][2]*pe->s1 + d[1][2]*pe->s2 + d[2][2]*pe->s3);

  rod->x0 = rx;  rod->y0 = ry;  rod->z0 = rz;

  // Calculate coordinates at L tip
  pn = &rod[0].nodes[rod->n_nodes-2];
  rx = pn->rx;  ry = pn->ry;  rz = pn->rz;

  q0 = 0.5*(rod->q0L + pn->q0);
  qx = 0.5*(rod->qxL + pn->qx);	 
  qy = 0.5*(rod->qyL + pn->qy);
  qz = 0.5*(rod->qzL + pn->qz);

  q  = sqrt( q0*q0 + qx*qx + qy*qy + qz*qz);
  q0 /= q; qx /= q; qy /= q;  qz /= q; 

  d[0][0] = -qy*qy+qx*qx-qz*qz+q0*q0; d[0][1] = 2*qy*qx+2*qz*q0; d[0][2] = -2*qy*q0+2*qx*qz;
  d[1][0] = 2*qy*qx-2*qz*q0; d[1][1] =  qy*qy-qx*qx-qz*qz+q0*q0; d[1][2] =  2*qy*qz+2*qx*q0;
  d[2][0] = 2*qy*q0+2*qx*qz; d[2][1] =  2*qy*qz-2*qx*q0; d[2][2] = -qy*qy-qx*qx+qz*qz+q0*q0;

  pe = &rod[0].elems[rod->n_elems-1];

  rx += rod[0].le*0.5*(d[0][0]*pe->s1 + d[1][0]*pe->s2 + d[2][0]*pe->s3); //half le
  ry += rod[0].le*0.5*(d[0][1]*pe->s1 + d[1][1]*pe->s2 + d[2][1]*pe->s3);
  rz += rod[0].le*0.5*(d[0][2]*pe->s1 + d[1][2]*pe->s2 + d[2][2]*pe->s3);

  rod->xL = rx;  rod->yL = ry;  rod->zL = rz;
}


void normal_modes(struct rod *rod, real dt, int t, int tprint, int n_rods, char *datadir)
{
  struct node *pn, *pn1, *pn2;
  int i, k; 

  real sxh[N], syh[N], szh[N], bxh[N], byh[N], bzh[N];                   // Strains in fourier space  
  real pxh[N+1], pyh[N+1], pzh[N+1], lxh[N+1], lyh[N+1], lzh[N+1];       // Momenta in fourier space  

  FILE *filep=0;
  char file1[64], file2[64], file3[64], file4[64], file5[64], file6[64], file7[64], file8[64], file9[64] ;
  real ds= rod->le;

  for(k=1; k<N+1; k++)  // k is the wavenumber, small k means large wavelength waves. k=0 is CoM motion.
    {
      pxh[k] =0;
      pyh[k] =0;
      pzh[k] =0;
      lxh[k] =0;
      lyh[k] =0;
      lzh[k] =0;
      
      for(i=1; i<N+1; i++)
	{
	  pn = &rod->nodes[i];
	  pxh[k] += pn->px*cos(k*i*Pi/N)*ds; 
	  pyh[k] += pn->py*cos(k*i*Pi/N)*ds;
	  pzh[k] += pn->pz*cos(k*i*Pi/N)*ds;
	  lxh[k] += pn->lx*cos(k*i*Pi/N)*ds; 
	  lyh[k] += pn->ly*cos(k*i*Pi/N)*ds;
	  lzh[k] += pn->lz*cos(k*i*Pi/N)*ds;
	}
    }

  for(k=0; k<N; k++)  // k=0 mode is overall extention of rod
    {
      sxh[k] =0;
      syh[k] =0;
      szh[k] =0;
      bxh[k] =0;
      byh[k] =0;
      bzh[k] =0;

      for(i=0; i<N; i++)
	{
	  pn1 = &rod->nodes[i];
	  pn2 = &rod->nodes[i+1];
	  if(N == 1)  // This is the only mode, corresponding to k=0, that can be observed for N=1
	    {
	      sxh[k] += ( (pn2->rx - pn1->rx)/ds );
	      syh[k] += ( (pn2->ry - pn1->ry)/ds );
	      szh[k] += ( (pn2->rz - pn1->rz)/ds - 1 );
	    } 
	  else
	    {
	      sxh[k] += ( (pn2->rx - pn1->rx)/ds )*cos(k*i*Pi/(N-1));
	      syh[k] += ( (pn2->ry - pn1->ry)/ds )*cos(k*i*Pi/(N-1));
	      szh[k] += ( (pn2->rz - pn1->rz)/ds - 1 )*cos(k*i*Pi/(N-1));
	    }
	}
    }
  
  if (t%tprint == 0)
    {
      strcpy(file1,datadir);
      strcat(file1,"/normalpx.ods");
      filep=fopen(file1,"a");
      fprintf (filep, "%f", t*dt );
      for(k=0; k<N+1; k++)
	fprintf (filep, "\t %.12e", pxh[k]*pxh[k] );
      fprintf (filep, "\n" );
      fclose(filep);

      strcpy(file2,datadir);
      strcat(file2,"/normalpy.ods");
      filep=fopen(file2,"a");
      fprintf (filep, "%f", t*dt );
      for(k=0; k<N+1; k++)
	fprintf (filep, "\t %.12e", pyh[k]*pyh[k] );
      fprintf (filep, "\n" );
      fclose(filep);

      strcpy(file3,datadir);
      strcat(file3,"/normalpz.ods");
      filep=fopen(file3,"a");
      fprintf (filep, "%f", t*dt );
      for(k=0; k<N+1; k++)
	fprintf (filep, "\t %.12e", pzh[k]*pzh[k] );
      fprintf (filep, "\n" );
      fclose(filep);

      strcpy(file4,datadir);
      strcat(file4,"/normallx.ods");
      filep=fopen(file4,"a");
      fprintf (filep, "%f", t*dt );
      for(k=0; k<N+1; k++)
	fprintf (filep, "\t %.12e", lxh[k]*lxh[k] );
      fprintf (filep, "\n" );
      fclose(filep);

      strcpy(file5,datadir);
      strcat(file5,"/normally.ods");
      filep=fopen(file5,"a");
      fprintf (filep, "%f", t*dt );
      for(k=0; k<N+1; k++)
	fprintf (filep, "\t %.12e", lyh[k]*lyh[k] );
      fprintf (filep, "\n" );
      fclose(filep);

      strcpy(file6,datadir);
      strcat(file6,"/normallz.ods");
      filep=fopen(file6,"a");
      fprintf (filep, "%f", t*dt );
      for(k=0; k<N+1; k++)
	fprintf (filep, "\t %.12e", lzh[k]*lzh[k] );
      fprintf (filep, "\n" );
      fclose(filep);

      /////////////////////////// PE /////////////
      strcpy(file7,datadir);
      strcat(file7,"/normalsx.ods");
      filep=fopen(file7,"a");
      fprintf (filep, "%f", t*dt );
      for(k=0; k<N; k++)
	fprintf (filep, "\t %.12e", sxh[k]*sxh[k]*ds/N);
      fprintf (filep, "\n" );
      fclose(filep);

      strcpy(file8,datadir);
      strcat(file8,"/normalsy.ods");
      filep=fopen(file8,"a");
      fprintf (filep, "%f", t*dt );
      for(k=0; k<N; k++)
	fprintf (filep, "\t %.12e", syh[k]*syh[k]*ds/N);
      fprintf (filep, "\n" );
      fclose(filep);

      strcpy(file9,datadir);
      strcat(file9,"/normalsz.ods");
      filep=fopen(file9,"a");
      fprintf (filep, "%f", t*dt );
      for(k=0; k<N; k++)
	fprintf (filep, "\t %.12e", szh[k]*szh[k]*ds/N);
      fprintf (filep, "\n" );
      fclose(filep);

    }
}
