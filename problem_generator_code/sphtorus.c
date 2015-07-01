#include "copyright.h"
/*============================================================================*/
/*! \file sphtorus.c
 *  \brief Problem generator for the torus problem (Stone et al. 1999)
 *
 * PURPOSE: Problem generator for the torus problem (Stone et al. 1999)
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef ISOTHERMAL
#error "Isothermal EOS cannot be used."
#endif

/*----------------------------------------------------------------------------*/
/* problem:  */
Real grav_pot(const Real x1, const Real x2, const Real x3);
Real grav_acc(const Real x1, const Real x2, const Real x3);
void sphoutflow_ix1(GridS *pG);
void sphoutflow_ox1(GridS *pG);
void sphoutflow_ix2(GridS *pG);
void sphoutflow_ox2(GridS *pG);

static Real gm, g, ptmass, en, cprime, beta, w0, rg, dist, acons, d0, rie,roe;

#ifdef MHD
Real A1(const Real x1, const Real x2, const Real x3)
{
  return 0.0;
}

Real A2(const Real x1, const Real x2, const Real x3)
{
  return 0.0;
}

Real A3(const Real x1, const Real x2, const Real x3)
{
  Real w,a=0.0,eq29,dens;
  w=x1*sin(x2);
  eq29 = (g*ptmass)/(w0*(en + 1.))*(w0/x1-0.5*SQR(w0/w) - cprime);
  if (eq29 > 0.0) {
    dens  = pow(eq29/acons,en);
    if (dens > 100.0*d0)
      a = SQR(dens)/(beta);
  }
  return a;
}

#define ND 100

Real magr(const GridS *pG, const int i, const int j, const int k)
{
  Real r,t,p,s,a,d,rd;
  int n;
  r=pG->px1i[i];
  t=pG->px2i[j];
  p=pG->px3i[k];
  s=2.0*SQR(r)*sin(pG->px2[j])*sin(0.5*pG->dx2)*pG->dx3;
/*  d=pG->dx2/(Real)ND;
  rd=r*d;
  a=0.5*(A2(r,t,p)+A2(r,t+pG->dx2,p))*rd-0.5*(A2(r,t,p+pG->dx3)+A2(r,t+pG->dx2,p+pG->dx3))*rd;
  for(n=1;n<ND;n++)
    a+=A2(r,t+((Real)n)*d,p)*rd-A2(r,t+d,p+pG->dx3)*rd;*/
  d=pG->dx3/(Real)ND;
  rd=r*d;
  a=0.5*(A3(r,t+pG->dx2,p)+A3(r,t+pG->dx2,p+pG->dx3))*rd*sin(t+pG->dx2)-0.5*(A3(r,t,p)+A3(r,t,p+pG->dx3))*rd*sin(t);
  for(n=1;n<ND;n++)
    a+=A3(r,t+pG->dx2,p+((Real)n)*d)*rd*sin(t+pG->dx2)-A3(r,t,p+((Real)n)*d)*rd*sin(t);
  return a/s;
}

Real magt(const GridS *pG, const int i, const int j, const int k)
{
  Real r,t,p,s,a,d,rd;
  int n;
  r=pG->px1i[i];
  t=pG->px2i[j];
  p=pG->px3i[k];
  s=pG->px1[i]*pG->dx1*sin(t)*pG->dx3;
/*  d=pG->dx1/(Real)ND;
  a=0.5*(A1(r,t,p)+A1(r+pG->dx1,t,p))*d-0.5*(A1(r,t,p+pG->dx3)+A1(r+pG->dx1,t,p+pG->dx3))*d;
  for(n=1;n<ND;n++)
    a+=A1(r+((Real)n)*d,t,p)*d-A1(r+((Real)n)*d,t,p+pG->dx3)*d;*/
  d=pG->dx3/(Real)ND;
  rd=sin(t)*d;
  a=0.5*(A3(r+pG->dx1,t,p)+A3(r+pG->dx1,t,p+pG->dx3))*rd*(r+pG->dx1)-0.5*(A3(r,t,p)+A3(r,t,p+pG->dx3))*rd*r;
  for(n=1;n<ND;n++)
    a+=A3(r+pG->dx1,t,p+((Real)n)*d)*rd*(r+pG->dx1)-A3(r,t,p+((Real)n)*d)*rd*r;
  return -a/s;
}

Real magp(const GridS *pG, const int i, const int j, const int k)
{
  return 0.0;
}
#endif

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real ***pr;
  Real f1l,f1r,f2l,f2r;
  Real eq29,dens;
  Real r, w, z;
  Real divb, divbmax;
  Real ld, lm, lv, vv, rp, rm, tp, tm, tv, rv, vt, pp, amp;
  int idb, jdb, kdb, ierr;
  int ii, jj, kk, ftorus;

  srand(myID_Comm_world);
  
  gm = par_getd("problem","gamma");
  //beta = par_getd("problem","beta");
  dist = par_getd("problem","dist");
  rg = par_getd("problem","rg");
  d0=par_getd("problem","d0");
  amp=par_getd("problem","amp");
  g = 1.0;
  ptmass = 1.0;
  w0 = 1.0;
  cprime = 0.5/dist;
  en = 1.0/(gm-1.0);
  acons=0.5*(dist-1.0)/dist/(en+1.0);
  rie=pDomain->RootMinX[0]+0.5*pDomain->dx[0];
  roe=pDomain->RootMaxX[0]-0.5*pDomain->dx[0];

  bvals_mhd_fun(pDomain, left_x1,  sphoutflow_ix1);
  bvals_mhd_fun(pDomain, right_x1, sphoutflow_ox1);
  bvals_mhd_fun(pDomain, left_x2,  sphoutflow_ix2);
  bvals_mhd_fun(pDomain, right_x2, sphoutflow_ox2);

//  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  pr=(Real***)calloc_3d_array(pGrid->Nx[2]+2*nghost, pGrid->Nx[1]+2*nghost, pGrid->Nx[0]+2*nghost, sizeof(Real));

  /* Background */ 
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].d  = d0;
        pGrid->U[k][j][i].M1 = 0.0;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;
        pr[k][j][i]=d0/pGrid->px1v[i];
      }
    }
  }
  /* Torus */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        ld=0.0; lm=0.0;
        ftorus=0;
        vt=0.0;
        for(jj=0;jj<10;jj++) {
          tp=pGrid->px2i[j]+pGrid->dx2*0.1*(jj+1);
          tm=pGrid->px2i[j]+pGrid->dx2*0.1*jj;
          tv=0.5*(tp+tm)+(1.0-0.5*(tp-tm)/tan(0.5*(tp-tm)))/tan(0.5*(tp+tm));
          for(ii=0;ii<10;ii++) {
            rp = pGrid->px1i[i]+pGrid->dx1*0.1*(ii+1);
            rm = pGrid->px1i[i]+pGrid->dx1*0.1*ii;
            rv= ((SQR(SQR(rp))-SQR(SQR(rm)))/4.0)/((CUBE(rp)-CUBE(rm))/3.0);
            lv= 1.0/3.0*(CUBE(rp)-CUBE(rm))*(cos(tm)-cos(tp));
            vt+=lv;
            w = rv*sin(tv);
            eq29 = (g*ptmass)/(w0*(en + 1.))*(w0/rv-0.5*SQR(w0/w) - cprime);
            if (eq29 > 0.0) {
              dens = pow(eq29/acons,en);
              pp=dens*eq29;
              if (pp > d0/rv) {
                ftorus=1;
                ld+=dens*lv;
                lm+=dens*lv*sqrt(g*ptmass*w0)/w;
              }
              else
                ld+=d0*lv;
            }
            else
              ld+=d0*lv;
          }
        }
        if(ftorus==1)
        {
          vv= 1.0/3.0*(CUBE(pGrid->px1i[i+1])-CUBE(pGrid->px1i[i]))*(cos(pGrid->px2i[j])-cos(pGrid->px2i[j+1]));
          ld/=vv; lm/=vv;
          pGrid->U[k][j][i].d = ld;
          pGrid->U[k][j][i].M3 = lm;
//          w=pGrid->px1v[i]*sin(pGrid->px2v[j]);
//          pGrid->U[k][j][i].M3 = ld*sqrt(g*ptmass*w0)/w;
          pr[k][j][i]=MAX(acons*pow(ld,gm),d0/pGrid->px1v[i])*(1+amp*((double)rand()/(double)RAND_MAX-0.5));
        }
      }
    }
  }


#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
        if(i!=pGrid->sil)
          pGrid->B1i[k][j][i] = magr(pGrid,i,j,k);
      }
      if(pGrid->sil>=is && pGrid->sil <= ie+1)
        pGrid->B1i[k][j][pGrid->sil] = pGrid->B1i[k][j][pGrid->sil+1];
    }
  }
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      if(j!=pGrid->sjl && j!=pGrid->sju)
      {
        for (i=is; i<=ie; i++) {
          pGrid->B2i[k][j][i] = magt(pGrid,i,j,k);
        }
      }
    }
    if(pGrid->sjl>=js && pGrid->sjl<=je+1)
    {
      for (i=is; i<=ie; i++) {
        pGrid->B2i[k][pGrid->sjl][i] = pGrid->B2i[k][pGrid->sjl+1][i];
      }
    }
    if(pGrid->sju>=js && pGrid->sju<=je+1)
    {
      for (i=is; i<=ie; i++) {
        pGrid->B2i[k][pGrid->sju][i] = pGrid->B2i[k][pGrid->sju-1][i];
      }
    }
  }
  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
      	pGrid->B3i[k][j][i] = 0.0;
      }
    }
  }
#endif

  divbmax=0.0; idb=is; jdb=js; kdb=ks;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].E  = pr[k][j][i]*en+0.5*SQR(pGrid->U[k][j][i].M3)/pGrid->U[k][j][i].d;
#ifdef MHD
        divb=fabs(((SQR(pGrid->px1i[i+1])*pGrid->B1i[k][j][i+1]-SQR(pGrid->px1i[i])*pGrid->B1i[k][j][i])*2.0*sin(pGrid->px2[j])*sin(0.5*pGrid->dx2)*pGrid->dx3
            +(sin(pGrid->px2i[j+1])*pGrid->B2i[k][j+1][i]-sin(pGrid->px2i[j])*pGrid->B2i[k][j][i])*pGrid->px1[i]*pGrid->dx1*pGrid->dx3
            +(pGrid->B3i[k+1][j][i]-pGrid->B3i[k][j][i])*pGrid->px1[i]*pGrid->dx1*pGrid->dx2)
            /(2.0/3.0*(CUBE(pGrid->px1i[i+1])-CUBE(pGrid->px1i[i]))*sin(pGrid->px2[j])*sin(0.5*pGrid->dx2)*pGrid->dx3));
        if(divb>divbmax)
        {
          idb=i, jdb=j, kdb=k;
          divbmax=divb;
        }
        pGrid->U[k][j][i].B1c = ((pGrid->px1i[i+1]-pGrid->px1v[i])*pGrid->B1i[k][j][i] + (pGrid->px1v[i]-pGrid->px1i[i])*pGrid->B1i[k][j][i+1])/pGrid->dx1;
        pGrid->U[k][j][i].B2c = ((pGrid->px2i[j+1]-pGrid->px2v[j])*pGrid->B2i[k][j][i] + (pGrid->px2v[j]-pGrid->px2i[j])*pGrid->B2i[k][j+1][i])/pGrid->dx2;
        pGrid->U[k][j][i].B3c = (pGrid->B3i[k][j][i] + pGrid->B3i[k+1][j][i])*0.5;
        pGrid->U[k][j][i].E += 0.5*(SQR(pGrid->U[k][j][i].B1c)+SQR(pGrid->U[k][j][i].B2c)+SQR(pGrid->U[k][j][i].B3c));
#endif /* MHD */
      }
    }
  }
#ifdef MHD
//  printf("divbmax=%g at i=%d j=%d k=%d, r=%g theta=%g phi=%g\n", divbmax, idb, jdb, kdb, pGrid->px1v[idb], pGrid->px2v[jdb], pGrid->px3v[kdb]);
#endif /* MHD */
  free_3d_array(pr);
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  DomainS *pDomain=(DomainS*)&(pM->Domain[0][0]);

  gm = par_getd("problem","gamma");
  rg = par_getd("problem","rg");
  g = 1.0;
  ptmass = 1.0;
  rie=pDomain->RootMinX[0]+0.5*pDomain->dx[0];
  roe=pDomain->RootMaxX[0]-0.5*pDomain->dx[0];

  bvals_mhd_fun(pDomain, left_x1,  sphoutflow_ix1);
  bvals_mhd_fun(pDomain, right_x1, sphoutflow_ox1);
  bvals_mhd_fun(pDomain, left_x2,  sphoutflow_ix2);
  bvals_mhd_fun(pDomain, right_x2, sphoutflow_ox2);

//  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}


void Userwork_in_loop(MeshS *pM)
{
  int i, j, k, n, m, l;
  GridS *pG;
  DomainS *pD;
  Real em, ek, eg, cd;
#ifdef MHD
  pD = (DomainS*)&(pM->Domain[0][0]);
  pG = pD->Grid;
  for (k=pG->ks; k<=pG->ke; k++) {
    for (j=pG->js; j<=pG->je; j++) {
      for (i=pG->is; i<=pG->ie; i++) {
        if(pG->U[k][j][i].d < 1e-5)
        {
          cd=pG->U[k][j][i].d;
          em=0.5*(SQR(pG->U[k][j][i].B1c)+SQR(pG->U[k][j][i].B2c)+SQR(pG->U[k][j][i].B3c));
          ek=0.5*(SQR(pG->U[k][j][i].M1)+SQR(pG->U[k][j][i].M2)+SQR(pG->U[k][j][i].M3))/cd;
          eg=pG->U[k][j][i].E-ek-em;
          pG->U[k][j][i].M1 *= 1e-5/cd;
          pG->U[k][j][i].M2 *= 1e-5/cd;
          pG->U[k][j][i].M3 *= 1e-5/cd;
          eg *= 1e-5/cd;
          ek *= 1e-5/cd;
          pG->U[k][j][i].d = 1e-5;
          pG->U[k][j][i].E = em+eg+ek;
        }
      }
    }
  }
#endif
}

void Userwork_after_loop(MeshS *pM)
{
}

Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return -ptmass*g/(MAX(rie,MIN(x1,roe))-rg);
}


Real grav_acc(const Real x1, const Real x2, const Real x3) {
//  if(x1 < rie || x1 > roe) return 0.0;
  return -ptmass*g/SQR(x1-rg);
}

/*  Boundary Condtions, outflowing, ix1, ox1, ix2, ox2  */
void sphoutflow_ix1(GridS *pG)
{
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
  Real pg;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pG->U[k][j][is-i] = pG->U[k][j][is];
/*      if(pG->U[k][j][is-i].M1 > 0.0)
        {
          pG->U[k][j][is-i].E -= 0.5*SQR(pG->U[k][j][is-i].M1)/pG->U[k][j][is-i].d;
          pG->U[k][j][is-i].M1 = 0.0;
        }*/
        pG->U[k][j][is-i].M1 = 0.0;
        pG->U[k][j][is-i].M2 = 0.0;
        pG->U[k][j][is-i].M3 = 0.0;
#ifdef MHD
        pg = (pG->U[k][j][is-i+1].E-0.5*(SQR(pG->U[k][j][is-i+1].M1)+SQR(pG->U[k][j][is-i+1].M2)+SQR(pG->U[k][j][is-i+1].M3))/pG->U[k][j][is-i+1].d-0.5*(SQR(pG->U[k][j][is-i+1].B1c)+SQR(pG->U[k][j][is-i+1].B2c)+SQR(pG->U[k][j][is-i+1].B3c)))*(gm-1.0);
#else 
        pg = (pG->U[k][j][is-i+1].E-0.5*(SQR(pG->U[k][j][is-i+1].M1)+SQR(pG->U[k][j][is-i+1].M2)+SQR(pG->U[k][j][is-i+1].M3))/pG->U[k][j][is-i+1].d)*(gm-1.0);
#endif
        pg+=grav_acc(pG->px1i[is-i+1],pG->px2v[j],pG->px3v[k])*pG->U[k][j][is-i].d*pG->dx1;
        pG->U[k][j][is-i].E=+pg/(gm-1.0)+0.5*SQR(pG->U[k][j][is-i].M1)/pG->U[k][j][is-i].d;
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost-1; i++) {
        pG->B1i[k][j][is-i] = pG->B1i[k][j][is];
      }
    }
  }

  if (pG->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pG->B2i[k][j][is-i] = pG->B2i[k][j][is];
      }
    }
  }

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pG->B3i[k][j][is-i] = pG->B3i[k][j][is];
      }
    }
  }
#endif /* MHD */

  return;
}


void sphoutflow_ox1(GridS *pG)
{
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pG->U[k][j][ie+i] = pG->U[k][j][ie];
        if(pG->U[k][j][ie+i].M1 < 0.0)
        {
          pG->U[k][j][ie+i].E -= 0.5*SQR(pG->U[k][j][ie+i].M1)/pG->U[k][j][ie+i].d;
          pG->U[k][j][ie+i].M1 = 0.0;
        }
      }
    }
  }

#ifdef MHD
/* i=ie+1 is not a boundary condition for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=2; i<=nghost; i++) {
        pG->B1i[k][j][ie+i] = pG->B1i[k][j][ie];
      }
    }
  }

  if (pG->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pG->B2i[k][j][ie+i] = pG->B2i[k][j][ie];
      }
    }
  }

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pG->B3i[k][j][ie+i] = pG->B3i[k][j][ie];
      }
    }
  }
#endif /* MHD */

  return;
}


void sphoutflow_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[k][js-j][i] = pG->U[k][js][i];
        if(pG->U[k][js-j][i].M2 > 0.0)
        {
          pG->U[k][js-j][i].E -= 0.5*SQR(pG->U[k][js-j][i].M2)/pG->U[k][js-j][i].d;
          pG->U[k][js-j][i].M2 = 0.0;
        }
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pG->B1i[k][js-j][i] = pG->B1i[k][js][i];
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->B2i[k][js-j][i] = pG->B2i[k][js][i];
      }
    }
  }

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->B3i[k][js-j][i] = pG->B3i[k][js][i];
      }
    }
  }
#endif /* MHD */

  return;
}


void sphoutflow_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[k][je+j][i] = pG->U[k][je][i];
        if(pG->U[k][je+j][i].M2 < 0.0)
        {
          pG->U[k][je+j][i].E -= 0.5*SQR(pG->U[k][je+j][i].M2)/pG->U[k][je+j][i].d;
          pG->U[k][je+j][i].M2 = 0.0;
        }
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        pG->B1i[k][je+j][i] = pG->B1i[k][je][i];
      }
    }
  }

/* j=je+1 is not a boundary condition for the interface field B2i */
  for (k=ks; k<=ke; k++) {
    for (j=2; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->B2i[k][je+j][i] = pG->B2i[k][je][i];
      }
    }
  }

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->B3i[k][je+j][i] = pG->B3i[k][je][i];
      }
    }
  }
#endif /* MHD */

  return;
}


