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
/* function prototypes and global variables*/
Real grav_acc(const Real x1, const Real x2, const Real x3);
void sphoutflow_ix1(GridS *pG);
void sphoutflow_ox1(GridS *pG);
void sphoutflow_ix2(GridS *pG);
void sphoutflow_ox2(GridS *pG);
static Real gm, g, ptmass, en, cprime, beta, w0, rg, dist, acons, d0;

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real ***pr;
  Real f1l,f1r,f2l,f2r;
  Real amp;

  /* read parameters */
  gm = par_getd("problem","gamma");
  beta = par_getd("problem","beta");
  dist = par_getd("problem","dist");
  rg = par_getd("problem","rg");
  d0=par_getd("problem","d0");
  amp=par_getd("problem","amp");
  /* calculate some constant values */
  g = 1.0;
  ptmass = 1.0;
  w0 = 1.0;
  cprime = 0.5/dist;
  en = 1.0/(gm-1.0);
  acons=0.5*(dist-1.0)/dist/(en+1.0);

  /* assign boundary conditions and gravitational force function*/
  bvals_mhd_fun(pDomain, left_x1,  sphoutflow_ix1);
  bvals_mhd_fun(pDomain, right_x1, sphoutflow_ox1);
  bvals_mhd_fun(pDomain, left_x2,  sphoutflow_ix2);
  bvals_mhd_fun(pDomain, right_x2, sphoutflow_ox2);
  x1GravAcc = grav_acc;

  /* allocate memory for the gas pressure */
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
        /*** Initial Condition: replace this with a torus in equilibrium ***/
        /* <r>[i] = pGrid->px1v[i], <theta>[j] = pGrid->px2v[j], <phi>[k] = pGrid->px3[k] */
      	/* U[k][j][i].d : density, U.M1,2,3 : r,theta,phi-momentum, U.E: total energy, pr : pressure */
        pGrid->U[k][j][i].d  = d0; //fix this
        pGrid->U[k][j][i].M1 = 0.0; // OK
        pGrid->U[k][j][i].M2 = 0.0; // OK
      	pGrid->U[k][j][i].M3 = pGrid->U[k][j][i].d / (pGrid->px1v[i]*sin(pGrid->px2v[j])); // OK
        pr[k][j][i]=d0/pGrid->px1v[i]; // fix this
        /*** basically, that's all! ***/
      }
    }
  }

  /* Calculate the total energy including kinetic and magnetic energies*/
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].E  = pr[k][j][i]*en+0.5*SQR(pGrid->U[k][j][i].M3)/pGrid->U[k][j][i].d;
#ifdef MHD
        pGrid->U[k][j][i].B1c = ((pGrid->px1i[i+1]-pGrid->px1v[i])*pGrid->B1i[k][j][i] + (pGrid->px1v[i]-pGrid->px1i[i])*pGrid->B1i[k][j][i+1])/pGrid->dx1;
        pGrid->U[k][j][i].B2c = ((pGrid->px2i[j+1]-pGrid->px2v[j])*pGrid->B2i[k][j][i] + (pGrid->px2v[j]-pGrid->px2i[j])*pGrid->B2i[k][j+1][i])/pGrid->dx2;
        pGrid->U[k][j][i].B3c = (pGrid->B3i[k][j][i] + pGrid->B3i[k+1][j][i])*0.5;
        pGrid->U[k][j][i].E += 0.5*(SQR(pGrid->U[k][j][i].B1c)+SQR(pGrid->U[k][j][i].B2c)+SQR(pGrid->U[k][j][i].B3c));
#endif /* MHD */
      }
    }
  }
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
}

void Userwork_after_loop(MeshS *pM)
{
}


Real grav_acc(const Real x1, const Real x2, const Real x3) {
  return -ptmass*g/SQR(x1-rg);
}

/*  Boundary Condtions, outflowing, ix1, ox1, ix2, ox2  */
void sphoutflow_ix1(GridS *pG)
{
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
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
        pg = (pG->U[k][j][is-i+1].E-0.5*(SQR(pG->U[k][j][is-i+1].M1)+SQR(pG->U[k][j][is-i+1].M2)+SQR(pG->U[k][j][is-i+1].M3))/pG->U[k][j][is-i+1].d)*(gm-1.0);
        pg+=grav_acc(pG->px1i[is-i+1],pG->px2v[j],pG->px3v[k])*pG->U[k][j][is-i].d*pG->dx1;
        pG->U[k][j][is-i].E=+pg/(gm-1.0)+0.5*SQR(pG->U[k][j][is-i].M1)/pG->U[k][j][is-i].d;
      }
    }
  }
  return;
}


void sphoutflow_ox1(GridS *pG)
{
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

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

  return;
}


void sphoutflow_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

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

  return;
}


void sphoutflow_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

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
  return;
}


