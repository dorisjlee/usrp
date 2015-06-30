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

#ifdef MHD
/* magnetic field functions */
Real magr(const GridS *pG, const int i, const int j, const int k)
{
  /* B1: r-component of the magnetic field */
  return 0.0;
}

Real magt(const GridS *pG, const int i, const int j, const int k)
{
  /* B2: theta-component of the magnetic field */
  return 0.0;
}

Real magp(const GridS *pG, const int i, const int j, const int k)
{
  /* B3: phi-component of the magnetic field */
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
  Real amp;

  /* initialize a random seed */
  srand(myID_Comm_world);
  
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


#ifdef MHD
  /* Calculate magnetic fields */
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
  DomainS *pDomain=(DomainS*)&(pM->Domain[0][0]);

  gm = par_getd("problem","gamma");
  rg = par_getd("problem","rg");
  g = 1.0;
  ptmass = 1.0;

  bvals_mhd_fun(pDomain, left_x1,  sphoutflow_ix1);
  bvals_mhd_fun(pDomain, right_x1, sphoutflow_ox1);
  bvals_mhd_fun(pDomain, left_x2,  sphoutflow_ix2);
  bvals_mhd_fun(pDomain, right_x2, sphoutflow_ox2);

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
  /* Density Floor */
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


