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
static Real G,M,gm,n,q,w0,omega0,d,d0, A,rg; //constants

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real ***pr;
  Real f1l,f1r,f2l,f2r;
  Real vf=1.0; //Effective volume (set to 1 for now)
  Real r,theta,phi,v3,rho,omega,p,args; //Positional variables and physical quantities
  PrimS  W;//Vector of primitives 
  ConsS  U;//Vector of Conservatives 
  //Constants and defined parameters
  n  =  par_getd("problem","n");
  q  =  par_getd("problem","q");
  d0 = par_getd("problem","d0");
  omega0  =  1.0;//par_getd("problem","omega0");
  d= par_getd("problem","d"); 
  gm  = 1.+1./n;
  rg = par_getd("problem","rg");
  // Normalized to natural unit system
  G = 1.0;
  M = 1.0; //pointmass
  w0 = 1.0; //radius of max density (and P) 
  A =pow((d-1)/(2*d),n)/(n+1); //Choose rho_max = 1
  //Checking Constants
  printf("G,M,gm,n,q,w0,omega0,d,d0, A:\n  %2f,%2f,%2f,%2f,%2f,%2f,%2f,%2f,%2f,%2f",G,M,gm,n,q,w0,omega0,d,d0, A);
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	r= pGrid->px1v[i];
	theta = pGrid->px2v[j]; 
	phi = pGrid->px3[k];
        // set default values for ambient gas
        //W.d = 46.0;//0.0001;
        W.V1=0.0;
        W.V2=0.0;
        W.V3=0.0;
        //W.P=46.0/r;//0.0001/r;
        //p = d0/r;//1.0e-4/r;
	//d = d0;//1.0e-4;
	printf("d,P: %2f,%2f",d, p);
	W.d =d0;
	//W.P =43.0;//okay
	//W.P =0.0001;//Not Okay
	//W.P = p;
	W.P = d0/r;
  	//Axissymmetric initial condition that satisfies hydrostatic equilibrium
  	//W.V3 = omega0*pow(r,1-q); // velocity (in phi direction)
	v3 = omega0*pow(r,1-q);
	//W.V3 = v3;
	args = G*M/((n+1)*w0*A)*(w0/r-0.5*pow(w0/(r*sin(theta)),2)-0.5/d);
	if (args>1e-5){//Keeping only real solutions
	    //Inside the disk?
	    //printf("(r,theta,phi): %2f,%2f,%2f \n",r,theta,phi);
	    //printf("v3: %2f \n",v3); 
	    rho= pow(args,1.5);
            if(rho>W.d) {
   	      //printf("rho: %2f \n",rho);
              p  =A*pow(rho,gm);
              //printf("p : %2f \n",p);
	      W.V3 = v3;
	      W.d = rho;
              if(p>W.P)
	        W.P = p;
            }
	}/*else{
            W.d = 0.0001;//d0;
            W.V1=0.0;
            W.V2=0.0;
            W.V3=0.0;
            W.P=0.0001/r;//d0/r;
	}*/
	/*W.V3 = 1.0;
	W.d = 1.0;
	W.P = 1.0;*/
 	U = Prim_to_Cons(&W); //Convert to vector of conservatives
        pGrid->U[k][j][i].d  = vf*U.d;
        pGrid->U[k][j][i].M1 = 0.0;// momentum_r 0 since initial vr =0
        pGrid->U[k][j][i].M2 = 0.0;// p_theta =0 
        pGrid->U[k][j][i].M3 = vf*U.M3;
        pGrid->U[k][j][i].E  = vf*U.E;
      }
    }
  }
  bvals_mhd_fun(pDomain, left_x1,  sphoutflow_ix1);
  bvals_mhd_fun(pDomain, right_x1, sphoutflow_ox1);
  bvals_mhd_fun(pDomain, left_x2,  sphoutflow_ix2);
  bvals_mhd_fun(pDomain, right_x2, sphoutflow_ox2);
  x1GravAcc = grav_acc;
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
  return -G*M/SQR(x1-rg);
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


