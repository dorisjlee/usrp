/*============================================================================*/
/*! \file ppi.c
 * Sets up a 2D axisymmetric torus to look at papaloizou pringle instability
 * Will be further extended to 3D to look at non-axisymmetric, global modes
 * PRIVATE FUNCTION PROTOTYPES:
 * - shkset2d_iib() - sets BCs on L-x1 (left edge) of grid.
 * - shkset2d_oib() - sets BCs on R-x1 (right edge) of grid.
 * - shkset2d_ijb() - sets BCs on L-x2 (bottom edge) of grid.
 * - shkset2d_ojb() - sets BCs on R-x2 (top edge) of grid.		      */
/*============================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

void shkset2d_iib(GridS *pGrid);
void shkset2d_oib(GridS *pGrid);
void shkset2d_ijb(GridS *pGrid);
void shkset2d_ojb(GridS *pGrid);

/* Make size of box and dimension of unit cell (r1 x r2) static globals so they
 * can be accessed by boundary value functions */
static Real Lx,Ly;
static int r1,r2;

/* Analytic solution at stopping time, shared with Userwork_after_loop to
 * compute L1 error */
static ConsS ***RootSoln=NULL;

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int kl,ku,irefine,ir,ix1,ix2,nx1,nx2,nx3,gcd;
  Real sin_a2, cos_a2, sin_a3, cos_a3;

  Real tlim;
  int err_test;
  Real x1,x2,x3,r,xs,xc,xf,xh,vs,vc,vf,vh;
  Real xfp,xrp,xsp,xsm,xrm,xfm,vfp,vrp,vsp,vsm,vrm,vfm;
  Real d0,v0,Mx,My,Mz,E0,r0;

  Real rootdx1, rootdx2;
  Prim1DS Wl, Wr;//Vector of primitives (left and right states)
  Cons1DS Ul, Ur;//Vector of Conservatives 
  ConsS ql, qr;//Not too sure what this is? 
  Real Bxl=0.0,Bxr=0.0;
  div_t id;   /* structure containing remainder and quotient */

/* Following are used to compute volume of cell crossed by initial interface
 * that is assigned to left/right states */
  int dll, dlr, drr, drl;
  Real afl_lx, afr_lx, afl_rx, afr_rx;
  Real afl_ly, afr_ly, afl_ry, afr_ry;
  Real vfl, vfr, B1r, B2r;

/* Parse left state read from input file: dl,pl,ul,vl,wl,bxl,byl,bzl */

  Wl.d = par_getd("problem","dl");
  Wl.P = par_getd("problem","pl"); //Adiabatic
  Wl.Vx = par_getd("problem","v1l");
  Wl.Vy = par_getd("problem","v2l");
  Wl.Vz = par_getd("problem","v3l");

/* Parse right state read from input file: dr,pr,ur,vr,wr,bxr,byr,bzr */

  Wr.d = par_getd("problem","dr");
  Wr.P = par_getd("problem","pr");
  Wr.Vx = par_getd("problem","v1r");
  Wr.Vy = par_getd("problem","v2r");
  Wr.Vz = par_getd("problem","v3r");

  Ul = Prim1D_to_Cons1D(&Wl, &Bxl);
  Ur = Prim1D_to_Cons1D(&Wr, &Bxr);

/* Initialize ql rotated to the (x1,x2,x3) coordinate system */
  ql.d   = Ul.d;
  ql.M1  = Ul.Mx*cos_a3 - Ul.My*sin_a3;
  ql.M2  = Ul.Mx*sin_a3 + Ul.My*cos_a3;
  ql.M3  = Ul.Mz;
  ql.E   = Ul.E;

/* Initialize qr rotated to the (x1,x2,x3) coordinate system */
  qr.d   = Ur.d;
  qr.M1  = Ur.Mx*cos_a3 - Ur.My*sin_a3;
  qr.M2  = Ur.Mx*sin_a3 + Ur.My*cos_a3;
  qr.M3  = Ur.Mz;
  qr.E   = Ur.E;

/* Initialize the grid */

  for (k=kl; k<=ku; k++) {
    for (j=0; j<=je+nghost; j++) {
      ix2 = j + pGrid->Disp[1];
      for (i=0; i<=ie+nghost; i++) {
	ix1 = i + pGrid->Disp[0];

/* cell is completely in the left state */
	if((drr = r2*(ix1) + r1*(ix2) - gcd*r1*r2) <= 0){
	  pGrid->U[k][j][i] = ql;
	}
/* cell is completely in the right state */
	else if((dll = r2*(ix1-1) + r1*(ix2-1) - gcd*r1*r2) >= 0){
	  pGrid->U[k][j][i] = qr;
	}
/* The more complicated case of a cell  split by the interface boundary */
	else{
	  dlr = r2*(ix1-1) + r1*(ix2) - gcd*r1*r2;

	  if(dlr < 0){ /* The boundary hits the right y-face */
	    afl_lx = 1.0;
	    afr_lx = 0.0;
	    afl_ry = (Real)(-dlr)/(Real)(r2);
	    afr_ry = 1.0 - afl_ry;
	  }
	  else if(dlr > 0){ /* The boundary hits the left x-face */
	    afl_lx = (Real)(-dll)/(Real)(r1);
	    afr_lx = 1.0 - afl_lx;
	    afl_ry = 0.0;
	    afr_ry = 1.0;
	  }
	  else{ /* dlr == 0.0, The boundary hits the grid cell corner */
	    afl_lx = 1.0;
	    afr_lx = 0.0;
	    afl_ry = 0.0;
	    afr_ry = 1.0;
	  }

	  drl = r2*(ix1) + r1*(ix2-1) - gcd*r1*r2;

	  if(drl < 0){ /* The boundary hits the right x-face */
	    afl_rx = (Real)(-drl)/(Real)(r1);
	    afr_rx = 1.0 - afl_rx;
	    afl_ly = 1.0;
	    afr_ly = 0.0;
	  }
	  else if(drl > 0){ /* The boundary hits the left y-face */
	    afl_rx = 0.0;
	    afr_rx = 1.0;
	    afl_ly = (Real)(-dll)/(Real)(r2);
	    afr_ly = 1.0 - afl_ly;
	  }
	  else{ /* drl == 0.0, The boundary hits the grid cell corner */
	    afl_rx = 0.0;
	    afr_rx = 1.0;
	    afl_ly = 1.0;
	    afr_ly = 0.0;
	  }

/* The boundary hits both x-interfaces */
	  if(dlr > 0 && drl < 0){ 
	    vfl = 0.5*(afl_lx + afl_rx);
	    vfr = 1.0 - vfl;
	  }
/* The boundary hits both y-interfaces */
	  else if(dlr < 0 && drl > 0){ 
	    vfl = 0.5*(afl_ly + afl_ry);
	    vfr = 1.0 - vfl;
	  }
/* The boundary hits both grid cell corners */
	  else if(dlr == 0 && drl == 0){ 
	    vfl = vfr = 0.5;
	  }
/* The boundary hits the left x- and left y-interface */
	  else if(dlr > 0 && drl > 0){
	    vfl = 0.5*afl_lx*afl_ly;
	    vfr = 1.0 - vfl;
	  }
/* dlr<0 && drl<0:  The boundary hits the right x- and right y-interface */
	  else{ 
	    vfr = 0.5*afr_rx*afr_ry;
	    vfl = 1.0 - vfr;
	  }


/* Initialize the volume averaged quantities */
	  pGrid->U[k][j][i].d  = vfl*ql.d + vfr*qr.d;
	  pGrid->U[k][j][i].M1 = vfl*ql.M1 + vfr*qr.M1;
	  pGrid->U[k][j][i].M2 = vfl*ql.M2 + vfr*qr.M2;
	  pGrid->U[k][j][i].M3 = vfl*ql.M3 + vfr*qr.M3;
	  pGrid->U[k][j][i].E  = vfl*ql.E + vfr*qr.E;
	}
      }
    }
  }

/* Set boundary value function pointers */

 /* bvals_mhd_fun(pDomain, left_x1,  shkset2d_iib);
  bvals_mhd_fun(pDomain, left_x2,  shkset2d_ijb);
  bvals_mhd_fun(pDomain, right_x1, shkset2d_oib);
  bvals_mhd_fun(pDomain, right_x2, shkset2d_ojb);
*/
  return;
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
void problem_read_restart(MeshS *pM, FILE *fp)
{//Default restart procedures.
  Omega_0 = par_getd_def("problem","omega",1.0e-3);
  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
  dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
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
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

