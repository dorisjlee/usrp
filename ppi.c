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

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int kl,ku,irefine,ir,ix1,ix2;

  Real x1,x2,x3,r,xs,xc,xf,xh,vs,vc,vf,vh;
  Real xfp,xrp,xsp,xsm,xrm,xfm,vfp,vrp,vsp,vsm,vrm,vfm;
  Real d0,v0,Mx,My,Mz,E0,r0;

  Prim1DS  W;//Vector of primitives (left and right states)
  Cons1DS  U;//Vector of Conservatives 
  ConsS  q;//Not too sure what this is? 

/* Following are used to compute volume of cell crossed by initial interface
 * that is assigned to left/right states */
/* int dll, dlr, drr, drl;
 * Real afl_lx, afr_lx, afl_rx, afr_rx;
 * Real afl_ly, afr_ly, afl_ry, afr_ry;
 * Real vfl, vfr, B1r, B2r; */
  Real vf;
  //Let the total volume be 1.0 for now (Boundary not determined)
  vf=1.0;

  // Reading values from input file: d,p,v1,v2,v3
  W.d = par_getd("problem","d");
  W.P = par_getd("problem","p");
  W.Vx = par_getd("problem","v1");
  W.Vy = par_getd("problem","v2");
  W.Vz = par_getd("problem","v3");

  U = Prim1D_to_Cons1D(&W);

  q.d   = U.d;
  q.M1  = U.Mx;
  q.M2  = U.My;
  q.M3  = U.Mz;
  q.E   = U.E;

  //Initializes the 2D grid (k not really considered here?)
  for (k=kl; k<=ku; k++) {
    for (j=0; j<=je+nghost; j++) {
      ix2 = j + pGrid->Disp[1];
      for (i=0; i<=ie+nghost; i++) {
	  ix1 = i + pGrid->Disp[0];
	  pGrid->U[k][j][i].d  = vf*q.d;
	  pGrid->U[k][j][i].M1 = vf*q.M1;
	  pGrid->U[k][j][i].M2 = vf*q.M2;
	  pGrid->U[k][j][i].M3 = vf*q.M3;
	  pGrid->U[k][j][i].E  = vf*q.E;
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

