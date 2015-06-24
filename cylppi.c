#include "copyright.h"
/*============================================================================*/
/*! \file ppi.c
 * \brief Sets up a 2D axisymmetric torus to look at papaloizou pringle instability. Will be further extended to 3D to look at non-axisymmetric, global modes */
/*============================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
/*============================================================================* 
 *  * PRIVATE FUNCTION PROTOTYPES:
 *   *============================================================================*/

static Real Lx,Ly;
static int r1,r2;

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int kl,ku,irefine,ir,ix1,ix2;
  int nx,ny,nz; //Box dimensions
  Real dt,pt; //Values for the torus
  Real d0,v0,Mx,My,Mz,E0,r0,drat,p0; //Ambient Conditions
  Real c , a ,r ,theta,p,phi;
  double x,y,z;
  PrimS  W;//Vector of primitives 
  ConsS  U;//Vector of Conservatives 
  int epsilon=5; //extra gridding near boundary
/* Following are used to compute volume of cell crossed by initial interface
 * that is assigned to left/right states */
  Real vf;
  vf=1.0; //Let the total volume be 1.0 for now (Boundary not determined)
  // Reading values from input file: d,p,v1,v2,v3
  dt = par_getd("problem","d");
  pt = par_getd("problem","p");
  W.V1 = par_getd("problem","v1");
  W.V2 = par_getd("problem","v2");
  W.V3 = par_getd("problem","v3");
  drat = par_getd("problem","drat"); //Density ratio of the cloud (used to determine ambient d and p)
  c = par_getd("problem","c"); //Distance from torus center to rotation axis
  a = par_getd("problem","a"); //Radius of the Torus cross section
  nx= par_getd("domain1","Nx1");
  ny= par_getd("domain1","Nx2");
  nz= par_getd("domain1","Nx3");
  d0=dt/drat;
  p0=pt/drat;
  printf("c: %2f",c);
  printf("a: %2f",a);  

  //Initializes the 2D grid (k not really considered here?)
  for (k=0; k<=ke+nghost; k++) {//from the edge (ke) plus boundary padding ghost cells
    for (j=0; j<=je+nghost; j++) {
      ix2 = j + pGrid->Disp[1];
      for (i=0; i<=ie+nghost; i++) {
	  ix1 = i + pGrid->Disp[0];
 	  /*printf("(i,j) : (%d,%d)   ",i,j);
	  printf("(i1,ix2) : (%d,%d)   ",ix1,ix2);*/
	  x = (double)(i-((double)nx/2.));
	  y = (double)(j-((double)ny/2.));
	  //z = (double)(k-32); //Maybe 80/2 = 40 shift is more appropriate
	  z = (double)(k-((double)nz/2.)); 
	  //printf("(x,y,z) : (%2f,%2f,%2f) ",x,y,z);
	  r = sqrt(pow(x,2)+pow(y,2)+pow(z,2)); //Pythagorean radii
	  theta = atan(y/x); //polar angle
	  phi = acos(z/r); //azimuthal angle
	  p = c*c - 2*c*sqrt(x*x+y*y) + r*r;
          /*printf("r: %2f		",r);
	  printf("phi: %2f    ",phi);	
	  printf("theta: %2f	",theta);
	  printf("p: %2f	 \n",p);*/
	  
	  if (p<=a*a){
	    //Inside the torus
	    //printf("inside!");
	    W.d = dt;
	    W.P = pt;
	  }
	  else{
	    //Outside the torus
	    //printf("outside!");
	    W.d = d0;
	    W.P = p0;
	  }
	  //W.d = dt;
	  //W.P =pt;
	  U = Prim_to_Cons(&W); //Convert to vector of conservatives
	  pGrid->U[k][j][i].d  = vf*U.d;
	  pGrid->U[k][j][i].M1 = vf*U.M1;
	  pGrid->U[k][j][i].M2 = vf*U.M2;
	  pGrid->U[k][j][i].M3 = vf*U.M3;
	  pGrid->U[k][j][i].E  = vf*U.E;
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
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

