#include "copyright.h"
/*============================================================================*/
/*! \file cylblast.c
 *  \brief Problem generator for blast wave in cylindrical coords.
 *
 * PURPOSE: Problem generator for blast wave in cylindrical coords.  Can only
 *   be run in 2D or 3D.  Input parameters are:
 *   -  problem/radius = radius of field initial overpressured region
 *   -  problem/pamb   = ambient pressure
 *   -  problem/prat   = ratio of interior to ambient pressure
 *   -  problem/b0     = initial azimuthal magnetic field (units sqrt(Pamb))
 *   -  problem/rho0   = background density
 *   -  problem/omega0 = initial azimuthal flow angular velocity
 *
 * REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
 *   multidimensional MHD", ApJ, 530, 508 (2000), and references therein. */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * grav_pot() - gravitational potential
 *============================================================================*/

static Real omega0,rho0;
static Real grav_pot(const Real x1, const Real x2, const Real x3);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pG = pDomain->Grid;
  int i,j,k;
  int is,ie,js,je,ks,ke,nx1,nx2,nx3;
  int il,iu,jl,ju,kl,ku;
  int nx,ny,nz; //Box dimensions
  Real dt,pt; //Values for the torus
  Real d0,v0,Mx,My,Mz,E0,r0,drat,p0; //Ambient Conditions
  Real c , a ,r ,theta,p,phi;
  Real phi0,x0,y0,z0,angle,radius,prat,pamb,b0;
  Real x1,x2,x3,x2i;
  Real x,y,z,Eint,Emag,Ekin;
  PrimS  W;//Vector of primitives 
  ConsS  U;//Vector of Conservatives 
//static Real grav_pot(const Real x1, const Real x2, const Real x3) i  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;


  int epsilon=5; //extra gridding near boundary
  Real vf =1.0; //Let the total volume be 1.0 for now (Boundary not determined)
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

  /* Read in initial conditions */
  radius = par_getd("problem","radius");
  pamb   = par_getd("problem","pamb");
  prat   = par_getd("problem","prat");
  rho0   = par_getd("problem","rho0");
  omega0 = par_getd("problem","omega0");
  b0     = par_getd("problem","b0");

  /* Placement of center of blast */
  r0   = par_getd("problem","r0");
  phi0 = par_getd("problem","phi0");
  z0   = par_getd("problem","z0");

  /* Orientation of field w.r.t. pos. x-axis */
  angle = (PI/180.0)*par_getd("problem","angle");

  x0 = r0*cos(phi0);
  y0 = r0*sin(phi0);

  /* Set up uniform ambient medium with circular over-pressured region */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3); //Computes the (r,theta,phi) position for a given (i,j,k) cell
        p = pow((c-x1),2)+pow(x3,2); // x1 = r ; x3 = z
	//x2i = x2 - 0.5*pG->dx2;
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
        U = Prim_to_Cons(&W); //Convert to vector of conservatives
        pG->U[k][j][i].d  = vf*U.d;
        pG->U[k][j][i].M1 = vf*U.M1;
        pG->U[k][j][i].M2 = vf*U.M2;
        pG->U[k][j][i].M3 = vf*U.M3;
        pG->U[k][j][i].E  = vf*U.E;

      }
    }
  }

  /* Enroll the gravitational function and radial BC */
  StaticGravPot = grav_pot;
  bvals_mhd_fun(pDomain,left_x1, do_nothing_bc);
  bvals_mhd_fun(pDomain,right_x1,do_nothing_bc);
  bvals_mhd_fun(pDomain,left_x2, do_nothing_bc);
  bvals_mhd_fun(pDomain,right_x2,do_nothing_bc);

  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
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

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3) {
 *  \brief  Gravitational potential*/
static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return 0.5*SQR(x1*omega0);
}

