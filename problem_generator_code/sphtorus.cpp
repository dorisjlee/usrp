/*============================================================================*/
/*! \file sphtorus.cpp
 *  \brief Problem generator for the torus problem (Stone et al. 1999)
 *
 * PURPOSE: Problem generator for the torus problem (Stone et al. 1999)
/*============================================================================*/
// Primary header
#include "../fluid/fluid.hpp"

// C++ headers
#include <iostream>   // endl
#include <cmath>      // sqrt
#include <cstdlib>    // srand

// Athena headers
#include "../athena.hpp"           // enums, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../mesh.hpp"             // MeshBlock
#include "../parameter_input.hpp"  // ParameterInput
#include "../fluid/eos/eos.hpp"    // ParameterInput
#include "../bvals/bvals.hpp" // EnrollFluidBValFunction
#include "../coordinates/coordinates.hpp" // Coordinates


#ifdef ISOTHERMAL
#error "Isothermal EOS cannot be used."
#endif

static Real gm,gmgas;


inline Real CUBE(Real x){
  return ((x)*(x)*(x));
}

inline Real MAX(Real x, Real y){
	return (x>y?x:y);
}

/*----------------------------------------------------------------------------*/
/* function prototypes and global variables*/
void stbv_iib(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke); //sets BCs on inner-x1 (left edge) of grid.
void stbv_ijb(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke); //sets BCs on inner-x2 (bottom edge) of grid.
void stbv_oib(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke); //sets BCs on outer-x1 (right edge) of grid.
void stbv_ojb(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke); //sets BCs on outer-x2 (top edge) of grid.
#ifdef MHD
  Real A1(  Real x1,   Real x2,   Real x3);
  Real A2(  Real x1,   Real x2,   Real x3);
  Real A3(  Real x1,   Real x2,   Real x3);
  Real magr(MeshBlock *pmb,   int i,   int j,   int k);
  Real magt(MeshBlock *pmb,   int i,   int j,   int k);
  Real magp(MeshBlock *pmb,   int i,   int j,   int k);
  Real A1(Real x1, Real x2,Real x3)
  {
    return 0.0;
  }

  Real A2(Real x1, Real x2,Real x3) 
  {
    return 0.0;
  }

  Real A3(Real x1, Real x2,Real x3)
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

  Real magr(MeshBlock *pmb,   int i,   int j,   int k)
  {
    Real r,t,p,s,a,d,rd;
    Coordinates *pco = pmb->pcoord;
    int n;
    r = pco->x1f(i);
    t = pco->x2f(j);
    p = pco->x3f(k);
    s=2.0*SQR(r)*sin(t)*sin(0.5*t)*pco->dx3(k);
  /*  d=pmb->dx2/(Real)ND;
    rd=r*d;
    a=0.5*(A2(r,t,p)+A2(r,t+pmb->dx2,p))*rd-0.5*(A2(r,t,p+pmb->dx3)+A2(r,t+pmb->dx2,p+pmb->dx3))*rd;
    for(n=1;n<ND;n++)
      a+=A2(r,t+((Real)n)*d,p)*rd-A2(r,t+d,p+pmb->dx3)*rd;*/
    d=pco->dx3(k)/(Real)ND;
    rd=r*d;
    a=0.5*(A3(r,t+pco->dx2(j),p)+A3(r,t+pco->dx2(j),p+pco->dx3(k)))*rd*sin(t+pco->dx2(j))-0.5*(A3(r,t,p)+A3(r,t,p+pco->dx3(k)))*rd*sin(t);
    for(n=1;n<ND;n++)
      a+=A3(r,t+pco->dx2(i),p+((Real)n)*d)*rd*sin(t+pco->dx2(j))-A3(r,t,p+((Real)n)*d)*rd*sin(t);
    return a/s;
  }

  Real magt(MeshBlock *pmb, int i, int j, int k)
  {
    Coordinates *pco = pmb->pcoord;
    Real r,t,p,s,a,d,rd;
    int n;
    r = pco->x1f(i);
    t = pco->x2f(j);
    p = pco->x3f(k);
    s=r*pco->dx1(i)*sin(t)*pco->dx3(k);
  /*  d=pmb->dx1/(Real)ND;
    a=0.5*(A1(r,t,p)+A1(r+pmb->dx1,t,p))*d-0.5*(A1(r,t,p+pmb->dx3)+A1(r+pmb->dx1,t,p+pmb->dx3))*d;
    for(n=1;n<ND;n++)
      a+=A1(r+((Real)n)*d,t,p)*d-A1(r+((Real)n)*d,t,p+pmb->dx3)*d;*/
    d=pco->dx3(k)/(Real)ND;
    rd=sin(t)*d;
    a=0.5*(A3(r+pco->dx1(i),t,p)+A3(r+pco->dx1(i),t,p+pco->dx3(k)))*rd*(r+pco->dx1(i))-0.5*(A3(r,t,p)+A3(r,t,p+pco->dx3(k)))*rd*r;
    for(n=1;n<ND;n++)
      a+=A3(r+pco->dx1(i),t,p+((Real)n)*d)*rd*(r+pco->dx1(i))-A3(r,t,p+((Real)n)*d)*rd*r;
    return -a/s;
  }

  Real magp(MeshBlock *pmb, int i, int j, int k)
  {
    return 0.0;
  }
#endif

//problem generator
void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  MeshBlock *pmb = pfl->pmy_block;
  Coordinates *pco = pmb->pcoord;
  int i, is = pmb->is, ie = pmb->ie;
  int j, js = pmb->js, je = pmb->je;
  int k, ks = pmb->ks, ke = pmb->ke;
  int ii,jj,ftorus;
  Real en, cprime, w0, rg, dist, acons, d0, amp,beta,divbmax,idb,jdb,kdb;
  Real ld, lm, lv, vv, rp, rm, tp, tm, tv, rv, vt, pp,eq29,dens,wt;
  AthenaArray<Real> pr;
  // ID_t rawid[10];
  /* initialize a random seed */
  //pmb->uid.GetRawUID(rawid);
  //std::cout << "MeshBlockID: "<< rawid[2];//pmb->uid.uid[0];
  std::srand(pmb->gid);             
  
  /* read parameters */
  gm = pin->GetReal("problem","GM");
  gmgas = pin->GetReal("fluid","gamma");
  dist = pin->GetReal("problem","dist");
  rg = pin->GetReal("problem","rg");
  d0= pin->GetReal("problem","d0");
  amp= pin->GetReal("problem","amp");
  beta = pin->GetReal("problem","beta");
  
  /* calculate some  ant values */
  w0 = 1.0;
  cprime = 0.5/dist;
  en = 1.0/(gmgas-1.0);
  acons=0.5*(dist-1.0)/dist/(en+1.0);

  /* assign boundary conditions and gravitational force function*/
  pmb->pbval->EnrollFluidBoundaryFunction(inner_x1, stbv_iib);
  pmb->pbval->EnrollFluidBoundaryFunction(inner_x2, stbv_ijb);
  pmb->pbval->EnrollFluidBoundaryFunction(outer_x1, stbv_oib);
  pmb->pbval->EnrollFluidBoundaryFunction(outer_x2, stbv_ojb);

  /* allocate memory for the gas pressure */
  pr.NewAthenaArray(pmb->block_size.nx3+2*(NGHOST),pmb->block_size.nx2+2*(NGHOST),pmb->block_size.nx1+2*(NGHOST));
  /* Background */ 
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pfl->u(IDN,k,j,i)  = d0;
        pfl->u(IM1,k,j,i) = 0.0;
        pfl->u(IM2,k,j,i) = 0.0;
        pfl->u(IM3,k,j,i) = 0.0;
        pr(k,j,i)=d0/pco->x1v(i);
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
          tp=pco->x2f(j)+pco->dx2f(j)*0.1*(jj+1);
          tm=pco->x2f(j)+pco->dx2f(j)*0.1*jj;
          tv=0.5*(tp+tm)+(1.0-0.5*(tp-tm)/tan(0.5*(tp-tm)))/tan(0.5*(tp+tm));
          for(ii=0;ii<10;ii++) {
            rp = pco->x1f(i)+pco->dx1f(i)*0.1*(ii+1);
            rm = pco->x1f(i)+pco->dx1f(i)*0.1*ii;
            rv= ((SQR(SQR(rp))-SQR(SQR(rm)))/4.0)/((CUBE(rp)-CUBE(rm))/3.0);
            lv= 1.0/3.0*(CUBE(rp)-CUBE(rm))*(cos(tm)-cos(tp));
            vt+=lv;
            wt = rv*sin(tv);
            eq29 = gm/(w0*(en + 1.))*(w0/rv-0.5*SQR(w0/wt) - cprime);
            if (eq29 > 0.0) {
              dens = pow(eq29/acons,en);
              pp=dens*eq29;
              if (pp > d0/rv) {
                ftorus=1;
                ld+=dens*lv;
                lm+=dens*lv*sqrt(gm*w0)/wt;
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
          vv= 1.0/3.0*(CUBE(pco->x1f(i+1))-CUBE(pco->x1f(i)))*(cos(pco->x2f(j))-cos(pco->x2f(j+1)));
          ld/=vv; lm/=vv;
          pfl->u(IDN,k,j,i) = ld;
          pfl->u(IM3,k,j,i) = lm;
          pr(k,j,i)=MAX(acons*pow(ld,gmgas),d0/pco->x1v(i))*(1+amp*((double)rand()/(double)RAND_MAX-0.5));
        }
        //pfl->u(IEN,k,j,i)=pr(k,j,i)*en+0.5*SQR(pfl->u(IM3,k,j,i))/pfl->u(IDN,k,j,i);
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
        if(i!=pGrid->sil)
          pfl->u(IB1,k,j,i)  = magr(pGrid,i,j,k);
      }
      //not too sure what this "sil" is referring to
      if(pGrid->sil>=is && pGrid->sil <= ie+1)
        pfl->u(IB1,k,j,pGrid->sil) = pfl->u(IB1,k,j,pGrid->sil+1);
    }
  }
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      if(j!=pGrid->sjl && j!=pGrid->sju)
      {
        for (i=is; i<=ie; i++) {
          pfl->u(IB2,k,j,i)  = magt(pGrid,i,j,k);
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
        pfl->u(IB3,k,j,i) = 0.0;
      }
    }
  }
#endif


divbmax=0.0; idb=is; jdb=js; kdb=ks;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pfl->u(IEN,k,j,i)=pr(k,j,i)*en+0.5*SQR(pfl->u(IM3,k,j,i))/pfl->u(IDN,k,j,i);
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


  pr.DeleteAthenaArray();
}


/*  Boundary Condtions, outflowing, ix1, ox1, ix2, ox2  */
void stbv_iib(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke)
{
  Coordinates *pco = pmb->pcoord;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
  Real pg;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=(NGHOST); i++) {
        a(IDN,k,j,is-i) = a(IDN,k,j,is);
        a(IM1,k,j,is-i) = 0.0;
        a(IM2,k,j,is-i) = a(IM2,k,j,is); //corotating ghost region
        a(IM3,k,j,is-i) = a(IM3,k,j,is);
        pg = (a(IEN,k,j,is-i+1)-0.5*(SQR(a(IM1,k,j,is-i+1))+SQR(a(IM2,k,j,is-i+1))+SQR(a(IM3,k,j,is-i+1)))/a(IDN,k,j,is-i+1))*(gmgas-1.0);
        pg-=gm/SQR(pco->x1f(is-i+1))*a(IDN,k,j,is-i+1)*pco->dx1v(is-i);
        a(IEN,k,j,is-i)=pg/(gmgas-1.0)+0.5*SQR(a(IM1,k,j,is-i))/a(IDN,k,j,is-i);
      }
    }
  }

  #ifdef MHD
  /* B1i is not set at i=is-nghost */
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=1; i<=nghost-1; i++) {
          a(IB1,k,j,is-i) = a(IB1,k,j,is);  
        }
      }
    }

    if (pmb->Nx[1] > 1) ju=je+1; else ju=je;
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=ju; j++) {
        for (i=1; i<=nghost; i++) {
          a(IB2,k,j,is-i) = a(IB2,k,j,is);  
        }
      }
    }

    if (pmb->Nx[2] > 1) ku=ke+1; else ku=ke;
    for (k=ks; k<=ku; k++) {
      for (j=js; j<=je; j++) {
        for (i=1; i<=nghost; i++) {
          a(IB3,k,j,is-i) = a(IB3,k,j,is);  
        }
      }
    }
  #endif /* MHD */

  return;
}


void stbv_oib(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke)
{
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=(NGHOST); i++) {
        a(IDN,k,j,ie+i) = a(IDN,k,j,ie);
        a(IM1,k,j,ie+i) = a(IM1,k,j,ie);
        a(IM2,k,j,ie+i) = a(IM2,k,j,ie);
        a(IM3,k,j,ie+i) = a(IM3,k,j,ie);
        a(IEN,k,j,ie+i) = a(IEN,k,j,ie);
        if(a(IM1,k,j,ie+i) < 0.0)
        {
          a(IEN,k,j,ie+i) -= 0.5*SQR(a(IM1,k,j,ie+i))/a(IDN,k,j,ie+i);
          a(IM1,k,j,ie+i) = 0.0;
        }
      }
    }
  }
  #ifdef MHD
  /* i=ie+1 is not a boundary condition for the interface field B1i */
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=2; i<=nghost; i++) {
          a(IB1,k,j,ie+i)=a(IB1,k,j,ie);
        }
      }
    }

    if (pmb->Nx[1] > 1) ju=je+1; else ju=je;
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=ju; j++) {
        for (i=1; i<=nghost; i++) {
          a(IB2,k,j,ie+i)=a(IB2,k,j,ie);
        }
      }
    }

    if (pmb->Nx[2] > 1) ku=ke+1; else ku=ke;
    for (k=ks; k<=ku; k++) {
      for (j=js; j<=je; j++) {
        for (i=1; i<=nghost; i++) {
          a(IB3,k,j,ie+i)=a(IB3,k,j,ie);
        }
      }
    }
  #endif /* MHD */

  return;
}


void stbv_ijb(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke)
{
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=(NGHOST); j++) {
      for (i=is; i<=ie; i++) {
        a(IDN,k,js-j,i) = a(IDN,k,js,i);
        a(IM1,k,js-j,i) = a(IM1,k,js,i);
        a(IM2,k,js-j,i) = a(IM2,k,js,i);
        a(IM3,k,js-j,i) = a(IM3,k,js,i);
        a(IEN,k,js-j,i) = a(IEN,k,js,i);
        if(a(IM2,k,js-j,i) > 0.0)
        {
          a(IEN,k,js-j,i) -= 0.5*SQR(a(IM2,k,js-j,i))/a(IDN,k,js-j,i);
          a(IM2,k,js-j,i) = 0.0;
        }
      }
    }
  }
#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        a(IB1,k,js-j,i)=a(IB1,k,js,i);
      }
    }
  }

/* B2i is not set at j=js-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        a(IB2,k,js-j,i)=a(IB2,k,js,i);
      }
    }
  }

  if (pmb->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        a(IB3,k,js-j,i)=a(IB3,k,js,i);
      }
    }
  }
#endif /* MHD */
  return;
}


void stbv_ojb(MeshBlock *pmb, AthenaArray<Real> &a,
              int is, int ie, int js, int je, int ks, int ke)
{
  int i,j,k;
#ifdef MHD
  int ku; /* k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=(NGHOST); j++) {
      for (i=is; i<=ie; i++) {
        a(IDN,k,je+j,i) = a(IDN,k,je,i);
        a(IM1,k,je+j,i) = a(IM1,k,je,i);
        a(IM2,k,je+j,i) = a(IM2,k,je,i);
        a(IM3,k,je+j,i) = a(IM3,k,je,i);
        a(IEN,k,je+j,i) = a(IEN,k,je,i);
        if(a(IM2,k,je+j,i) < 0.0)
        {
          a(IEN,k,je+j,i) -= 0.5*SQR(a(IM2,k,je+j,i))/a(IDN,k,je+j,i);
          a(IM2,k,je+j,i) = 0.0;
        }
      }
    }
  }
#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
        a(IB1,k,je+j,i)=a(IB1,k,je,i);
      }
    }
  }

/* j=je+1 is not a boundary condition for the interface field B2i */
  for (k=ks; k<=ke; k++) {
    for (j=2; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        a(IB2,k,je+j,i)=a(IB2,k,je,i);
      }
    s
  }

  if (pmb->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        a(IB3,k,je+j,i)=a(IB3,k,je,i);
      }
    }
  }
#endif /* MHD */
  return;
}

