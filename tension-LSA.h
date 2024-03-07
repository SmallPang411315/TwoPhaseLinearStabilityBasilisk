#include "iforce-LSA.h" 
#include "curvature.h"
/**
The surface tension coefficient is associated to each VOF tracer. */

attribute {
  double sigma;
}

event stability (i++) {

  /**
  We first compute the minimum and maximum values of $\alpha/f_m =
  1/\rho$, as well as $\Delta_{min}$. */

  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin)) {
    if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
    if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
    if (Delta < dmin) dmin = Delta;
  }
  double rhom = (1./amin + 1./amax)/2.;

  /**
  The maximum timestep is set using the sum of surface tension
  coefficients. */

  double sigma = 0.;
  for (scalar c in interfaces)
    sigma += c.sigma;
  if (sigma) {
    double dt = sqrt (rhom*cube(dmin)/(pi*sigma));
    if (dt < dtmax)
      dtmax = dt;
  }
}


scalar fPhi[],gPhi[],fx[],fy[],gx[],gy[],phiPotentialTemp[];
vector  dPhibase[],dPhibaseTemp[],dphi[],dphiTemp[],dPhibasex[],dPhibasexTemp[],dPhibasey[],dPhibaseyTemp[],dphix[],dphixTemp[],dphiy[],dphiyTemp[],df[],NormBase[],KNphi[];
vector nabla_KNxphi[];
vector nabla_KNyphi[];
event acceleration (i++)
{
  for (scalar f in interfaces)
    if (f.sigma) {
      boundary({Phibase});
      gradients ({Phibase}, {dPhibase}); // \par{Phi}/\par{}
      foreach(){
        double dPhix = dPhibase.x[];
        double dPhiy = dPhibase.y[];
        point = locate (x, -y);
        dPhibase.x[] = dPhix;
        dPhibase.y[] = -dPhiy;
      }
      boundary({Phibase});
//       dPhibaseTemp = dPhibase;
//       foreach(){
//         dPhibase.x[] = (4.*dPhibaseTemp.x[] +
// 	    2.*(dPhibaseTemp.x[0,1] + dPhibaseTemp.x[0,-1] + dPhibaseTemp.x[1,0] + dPhibaseTemp.x[-1,0]) +
// 	    dPhibaseTemp.x[-1,-1] + dPhibaseTemp.x[1,-1] + dPhibaseTemp.x[1,1] + dPhibaseTemp.x[-1,1])/16.;
//         dPhibase.y[] = (4.*dPhibaseTemp.y[] +
// 	    2.*(dPhibaseTemp.y[0,1] + dPhibaseTemp.y[0,-1] + dPhibaseTemp.y[1,0] + dPhibaseTemp.y[-1,0]) +
// 	    dPhibaseTemp.y[-1,-1] + dPhibaseTemp.y[1,-1] + dPhibaseTemp.y[1,1] + dPhibaseTemp.y[-1,1])/16.;
//       }

      phi.gradient = minmod2;
      boundary({phi});
      gradients ({phi}, {dphi});         // \par{phi}/\par{}
//       dphiTemp = dphi;
//       foreach(){
//         dphi.x[] = (4.*dphiTemp.x[] +
// 	    2.*(dphiTemp.x[0,1] + dphiTemp.x[0,-1] + dphiTemp.x[1,0] + dphiTemp.x[-1,0]) +
// 	    dphiTemp.x[-1,-1] + dphiTemp.x[1,-1] + dphiTemp.x[1,1] + dphiTemp.x[-1,1])/16.;
//         dphi.y[] = (4.*dphiTemp.y[] +
// 	    2.*(dphiTemp.y[0,1] + dphiTemp.y[0,-1] + dphiTemp.y[1,0] + dphiTemp.y[-1,0]) +
// 	    dphiTemp.y[-1,-1] + dphiTemp.y[1,-1] + dphiTemp.y[1,1] + dphiTemp.y[-1,1])/16.;
//       }

      boundary({dPhibase.x,dPhibase.y});
      gradients ({dPhibase.x}, {dPhibasex}); // \par{Phi_x}/\par{}
      gradients ({dPhibase.y}, {dPhibasey}); // \par{Phi_y}/\par{}
//       dPhibasexTemp = dPhibasex;
//       dPhibaseyTemp = dPhibasey;
//       foreach(){
//         dPhibasex.x[] = (4.*dPhibasexTemp.x[] +
// 	    2.*(dPhibasexTemp.x[0,1] + dPhibasexTemp.x[0,-1] + dPhibasexTemp.x[1,0] + dPhibasexTemp.x[-1,0]) +
// 	    dPhibasexTemp.x[-1,-1] + dPhibasexTemp.x[1,-1] + dPhibasexTemp.x[1,1] + dPhibasexTemp.x[-1,1])/16.;
//         dPhibasex.y[] = (4.*dPhibasexTemp.y[] +
// 	    2.*(dPhibasexTemp.y[0,1] + dPhibasexTemp.y[0,-1] + dPhibasexTemp.y[1,0] + dPhibasexTemp.y[-1,0]) +
// 	    dPhibasexTemp.y[-1,-1] + dPhibasexTemp.y[1,-1] + dPhibasexTemp.y[1,1] + dPhibasexTemp.y[-1,1])/16.;
//
//         dPhibasey.x[] = dPhibasex.y[];
//         dPhibasey.y[] = (4.*dPhibaseyTemp.y[] +
// 	    2.*(dPhibaseyTemp.y[0,1] + dPhibaseyTemp.y[0,-1] + dPhibaseyTemp.y[1,0] + dPhibaseyTemp.y[-1,0]) +
// 	    dPhibaseyTemp.y[-1,-1] + dPhibaseyTemp.y[1,-1] + dPhibaseyTemp.y[1,1] + dPhibaseyTemp.y[-1,1])/16.;
//       }

      boundary({dphi.x,dphi.y});
//       dphi.x.gradient = minmod2;
//       dphi.y.gradient = minmod2;
      gradients ({dphi.x}, {dphix});         // \par{phi_x}/\par{}
      gradients ({dphi.y}, {dphiy});         // \par{phi_x}/\par{}
//       double filter = 0.4;
//     foreach() {
//         dphix.x[] = (filter*(phi[1,1] + phi[-1,1] - 2.*phi[0,1]) +
// 		(phi[1] + phi[-1] - 2.*phi[]) +
// 		filter*(phi[1,-1] + phi[-1,-1] - 2.*phi[0,-1]))/
//     ((1. + 2.*filter)*sq(Delta));
//         dphiy.y[] = (filter*(phi[1,1] + phi[1,-1] - 2.*phi[1]) +
// 		(phi[0,1] + phi[0,-1] - 2.*phi[]) +
// 		filter*(phi[-1,1] + phi[-1,-1] - 2.*phi[-1]))/
//     ((1. + 2.*filter)*sq(Delta));
//       }
//       foreach() {
//         dphix.y[] = (phi[1,1] + phi[-1,-1] - phi[1,-1] - phi[-1,1])/(4.*sq(Delta));
//         dphiy.x[] = dphix.y[];
//       }

      foreach(){
        fPhi[] = 1./sqrt(sq(dPhibase.x[])+sq(dPhibase.y[]));
        gPhi[] = -(dphi.x[]*dPhibase.x[]+dphi.y[]*dPhibase.y[]);
      }

        foreach(){
          fx[] = -cube(fPhi[])*(dPhibase.x[]*dPhibasex.x[]+dPhibase.y[]*dPhibasex.y[]);
          fy[] = -cube(fPhi[])*(dPhibase.x[]*dPhibasex.y[]+dPhibase.y[]*dPhibasey.y[]);
          gx[] = -(dphi.x[]*dPhibasex.x[]+dphix.x[]*dPhibase.x[]+dphi.y[]*dPhibasex.y[]  +dphix.y[]*dPhibase.y[]);
          gy[] = -(dphi.y[]*dPhibasey.y[]+dphiy.y[]*dPhibase.y[]+dphi.x[]*dPhibasey.x[]  +dphiy.x[]*dPhibase.x[]);
        }

        scalar phiPotential = f.phiPotential;
        if (phiPotential.i)
        foreach(){
           phiPotential[] += (fx[]*dphi.x[]+fPhi[]*dphix.x[]+fy[]*dphi.y[]+fPhi[]*dphiy.y[]+gx[]*cube(fPhi[])*dPhibase.x[]+gPhi[]*3.*sq(fPhi[])*fx[]*dPhibase.x[]+gPhi[]*cube(fPhi[])*dPhibasex.x[]+gy[]*cube(fPhi[])*dPhibase.y[]+gPhi[]*3.*sq(fPhi[])*fy[]*dPhibase.y[]+gPhi[]*cube(fPhi[])*dPhibasey.y[]);
        }
        else{
          phiPotential = new scalar;
          foreach(){
           phiPotential[] = (fx[]*dphi.x[]+fPhi[]*dphix.x[]+fy[]*dphi.y[]+fPhi[]*dphiy.y[]+gx[]*cube(fPhi[])*dPhibase.x[]+gPhi[]*3.*sq(fPhi[])*fx[]*dPhibase.x[]+gPhi[]*cube(fPhi[])*dPhibasex.x[]+gy[]*cube(fPhi[])*dPhibase.y[]+gPhi[]*3.*sq(fPhi[])*fy[]*dPhibase.y[]+gPhi[]*cube(fPhi[])*dPhibasey.y[]);
        }
          f.phiPotential = phiPotential;
        }

//         phiPotentialTemp = phiPotential;
//         foreach(){
//         phiPotential[] = (4.*phiPotentialTemp[] +
// 	    2.*(phiPotentialTemp[0,1] + phiPotentialTemp[0,-1] + phiPotentialTemp[1,0] + phiPotentialTemp[-1,0]) +  phiPotentialTemp[-1,-1] + phiPotentialTemp[1,-1] + phiPotentialTemp[1,1] + phiPotentialTemp[-1,1])/16.;
//       }


//         foreach(){
//           if (t>95&&t<100&&f[] != f[-1] && fm.x[] > 0.)
//           printf("%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
//                  t,x,y,
//                  fx[]*dphi.x[],
//                  fPhi[]*dphix.x[],
//                  fy[]*dphi.y[],
//                  fPhi[]*dphiy.y[],
//                  gx[]*cube(fPhi[])*dPhibase.x[],
//                  gPhi[]*3.*sq(fPhi[])*fx[]*dPhibase.x[],
//                  gPhi[]*cube(fPhi[])*dPhibasex.x[],
//                  gy[]*cube(fPhi[])*dPhibase.y[],
//                  gPhi[]*3.*sq(fPhi[])*fy[]*dPhibase.y[],
//                  gPhi[]*cube(fPhi[])*dPhibasey.y[]);
//         }

/////////////////// Surface tension term 2 ///////////////////////////
        scalar KappaBase = f.KappaBase;
        if (KappaBase.i){
        curvature (f, KappaBase, f.sigma, add = true);
        }
        else{
          KappaBase = new scalar;
          curvature (f, KappaBase, f.sigma, add = false);
          f.KappaBase = KappaBase;
        }

      scalar KappaBasePhi = f.KappaBasePhi;
      if (KappaBasePhi.i){
        foreach()KappaBasePhi[] += (dPhibasex.x[]*sq(dPhibase.y[])-2.*dPhibasex.y[]*dPhibase.x[]*dPhibase.y[]+dPhibasey.y[]*sq(dPhibase.x[]))*cube(fPhi[]);
      }
      else{
        KappaBasePhi = new scalar;
        foreach()KappaBasePhi[] = (dPhibasex.x[]*sq(dPhibase.y[])-2.*dPhibasex.y[]*dPhibase.x[]*dPhibase.y[]+dPhibasey.y[]*sq(dPhibase.x[]))*cube(fPhi[]);
        f.KappaBasePhi = KappaBasePhi;
      }

      vector nPerturbed = f.nPerturbed;
      if (nPerturbed.x.i && nPerturbed.y.i){
        foreach() {
          nPerturbed.x[] += dphi.x[]*fPhi[]+gPhi[]*cube(fPhi[])*dPhibase.x[];
          nPerturbed.y[] += dphi.y[]*fPhi[]+gPhi[]*cube(fPhi[])*dPhibase.y[];
        }
      }
      else{
        nPerturbed = new vector;
        foreach() {
          nPerturbed.x[] = dphi.x[]*fPhi[]+gPhi[]*cube(fPhi[])*dPhibase.x[];
          nPerturbed.y[] = dphi.y[]*fPhi[]+gPhi[]*cube(fPhi[])*dPhibase.y[];
        }
        f.nPerturbed = nPerturbed;
      }
/////////////////// Surface tension term 3 ///////////////////////////
      foreach(){
        NormBase.x[] = dPhibase.x[]*fPhi[];
        NormBase.y[] = dPhibase.y[]*fPhi[];
      }
      foreach() {
        KNphi.x[] = KappaBasePhi[]*NormBase.x[]*phi[];
        KNphi.y[] = KappaBasePhi[]*NormBase.y[]*phi[];
      }

      boundary({KNphi.x,KNphi.y});
      gradients({KNphi.x},{nabla_KNxphi});
      gradients({KNphi.y},{nabla_KNyphi});

      vector sigmaForcing3 = f.sigmaForcing3;
      if (sigmaForcing3.x.i){
        foreach() {
          sigmaForcing3.x[] += (NormBase.x[]*nabla_KNxphi.x[]+NormBase.y[]*nabla_KNxphi.y[]);
          sigmaForcing3.y[] += (NormBase.x[]*nabla_KNyphi.x[]+NormBase.y[]*nabla_KNyphi.y[]);
        }
      }
      else{
        sigmaForcing3 = new vector;
        foreach() {
          sigmaForcing3.x[] = (NormBase.x[]*nabla_KNxphi.x[]+NormBase.y[]*nabla_KNxphi.y[]);
          sigmaForcing3.y[] = (NormBase.x[]*nabla_KNyphi.x[]+NormBase.y[]*nabla_KNyphi.y[]);
        }
        f.sigmaForcing3 = sigmaForcing3;
      }

    }
}
