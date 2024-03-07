/**
# Two-phase interfacial flows
	2022.03.09 (GONGMINJIANG)
		if i wanna do a LSA solver, no volume fraction is considered. 
		So $f[]$ in this header file need to be deleted. But now the 
		stage, the averaged density and viscosity are testd based on
		Heaviside function.
 */

#include "heaviside.h"
#include "fractions.h"

scalar f[], Phibase[], phi[];
scalar * interfaces = {f};
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;

/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;

  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */
  
  if (mu1 || mu2){
    mu      = new face vector;
    mu_base = new face vector;
  }
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sphi[];
scalar sf[];
#else
# define sphi phi
# define sf f
#endif

event tracer_advection (i++)
{
#ifndef sf
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] +
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf[] = (8.*f[] +
	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif // !sf

#ifndef sphi
#if dimension <= 2
  foreach(){
    sphi[] = (4.*phi[] + 
	    2.*(phi[0,1] + phi[0,-1] + phi[1,0] + phi[-1,0]) +
		phi[-1,-1] + phi[1,-1] + phi[1,1] + phi[-1,1])/16.;
  }
#else // dimension == 3
  foreach()
    sphi[] = (8.*phi[] +
	    4.*(phi[-1] + phi[1] + phi[0,1] + phi[0,-1] + phi[0,0,1] + phi[0,0,-1]) +
	    2.*(phi[-1,1] + phi[-1,0,1] + phi[-1,0,-1] + phi[-1,-1] + 
		phi[0,1,1] + phi[0,1,-1] + phi[0,-1,1] + phi[0,-1,-1] +
		phi[1,1] + phi[1,0,1] + phi[1,-1] + phi[1,0,-1]) +
	    phi[1,-1,1] + phi[-1,1,1] + phi[-1,1,-1] + phi[1,1,1] +
	    phi[1,1,-1] + phi[-1,-1,-1] + phi[1,-1,-1] + phi[-1,-1,1])/64.;
#endif
#endif // !sphi

#if TREE
  sf.prolongation = refine_bilinear;
  sphi.prolongation = refine_bilinear;
  boundary ({sf, sphi});
  boundary ({sphi});
#endif
}

event properties (i++)
{
  foreach_face() {
    // double ff = (sf[] + sf[-1])/2.;
    // alphav.x[] = fm.x[]/rho(ff);
    double fphi = (Phibase[] + Phibase[-1])/2.;
	alphav.x[] = fm.x[]/( value_Heaviside(-fphi, Delta)*(rho1 - rho2) + rho2 );
    if (mu1 || mu2) {
      face vector muv = mu;
      face vector mu_basev = mu_base;
      muv.x[]      = fm.x[]*( value_Heaviside(-fphi, Delta)*(mu1 - mu2) + mu2  );
      mu_basev.x[] = fm.x[]*( value_dHeaviside(-fphi, Delta)*phi[]*(mu1 - mu2) );
    }
  }
  foreach(){
    rhov[] = cm[]*( value_Heaviside(-Phibase[], Delta)*(rho1 - rho2) + rho2 );//cm=1 in cartesian grid
  }
  boundary({f});

#if TREE  
  sf.prolongation = fraction_refine;
  sphi.prolongation = fraction_refine;
  boundary ({sphi,sf});
#endif
}
