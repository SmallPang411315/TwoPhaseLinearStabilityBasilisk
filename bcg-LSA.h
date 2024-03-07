/**
	2022.03.09 (GONG MINJIANG): 
		u'\dot \nabla u' needs bcg scheme  
*/

#include"heaviside.h"
extern scalar Phibase, phi;
extern double rho1, rho2;
extern scalar rhov;
// extern face vector alpha;
// extern scalar rho, dhea;        
// extern double rho1, rho2, mu1, mu2, rho_til, mu_til;//all declared in two-phase-LS.h

struct Advection {
  scalar * tracers;
  face vector uf;
  double dt;
  scalar * src; // optional
};

struct Advection_Linear {
  scalar * tracers;
  scalar * Ubase;
  face vector uf;
  double dt;
  scalar * src; // optional
};

struct Advection_Redistance {
  scalar * interfaces;
  scalar * tracers;
  face vector uf;
  double dt;
  scalar * src; // optional
};

void tracer_fluxes (scalar f,
			face vector uf,
			face vector flux,
		    double dt,
		    (const) scalar src)
{
  vector g[];
  gradients ({f}, {g});

  foreach_face() {
    double un = dt*uf.x[]/(fm.x[]*Delta + SEPS), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = f[i] + (src[] + src[-1])*dt/4. + s*(1. - s*un)*g.x[i]*Delta/2.;//f2 is u^{n+1/2} and f is u

    #if dimension > 1
    if (fm.y[i] && fm.y[i,1]) {
      double vn = (uf.y[i] + uf.y[i,1])/(fm.y[i] + fm.y[i,1]);
      double fyy = vn < 0. ? f[i,1] - f[i] : f[i] - f[i,-1];
      f2 -= dt*vn*fyy/(2.*Delta);
    }
    #endif
    #if dimension > 2
    if (fm.z[i] && fm.z[i,0,1]) {
      double wn = (uf.z[i] + uf.z[i,0,1])/(fm.z[i] + fm.z[i,0,1]);
      double fzz = wn < 0. ? f[i,0,1] - f[i] : f[i] - f[i,0,-1];
      f2 -= dt*wn*fzz/(2.*Delta);
    }
    #endif

    flux.x[] = f2*uf.x[];
  }

  boundary_flux ({flux});
}
/*
void tracer_fluxes2 (scalar f,
			face vector uf,
			face vector flux,
		    double dt,
		    (const) scalar src)
{
  foreach_face() {
    double un = dt*uf.x[]/(fm.x[]*Delta + SEPS), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = f[i];

    flux.x[] = f2*uf.x[];
  }

  boundary_flux ({flux});
}*/

void advection (struct Advection p)
{
  /**
  If *src* is not provided we set all the source terms to zero. */
  
  scalar * lsrc   = p.src;
  if (!lsrc) {
    const scalar zero[] = 0.;
    for (scalar s in p.tracers)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
  scalar f, src;
  
  face vector flux[];
  for (f,src in p.tracers,lsrc) {
	  tracer_fluxes(f, p.uf, flux, p.dt, src);
      foreach()
        foreach_dimension(){
		  f[] += p.dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
	    }  	  
  }
  
  boundary (p.tracers);
  
  if (!p.src)
    free (lsrc);
}

void advection_linear (struct Advection_Linear p)
{
  /**
  If *src* is not provided we set all the source terms to zero. */
  
  scalar * lsrc     = p.src;
  scalar * lUbase   = p.Ubase;
  if (!lsrc) {
    const scalar zero[] = 0.;
    for (scalar s in p.tracers)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
  scalar f, Ubase, src;
  
  face vector flux[];
  for (f, Ubase, src in p.tracers, lUbase, lsrc) {
	  tracer_fluxes(Ubase, p.uf, flux, p.dt, src);
      foreach()
        foreach_dimension(){
		  f[] += p.dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
	    }  	  
  }
  
  boundary (p.tracers);
  
  if (!p.src)
    free (lsrc);
}


void tracer_fluxes_Ubase (scalar f,
			face vector uf,
			face vector flux,
		    double dt,
		    (const) scalar src)
{
  foreach_face() {
    double un = dt*uf.x[]/(fm.x[]*Delta + SEPS), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = f[i];

    flux.x[] = f2*uf.x[];
  }

  boundary_flux ({flux});
}

void advection_linear_base (struct Advection_Linear p)
{
  scalar * lsrc     = p.src;
  scalar * lUbase   = p.Ubase;
  if (!lsrc) {
    const scalar zero[] = 0.;
    for (scalar s in p.tracers)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
  scalar f, Ubase, src;

  face vector flux[];
  for (f, Ubase, src in p.tracers, lUbase, lsrc) {
	  tracer_fluxes_Ubase(Ubase, p.uf, flux, p.dt, src);
      foreach()
        foreach_dimension(){
		  f[] += 2./rhov[]*value_dHeaviside(Phibase[], Delta)*phi[]*(rho1-rho2)*p.dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
//           f[] += 2./rhov[]*(value_dHeaviside(-Phibase[], Delta)+value_dHeaviside(-Phibase[-1], Delta))/2.*(phi[]+phi[-1])/2.*(rho1-rho2)*p.dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
	    }
  }

  boundary (p.tracers);

  if (!p.src)
    free (lsrc);
}

/**
	The below is solving 
	\partial_t\phi^{\prime} 
	+ \mathbf{u_f}^{\prime}\cdot \nabla\varPhi   
	+ \mathbf{U}\cdot \nabla \phi^{\prime}     = 0
*/

/**
	\nabla \cdot (\phi^{\prime} U)
*/
void phi_tracer_fluxes (scalar f,
			face vector uf,
			face vector flux,
		    double dt,
		    (const) scalar src)
{
  vector g[];
  gradients ({f}, {g});

  foreach_face() {
    double un = dt*uf.x[]/(fm.x[]*Delta + SEPS), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = f[i] + (src[] + src[-1])*dt/4. + s*(1. - s*un)*g.x[i]*Delta/2.;//f2 is u^{n+1/2} and f is u

    #if dimension > 1
    if (fm.y[i] && fm.y[i,1]) {
      double vn = (uf.y[i] + uf.y[i,1])/(fm.y[i] + fm.y[i,1]);
      double fyy = vn < 0. ? f[i,1] - f[i] : f[i] - f[i,-1];
      f2 -= dt*vn*fyy/(2.*Delta);
    }
    #endif
    #if dimension > 2
    if (fm.z[i] && fm.z[i,0,1]) {
      double wn = (uf.z[i] + uf.z[i,0,1])/(fm.z[i] + fm.z[i,0,1]);
      double fzz = wn < 0. ? f[i,0,1] - f[i] : f[i] - f[i,0,-1];
      f2 -= dt*wn*fzz/(2.*Delta);
    }
    #endif

    flux.x[] = f2*uf.x[];
  }

  boundary_flux ({flux});
}

void phi_tracer_fluxes1 (scalar f,
		    face vector uf,
		    face vector flux,
		    double dt,
		    (const) scalar src)
{
   vector g[];
   gradients ({f}, {g});

  foreach_face() {
     double un = dt*uf.x[]/(fm.x[]*Delta + SEPS), s = sign(un);
     int i = -(s + 1.)/2.;
     double f2 = f[i];

     flux.x[] = f2*uf.x[];
  }
  boundary_flux ({flux});
}

/**
	\nabla \cdot (\varPhi u^{\prime})
*/
void phi_tracer_fluxes2 (scalar f, // f is Phibase
		    face vector uf,
		    face vector flux,
		    double dt,
		    (const) scalar src)
{
   vector g[];
   gradients ({f}, {g});

  foreach_face() {
     double un = dt*uf.x[]/(fm.x[]*Delta + SEPS), s = sign(un);
     int i = -(s + 1.)/2.;
     double f2 = f[i];

     flux.x[] = f2*uf.x[];
  }
  boundary_flux ({flux});
}


struct Advection_LinearPhi1 {
  scalar * tracers;
  face vector Ufbase;
  double dt;
  scalar * src; // optional
};

struct Advection_LinearPhi2 {
  scalar * tracers;
  scalar * Phibase;
  face vector uf;
  double dt;
  scalar * src; // optional
};

void phi_linearadvection1 (struct Advection_LinearPhi1 p)
{ 
  scalar * lsrc   = p.src;
  if (!lsrc) {
    const scalar zero[] = 0.;
    for (scalar s in p.tracers)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
  scalar f,src;
  
  face vector flux1[];
  
  for (f,src in p.tracers,lsrc) {
    phi_tracer_fluxes (f, p.Ufbase, flux1, p.dt, src);
    foreach()
      foreach_dimension(){
        f[] += p.dt*((flux1.x[] - flux1.x[1])/(Delta*cm[]));
	  }
  }
  
  boundary (p.tracers);

  if (!p.src)
    free (lsrc);
}

void phi_linearadvection2 (struct Advection_LinearPhi2 p)
{
  scalar * lsrc       = p.src;
//   scalar * lPhibase   = p.Phibase;
  if (!lsrc) {
    const scalar zero[] = 0.;
    for (scalar s in p.tracers)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
//   scalar f, Phibase, src;
  scalar f, src;

  face vector flux2[];
//   for (f, Phibase, src in p.tracers, lPhibase, lsrc) {
    for (f, src in p.tracers, lsrc) {
	  phi_tracer_fluxes(Phibase, p.uf, flux2, p.dt, src);
      foreach()
        foreach_dimension(){
		  f[] += p.dt*(flux2.x[] - flux2.x[1])/(Delta*cm[]);
	    }
  }

  boundary (p.tracers);

  if (!p.src)
    free (lsrc);
}
