/**
# Tracer advection event

This event integrates advection equations of the form
$$
\partial_t\phi^{\prime} + \mathbf{u_f}^{\prime}\cdot \nabla\varPhi   
                        + \mathbf{U}\cdot \nabla \phi^{\prime} 
						= 0
$$
where $\mathbf{u_f}$ is the velocity field and $f_i$ are a list of
passive tracers.

The `tracers` list is defined elsewhere (typically by the user), the
face vector field `uf` and the timestep `dt` are defined by a
solver. */

extern scalar * tracers;
extern scalar Phibase;
extern vector u, Ubase;
extern face vector uf, Ufbase;
extern double dt; 

/**
On adaptive meshes, tracers need to use linear interpolation (rather
than the default bilinear interpolation) to ensure conservation when
refining cells. */

#if TREE
event defaults (i = 0) {
  for (scalar s in tracers) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
  }
}
#endif

/**
The integration is performed using the Bell-Collela-Glaz scheme. */

#include "bcg-LSA.h"
// #include "level_set_reinitialization-LSA.h"

event tracer_advection (i++,last) {
  phi_linearadvection1 (tracers, Ufbase, dt);
  phi_linearadvection2 (tracers, (scalar *){Phibase}, uf, dt);
//   redistance_Linearphi(tracers, Phibase, dt);
}

/**
Diffusion can be added by overloading this hook. */

event tracer_diffusion (i++,last);
