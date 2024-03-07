/**
# Interfacial forces

We assume that the interfacial acceleration can be expressed as
$$
\phi\mathbf{n}\delta_s/\rho
$$
with $\mathbf{n}$ the interface normal, $\delta_s$ the interface Dirac
function, $\rho$ the density and $\phi$ a generic scalar field. Using
a CSF/Peskin-like approximation, this can be expressed as
$$
\phi\nabla f/\rho
$$
with $f$ the volume fraction field describing the interface.

The interfacial force potential $\phi$ is associated to each VOF
tracer. This is done easily by adding the following [field
attributes](/Basilisk C#field-attributes). */

attribute {
  scalar phiPotential;
  scalar KappaBase;
  scalar KappaBasePhi;
  vector nPerturbed;
  vector sigmaForcing3;
}

/**
Interfacial forces are a source term in the right-hand-side of the
evolution equation for the velocity of the [centered Navier--Stokes
solver](navier-stokes/centered.h) i.e. it is an acceleration. If
necessary, we allocate a new vector field to store it. */

event defaults (i = 0) {  
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *){a});
  }
}

/**
The calculation of the acceleration is done by this event, overloaded
from [its definition](navier-stokes/centered.h#acceleration-term) in
the centered Navier--Stokes solver. */

event acceleration (i++)
{
  /**
  We check for all VOF interfaces for which $\phi$ is allocated. The
  corresponding volume fraction fields will be stored in *list*. */

  scalar * list = NULL;
  for (scalar f in interfaces)
    if (f.phiPotential.i) {
      list = list_add (list, f);
      /**
      To avoid undeterminations due to round-off errors, we remove
      values of the volume fraction larger than one or smaller than
      zero. */

      foreach()
        f[] = clamp (f[], 0., 1.);
      boundary ({f});
    }
  /**
  On trees we need to make sure that the volume fraction gradient
  is computed exactly like the pressure gradient. This is necessary to
  ensure well-balancing of the pressure gradient and interfacial force
  term. To do so, we apply the same prolongation to the volume
  fraction field as applied to the pressure field. */
  
#if TREE
  for (scalar f in list)
    f.prolongation = p.prolongation;
  boundary (list);
#endif

  /**
  Finally, for each interface for which $\phi$ is allocated, we
  compute the interfacial force acceleration
  $$
  \phi\mathbf{n}\delta_s/\rho \approx \alpha\phi\nabla f
  $$ 
  */


  face vector ia = a;
  foreach_face()
    for (scalar f in list)
      if (f[] != f[-1]  && fm.x[] > 0.) {
        scalar  phiPotential = f.phiPotential;
        scalar    KappaBase  = f.KappaBase;
        scalar KappaBasePhi  = f.KappaBasePhi;
        vector nPerturbed    = f.nPerturbed;
        vector sigmaForcing3 = f.sigmaForcing3;

      double phiPotentialf  = (phiPotential[] + phiPotential[-1])/2.;
      double KappaBasef     = (KappaBase[] < nodata && KappaBase[-1] < nodata) ?
	  (KappaBase[] + KappaBase[-1])/2. :
	  KappaBase[] < nodata ? KappaBase[] :
	  KappaBase[-1] < nodata ? KappaBase[-1] :
	  0.;
      double nPerturbedf    =
      (nPerturbed.x[] < nodata && nPerturbed.x[-1] < nodata) ?
	  (nPerturbed.x[] + nPerturbed.x[-1])/2. :
	  nPerturbed.x[] < nodata ? nPerturbed.x[] :
	  nPerturbed.x[-1] < nodata ? nPerturbed.x[-1] :
	  0.;
      double KappaBasePhif = (KappaBasePhi[] + KappaBasePhi[-1])/2.;

      double sigmaForcing3f = (sigmaForcing3.x[]+sigmaForcing3.x[-1])/2.;

      ia.x[] += alpha.x[]/fm.x[]*f.sigma*phiPotentialf*(f[] - f[-1])/Delta;
      ia.x[] += alpha.x[]/fm.x[]*f.sigma*KappaBasePhif*nPerturbedf*fabs((f[] - f[-1])/Delta);
//       ia.x[] += KappaBasef*nPerturbed.x[]*sqrt(sq((f[] - f[-1])/Delta)+sq((f[0,1]+f[-1,1]-f[-1,-1]-f[0,-1])/(4*Delta)));
      ia.x[] += alpha.x[]/fm.x[]*f.sigma*sigmaForcing3f*fabs((f[] - f[-1])/Delta);
      }

  /**
  On trees, we need to restore the prolongation values for the
  volume fraction field. */
  
#if TREE
  for (scalar f in list)
    f.prolongation = fraction_refine;
  boundary (list);
#endif

  /**
  Finally we free the potential fields and the list of volume
  fractions. */

  for (scalar f in list) {
    scalar phiPotential = f.phiPotential;
    scalar KappaBase  = f.KappaBase;
    scalar KappaBasePhi  = f.KappaBasePhi;
    vector nPerturbed = f.nPerturbed;
    vector sigmaForcing3   = f.sigmaForcing3;

    delete ({phiPotential,KappaBase,KappaBasePhi,nPerturbed.x,nPerturbed.y,sigmaForcing3.x,sigmaForcing3.y});
//     delete ({phiPotential,KappaBasePhi,nPerturbed.x,nPerturbed.y,sigmaForcing3.x,sigmaForcing3.y});
    f.phiPotential.i = 0;
    f.KappaBase.i = 0;
    f.KappaBasePhi.i = 0;
    f.nPerturbed.x.i=0;
    f.nPerturbed.y.i=0;
    f.sigmaForcing3.x.i=0;
    f.sigmaForcing3.y.i=0;
  }
  free (list);
}

/**
## References

See Section 3, pages 8-9 of:

~~~bib
@hal{popinet2018, hal-01528255}
~~~
*/
