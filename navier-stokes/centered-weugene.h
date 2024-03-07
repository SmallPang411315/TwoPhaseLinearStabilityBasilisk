#include "run.h"
#include "timestep.h"
#include "bcg.h"
#if EMBED
#include "viscosity-embed.h"
#else
#ifndef BRINKMAN_PENALIZATION
#include "viscosity.h"
#else
#include "viscosity-weugene.h"
#endif
#endif

scalar p[];
vector u[], g[];
scalar pf[];
face vector uf[];

(const) face vector mu = zerof, a = zerof, alpha = unityf, kappa = zerof;
(const) scalar rho = unity;
mgstats mgp, mgpf, mgu;
bool stokes = false;

#if EMBED
# define neumann_pressure(i) (alpha.n[i] ? a.n[i]*fm.n[i]/alpha.n[i] :	\
			      a.n[i]*rho[]/(cm[] + SEPS))
#else
# define neumann_pressure(i) (a.n[i]*fm.n[i]/alpha.n[i])
#endif

p[right] = neumann (neumann_pressure(ghost));
p[left]  = neumann (- neumann_pressure(0));

#if AXI
uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which
                             // is zero on the axis of symmetry
p[top]    = neumann (neumann_pressure(ghost));
#else // !AXI
#  if dimension > 1
p[top]    = neumann (neumann_pressure(ghost));
p[bottom] = neumann (- neumann_pressure(0));
#  endif
#  if dimension > 2
p[front]  = neumann (neumann_pressure(ghost));
p[back]   = neumann (- neumann_pressure(0));
#  endif
#endif // !AXI

#if TREE && EMBED
void pressure_embed_gradient (Point point, scalar p, coord * g)
{
  foreach_dimension()
    g->x = rho[]/(cm[] + SEPS)*(a.x[] + a.x[1])/2.;
}
#endif // TREE && EMBED

event defaults (i = 0)
{

  CFL = 0.8;

  p.nodump = pf.nodump = true;
  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    face vector alphav = alpha;
    foreach_face()
      alphav.x[] = fm.x[];
    boundary ((scalar *){alpha});
}

#if TREE
  uf.x.refine = refine_face_solenoidal;
#if EMBED
  uf.x.refine = refine_face;
  foreach_dimension()
    uf.x.prolongation = refine_embed_face_x;
  for (scalar s in {p, pf, u, g}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
  }
  for (scalar s in {p, pf})
    s.embed_gradient = pressure_embed_gradient;
#endif // EMBED
#endif // TREE
}

double dtmax;

event init (i = 0)
{
  boundary ((scalar *){u});
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*face_value (u.x, 0);
  boundary ((scalar *){uf});
  event ("properties");
  dtmax = DT;
  event ("stability");
}

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
}

event vof (i++,last);
event tracer_advection (i++,last);

event properties (i++,last) {
  boundary ({alpha, mu, rho, kappa}); // Weugene: kappa added
}

event tracer_diffusion (i++,last);

void prediction()
{
  vector du;
  foreach_dimension() {
    scalar s = new scalar;
    du.x = s;
  }

  if (u.x.gradient)
    foreach()
      foreach_dimension() {
#if EMBED
        if (!fs.x[] || !fs.x[1])
	  du.x[] = 0.;
	else
#endif
	  du.x[] = u.x.gradient (u.x[-1], u.x[], u.x[1])/Delta;
      }
  else
    foreach()
      foreach_dimension() {
#if EMBED
        if (!fs.x[] || !fs.x[1])
	  du.x[] = 0.;
	else
#endif
	  du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
    }
  boundary ((scalar *){du});

  trash ({uf});
  foreach_face() {
    double un = dt*(u.x[] + u.x[-1])/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    uf.x[] = u.x[i] + (g.x[] + g.x[-1])*dt/4. + s*(1. - s*un)*du.x[i]*Delta/2.;
    #if dimension > 1
    if (fm.y[i,0] && fm.y[i,1]) {
      double fyy = u.y[i] < 0. ? u.x[i,1] - u.x[i] : u.x[i] - u.x[i,-1];
      uf.x[] -= dt*u.y[i]*fyy/(2.*Delta);
    }
    #endif
    #if dimension > 2
    if (fm.z[i,0,0] && fm.z[i,0,1]) {
      double fzz = u.z[i] < 0. ? u.x[i,0,1] - u.x[i] : u.x[i] - u.x[i,0,-1];
      uf.x[] -= dt*u.z[i]*fzz/(2.*Delta);
    }
    #endif
    uf.x[] *= fm.x[];
  }
  boundary ((scalar *){uf});

  delete ((scalar *){du});
}

event advection_term (i++,last)
{
  if (!stokes) {
    prediction();
    mgpf = project (uf, pf, alpha, dt/2., mgpf.nrelax);
    advection ((scalar *){u}, uf, dt, (scalar *){g});
  }
}

static void correction (double dt)
{
  foreach()
    foreach_dimension()
      u.x[] += dt*g.x[];
  boundary ((scalar *){u});
}

event viscous_term (i++,last)
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity (u, mu, rho, dt, mgu.nrelax);
    correction (-dt);
  }

  if (!is_constant(a.x)) {
    face vector af = a;
    trash ({af});
    foreach_face()
      af.x[] = 0.;
  }
}

event acceleration (i++,last)
{
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*(face_value (u.x, 0) + dt*a.x[]);
  boundary ((scalar *){uf, a});
}

void centered_gradient (scalar p, vector g)
{
  face vector gf[];
  foreach_face()
    gf.x[] = fm.x[]*a.x[] - alpha.x[]*(p[] - p[-1])/Delta;
  boundary_flux ({gf});

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);
  boundary ((scalar *){g});
}

event projection (i++,last)
{
  mgp = project (uf, p, alpha, dt, mgp.nrelax);
  centered_gradient (p, g);

  correction (dt);
}

#if BRINKMAN_PENALIZATION
event brinkman_penalization(i++, last){
    brinkman_correction(u, uf, rho, dt);
}
#endif

#if TREE
event adapt (i++,last) {
#if EMBED
  fractions_cleanup (cs, fs);
  foreach_face()
    if (uf.x[] && !fs.x[])
      uf.x[] = 0.;
  boundary ((scalar *){uf});
#endif
  event ("properties");
}
#endif

void MinMaxValues(scalar * list, double * arr_eps) {// for each scalar min and max
  double arr[10][2];
  int ilist = 0;
  for (scalar s in list) {
    double mina= HUGE, maxa= -HUGE;
    foreach( reduction(min:mina) reduction(max:maxa) ){
      if (fabs(s[]) < mina) mina = fabs(s[]);
      if (fabs(s[]) > maxa) maxa = fabs(s[]);
    }
    arr[ilist][0] = mina;
    arr[ilist][1] = maxa;
    ilist++;
//        fprintf(stderr, "arr for i=%d", ilist);
  }

  for (int i = 0; i < ilist; i++){
#if EPS_MAXA == 1
    arr_eps[i] *=arr[i][1];
#else
    arr_eps[i] *= 0.5*(arr[i][0] + arr[i][1]);
#endif
#ifdef DEBUG_MINMAXVALUES
    fprintf(stderr, "MinMaxValues: i=%d, min=%g, max=%g, eps=%g\n", i, arr[i][0], arr[i][1], arr_eps[i]);
#endif
  }
}

stats statsf_weugene (scalar f, scalar fs)
{
    double dvr, min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
    foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
    reduction(max:max) reduction(min:min))
    if (fs[] < 1. && f[] != nodata) {
        dvr = dv()*(1. - fs[]);
        volume += dvr;
        sum    += dvr*f[];
        sum2   += dvr*sq(f[]);
        if (f[] > max) max = f[];
        if (f[] < min) min = f[];
    }
    stats s;
    s.min = min, s.max = max, s.sum = sum, s.volume = volume;
    if (volume > 0.)
        sum2 -= sum*sum/volume;
    s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
    return s;
}

stats statsf_weugene2 (scalar f, scalar fs)
{
    double dvr, min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
    foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
    reduction(max:max) reduction(min:min))
    if (fs[] == 0. && f[] != nodata) {
        dvr = dv()*(1. - fs[]);
        volume += dvr;
        sum    += dvr*f[];
        sum2   += dvr*sq(f[]);
        if (f[] > max) max = f[];
        if (f[] < min) min = f[];
    }
    stats s;
    s.min = min, s.max = max, s.sum = sum, s.volume = volume;
    if (volume > 0.)
        sum2 -= sum*sum/volume;
    s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
    return s;
}

norm normf_weugene (scalar f, scalar fs)
{
  double dvr, avg = 0., rms = 0., max = 0., volume = 0.;
  foreach(reduction(max:max) reduction(+:avg)
  reduction(+:rms) reduction(+:volume))
  if (fs[] < 1. && f[] != nodata) {
    dvr = dv()*(1. - fs[]);
    double v = fabs(f[]);
    if (v > max) max = v;
    volume += dvr;
    avg    += dvr*v;
    rms    += dvr*sq(v);
  }
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}

double change_weugene (scalar s, scalar sn, scalar fs)
{
  double max = 0.;
  foreach(reduction(max:max)) {
    if (fs[] < 1) {
      double ds = fabs (s[] - sn[]);
      if (ds > max)
        max = ds;
    }
    sn[] = s[];
  }
  return max;
}
