/**
	2022.03.09 (GONGMINJIANG)
		A re-initialization precudure of Level set function \phi is doing
		in this header file.
*/

/**
	Compute Re-distanced function \phi
*/

void tracer_RDF (scalar f,
		    scalar sgnphi,
		    scalar rhs,
		    scalar s,
		    double dt,
		    (const) scalar src)
{
  vector g[];
  scalar a[], ap[], am[], 
         b[], bp[], bm[],
         c[], cp[], cm[],
         d[], dp[], dm[];
  boundary({f});
  foreach() f.gradient = minmod2;
  gradients ({f}, {g});
//   foreach_boundary (left)   g.x[] = (f[1,0]-f[0,0])/Delta;
//   foreach_boundary (top)    g.y[] = (f[0,0]-f[0,-1])/Delta;
//   foreach_boundary (right)  g.x[] = (f[0,0]-f[-1,0])/Delta;
//   foreach_boundary (bottom) g.y[] = (f[0,1]-f[0,0])/Delta;

  foreach(){
    a[]     = (f[0,0]-f[-1,0])/Delta;
    ap[]    = max(a[],0.0);           //a+
    am[]    = min(a[],0.0);           //a-
	b[]     = (f[1,0]-f[0,0])/Delta;
	bp[]    = max(b[],0.0);
	bm[]    = min(b[],0.0);
	c[]     = (f[0,0]-f[0,-1])/Delta;
	cp[]    = max(c[],0.0);
	cm[]    = min(c[],0.0);
	d[]     = (f[0,1]-f[0,0])/Delta;
	dp[]    = max(d[],0.0);
	dm[]    = min(d[],0.0);
	s[]     = f[]/sqrt( sq(f[]) + sq(Delta) );
  }
//   foreach_boundary (left){
// 	a[]     = (f[1,0]-f[0,0])/Delta;
//     ap[]    = max(a[],0.0);           //a+
//     am[]    = min(a[],0.0);           //a-
//   }
//
//   foreach_boundary (top){
// 	d[]     = (f[0,0]-f[0,-1])/Delta;
// 	dp[]    = max(d[],0.0);
// 	dm[]    = min(d[],0.0);
//   }
//   foreach_boundary (right){
// 	b[]     = (f[0,0]-f[-1,0])/Delta;
// 	bp[]    = max(b[],0.0);
// 	bm[]    = min(b[],0.0);
//   }
//   foreach_boundary (bottom){
// 	c[]     = (f[0,1]-f[0,0])/Delta;
// 	cp[]    = max(c[],0.0);
// 	cm[]    = min(c[],0.0);
//   }
  foreach(){
	if(sgnphi[]>Delta)      rhs[] = sqrt( max(sq(ap[]),sq(bm[])) + max(sq(cp[]),sq(dm[])) ) - 1 ;
// 	if(sgnphi[]>1e-9)      rhs[] = sqrt( max(sq(ap[]),sq(bm[])) + max(sq(cp[]),sq(dm[])) ) - 1 ;
    else
		if(sgnphi[]<-Delta) rhs[] = sqrt( max(sq(am[]),sq(bp[])) + max(sq(cm[]),sq(dp[])) ) - 1 ;
// 		if(sgnphi[]<-1e-9) rhs[] = sqrt( max(sq(am[]),sq(bp[])) + max(sq(cm[]),sq(dp[])) ) - 1 ;
		else               rhs[] = 0.0 ;
  }
  boundary({rhs});
}  

/**
	Compute Re-distanced function \phi
*/
void redistance (struct Advection p)
{  
	scalar * lsrc = p.src;
	if (!lsrc) {
		const scalar zero[] = 0.;
		for (scalar s in p.tracers)
		  lsrc = list_append (lsrc, zero);
	  }

	  assert (list_len(p.tracers) == list_len(lsrc));
	  scalar f, src;
	  for (f,src in p.tracers,lsrc) {
		vector g[];
		scalar phiPre[], sgnphi[], rhs[], s[], dtp[];
		double tpStep     = p.dt/10;
		double RedistErr  = 100.;
		
		foreach() {
			if ( f[] < -Delta) sgnphi[] = -1.0;
				else
					if(f[] > Delta)  sgnphi[] = 1.0;
					else  sgnphi[] = 0.0;  
		}
        while(RedistErr > 5e-2){
			RedistErr     = 0.0;
			tracer_RDF (f, sgnphi, rhs, s, p.dt, src);
			foreach() {
				phiPre[]  = f[];
				f[]      -= tpStep*sgnphi[]*rhs[];
			}
			foreach(reduction(+:RedistErr)){
				RedistErr += sq(f[]-phiPre[]);
			}
			RedistErr      = sqrt(RedistErr);
// 			fprintf(stderr, "RedistErr=%g dtp=%g\n", RedistErr, tpStep);
		}
	}
	  boundary (p.tracers);

	  if (!p.src)
		free (lsrc);
} 
