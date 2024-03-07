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
// 	vector dPhibase[];
// 	gradients({f},{dPhibase});
// 	foreach(){
// 		rhs[] = sqrt(sq(dPhibase.x[])+sq(dPhibase.y[]))-1;
// 		printf("rhs=%g\n",rhs[]);
// 	}

  scalar a[], ap[], am[], 
         b[], bp[], bm[],
         c[], cp[], cm[],
         d[], dp[], dm[];

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

  foreach_boundary (top){
	f[ghost]= 2*f[0,0]-f[0,-1];
	a[]     = (f[0,0]-f[-1,0])/Delta;
    ap[]    = max(a[],0.0);           //a+
    am[]    = min(a[],0.0);           //a-
	b[]     = (f[1,0]-f[0,0])/Delta;
	bp[]    = max(b[],0.0);
	bm[]    = min(b[],0.0);
	c[]     = (f[0,0]-f[0,-1])/Delta;
	cp[]    = max(c[],0.0);
	cm[]    = min(c[],0.0);
	d[]     = (f[ghost]-f[0,0])/Delta;
	dp[]    = max(d[],0.0);
	dm[]    = min(d[],0.0);
  }

  foreach_boundary (bottom){
	f[ghost]= 2*f[0,0]-f[0,1];
	a[]     = (f[0,0]-f[-1,0])/Delta;
    ap[]    = max(a[],0.0);           //a+
    am[]    = min(a[],0.0);           //a-
	b[]     = (f[1,0]-f[0,0])/Delta;
	bp[]    = max(b[],0.0);
	bm[]    = min(b[],0.0);
	c[]     = (f[0,0]-f[ghost])/Delta;
	cp[]    = max(c[],0.0);
	cm[]    = min(c[],0.0);
	d[]     = (f[0,1]-f[0,0])/Delta;
	dp[]    = max(d[],0.0);
	dm[]    = min(d[],0.0);
  }

  foreach_boundary (left){
	f[ghost]= 2*f[0,0]-f[1,0];
	a[]     = (f[0,0]-f[ghost])/Delta;
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
  }

/*
  foreach_boundary (right){
	b[]     = (f[0,0]-f[-1,0])/Delta;
	bp[]    = max(b[],0.0);
	bm[]    = min(b[],0.0);
  }
*/
  foreach(){
	  if(y>0.){
		  rhs[] = 0.0 ;

		  if(sgnphi[]>Delta) rhs[] = sqrt( max(sq(ap[]),sq(bm[])) + max(sq(cp[]),sq(dm[])) ) - 1 ;

		  if(sgnphi[]<-Delta) rhs[] = sqrt( max(sq(am[]),sq(bp[])) + max(sq(cm[]),sq(dp[])) ) - 1 ;
	}
	else{
		rhs[] = 0.0 ;
	}
  }
  foreach_boundary(left) rhs[ghost] = rhs[];
//   boundary({rhs});
}  

/**
	Compute Re-distanced function \phi
*/
void redistance (struct Advection_Redistance p)
{  
	scalar * lsrc = p.src;
	if (!lsrc) {
		const scalar zero[] = 0.;
		for (scalar s in p.tracers)
		  lsrc = list_append (lsrc, zero);
	  }

	  assert (list_len(p.tracers) == list_len(lsrc));
	  scalar frac, f, src;
	  for (frac, f, src in p.interfaces,p.tracers,lsrc) {
		vector g[];
		scalar phiPre[], sgnphi[], rhs[], s[], dtp[];
		double tpStep     = p.dt/100.;
		double RedistErr  = 100.;
		double maxGradf   = -100.;
// 		int itenum        = 0;
		
		foreach(){
			sgnphi[] = -2.*frac[]+1;//好一些
// 			sgnphi[] = value_Heaviside(f[], Delta);
// 			printf("x=%g,y=%g,sgnphi=%g\n",x,y,sgnphi[]);
		}

        while(RedistErr > 1e-4){
			RedistErr     = 0.0;
			boundary({f});
			tracer_RDF (f, sgnphi, rhs, s, p.dt, src);
			gradients({f},{g});
			foreach() {
				phiPre[]  = f[];
				if(y>0.) {
					f[]      -= tpStep*sgnphi[]*rhs[];
// 					printf("x=%g,y=%g,gradPhi=%g,phiPre=%g,f=%g\n",x,y,sqrt(sq(g.x[])+sq(g.y[])),phiPre[],f[]);
				}
			}

// 			itenum++;
// 			char name[80];
// 			sprintf (name, "Phibase-%d", itenum);
// 			output_ppm (Phibase, file = name);

			foreach(reduction(max:maxGradf)){
				if (sqrt(sq(g.x[])+sq(g.y[])) > maxGradf) maxGradf = sqrt(sq(g.x[])+sq(g.y[]));
			}

			foreach(reduction(+:RedistErr)){
				RedistErr += sq(f[]-phiPre[]);
			}

			RedistErr      = sqrt(RedistErr);
// 			fprintf(stderr, "RedistErr=%g,maxGradf=%g,dtp=%g\n", RedistErr, maxGradf, tpStep);
		}
	}
	  boundary (p.tracers);

	  if (!p.src)
		free (lsrc);
} 
