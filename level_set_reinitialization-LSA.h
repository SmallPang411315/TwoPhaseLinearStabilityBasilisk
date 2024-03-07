/**
	2022.03.11 (GONGMINJIANG)
		A re-initialization precudure of Level set function \phi' in 
		LSA frame is doing in this header file.
*/

/**
	Compute Re-distanced function \phi'
*/

#include "tecplot.h"

struct ReInitial {
  scalar * tracers;
  scalar Phibase;
  double dt;
  scalar * src; // optional
};


void tracer_RDF_linearphi (scalar f,
		    scalar Phibase,
		    scalar rhs,
		    (const) scalar src)
{
  vector g[], gBase[];
  scalar a[], ap[], am[], 
         b[], bp[], bm[],
         c[], cp[], cm[],
         d[], dp[], dm[],
		 a_base[], ap_base[], am_base[], 
         b_base[], bp_base[], bm_base[],
         c_base[], cp_base[], cm_base[],
         d_base[], dp_base[], dm_base[];
  boundary({f, Phibase});
  foreach(){
	  f.gradient       = minmod2;
	  Phibase.gradient = minmod2;
  }
  gradients ({f}, {g});
  gradients ({Phibase}, {gBase});
  
  foreach_boundary (left)   gBase.x[] = (Phibase[1,0]-Phibase[0,0])/Delta;
  foreach_boundary (top)    gBase.y[] = (Phibase[0,0]-Phibase[0,-1])/Delta;
  foreach_boundary (right)  gBase.x[] = (Phibase[0,0]-Phibase[-1,0])/Delta;
  foreach_boundary (bottom) gBase.y[] = (Phibase[0,1]-Phibase[0,0])/Delta;

  foreach() rhs[] = Phibase[]*( g.x[]*gBase.x[] + g.y[]*gBase.y[] );
  
//   foreach() fprintf(stdout, "rhs=%g\n", rhs[]);
//   exit(1);
//   boundary({rhs});
}  

/**
	Compute Re-distanced function \phi
*/
void redistance_Linearphi (struct ReInitial p)
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
		scalar phiPre[], rhs[], dtp[];
		double tpStep     = p.dt;
		double RedistErr  = 100.;
		
		int numiter= 0 ;
        while(RedistErr > 1e-4){
			RedistErr     = 0.0;
			tracer_RDF_linearphi (f, p.Phibase, rhs, src);

			foreach() {
				phiPre[]  = f[];
				f[]      -= tpStep*rhs[];
			}

			foreach(reduction(+:RedistErr)){
				RedistErr += sq(f[]-phiPre[]);
			}
//
			RedistErr      = sqrt(RedistErr);
			fprintf(stdout, "RedistErr=%g dtp=%g\n", RedistErr, tpStep);
// 			exit(1);
// 			if(++numiter%50==0){
// 				struct OutputTec tec;
// 				//Give primary variables to the tec_cc list. "cc" means "cell center".
// 				  tec.tec_cc = list_copy({f});
// 				//Give the name of variables to varname.
// 				#if dimension>1
// 				  sprintf(tec.varname,"x y f");
// 				#endif
// 				  output_tec(tec,numiter);
// // 				exit(1);
// 			}
		}
	}
	  boundary (p.tracers);

	  if (!p.src)
		free (lsrc);
} 
