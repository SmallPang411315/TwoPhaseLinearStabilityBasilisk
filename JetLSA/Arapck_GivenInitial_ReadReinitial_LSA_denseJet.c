#include <complex.h>  // creal, cimag.
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "arpack.h"
#include "debug_c.h"  // debug arpack.
#include "stat_c.h"   // arpack statistics.

#include "grid/quadtree.h"
#include "navier-stokes/centered-LSA.h"
#include "two-phase-LS.h"
#include "tracer-LSA.h"
#include "tension-LSA.h"
#include "view.h"
#include "tecplot.h"

#define Re 316.
#define SIGMA 0.4
#define L0 32.
#define Nele 32768

double t_LSA = 0.;
double tend  = 1.5;
double delta_t = 0.05;

scalar phi_stored[];
vector u_stored[];
scalar omega[], f0[];
scalar * tracers = {phi};

double  invec[Nele* 3];
double  outvec[Nele* 3];

double Xc[Nele], Yc[Nele];

void set_invec(a_int  ntot, double * work_in)
{
  for (int i = 0; i < (int) ntot; i++){
    invec[i] = work_in[i];
  }
}

void put_outvec(a_int  ntot, double * work_out)
{
  for (int i = 0; i < (int) ntot; i++){
    work_out[i] = outvec[i];
  }
}

int main() { 
  debug_c(6, -6, 1,
          1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1); // set debug flags.

  X0 = 0.0;
  Y0 = -2.;
  N = 512;
  init_grid (N);
  origin (X0, Y0);
  size (L0);

  /**
  We set the density and viscosity of each phase as well as the
  surface tension coefficient and start the simulation. */

  rho1 = 1., rho2 = 0.5;
  mu1 = 1./Re, mu2 = 1./Re;
  f.sigma = SIGMA;

  mask(y>2?top:none);
  mask(y<-2?bottom:none);

  input_grd_U(u_stored.x, file = "/mnt/d/basilisk/gongminjiang/paper2/NewJetLSA/r0.5/We-1_0.1/initial/ux.asc");
  input_grd_V(u_stored.y, file = "/mnt/d/basilisk/gongminjiang/paper2/NewJetLSA/r0.5/We-1_0.1/initial/uy.asc");
  input_grd_Phi(phi_stored, file = "/mnt/d/basilisk/gongminjiang/paper2/NewJetLSA/r0.5/We-1_0.1/initial/phi.asc");

  output_ppm (u_stored.x, file = "u_storedx.png");
  output_ppm (u_stored.y, file = "u_storedy.png");
  output_ppm (phi_stored, file = "phi_stored.png");

  double u_storedxd[Nele],u_storedyd[Nele],phi_storedd[Nele],uvphi_storedd[Nele*3];

  int id_stored = 0;
  scalar idc_stored[];
  foreach(){
    id_stored++;
    idc_stored[] = id_stored;
  }

  foreach(){
    int idcc_stored          = (int)idc_stored[] - 1;
    u_storedxd[idcc_stored]  = u_stored.x[];
    u_storedyd[idcc_stored]  = u_stored.y[];
    phi_storedd[idcc_stored] = phi_stored[];
  }

  for (int ii = 0; ii<Nele; ii++)          uvphi_storedd[ii] = u_storedxd[ii];
  for (int ii = Nele; ii<2*Nele; ii++)     uvphi_storedd[ii] = u_storedyd[ii-Nele];
  for (int ii = 2*Nele; ii<3*Nele; ii++)   uvphi_storedd[ii] = phi_storedd[ii-2*Nele];

  a_int ido = 0;
  char bmat[] = "I";
  a_int ntot = Nele*3;
  char which[] = "LM";
  a_int nvec = 2;
  double evtol = 1e-6; // small evtol => more stable checks after EV computation.
  double resid[ntot];
  a_int kdim = 32;//4 * nvec;
  //double V[kdim * ntot];
  double * V = (double *)malloc((kdim*ntot) * sizeof(double));//need memory free
  a_int iparam[11];
  a_int ipntr[14];
  double workd[3 * ntot]; 
  a_int rvec = 1;
  char howmny[] = "A";
  double * dr =
      (double *)malloc((nvec + 1) * sizeof(double));
  double * di =
      (double *)malloc((nvec + 1) * sizeof(double)); 
  a_int select[kdim];
  int i; // C99 compliant.
  //for (i = 0; i < kdim; i++) select[i] = 1;
  //double z[(ntot + 1) * (nvec + 1)];
  double * z =  (double *)malloc((ntot + 1) * (nvec + 1) * sizeof(double));//need memory free
  double sigmar = 0.;
  double sigmai = 0.;
  int k;
  for (k = 0; k < 3 *  ntot; ++k) workd[k] = 0;
  double workl[3 * (kdim * kdim) + 6 * kdim];
  for (k = 0; k < 3 * (kdim * kdim) + 6 * kdim; ++k) workl[k] = 0;
  a_int lworkl = 3 * (kdim * kdim) + 6 * kdim;
  //double rwork[kdim];
  double workev[3 * kdim];
  a_int info = 1;//0 when no specific starting vector
                 //1 when setting the starting vector
  for (i = 0; i <  ntot; i++) resid[i] = uvphi_storedd[i];//starting vector

  iparam[0] = 1;
  iparam[2] = 1000;// max iter number
  iparam[3] = 1;
  //iparam[4] = 0;  // number of ev found by arpack.
  iparam[6] = 1;

  int numiter = 0;
  char sname[40];
  
  dnaupd_c(&ido, bmat, ntot, which, nvec, evtol, resid, kdim, V, ntot, iparam, ipntr, workd, workl, lworkl, &info);
  
  while (ido == -1 || ido == 1) {
//     FILE * FixPoiHis;
//     FixPoiHis = fopen ("FixPoiHis", "w+");

    set_invec(ntot, &(workd[ipntr[0] - 1]));
    run();
//     exit(1);
    put_outvec(ntot, &(workd[ipntr[1] - 1]));
    
//     fclose(FixPoiHis);
        
    dnaupd_c(&ido, bmat, ntot, which, nvec, evtol, resid, kdim, V, ntot, iparam, ipntr, workd, workl, lworkl, &info);

    numiter++;
    fprintf (stderr, "\n numiter = %d, info = %d, ido = %d\n\n", numiter, info, ido);

    if(numiter%20==0){
      FILE * EigenVec;
      sprintf(sname, "EigenVec_%d", numiter);
      EigenVec= fopen (sname, "w+");
      fprintf(EigenVec,"eigennum x y u v phi\n");
      for (i = 1; i < kdim; i++){
        for (int j = 0; j < ntot/3; j++){
          fprintf(EigenVec, "%d %g %g %g %g %g\n", i, Xc[j], Yc[j], V[i* ntot + j], V[i*ntot + j + ntot/3], V[i*ntot + j + 2*ntot/3]);
        }
        fprintf(EigenVec, "\n");
      }
      fclose(EigenVec);
    }
  }
  //if (iparam[4] != nvec) {printf("Error: iparam[4] %ld, nvec %ld\n", iparam[4], nvec); return 1;}//check number of ev found by arpack.
        
  dneupd_c(rvec, howmny, select, dr, di, z, ntot, sigmar, sigmai, workev, bmat, ntot, which, nvec, evtol, resid, kdim, V, ntot, iparam, ipntr, workd, workl, lworkl, &info);      
  
  FILE * EigenVec;
  EigenVec= fopen ("EigenVec", "w+");
  fprintf(EigenVec,"eigennum x y u v phi\n");
  for (i = 0; i < kdim; i++){
      for (int j = 0; j < ntot/3; j++){
        fprintf(EigenVec, "%d %g %g %g %g %g\n", i, Xc[j], Yc[j], V[i* ntot + j], V[i*ntot + j + ntot/3], V[i*ntot + j + 2*ntot/3]);
    }
    fprintf(EigenVec, "\n");
  }
  fclose(EigenVec);
  free(dr);
  free(di);
  free(V);
  free(z);
  return(0);
}

u.n[left]  = dirichlet (0);
u.t[left]  = dirichlet (0);
p[left]    = neumann (0);
Ubase.n[left] = neumann(0.);//dirichlet(fabs(y)<=1.?(1+Lam)/(1-Lam):1.);
Ubase.t[left] = neumann(0.);
// Phibase[left] = neumann(0.0);
phi[left]     = dirichlet (0);
// f[left]    = f0[];//fabs(y)<1.?1.:0.;
f[left]       = neumann(0.);

u.n[top]   = dirichlet(0);
u.t[top]   = dirichlet(0);
p[top]     = neumann(0);
Ubase.n[top] = neumann(0.);
Ubase.t[top] = neumann(0.);
// Phibase[top] = neumann(1.);
phi[top]     = dirichlet (0);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
pf[right]  = dirichlet (0.);//右边的uy没有震荡了
p[right]   = dirichlet (0.);
Ubase.n[right] = neumann(0.);
Ubase.t[right] = neumann(0.);
Phibase[right] = neumann(0.);
phi[right]     = dirichlet (0);
f[right]       = neumann(0.);

u.n[bottom]   = dirichlet(0);
u.t[bottom]   = dirichlet(0);
p[bottom]     = neumann(0);
Ubase.n[bottom] = neumann(0);
Ubase.t[bottom] = neumann(0);
// Phibase[bottom] = neumann(1.);
phi[bottom]     = dirichlet (0);

event init (t = 0)
{ 
  mask(y>2?top:none);
  mask(y<-2?bottom:none);
      fraction (f0, -sq(y)+sq(1.));
      input_grd_U(Ubase.x, file = "/mnt/d/basilisk/gongminjiang/paper2/NewJet/jet/BaseFlow/r0.5/Lam1.3/mask/U.asc");
      input_grd_V(Ubase.y, file = "/mnt/d/basilisk/gongminjiang/paper2/NewJet/jet/BaseFlow/r0.5/Lam1.3/mask/V.asc");
      input_grd_Phi(Phibase, file = "/mnt/d/basilisk/gongminjiang/paper2/NewJet/jet/BaseFlow/r0.5/Lam1.3/mask/Phibase_reinitialization.asc");
      //input_grd_f(f, file = "/mnt/d/basilisk/gongminjiang/paper2/NewJet/jet/BaseFlow/r0.5/Lam1.3/mask/f.asc");
	  foreach() f[]= cm[]*( value_Heaviside(-Phibase[], Delta)*(1 - 0) + 0 )

  boundary ((scalar *){Ubase,Phibase,f});

  int id = 0;
  scalar idc[];
  foreach(){
    id++;
    idc[] = id;
  }
  
  foreach(){
    int idcc = (int)idc[] - 1;
    u.x[] = invec[idcc];
    u.y[] = invec[idcc + Nele];
    phi[] = invec[idcc + 2*Nele];
  }
  boundary((scalar *){u,phi});
}

event set_outvec (t = end)
{ 
  int id = 0;
  scalar idc[];
  foreach(){
    id++;
    idc[] = id;
  }
  
  foreach(){
    int idcc = (int)idc[] - 1;
    outvec[idcc]          = u.x[];
    outvec[idcc + Nele]   = u.y[];
    outvec[idcc + 2*Nele] = phi[];
  }
  
  foreach(){
    int idcc = (int)idc[] - 1;
    Xc[idcc] = x;
    Yc[idcc] = y;
  }
}

event logfile (i++){
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

//   foreach(){
//       if(fabs(x-10.)<0.5 && fabs(y-10)<0.5){
//         fprintf(FixPoiHis, "%g %g %g %g %g\n", t, x, y, u.x[], u.y[]);
//       }
//   }
}

/*event movie (t = 0.; t+=5e-2; t <= tend)
{
  output_ppm (Ubase.x, file = "Ux.mp4", spread = -1, linear = true, n = 1024, box = {{0.,-2},{L0,2}});
  output_ppm (Ubase.y, file = "Uy.mp4", spread = -1, linear = true, n = 1024, box = {{0.,-2},{L0,2}});
  output_ppm (u.x, file = "ux.mp4", spread = -1, linear = true, n = 1024,   box = {{0.,-2},{L0,2}});
  output_ppm (u.y, file = "uy.mp4", spread = -1, linear = true, n = 1024,   box = {{0.,-2},{L0,2}});
  output_ppm (phi, file = "phi.mp4", spread = -1, linear = true, n = 1024,   box = {{0.,-2},{L0,2}});
  output_ppm (f, file = "f.mp4", spread = -1, linear = true, n = 1024,	      box = {{0.,-2},{L0,2}});
}

event snapshot (t = 0.; t+=delta_t; t <= tend) {
 vorticity (u, omega);
 char name[80];
 sprintf (name, "snapshot-%g", t);
 scalar pid[];
 foreach()
 pid[] = fmod(pid()*(npe() + 37), npe());
 boundary ({pid});
 dump (name);
}*/

// event out_grd(i=end)
// {
//   FILE * fpux = fopen("ux.asc","w");
//   output_grd(u.x, fpux, 1./N, true, {{0,-2},{L0,2}});
//   fclose(fpux);
//   FILE * fpuy = fopen("uy.asc","w");
//   output_grd(u.y, fpuy, 1./N, true, {{0,-2},{L0,2}});
//   fclose(fpuy);
//   FILE * fpphi = fopen("phi.asc","w");
//   output_grd(phi, fpphi, 1./N, true, {{0,-2},{L0,2}});
//   fclose(fpphi);
// }

event output(t = 0; t += delta_t; t <= tend){
   struct OutputTec tec;
//Give primary variables to the tec_cc list. "cc" means "cell center".
  tec.tec_cc = list_copy({p, rho, phi, dphi, Phibase, dPhibase, dPhibasex, dPhibasey, u, Ubase, dphix, dphiy, fPhi, gPhi, fx, fy, gx, gy});
//Give the name of variables to varname.
#if dimension>1
  sprintf(tec.varname,"x y p rho phi dphix dphiy Phibase dPhibase.x dPhibase.y dPhibasexx dPhibasexy dPhibaseyx dPhibaseyy u.x u.y Ubase.x Ubase.y dphixx dphixy dphiyx dphiyy fPhi gPhi fx fy gx gy");
#endif

//Extra variables can be output by adding them to the list of scalar and name.
  tec.tec_cc = list_concat(tec.tec_cc,{f});
  sprintf(tec.extname," f");

  output_tec(tec,i);
}
/*
event output(t = 0.00; t += delta_t; t <= tend){
	double ymin=0.001,ymax=19.999,xmin=-4.999,xmax=14.999;
	int ny=N,nx=N;
	char name[40];
	sprintf (name, "tec_%d.plt", i);
	output_field_2D_Matlab_zslice ( (scalar*) {p,u.x,u.y,Ubase.x,Ubase.y}, fopen(name, "w") , nx, ny, xmin, xmax, ymin, ymax);
}*/
