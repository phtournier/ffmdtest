---
name: fem-bem
category: bem
layout: example
---

~~~freefem
load "bem"
load "msh3"
load "iovtk"
load "PETSc-complex"
load "medit"
include "getARGV.idp"

real k0 = getARGV("-k0", 8.);
int NpWL = getARGV("-npwl", 25);
int prec = getARGV("-prec", 0);

include "mesh.idp"

// =======| Discretization in FE spaces |=======
load "Element_P3"
fespace Uh(Th, P2);
fespace UhBem(ThL,P2);
fespace UhOut(ThOut,P1);
fespace Vh=Uh*UhBem; // composite space
 
Uh hFEM = hTriangle;
if (mpirank == 0){
  cout << "size of FEM's mesh = " << hFEM[].max << endl;
  cout << "FEM's ndof = " << Uh.ndof << endl;
  cout << "BEM's ndof = " << UhBem.ndof << endl;
}
// =======| Definitions of BIO' kernels 
BemKernel kerSL("SL",k=k0);
BemKernel kerDL("DL",k=k0);
BemKernel kerTDL("TDL",k=k0);
BemKernel kerHS("HS",k=k0);

// =======| Definition of potentials |=======
varf SLpot(u,v) = int1d(ThL)(POT(BemPotential("SL",k=k0),u,v));
HMatrix<complex> HPotSL = SLpot(UhBem, UhOut);

varf DLpot(u,v) = int1d(ThL)(POT(BemPotential("DL",k=k0),u,v)); 
HMatrix<complex> HPotDL = DLpot(UhBem, UhOut);

// =======| Définition of incident wave
func uinc = exp(1i*k0*(cos(alpha)*x + sin(alpha)*y));
func dxuinc = 1i*k0*cos(alpha)*uinc;
func dyuinc = 1i*k0*sin(alpha)*uinc;
Uh<complex> gd = -uinc;

macro Grad(u) [dx(u), dy(u)] //

Uh<complex> ufem,v1;
UhBem<complex> ubem,v2;
~~~

## Simple Layer FEM-BEM coupling (Bielack-MacCamy) - Ansatz_SL file

$$
\begin{cases}
\begin{aligned}
\displaystyle a_{\mathrm{int}}(u,v)  % Volume term         
+\int_{\Gamma} \left(\frac{1}{2} I_d+ K_k^{\prime}\right) \phi\ \gamma_D^{\mathrm{int}}(v)  % Top right block &=& \int_{\Omega^{\mathrm{int}}} f\ v - \int_{\Gamma} \gamma_N^{\mathrm{ext}}(u_0) \gamma_D^{\mathrm{int}}(v) \\
\displaystyle \int_{\Gamma} \gamma_D^{\mathrm{int}}(u) \psi  % bottom left term
-\int_\Gamma V_k (\phi) \ \psi  % bottom right term &=& \int_{\Gamma} \gamma_D^{\mathrm{ext}}(u_0) \psi
\end{aligned}
\end{cases}
$$

~~~freefem
varf vfFemBemAnsatzSL(<[ufem],[ubem]>, <[v1],[v2]>) =
    int2d(Th)(-k0^2*ufem*v1 + Grad(ufem)'*Grad(v1)) // Top left term
  + int1dx1d(ThL)(ThL)(BEM(kerTDL,ubem,v1)) // Top right term
  + int1d(ThL)(0.5*ubem*v1) // Top right term
  - int1d(ThL)(ufem*v2)  // Bottom left term
  + int1dx1d(ThL)(ThL)(BEM(kerSL,ubem,v2))  // Bottom right term
  + on(labObstacle,ufem = gd) ; // Boundary condition

tgv = -1;

include "prec.idp"

string petscparams = "-ksp_rtol 1e-5 -ksp_type fgmres -ksp_side right "
+ "-pc_type fieldsplit -fieldsplit_0_pc_type lu -fieldsplit_0_pc_mat_solver_type mumps "
+ "-fieldsplit_1_ksp_type gmres -fieldsplit_1_pc_type none -fieldsplit_1_ksp_rtol 1e-2 "
+ "-fieldsplit_1_ksp_gmres_restart 1000 -ksp_converged_reason -ksp_view_final_residual";

Mat<complex> HCAnsatzSL = vfFemBemAnsatzSL(Vh,Vh,sparams=petscparams);

complex[int] AnsatzSLRHS = vfFemBemAnsatzSL(0,Vh);

complex[int] u = HCAnsatzSL^-1*AnsatzSLRHS;
[ufem[],ubem[]] = u;

// =======| Reconstruct u in unbounded domain |========
UhOut<complex> uext, uFstPart, uSndPart;

uext[] = HPotSL*ubem[];
uext = uext + uinc; // Compute the total wave
ufem = ufem + uinc;

UhOut<real> uextAbs = abs(uext);
Uh<real> ufemAbs = abs(ufem);
plot(ufemAbs,uextAbs,fill=1,value=1);
~~~

## Double Layer FEM-BEM coupling - Ansatz_DL file

$$
\begin{cases}
\begin{aligned}
\displaystyle a_{\mathrm{int}}(u,v)  % Volume term        
+\int_{\Gamma} W_k\ \phi\ \gamma_D^{\mathrm{int}}(v)  % Top right block
&=& \int_{\Omega^{\mathrm{int}}} f\ v - \int_{\Gamma} \gamma_N^{\mathrm{ext}}(u_0) \gamma_D^{\mathrm{int}}(v) \\
\displaystyle \int_{\Gamma} \gamma_D^{\mathrm{int}}(u) \psi  % bottom left term
-\int_\Gamma \left(\frac{1}{2} I_d + K_k\right) \phi \ \psi  % bottom right term 
&=& \int_{\Gamma} \gamma_D^{\mathrm{ext}}(u_0) \psi
\end{aligned}
\end{cases}
$$

~~~freefem
varf vfFemBemAnsatzDL(<[ufem],[phi]>,<[v1],[v2]>) =
    int2d(Th)(-k0^2 * ufem*v1 + Grad(ufem)'*Grad(v1)) // Volume term
  + int1dx1d(ThL)(ThL)(BEM(kerHS,phi,v1)) // Top right term
  - int1d(ThL)(ufem*v2) // Bottom left term
  + int1d(ThL)(0.5*phi*v2)                        // Bottom right term
  + int1dx1d(ThL)(ThL)(BEM(kerDL,phi,v2))  // Bottom right term
  + on(labObstacle,ufem = gd) ; // RHS

Mat<complex> HCAnsatzDL = vfFemBemAnsatzDL(Vh,Vh,sparams=petscparams);

complex[int] AnsatzDLRHS = vfFemBemAnsatzDL(0,Vh);

u = 0;
u = HCAnsatzDL^-1*AnsatzDLRHS;
[ufem[],ubem[]] = u;

// =======| Reconstruct u in unbounded domain |========
uext[] = HPotDL*ubem[];

uext = uext + uinc; // Compute the total wave
ufem = ufem + uinc;

uextAbs = abs(uext);
ufemAbs = abs(ufem);
plot(ufemAbs,uextAbs,fill=1,value=1);
~~~

## Ansatz combined FEM-BEM coupling - Ansatz_Combined file

$$
\begin{cases}
\begin{aligned}
\displaystyle a_{\mathrm{int}}(u,v)  % Volume term
+\int_{\Gamma} \left[W_k - i Z \left(\frac{1}{2} I_d+ K_k^{\prime}\right) \right] \phi\ \gamma_D^{\mathrm{int}}(v)  % Top right block
&=& \int_{\Omega^{\mathrm{int}}} f\ v - \int_{\Gamma} \gamma_N^{\mathrm{ext}}(u_0) \gamma_D^{\mathrm{int}}(v) \\
\displaystyle\int_{\Gamma} \gamma_D^{\mathrm{int}}(u) \psi  % bottom left term
-\int_\Gamma \left[\left(\frac{1}{2}I_d + K_k\right) - i Z V_k\right] \phi \ \psi  % bottom right term
&=& \int_{\Gamma} \gamma_D^{\mathrm{ext}}(u_0) \psi
\end{aligned}
\end{cases}
$$

~~~freefem
real eta = k0;
BemKernel KerAstar = -1i * eta * kerSL + kerDL;
BemKernel KerBstar = -1i * eta * kerTDL + kerHS;

varf vfFemBemAnsatzCombined(<[ufem],[phi]>,<[v1],[v2]>) =
    int2d(Th)(-k0^2 * ufem*v1 + Grad(ufem)'*Grad(v1)) // Volume term
  - int1d(ThL)(1i*eta*0.5 * phi*v1) // Top right term
  + int1dx1d(ThL)(ThL)(BEM(KerBstar,phi,v1)) // Top right term
  - int1d(ThL)(ufem*v2) // Bottom left term
  + int1d(ThL)(0.5*phi*v2) // Bottom right term
  + int1dx1d(ThL)(ThL)(BEM(KerAstar,phi,v2)) // Bottom right term
  + on(labObstacle, ufem = gd) ; // RHS

Mat<complex> HCAnsatzCombined = vfFemBemAnsatzCombined(Vh,Vh,sparams=petscparams);

complex[int] AnsatzCombinedRHS = vfFemBemAnsatzCombined(0,Vh);

u = 0;
u = HCAnsatzCombined^-1*AnsatzCombinedRHS;
[ufem[],ubem[]] = u; 

// =======| Reconstruct u in unbounded domain |========
uFstPart[] = HPotDL * ubem[];
uSndPart[] = HPotSL * ubem[];
uext[] = uFstPart[] - 1i * eta * uSndPart[]; 

uext = uext + uinc; // Compute the total wave
ufem = ufem + uinc;

uextAbs = abs(uext);
ufemAbs = abs(ufem);
plot(ufemAbs,uextAbs,fill=1,value=1);
~~~

## Direct JN (Johnson-Nédélec) - Direct_JN file
$\varphi \coloneqq \gamma_N^{\mathrm{ext}}(u)$

$$
\begin{cases}
\begin{aligned}
\displaystyle a_{\mathrm{int}}(u,v)  % Volume term
+\int_{\Gamma} \varphi \ \gamma_D^{\mathrm{int}}(v)  % Top right block
&=& \int_{\Omega^{\mathrm{int}}} f\ v - \int_{\Gamma} \gamma_N^{\mathrm{ext}}(u_0) \gamma_D^{\mathrm{int}}(v) \\
\displaystyle \int_{\Gamma} \left(\frac{1}{2}I_d - K_k\right)\gamma_D^{\mathrm{int}}(u) \psi  % bottom left term
-\int_\Gamma V_k (\varphi) \ \psi  % bottom right term
&=& \int_{\Gamma} \left(\frac{1}{2}I_d - K_k\right) \gamma_D^{\mathrm{ext}}(u_0) \psi
\end{aligned}
\end{cases}
$$

~~~freefem
varf vfFemBemDirectJN(<[ufem],[t]>,<[v1],[v2]>) =
    int2d(Th)(-k0^2 * ufem*v1 + Grad(ufem)'*Grad(v1)) // Volume term
  + int1d(ThL)(t*v1) // Top right term
  + int1d(ThL)((-0.5)*ufem*v2) // Bottom left term
  + int1dx1d(ThL)(ThL)(BEM(kerDL,ufem,v2)) // Bottom right term
  + int1dx1d(ThL)(ThL)(BEM(kerSL,t,v2))  // Bottom right term
  + on(labObstacle,ufem = gd) ; // RHS

Mat<complex> HCDirectJN = vfFemBemDirectJN(Vh,Vh,sparams=petscparams);

complex[int] DirectJNRHS = vfFemBemDirectJN(0,Vh);

u = 0;
u = HCDirectJN^-1*DirectJNRHS;
[ufem[],ubem[]] = u;

// =======| Reconstruct u in unbounded domain |========
UhBem<complex> ufemTrace;
ufemTrace = ufem;

uFstPart[] = HPotDL * ufemTrace[];
uSndPart[] = HPotSL * ubem[];
uext = uFstPart + uSndPart;

uext = uext + uinc; // Compute the total wave
ufem = ufem + uinc;

uextAbs = abs(uext);
ufemAbs = abs(ufem);
plot(ufemAbs,uextAbs,fill=1,value=1);
~~~

## Direct TraceNeu - Direct_TraceNeu file
$\varphi \coloneqq \gamma_N^{\mathrm{ext}}(u)$

$$
\begin{cases}
\begin{aligned}
\displaystyle a_{\mathrm{int}}(u,v)  % Volume term
+\int_{\Gamma} \varphi \ \gamma_D^{\mathrm{int}}(v)  % Top right block
&=& \int_{\Omega^{\mathrm{int}}} f\ v - \int_{\Gamma} \gamma_N^{\mathrm{ext}}(u_0) \gamma_D^{\mathrm{int}}(v) \\
\displaystyle \int_{\Gamma} W_k \gamma_D^{\mathrm{int}}(u) \psi  % bottom left term
-\int_\Gamma \left(\frac{1}{2}I_d - K_k^{\prime}\right) (\varphi) \ \psi  % bottom right term
&=& \int_{\Gamma} W_k \gamma_D^{\mathrm{ext}}(u_0) \psi
\end{aligned}
\end{cases}
$$

~~~freefem
varf vfFemBemDirectTraceNeu(<[ufem],[t]>,<[v1],[v2]>) =
    int2d(Th)(-k0^2 * ufem*v1 + Grad(ufem)'*Grad(v1)) // Volume term
  + int1d(ThL)(t*v1) // Top right term
  + int1dx1d(ThL)(ThL)(BEM(kerHS,ufem,v2)) // Bottom left term
  + int1d(ThL)((-0.5)*t*v2) // Bottom right term 
  + int1dx1d(ThL)(ThL)(BEM(kerTDL,t,v2))  // Bottom right term
  + on(labObstacle,ufem = gd) ; // RHS

Mat<complex> HCDirectTraceNeu = vfFemBemDirectTraceNeu(Vh,Vh,sparams=petscparams);

complex[int] DirectTraceNeuRHS = vfFemBemDirectTraceNeu(0,Vh);

u = 0;
u = HCDirectTraceNeu^-1*DirectTraceNeuRHS;
[ufem[],ubem[]] = u;

// =======| Reconstruct u in unbounded domain |========
ufemTrace = ufem;

uFstPart[] = HPotDL * ufemTrace[];
uSndPart[] = HPotSL * ubem[];
uext = uFstPart + uSndPart;

uext = uext + uinc; // Compute the total wave
ufem = ufem + uinc;

uextAbs = abs(uext);
ufemAbs = abs(ufem);
plot(ufemAbs,uextAbs,fill=1,value=1);
~~~

## Direct Combined - Direct_Combined file
$\varphi \coloneqq \gamma_N^{\mathrm{ext}}(u)$

$$
\begin{cases}
\begin{aligned}
\displaystyle a_{\mathrm{int}}(u,v)  % Volume term         
+\int_{\Gamma} \varphi \ \gamma_D^{\mathrm{int}}(v)  % Top right block
&=& \int_{\Omega^{\mathrm{int}}} f\ v - \int_{\Gamma} \gamma_N^{\mathrm{ext}}(u_0) \gamma_D^{\mathrm{int}}(v) \\
\displaystyle \int_{\Gamma} \left[W_k - i  \left(\frac{1}{2} I_d - K_k\right) Z \right]  \ \gamma_D^{\mathrm{int}}(u) \psi  % bottom left term
-\int_\Gamma \left[\left(\frac{1}{2}I_d - K_k^{ \prime}\right) - i \ V_k\ Z\right] (\varphi) \ \psi  % bottom right term                                                                  
&=& \int_{\Gamma} \left[W_k - i  \left(\frac{1}{2} I_d - K_k\right) Z \right] \gamma_D^{\mathrm{ext}}(u_0) \psi
\end{aligned}
\end{cases}
$$

~~~freefem
eta = k0;
BemKernel KerA = 1i * eta * kerSL + kerTDL;
BemKernel KerB = 1i * eta * kerDL + kerHS;

varf vfFemBemDirectCombined(<[ufem],[phi]>,<[v1],[v2]>) =
    int2d(Th)(-k0^2 * ufem*v1 + Grad(ufem)'*Grad(v1)) // Volume term
  + int1d(ThL)(phi*v1) // Top right term
  + int1d(ThL)((-1i*eta * 0.5)*ufem*v2)   // Bottom left term
  + int1dx1d(ThL)(ThL)(BEM(KerB,ufem,v2))  // Bottom left term
  + int1d(ThL)((-0.5)*phi*v2) // Bottom right term
  + int1dx1d(ThL)(ThL)(BEM(KerA,phi,v2))  // Bottom right term
  + on(labObstacle,ufem = gd) ; // RHS

Mat<complex> HCDirectCombined = vfFemBemDirectCombined(Vh,Vh,sparams=petscparams);

complex[int] DirectCombinedRHS = vfFemBemDirectCombined(0,Vh);

u = 0;
u = HCDirectCombined^-1*DirectCombinedRHS;
[ufem[],ubem[]] = u;

// =======| Reconstruct u in unbounded domain |========
ufemTrace = ufem;

uFstPart[] = HPotDL * ufemTrace[];
uSndPart[] = HPotSL * ubem[];
uext = uFstPart + uSndPart;

uext = uext + uinc; // Compute the total wave
ufem = ufem + uinc;

uextAbs = abs(uext);
ufemAbs = abs(ufem);
plot(ufemAbs,uextAbs,fill=1,value=1);
~~~
