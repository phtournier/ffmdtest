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
\delta
$$

$$
\begin{cases}
\begin{aligned}
\displaystyle a_{\mathrm{int}}(u,v)  % Volume term         
+ \int_{\Gamma} \left(\frac{1}{2} I_d+ K_k^{\prime}\right) \phi\ \gamma_D^{\mathrm{int}}(v)  % Top right block
&=& \int_{\Omega^{\mathrm{int}}} f\ v - \int_{\Gamma} \gamma_N^{\mathrm{ext}}(u_0) \gamma_D^{\mathrm{int}}(v) \\
\displaystyle \int_{\Gamma} \gamma_D^{\mathrm{int}}(u) \psi  % bottom left term
-  \int_\Gamma V_k (\phi) \ \psi  % bottom right term
&=& \int_{\Gamma} \gamma_D^{\mathrm{ext}}(u_0) \psi
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

\begin{equation}
\begin{cases}
\begin{aligned}
\displaystyle a_{\mathrm{int}}(u,v)  % Volume term        
+ \int_{\Gamma} W_k\ \phi\ \gamma_D^{\mathrm{int}}(v)  % Top right block
&=& \int_{\Omega^{\mathrm{int}}} f\ v - \int_{\Gamma} \gamma_N^{\mathrm{ext}}(u_0) \gamma_D^{\mathrm{int}}(v) \\
\displaystyle \int_{\Gamma} \gamma_D^{\mathrm{int}}(u) \psi  % bottom left term
- \int_\Gamma \left(\frac{1}{2} I_d + K_k\right) \phi \ \psi  % bottom right term 
&=& \int_{\Gamma} \gamma_D^{\mathrm{ext}}(u_0) \psi
\end{aligned}
\end{cases}
\end{equation}

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

\begin{equation}
\begin{cases}
\begin{aligned}
\displaystyle a_{\mathrm{int}}(u,v)  % Volume term
+ \int_{\Gamma} \left[W_k - i Z \left(\frac{1}{2} I_d+ K_k^{\prime}\right) \right] \phi\ \gamma_D^{\mathrm{int}}(v)  % Top right block
&=& \int_{\Omega^{\mathrm{int}}} f\ v - \int_{\Gamma} \gamma_N^{\mathrm{ext}}(u_0) \gamma_D^{\mathrm{int}}(v) \\
\displaystyle\int_{\Gamma} \gamma_D^{\mathrm{int}}(u) \psi  % bottom left term
-  \int_\Gamma \left[\left(\frac{1}{2}I_d + K_k\right) - i Z V_k\right] \phi \ \psi  % bottom right term
&=& \int_{\Gamma} \gamma_D^{\mathrm{ext}}(u_0) \psi
\end{aligned}
\end{cases}
\end{equation}

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

\begin{equation}
\begin{cases}
\begin{aligned}
\displaystyle a_{\mathrm{int}}(u,v)  % Volume term
+ \int_{\Gamma} \varphi \ \gamma_D^{\mathrm{int}}(v)  % Top right block
&=& \int_{\Omega^{\mathrm{int}}} f\ v - \int_{\Gamma} \gamma_N^{\mathrm{ext}}(u_0) \gamma_D^{\mathrm{int}}(v) \\
\displaystyle \int_{\Gamma} \left(\frac{1}{2}I_d - K_k\right)\gamma_D^{\mathrm{int}}(u) \psi  % bottom left term
-  \int_\Gamma V_k (\varphi) \ \psi  % bottom right term
&=& \int_{\Gamma} \left(\frac{1}{2}I_d - K_k\right) \gamma_D^{\mathrm{ext}}(u_0) \psi
\end{aligned}
\end{cases}
\end{equation}

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

\begin{equation}
\begin{cases}
\begin{aligned}
\displaystyle a_{\mathrm{int}}(u,v)  % Volume term
+ \int_{\Gamma} \varphi \ \gamma_D^{\mathrm{int}}(v)  % Top right block
&=& \int_{\Omega^{\mathrm{int}}} f\ v - \int_{\Gamma} \gamma_N^{\mathrm{ext}}(u_0) \gamma_D^{\mathrm{int}}(v) \\
\displaystyle \int_{\Gamma} W_k \gamma_D^{\mathrm{int}}(u) \psi  % bottom left term
-  \int_\Gamma \left(\frac{1}{2}I_d - K_k^{\prime}\right) (\varphi) \ \psi  % bottom right term
&=& \int_{\Gamma} W_k \gamma_D^{\mathrm{ext}}(u_0) \psi
\end{aligned}
\end{cases}
\end{equation}

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

\begin{equation}
\begin{cases}
\begin{aligned}
\displaystyle a_{\mathrm{int}}(u,v)  % Volume term         
+ \int_{\Gamma} \varphi \ \gamma_D^{\mathrm{int}}(v)  % Top right block
&=& \int_{\Omega^{\mathrm{int}}} f\ v - \int_{\Gamma} \gamma_N^{\mathrm{ext}}(u_0) \gamma_D^{\mathrm{int}}(v) \\
\displaystyle \int_{\Gamma} \left[W_k - i  \left(\frac{1}{2} I_d - K_k\right) Z \right]  \ \gamma_D^{\mathrm{int}}(u) \psi  % bottom left term
-  \int_\Gamma \left[\left(\frac{1}{2}I_d - K_k^{ \prime}\right) - i \ V_k\ Z\right] (\varphi) \ \psi  % bottom right term                                                                  
&=& \int_{\Gamma} \left[W_k - i  \left(\frac{1}{2} I_d - K_k\right) Z \right] \gamma_D^{\mathrm{ext}}(u_0) \psi
\end{aligned}
\end{cases}
\end{equation}

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


# 3d-Potential Flow

$\delta$

Reads a 2d mesh, construct a 3d mesh and solves the potential flow equation

\begin{equation*}
\begin{cases}
\begin{aligned}
	-\Delta u  &= 0 \hbox{in $\Omega_3$}\\
	\frac{\partial u}{\partial n}|_{\partial\Omega_3} &= q
\end{aligned}
\end{cases}
\end{equation*}

The solution exists only if $\int_{\partial\Omega_3}q=0$ and uniqueness holds up to a constant.
### Variational form
Given a small positive parameter $\epsilon$, find $u\in H^1(\Omega_3)$ such that
$$
\displaystyle{
	\mu\int_{\Omega_3}{\nabla u\cdot v +\epsilon u v}  = \int_{\partial\Omega_3}{q v},
  \quad \forall v\in H^1(\Omega_3)
}
$$
## Step 1: Build a 3d mesh from a 2d mesh and and a depth function

To build the 3D mesh we need to load  module msh3 and to display the results we will call Medit so we need to load the corresponding module.
~~~freefem
load "msh3"
load "medit"
~~~~
Verbosity adjusts the messages sent by FreeFem++.
~~~freefem
verbosity=3;
~~~
Let us read the 2d mesh contained in file ("lac-leman-v4.msh". The format is
- On the first line : the number of vertices, the number of triangles and the number of boundary edges
- Then the coordinates of the vertices and their region number
- Then the identification number (i.e; the position in the above list) of the 3 vertices of each triangle
- Finally the connectivity of the boundary edges: the 2 vertices and the label number of the edge
~~~freefem
mesh Th2("lac-leman-v4.msh");
fespace Vh2(Th2,P1);
Vh2 d;
~~~
The 2 lines above have defined a Finite Element Space on the 2d mesh with $P^1$ triangular elements.

To construct a 3d mesh we need the bathymetry of the lake, i.e. a function for the lake's depth at each point of the surface, $d(x,y),~(x,y)\in \Omega$.
This is done in a very approximate way as the result of a Laplace equation
$$
\left\{
\begin{align*}
	-\Delta d  &= 1 \hbox{ in }\Omega
  \\
	d|_{\partial\Omega} &= 0 
\end{align*}
\right.
$$
In variational form, this is
~~~freefem
{  Vh2 v; 
	macro Grad(u) [dx(u),dy(u)] //
	solve P(d,v)= int2d(Th2)(Grad(d)'*Grad(v))+int2d(Th2)(v)
	+on(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,d=-1);
	d = d*5/abs(d[].min);
	plot(d,wait=1);
}
~~~~
The macro is here to make the formula easier to read.
There are many labels for $\partial\Omega$ each corresponding to a portion of the lake boundary betwen two rivers.
The code has been surrounded by braces so as to make all variables local. This is by no means compulsory but this way one can reuse the names of the variables.

#### Build the 3d mesh using "buildlayers"
The mesh will have nn=10 layers:
~~~freefem
Vh2 ux,uz,p2;
int[int] rup=[0,200],  rdown=[0,100],rmid(17*2);
for(int i=0;i<rmid.n;++i)
  rmid[i]=1+i/2;
cout << rmid << endl;
real maxdeep = d[].min;
int nn=10;
mesh3 Th=buildlayers(Th2,nn,
  coef= d/maxdeep,
  zbound=[d,0],
  labelmid=rmid, 
  reffaceup = rup,
  reffacelow = rdown);
~~~~
The mesh is displyed with "medit" (use Command+Q to exit medit)
~~~freefem
if(!NoUseOfWait)medit("Leman",Th);
~~~
### Step2: Solve the problem
As usual we need to work witha finite element space on the 3d mesh and define some functions
~~~freefem
fespace Vh(Th,P13d);
Vh p,q;
macro Grad(u) [dx(u),dy(u),dz(u)] //

~~~~
In the following $\texttt{din}$ and $\texttt{dout}$ are the debits from borders 1 and borders 2, in and out of the lake.  Some rescaling is needed to satisfy the necessary condition for existence of a solution.
~~~freefem
real ain=int2d(Th,1)(1.);
real aout=int2d(Th,2)(1.);
cout << " area " << ain << " " << aout << endl;
real din=1./ain;
real dout=-1./aout;
~~~~
At last the problem is solved and the solution is displayed statically and inteeractively:
~~~freefem
solve P(p,q)= int3d(Th)(Grad(p)'*Grad(q)+1e-5*p*q)-int2d(Th,1)(q*din)+int2d(Th,2)(q*dout);

plot(p,wait=1,nbiso=30,value=1);
if(!NoUseOfWait) medit("potentiel",Th,p,wait=1);
~~~
The depth function d:
![](bathymetry.png)
The flow function p:
![](flow.png)
