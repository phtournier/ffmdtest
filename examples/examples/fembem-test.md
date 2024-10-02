---
name: fem-bem
category: bem
layout: example
---

~~~freefem
load "bem"
load "msh3"
load "PETSc-complex"
include "getARGV.idp"

real k0 = getARGV("-k0", 16.);
int NpWL = getARGV("-npwl", 8);
int prec = getARGV("-prec", 0);

string obstacleType = getARGV("-obstacle", "rectangular_cavity");
string inclination = getARGV("-inclination", "inclined");

real alpha; // Incident angle
if (inclination == "inclined"){
    alpha = 0;
}
else if (inclination == "behind"){
    alpha = -pi;
}

real lambda0 = 2*pi/k0; // wavelength
int NpLength = NpWL * 1./lambda0;
if (mpirank == 0) cout << "wavelength = " << lambda0 << endl;

// ==================================================
// ===============| Label definition |===============
// ==================================================
int labFemBem = 1;
int labObstacle = 2;
int labOther = 3;

// ==================================================
// =================| Volume domain |================
// ==================================================
int np = NpLength*8;
int Rfic = 2;
int L = 4;
real eps = 0.01;

mesh Th, ThOut;
if (obstacleType == "elliptic_cavity"){
    // =======| Boundary of obstacle |=======
    real phi0 = 7*pi/10;
    real lphi0 = 3.616;
    real lphi1 = 4.272;
    real phi1 = acos(1./1.3*cos(phi0));
    border bphi0(t=phi0,-phi0) {x = cos(t); y = 0.5*sin(t); label=labObstacle;};
    border bd0(t=0.5*sin(-phi0),0.6*sin(-phi1)) {x = cos(-phi0);y=t;label=labObstacle;};
    border bphi1(t=-phi1,phi1) {x = 1.3*cos(t); y = 0.6*sin(t); label=labObstacle;};
    border bd1(t=0.6*sin(phi1),0.5*sin(phi0)) {x = cos(phi0);y=t;label=labObstacle;};
    real lseg = 0.6*sin(phi1)-0.5*sin(phi0);

    border bInterface(t = 0, 2*pi){x=Rfic*cos(t); y=Rfic*sin(t); label=labFemBem;} // Interface for FEM-BEM coupling
    border bInterfaceOut(t = 0, 2*pi){x=(1+eps)*Rfic*cos(t); y=(1+eps)*Rfic*sin(t); label=labOther;}
    // =======| Boundary of viewing area |=======
    border b1(t=-L, L){x=t; y=-L; label=labOther;}
    border b2(t=-L, L){x=L; y=t; label=labOther;}
    border b3(t=-L, L){x=-t; y=L; label=labOther;}
    border b4(t=-L, L){x=-L; y=-t; label=labOther;}

    Th = buildmesh(bphi0(-NpLength*lphi0)
                        + bd0(-NpLength*lseg)
                        + bphi1(-NpLength*lphi1)
                        + bd1(-NpLength*lseg)
                        + bInterface(NpLength * 2*pi*Rfic)); // elliptic cavity
    ThOut = buildmesh(b1(np)
                        + b2(np)
                        + b3(np)
                        + b4(np)
                        + bInterfaceOut(-NpLength * 2*pi*Rfic*(1+eps)));
}
else if (obstacleType == "rectangular_cavity"){
    real Lext = 1.5;
    real Lint = 1.3;
    real H = 0.6;
    real h = 0.1;
    // =======| Boundary of obstacle |=======
    border Obstacle1(t=-Lext/2, Lext/2){x=t; y=H/2; label=labObstacle;}
    border Obstacle2(t=0, H){x=Lext/2; y=H/2-t; label=labObstacle;}
    border Obstacle3(t=0, Lext){x=Lext/2-t; y=-H/2; label=labObstacle;}
    border Obstacle4(t=0, h){x=-Lext/2; y=-H/2+t; label=labObstacle;}
    border Obstacle5(t=0, Lint){x=-Lext/2+t; y=-H/2+h; label=labObstacle;}
    border Obstacle6(t=0, H-2*h){x=-Lext/2+Lint; y=-H/2+h+t; label=labObstacle;}
    border Obstacle7(t=0, Lint){x=-Lext/2+Lint-t; y=H/2-h; label=labObstacle;}
    border Obstacle8(t=0, h){x=-Lext/2; y=H/2-h+t; label=labObstacle;}

    border bInterface(t = 0, 2*pi){x=Rfic*cos(t); y=Rfic*sin(t); label=labFemBem;} // Interface for FEM-BEM coupling
    border bInterfaceOut(t = 0, 2*pi){x=(1+eps)*Rfic*cos(t); y=(1+eps)*Rfic*sin(t); label=labOther;}
    // =======| Boundary of viewing area |=======
    border b1(t=-L, L){x=t; y=-L; label=labOther;}
    border b2(t=-L, L){x=L; y=t; label=labOther;}
    border b3(t=-L, L){x=-t; y=L; label=labOther;}
    border b4(t=-L, L){x=-L; y=-t; label=labOther;}
    border circleOut(t = 0, 2*pi){x=(1+eps)*Rfic*cos(t); y=(1+eps)*Rfic*sin(t); label=labOther;}

    Th = buildmesh(Obstacle1(NpLength * Lext)
                    + Obstacle2(NpLength * H)
                    + Obstacle3(NpLength * Lext)
                    + Obstacle4(NpLength * h)
                    + Obstacle5(NpLength * Lint)
                    + Obstacle6(NpLength * (H-2*h))
                    + Obstacle7(NpLength * Lint)
                    + Obstacle8(NpLength * h)
                    + bInterface(NpLength * 2*pi*Rfic)); // rectangular cavity (closed right)
   
    ThOut = buildmesh(b1(np)+b2(np)+b3(np)+b4(np) + bInterfaceOut(-NpLength * 2*pi*Rfic*(1+eps)));    
}
else if (obstacleType == "curved_square"){
    real sqSize = 1;
    // =======| Boundary of obstacle |=======
    border LeftEdge(t=0, sqSize/2){x=-sqSize/2; y=-sqSize/4 + t; label=labObstacle;}
    border TopLeftCorner(t = pi/2, pi){x=-sqSize/4 + sqSize/4 * cos(t); y=sqSize/4 + sqSize/4 * sin(t); label=labObstacle;}
    border TopEdge(t=0, sqSize/2){x=-sqSize/4 + t; y=sqSize/2; label=labObstacle;}
    border TopRightCorner(t = 0, pi/2){x=sqSize/4 + sqSize/4 * cos(t); y=sqSize/4 + sqSize/4 * sin(t); label=labObstacle;}
    border RightEdge(t=0, sqSize/2){x=sqSize/2; y=sqSize/4 - t; label=labObstacle;}
    border BotRightCorner(t = 3*pi/2, 2*pi){x=sqSize/4 + sqSize/4 * cos(t); y=-sqSize/4 + sqSize/4 * sin(t); label=labObstacle;}
    border BotEdge(t=0, sqSize/2){x=sqSize/4 - t; y=-sqSize/2; label=labObstacle;}
    border BotLeftCorner(t = pi, 3*pi/2){x=-sqSize/4 + sqSize/4 * cos(t); y=-sqSize/4 + sqSize/4 * sin(t); label=labObstacle;}   

    border bInterface(t = 0, 2*pi){x=Rfic*cos(t); y=Rfic*sin(t); label=labFemBem;} // Interface for FEM-BEM coupling
    border bInterfaceOut(t = 0, 2*pi){x=(1+eps)*Rfic*cos(t); y=(1+eps)*Rfic*sin(t); label=labOther;}
    // =======| Boundary of viewing area |=======
    border b1(t=-L, L){x=t; y=-L; label=labOther;}
    border b2(t=-L, L){x=L; y=t; label=labOther;}
    border b3(t=-L, L){x=-t; y=L; label=labOther;}
    border b4(t=-L, L){x=-L; y=-t; label=labOther;}
    border circleOut(t = 0, 2*pi){x=(1+eps)*Rfic*cos(t); y=(1+eps)*Rfic*sin(t); label=labOther;}

    Th = buildmesh(LeftEdge(NpLength * sqSize/2)
                    + TopLeftCorner(-NpLength * pi/2 * sqSize/4)
                    + TopEdge(NpLength * sqSize/2)
                    + TopRightCorner(-NpLength * pi/2 * sqSize/4)
                    + RightEdge(NpLength * sqSize/2)
                    + BotRightCorner(-NpLength * pi/2 * sqSize/4)
                    + BotEdge(NpLength * sqSize/2)
                    + BotLeftCorner(-NpLength * pi/2 * sqSize/4)
                    + bInterface(NpLength * 2*pi*Rfic));
   
    ThOut = buildmesh(b1(np)+b2(np)+b3(np)+b4(np) + bInterfaceOut(-NpLength * 2*pi*Rfic*(1+eps)));
}

// =======| Boundary extraction |=======
int[int] labInterface = [labFemBem];
meshL ThL = extract(Th, label=labInterface); // BEM mesh
labInterface = [labObstacle];
meshL ThLObstacle = extract(Th, label=labInterface); 
ThL = OrientNormal(ThL,unbounded=1); // Normal orientation

// =======| Discretization in FE spaces |=======
load "Element_P3"
fespace Uh(Th, P2);
fespace UhBem(ThL,P2);
fespace UhOut(ThOut,P1);
fespace Vh=Uh*UhBem; // composite space
 
Uh hFEM = hTriangle;
if (mpirank == 0){
  cout << "size of FEM mesh = " << hFEM[].max << endl;
  cout << "FEM ndof = " << Uh.ndof << endl;
  cout << "BEM ndof = " << UhBem.ndof << endl;
}
// =======| Definitions of BIO kernels 
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
$\varphi := \gamma_N^{\mathrm{ext}}(u)$

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
$\varphi := \gamma_N^{\mathrm{ext}}(u)$

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
$\varphi := \gamma_N^{\mathrm{ext}}(u)$

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
