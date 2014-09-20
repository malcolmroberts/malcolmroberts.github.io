// Animation of solutions of the 1D heat equation.
// use asymptote (asymptote.sf.net) to compile using "asy heat1.asy" in Linux
// or some other such things in Windows/Mac
// Change the initial conditions by hand in the function.
// Written by Malcolm Roberts, 2010-03-31

import graph;
import animation;

size(0,300);

defaultpen(3);
dotfactor=4;

real t1=0; 
real t2=1;
int Nx=256;
real beta=0.3;
real u0=0, uL=0;

xaxis(Label("$x$",align=3SW),0);
yaxis("$u(x,t)$",0,1.2);

animation a;
typedef real realfcn(real);

string BC=getstring("Boundary conditions (dirichelt or neumann)");
bool odd;
if(BC=="dirichlet") {
  odd=true;
  u0=getreal("u0");
  uL=getreal("uL");
}
if(BC=="neumann")
  odd=false;

// the initial conditions
int icn=getint("Initial conditions (1,2,3)");
real ic(real x) {
  if(icn==1)
    return x;
  if(icn==2)
    return sin(x)+2*sin(3*x);
  if(icn==3) {
  if (x>pi/3 && x< 2pi/3)
    return 1;
  return 0;
  }
  return 0;
}
// the steady-state solution
real steady(real x) {
  if(odd)
    return u0 + (uL-u0)*x/pi;
  return 0; // not implemented for Neumann BCs
}
// symetrize (depending on BC)
real f(real x) {
  if(odd) {
    if(x<0) return -ic(abs(x))+steady(abs(x));
    return ic(x)-steady(x);
  }
  return ic(abs(x)-steady(abs(x)));
}
// discretize
real[] xi, fi;
real dx=2*pi/Nx;
for (int i=0; i <= Nx; ++i) {
  xi[i]=i*dx-pi;
  fi[i]=f(xi[i]);
}
pair[] A=fft(fi);


int nmodes=20;
nmodes=getint("number of modes");
realfcn F(real t) {
  return new real(real x) {
    pair val=0;
    if(!odd)
      val=A[0].x/2;
    for (int i=1; i <= nmodes; ++i) {
      // (-1)^i from the shift theorem
      if(odd)
	val+=exp(-beta*i*i*t)*A[i].y*sin(i*x)*(-1)^i;
      else
	val+=exp(-beta*i*i*t)*A[i].x*cos(i*x)*(-1)^i;
    }
    return steady(x)+2*val.x/Nx;
  };
};

// animate
int n=4;
n=getint("number of frames");
real dt=(t2-t1)/n;
for(int i=0; i <= n; ++i) {
  save();
  
  real t=t1+dt*i;

  draw(graph(F(t),0,pi));

  draw(graph(steady,0,pi),blue+dashed);

  draw(graph(ic,0,pi),red+dashed);
  
  a.add(); // Add currentpicture to animation.
  restore();
}

// Produce the final merged gif.
a.movie(loops=10,delay=250);
