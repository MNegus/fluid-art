/**
# Many particles

![All particles](tpp/mov.mp4)

~~~gnuplot A single particle path
set size ratio -1
set key off
set xlabel 'x'
set ylabel 'y'
set key right outside
plot 'out' u 3:4:2 w l palette lw 2 t 'pid()'
~~~
 */
#include "navier-stokes/centered.h"
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"

int maxlevel = 10;
u.t[top] = dirichlet (0.);
#define RAD (sqrt(sq(x) + sq(y)))
#define ST (x/RAD)

Particles many;
const face vector nu[] = {0.005, 0.005};

int main() {
  mu = nu;
  L0 = 20.;
  X0 = Y0 = -L0/2.;
  N = 512;
  run();
}

event init (t = 0) {
  double k = 3.83170597;
  scalar psi[];
  foreach()
    psi[] = ((RAD > 1)*((1/RAD))*ST +
	     (RAD < 1)*((-2*j1(k*RAD)*ST/(k*j0(k))) + (RAD*ST)));
  boundary({psi});
  foreach() {
    u.x[] = -(psi[0, 1] - psi[0, -1])/(2*Delta);
    u.y[] = (psi[1] - psi[-1])/(2*Delta);
  }
  many = init_tp_cells();
}

event adapt (i++)
  adapt_wavelet ((scalar*){u}, (double[]){0.05, 0.05}, maxlevel, 5);
 
event tracking (t += 0.1) {
  long unsigned int ptaga =  4*49164;
  int pnr = 0;
  double xp = HUGE, yp = HUGE;
  int lpid = npe();
  foreach_particle_in(many reduction(+:pnr) reduction(min:yp)
		      reduction(min:xp) reduction(min:lpid)) {
    pnr++;
    if (p().tag == ptaga) {
      lpid = pid();
      xp = x;
      yp = y;
    }
  }
  if (pid() == 0)
    printf ("%g %d %g %g %d %ld\n",t, lpid, xp, yp, pnr, grid->tn);
}

event render (t += 0.5; t < 25) {
  scatter (many, s = 5, pc = {sin(pid()), cos(pid()), sin(pid()*2.4)});
  box();
  save ("mov.mp4");
}
