/**
# Bénard–von Kármán Vortex Street for flow around a cylinder at Re=160

An example of 2D viscous flow around a simple solid boundary. Fluid is
injected to the left of a channel bounded by solid walls with a slip
boundary condition. A passive tracer is injected in the bottom half of
the inlet.

![Animation of the vorticity field.](karman/vort.mp4)(loop)

![Animation of the tracer field.](karman/f.mp4)(loop)

We use the centered Navier-Stokes solver, with embedded boundaries and
advect the passive tracer *f*. */

#include "embed.h"
#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "tracer.h"
#include "view.h" // Creating movies using bview
#include <omp.h> // For openMP parallel

scalar f[];
scalar * tracers = {f};

double Reynolds = 160.;
int MAXLEVEL = 9;
face vector muv[];

int gfs_output_no = 0;
scalar omega[];

double END_TIME = 2.5;
double INJECT_POSITION;
double RADIUS = 0.25;
int gfs_output_no;

/**
The domain is eight units long, centered vertically. */

int main() {
  L0 = 2.;
  origin (-0.5, -L0/2.);
  N = 512;
  mu = muv;

  // INJECT_POSITION = -0.0;
  // f[left] = dirichlet(y < INJECT_POSITION);
  // run();

  // for (INJECT_POSITION = -0.002; INJECT_POSITION > -0.02; INJECT_POSITION -= 0.002) {
  //   gfs_output_no = 0;
  //   f[left]  = dirichlet(y < INJECT_POSITION);
  //   run();
  // }

  for (INJECT_POSITION = -0.50; INJECT_POSITION >= -1; INJECT_POSITION -= 0.02) {
    gfs_output_no = 0;
    f[left]  = dirichlet(y < INJECT_POSITION);
    run();
  }

  // for (INJECT_POSITION = 0.02; INJECT_POSITION <= 0.04; INJECT_POSITION += 0.02) {
  //   f[left] = dirichlet(y > INJECT_POSITION);
  //   run();
  // }
   
}

/**
We set a constant viscosity corresponding to a Reynolds number of 160,
based on the cylinder diameter (0.125) and the inflow velocity (1). */

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*RADIUS/Reynolds;
}

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The top and bottom walls are free-slip and the cylinder is no-slip. */

u.n[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);


event init (t = 0)
{

  /**
  The domain is the intersection of a channel of width unity and a
  circle of diameter 0.125. */
  fprintf(stderr, "Starting with INJECT_POSITION = %g\n", INJECT_POSITION);

  vertex scalar phi[];
  foreach_vertex() {
    phi[] = intersection (L0 - y, L0 + y);
    phi[] = intersection (phi[], sq(x) + sq(y) - sq(RADIUS/2.));
  }
  boundary ({phi});
  fractions (phi, cs, fs);

  // fraction(cs, -sq(x) - sq(y) + sq(0.125/2.));

  // foreach() {
  //   foreach_dimension()
  //     u.x[] = (1. - cs[])*u.x[];
  // }
  // boundary ((scalar *){u});

}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++) {
  
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}

/**
We produce animations of the vorticity and tracer fields... */

event movies (t += 1e-1)
{
        // scalar omega[], m[];
        // vorticity (u, omega);
        // foreach()
        //   m[] = cs[] - 0.5;
        // boundary ({m});
        // output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
        //       min = -10, max = 10, linear = true, mask = m);
        // output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
        //       linear = false, min = 0, max = 1, mask = m);

        char time_str[80];
        sprintf(time_str, "t = %g\n", t);

        /* Zoomed out view */
        // Set up bview box
        view (width = 1024, height = 1024, fov = 24.0, tx = -0.25);

        /* Movie of the volume fraction of the droplet */
        clear();
        // draw_vof("f", lw = 2);
        isoline("f", 0.5, 1, lw=2);
        squares("f", linear = true, spread = -1, map = cool_warm); 
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);

        char filename[80];
        sprintf(filename, "tracer_%g.mp4", -INJECT_POSITION);
        save (filename);

}

event gfs_output (t += 1e-1) {
/* Saves a gfs file */
    vorticity(u, omega);

    char gfs_filename[80];
    sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
    output_gfs(file = gfs_filename);

    gfs_output_no++;
}


/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

event adapt (i++) {
  adapt_wavelet ({cs,u,f}, (double[]){1e-3,1e-3,1e-3,1e-3}, maxlevel = MAXLEVEL);
}

event end(t = END_TIME) {
  char filename[80];
  sprintf(filename, "f_output_%g.txt", -INJECT_POSITION);
  FILE *endfile = fopen(filename, "w");
  output_field({f}, endfile);
  fclose(endfile);

  fprintf(stderr, "Reached end event\n");

}

/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/
