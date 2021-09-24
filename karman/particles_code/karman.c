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
// #include "tracer.h"
#include "view.h" // Creating movies using bview
#include <omp.h> // For openMP parallel
#include "../Antoonvh/tracer-particles.h"
#include "../Antoonvh/scatter2.h"


// scalar upper_f[];
// scalar lower_f[];
// scalar * tracers = {upper_f,lower_f};
Particles trackers;

double Reynolds = 160.;
int MAXLEVEL = 8;
face vector muv[];

int gfs_output_no = 0;
int field_output_no = 0;
scalar omega[];

double END_TIME = 20.0;
double DIAMETER = 0.125;
double INJECT_POSITION;
int gfs_output_no;

/**
The domain is eight units long, centered vertically. */

int main() {
  L0 = 8.;
  origin (-0.5, -L0/2.);
  N = 512;
  mu = muv;

  INJECT_POSITION = 0.001;

  // upper_f[left] = dirichlet(y > INJECT_POSITION);
  // lower_f[left] = dirichlet(y < - INJECT_POSITION);

  gfs_output_no = 0;
  
  run();
   
}

/**
We set a constant viscosity corresponding to a Reynolds number of 160,
based on the cylinder diameter and the inflow velocity (1). */

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*DIAMETER/Reynolds;
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
    phi[] = intersection (phi[], sq(x) + sq(y) - sq(DIAMETER/2.));
  }
  boundary ({phi});
  fractions (phi, cs, fs);

  /* Initialise tracer particles */
  trackers = new_tracer_particles(0);
  // foreach_particle_in(Pupper) {
  //   p().x = -0.5;
  //   p().y = INJECT_POSITION;
  // }

  // Pupper = new_tracer_particles(1);
  // foreach_particle_in(Pupper) {
  //   p().x = -0.5;
  //   p().y = -INJECT_POSITION;
  // }



}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++) {
  
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}

/* Add particles */
event add_particles(t += 1e-4) {

  // Add upper particle
  particle upper_p = {.x = -0.5, .y = INJECT_POSITION};
  set_a_particle_attributes(&upper_p);
  add_particle(upper_p, trackers);

  // Add lower particle
  particle lower_p = {.x = -0.5, .y = -INJECT_POSITION};
  set_a_particle_attributes(&lower_p);
  add_particle(lower_p, trackers);
}

/* Remove particles */
event remover (i++) {
  remove_particles(trackers, x > L0 - 0.5);
}

/**
We produce animations of the vorticity and tracer fields... */

event movies (t += 1e-2) {
    char time_str[80];
    sprintf(time_str, "t = %g\n", t);

    /* Zoomed out view */
    // Set up bview box
    view (width = 2048, height = 682, fov = 6, tx = -0.42);

    /* Movie of the volume fraction of the droplet */
    clear();
    draw_vof("cs", filled=1, lc = {1, 1, 1});
    scatter (trackers, s = 10, pc = {1, 1, 1});

    char filename[80];
    sprintf(filename, "tracer_%g.mp4", -INJECT_POSITION);
    save (filename);

}

event gfs_output (t += 5e-1) {
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
  // adapt_wavelet ({cs,u,upper_f,lower_f}, (double[]){1e-3,1e-3,1e-3,1e-3,1e-3}, maxlevel = MAXLEVEL);
  adapt_wavelet ({cs,u}, (double[]){1e-3,1e-3,1e-3}, maxlevel = MAXLEVEL);
}

event end(t = END_TIME) {
  

  fprintf(stderr, "Reached end event\n");

}

/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/
