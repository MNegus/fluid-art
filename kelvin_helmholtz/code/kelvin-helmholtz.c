

#include "navier-stokes/centered.h"
#include "two-phase.h" // Implements two-phase flow
#include "view.h" // Creating movies using bview
#include <omp.h> // For openMP parallel


/* Physical parameters */
double REYNOLDS = 5000.;
double RHO_R = 0.75;
double MU_R = 1;

int MINLEVEL = 6;
int MAXLEVEL = 11;
int BOX_WIDTH = 4;

int gfs_output_no = 0;
scalar omega[];

double END_TIME = 5;

/* Boundary conditions */

u.n[top] = neumann(0.); // Allows outflow through boundary
u.n[bottom] = neumann(0.); // Allows outflow through boundary



int main() {
    /* Create the computational domain */
    init_grid(1 << MINLEVEL); // Create grid according to the minimum level
    size(BOX_WIDTH); // Size of the domain
    origin(-BOX_WIDTH/2, -BOX_WIDTH/2);

    periodic(right);

    /* Set physical constants */
    rho1 = 1.; // Density of the top phase
    rho2 = RHO_R; // Density of the bottom phase
    mu1 = 1. / REYNOLDS; // Viscosity of top phase
    mu2 = mu1 * MU_R; // Viscosity of bottom phase

    /* Run simulation */
    run();
}

event init(t = 0) {

    /* Define volume fraction along middle of box */
    fraction(f, y - 0.001 * sin(100 * 2 * pi * x / BOX_WIDTH));
    // fraction(f, y);

    refine((y < 0.02) && (y > -0.02) && level < MAXLEVEL);

    /* Initialise the velocities */
    foreach() {
        u.x[] = f[];
    }
}

event adapt (i++) {
  adapt_wavelet ({u,f}, (double[]){1e-3,1e-3,1e-4}, maxlevel = MAXLEVEL);
}

event logfile (i++) {
  
  fprintf (stderr, "%d %g\n", i, t);
}


event gfs_output (t += 0.1) {
/* Saves a gfs file */
    vorticity(u, omega);

    char gfs_filename[80];
    sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
    output_gfs(file = gfs_filename);

    gfs_output_no++;
}

event movies (t += 1e-1) {

        char time_str[80];
        sprintf(time_str, "t = %g\n", t);

        /* Zoomed out view */
        // Set up bview box
        view (width = 1024, height = 1024, fov = 24.0);

        /* Movie of the volume fraction of the droplet */
        clear();
        draw_vof("f", lw = 2);
        squares("f", linear = true, spread = -1, map = cool_warm); 
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);

        save ("tracer.mp4");

}


event end(t = END_TIME) {
    fprintf(stderr, "Done!\n");
}