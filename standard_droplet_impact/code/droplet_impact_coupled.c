/* droplet_impact_coupled.c
    A 2D droplet falling towards an elastic membrane lying 
    along the boundary at y = 0. The solution for the membrane is given as a 
    function of pressure by the routines defined in wave-equation.h, and its 
    velocity is fed back into Basilisk by altering the boundary condition. 

    Runs until the turnover point approximately reaches the initial droplet 
    radius. 

    Author: Michael Negus
*/

#define FILTERED
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

#include "parameters.h" // Includes all defined parameters
#include "axi.h" // Axisymmetric coordinates
#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h" // Implements two-phase flow
#include "view.h" // Creating movies using bview
#include "tension.h" // Surface tension of droplet
#include "tag.h" // For removing small droplets
#include "heights.h"
#include "contact.h" // For imposing contact angle on the surface
#include <omp.h> // For openMP parallel
#include <stdlib.h>

/* Computational constants derived from parameters */
double MIN_CELL_SIZE; // Size of the smallest cell
double DROP_REFINED_WIDTH; // Width of the refined area around the droplet
double MEMBRANE_REFINE_NO; // Number of grid cells above the membrane to refine
double MEMBRANE_REFINED_HEIGHT; // Width of the refined area above the membrane 
double DROP_CENTRE; // Initial centre position of the droplet
double IMPACT_TIME; // Theoretical time of impact

/* Global variables */
double start_wall_time; // Time the simulation was started
double end_wall_time; // Time the simulation finished
int gfs_output_no = 0; // Records how many GFS files have been outputted
int log_output_no = 0; // Records how many plate data files there have been
int interface_output_no = 0; // Records how many interface files there have been
int membrane_output_no = 0; // Records how many membrane outputs there have been
int start_membrane = 0; // Boolean to indicate if membrane motion has started
double drop_thresh = 1e-4; // Remove droplets threshold
double pinch_off_time = 0.; // Time of pinch-off

/* Contact angle variables */ 
vector h[];  // Height function
double theta0 = 90;  // Contact angle in degrees

/* Boundary conditions */

// Conditions on surface
uf.n[left] = dirichlet(0.); // No flow through surface
uf.t[left] = dirichlet(0.); // No slip at surface
h.t[left] = contact_angle (theta0*pi/180.); // RC contact angle

// Conditions for entry from above
u.n[right] = neumann(0.); // Free flow condition
p[right] = dirichlet(0.); // 0 pressure far from surface

// Conditions far from the droplet in the radial direction
u.n[top] = neumann(0.); // Allows outflow through boundary
u.t[top] = dirichlet(0.); // Stationary vertical flow
p[top] = dirichlet(0.); // 0 pressure far from surface

int main() {
/* Main function for running the simulation */

    /* Create the computational domain */
    init_grid(1 << MINLEVEL); // Create grid according to the minimum level
    size(BOX_WIDTH); // Size of the domain

    /* Set physical constants */
    rho1 = 1.; // Density of water phase
    rho2 = RHO_R; // Density of air phase
    mu1 = 1. / REYNOLDS; // Viscosity of water phase
    mu2 = mu1 * MU_R; // Viscosity of air phase
    f.height = h; // For contact angle calculation
    f.sigma = 1. / WEBER; // Surface tension at interface

    /* Derived constants */
    MIN_CELL_SIZE = BOX_WIDTH / pow(2, MAXLEVEL); // Size of the smallest cell
    DROP_REFINED_WIDTH = 0.05; // Refined region around droplet
    DROP_CENTRE = INITIAL_DROP_HEIGHT + DROP_RADIUS; // Initial centre of drop
    IMPACT_TIME = INITIAL_DROP_HEIGHT / (-DROP_VEL); // Theoretical impact time
    MEMBRANE_REFINE_NO = 8; // Number of cells above membrane to refine by
    MEMBRANE_REFINED_HEIGHT = MEMBRANE_REFINE_NO * MIN_CELL_SIZE; 
    

    /* Creates log file */
    FILE *logfile = fopen("log", "w");
    fclose(logfile);

    /* Poisson solver constants */
    DT = 1.0e-4; // Minimum timestep
    NITERMIN = 1; // Min number of iterations (default 1)
    NITERMAX = 300; // Max number of iterations (default 100)
    TOLERANCE = 1e-6; // Possion solver tolerance (default 1e-3)

    /* Runs the simulation */
    run(); 
}


event init(t = 0) {
/* Initialises the flow as a spherical droplet falling downwards */

    // Records the wall time
    start_wall_time = omp_get_wtime();

    /* Refines around the droplet */
    refine((sq(x - DROP_CENTRE) + sq(y) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
        && (sq(x - DROP_CENTRE) + sq(y)  > sq(DROP_RADIUS - DROP_REFINED_WIDTH)) \
        && (level < MAXLEVEL));
    
    /* Initialises the droplet volume fraction */
    fraction(f, -sq(x - DROP_CENTRE) - sq(y) + sq(DROP_RADIUS));

    /* Initialise the droplet velocity downwards */
    foreach() {
        u.x[] = DROP_VEL * f[];
    }
    boundary ((scalar *){u});
}


event refinement (i++) {
/* Adaptive grid refinement */

    // Adapts with respect to velocities and volume fraction 
    adapt_wavelet ({u.x, u.y, f}, (double[]){1e-3, 1e-3, 1e-6},
        minlevel = MINLEVEL, maxlevel = MAXLEVEL);

    /* Attempts to refine above the membrane, doubling the refine height until
    successful */
    double refine_height = MEMBRANE_REFINED_HEIGHT;
    int adequate_refinement = 0;

    while (adequate_refinement == 0) {

        // Attempts to refine
        refine((y < MEMBRANE_RADIUS) && (x <= refine_height) \
            && level < MAXLEVEL);

        // Loops and check if refinement was successful
        adequate_refinement = 1;
        foreach_boundary(left) {
            if ((y < MEMBRANE_RADIUS) && (level < MAXLEVEL)) {
                adequate_refinement = 0;
                break;
            }
        }

        // If refinement was unsuccessful, then double the refined height
        if (adequate_refinement == 0) refine_height = 2 * refine_height;
    }
    
}


event gravity (i++) {
/* Adds acceleration due to gravity in the vertical direction */
    face vector av = a; // Acceleration at each face
    foreach_face(x) av.y[] -= 1./sq(FR); // Adds acceleration due to gravity
}


event small_droplet_removal (t += 1e-4) {
/* Removes any small droplets or bubbles that have formed, that are smaller than
    a specific size */
    // Removes droplets of diameter 5 cells or less
    int remove_droplet_radius = min(16, (int)(0.2 / MIN_CELL_SIZE));
    remove_droplets(f, remove_droplet_radius);

    // Also remove air bubbles
    remove_droplets(f, remove_droplet_radius, 1e-4, true);
    
}


event output_data (t += LOG_OUTPUT_TIMESTEP) {
/* Outputs data about the flow */
    /* Outputs data to log file */
    printf("t = %.5f, v = %.8f\n", t, 2 * pi * statsf(f).sum);
    
    FILE *logfile = fopen("log", "a");
    fprintf(logfile, "t = %.5f, v = %.8f\n", t, 2 * pi * statsf(f).sum);
    fclose(logfile);

    log_output_no++;
}


event output_interface (t += PLATE_OUTPUT_TIMESTEP) {
/* Outputs the interface locations of the droplet */
    // Creates text file to save output to
    char interface_filename[80];
    sprintf(interface_filename, "interface_%d.txt", interface_output_no);
    FILE *interface_file = fopen(interface_filename, "w");

    // Outputs the interface locations and closes the file
    output_facets(f, interface_file);
    fclose(interface_file);

    interface_output_no++;
}


event gfs_output (t += GFS_OUTPUT_TIMESTEP) {
/* Saves a gfs file */
    char gfs_filename[80];
    sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
    output_gfs(file = gfs_filename);

    gfs_output_no++;
}


event movies (t += 5e-3) {
/* Produces movies using bview */ 
    if (MOVIES) {
        // Creates a string with the time to put on the plots
        char time_str[80];
        sprintf(time_str, "t = %g\n", t);

        /* Zoomed out view */
        // Set up bview box
        view (width = 1024, height = 1024, fov = 20.0, ty = -0.5, tx = 0.5, quat = {0, 0, -0.707, 0.707});

        /* Movie of the volume fraction of the droplet */
        clear();
        draw_vof("f", lw = 2);
        squares("f", linear = true, spread = -1, map = cool_warm); // RC - minor changes here and beyond
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("tracer.mp4");

        /* Movie of the horiztonal velocity */
        clear();
        draw_vof("f", lw = 2);
        squares("u.x", spread = -1, linear = true, map = cool_warm);
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("horizontal_vel.mp4");


        /* Movie of the vertical velocity */
        clear();
        draw_vof("f", lw = 2);
        squares("u.y", min = -1.5, max = 1.5, linear = true, spread = -1, map = cool_warm);
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("vertical_vel.mp4");

        /* Movie of the pressure */
        clear();
        draw_vof("f", lw = 2);
        squares("p", spread = -1, linear = true, map = cool_warm);
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("pressure.mp4");

        /* Zoomed in view of pressure around entrapped bubble */
        // Set up bview box
        view (width = 1024, height = 1024, fov = 2.0, ty = -0.05, tx = -0.05);
        clear();
        draw_vof("f", lw = 2);
        squares("u.y", min = -1.5, max = 1.5, linear = true, spread = -1, map = cool_warm);
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("zoomed_vertical_vel.mp4");
    }
}


event end (t = MAX_TIME) {
/* Ends the simulation */ 

    end_wall_time = omp_get_wtime(); // Records the time of finish

    fprintf(stderr, "Finished after %g seconds\n", \
        end_wall_time - start_wall_time);
    FILE *logfile = fopen("log", "a");
    fprintf(logfile, "Finished after %g seconds\n", end_wall_time - start_wall_time);
    fclose(logfile);

}



// void output_arrays(double *w_arr, double *w_deriv_arr, double *p_arr) {
// /* output_membrane
// Outputs the x positions of the membrane into a text file
// */
//     char w_filename[40];
//     sprintf(w_filename, "w_%d.txt", membrane_output_no);
//     FILE *w_file = fopen(w_filename, "w");

//     char w_deriv_filename[40];
//     sprintf(w_deriv_filename, "w_deriv_%d.txt", membrane_output_no);
//     FILE *w_deriv_file = fopen(w_deriv_filename, "w");

//     char p_filename[40];
//     sprintf(p_filename, "p_%d.txt", membrane_output_no);
//     FILE *p_file = fopen(p_filename, "w");

//     // Outputs from x = 0 to L - dx
//     #pragma omp parallel for
//     for (int k = 0; k < M; k++) {
//         double x = k * DELTA_X;
//         // fprintf(w_file, "%.10f, %.10f\n", x, w_arr[k]);
//         // fprintf(w_deriv_file, "%.10f, %.10f\n", x, w_deriv_arr[k]);
//         // fprintf(p_file, "%.10f, %.10f\n", x, p_arr[k]);
//         fprintf(w_file, "%g, %g\n", x, w_arr[k]);
//         fprintf(w_deriv_file, "%g, %g\n", x, w_deriv_arr[k]);
//         fprintf(p_file, "%g, %g\n", x, p_arr[k]);
//     }

//     // Outputs x = L, where w and w_deriv = 0
//     double x = M * DELTA_X;
//     fprintf(w_file, "%.10f, %.10f\n", x, 0.0);
//     fprintf(p_file, "%.10f, %.10f\n", x, 0.0);
//     fprintf(w_deriv_file, "%.10f, %.10f", x, 0.0);

//     fclose(w_file);
//     fclose(p_file);
//     fclose(w_deriv_file);

//     membrane_output_no++;
// }


// void output_arrays_stationary(double *p_arr) {
// /* output_membrane_stationary
// Outputs the x positions of the pressure in a text file
// */
//     char p_filename[40];
//     sprintf(p_filename, "p_%d.txt", membrane_output_no);
//     FILE *p_file = fopen(p_filename, "w");

//     // Outputs from x = 0 to L - dx
//     #pragma omp parallel for
//     for (int k = 0; k < M; k++) {
//         double x = k * DELTA_X;
//         fprintf(p_file, "%g, %g\n", x, p_arr[k]);
//     }

//     // Outputs x = L, where w and w_deriv = 0
//     double x = M * DELTA_X;
//     fprintf(p_file, "%.10f, %.10f\n", x, 0.0);

//     fclose(p_file);

//     membrane_output_no++;
// }

