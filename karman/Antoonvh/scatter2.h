/**
# A Scattering bview function

the `scatter` function can be used to scatter
[`Particles`](particle.h) in bview. The formulation is inspired by
other bview functions in [draw.h](src/draw.h).
 */
#if dimension == 3 
void glPointParameterfv(GLenum pname, const GLfloat * params);
#endif

struct _scatter {
  Particles p;    // particles
  float s, pc[3], coefs[3]; // point size, colour and distance attenuation coefs.
};

trace
void scatter (struct _scatter p){
  bview * view = draw();
#if dimension == 2
  glTranslatef (0., 0., view->lc*view->fov/24.); //from draw_lines()
#else // Dimension == 3
  if (!p.coefs[0]){ // A guess:
    p.coefs[0] = 0.01;
    p.coefs[1] = 0.2;
    p.coefs[2] = 0.5;
  }
  glPointParameterfv(GL_POINT_DISTANCE_ATTENUATION, p.coefs);
#endif
  glEnable (GL_BLEND);
  glEnable (GL_POINT_SMOOTH);
  glColor3f(p.pc[0], p.pc[1], p.pc[2]);
  if (!p.s)
    p.s = 20;
  glPointSize(p.s);

  glBegin (GL_POINTS);
  foreach_particle_in(p.p) {
#if dimension == 2    
    glvertex2d (view, x, y);
#else // dimension == 3
    glvertex3d (view, x, y, z);
#endif
  }
  glEnd();

  view->ni++; 
}
/**
## Curve
From [A. Castillo's](../acastillo/filaments/filaments.h) sandbox.
 */
struct _curve {
  Particles p;    // particles
  float lw, lc[3]; // line width and color.
  bool noloop;
};

void draw_curve(struct _curve p){
  if (!p.lc[0]) {
    p.lc[0] = 0, p.lc[1] = 0, p.lc[2] = 0;
  } 
  glColor3f(p.lc[0], p.lc[1], p.lc[2]);
  if (!p.lw)
    p.lw = 5;
  glLineWidth(p.lw);
  glEnable (GL_BLEND);
  glEnable (GL_LINE_SMOOTH);
  particles pp = pl[(int)p.p];
  foreach_particle_in(p.p) {
    glBegin(GL_LINES);
    int in = j < pn[p.p] - 1 ? j + 1 : p.noloop ? j : 0;
    glVertex3f(x, y, z);
    glVertex3f(pp[in].x, pp[in].y, pp[in].z);
    glEnd();
  }
}
/**
## Usage

* [Visualizing flow tracer particles](tracer-particles.h)
* [A confetti cannon](confetti.c)

*/
