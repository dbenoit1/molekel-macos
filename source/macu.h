/*  MOLEKEL, Version 4.3, Date: 11.Nov.02
 *  Copyright (C) 2000-2002 Stefan Portmann (CSCS/ETHZ)
 *  (original IRIX GL implementation, concept and data structure
 *   by Peter F. Fluekiger, CSCS/UNI Geneva)
 *
 *  This software makes use of the 
 *  GLUT http://reality.sgi.com/mjk/glut3/glut3.html
 *  GLUI http://www.cs.unc.edu/~rademach/glui/
 *  libtiff http://www.libtiff.org/tiff-v3.5.5.tar.gz
 *  libjpeg ftp://ftp.uu.net/graphics/jpeg
 *  and in some versions of the 
 *  Mesa http://www.mesa3d.org/
 *  and the 
 *  libimage https://toolbox.sgi.com/toolbox/src/haeberli/libimage/index.html#dl
 *  libraries.
 *  An adapted version of the tr library by Brian Paul
 *  (http://www.mesa3d.org/brianp/TR.html)
 *  is part of the distribution.
 *
 *  The binary code is available free of charge but is not in the
 *  public domain. See license for details on conditions and restrictions.
 *  Source code is only available in the framework of a collaboration. 
 *
 *  Info: http://www.cscs.ch/molekel/
**/


/* structures for macu */

typedef struct { int ix, iy, iz;
		 int weight, index;
	       } Point; 


typedef struct { float  **value;
                 Point  **point; } Plane;
typedef struct { float xmin, xmax, ymin, ymax, zmin, zmax;
                 int   nx, ny, nz, len;
               } MacuHeader;
typedef struct { unsigned edge0 : 4;
                 unsigned edge1 : 4;
                 unsigned edge2 : 4;
               } tri_pattern;
typedef struct { float **value, vmin, vmax;
                 float p1[3], p2[3], p3[3];
                 short d1, d2;
               } _2D_grid;

#define MACU_SURFACE  110
#define DOT_SURFACE   111
/*----------------------------------------*/
int init_cubes(void);
void cubes(void);
int alloc_tridot(void);
void set_cutplane(Mol *mp);
void readmacu(char *filename, int key);
Plane *plane_alloc(int nx, int ny);
void  *alloc_3D(int nx, int ny, int nz, size_t size);
void  free_plane(Plane *p);
void free_3D(Mol *mp);
void shift_planes(void);
void rectangle(void);
Vector ***do_normals(void);
void surfaceNormals(Vector ***gc);
static float cubeInterpolGradient(float dx, float dy, float dz, Vector *grad,
       float p000, float p100, float p110, float p010,
       float p001, float p101, float p111, float p011);
static void cubeGradient(float dx, float dy, float dz, Vector *grad,
       Vector *p000, Vector *p100, Vector *p110, Vector *p010,
       Vector *p001, Vector *p101, Vector *p111, Vector *p011);
void do_points(void);
void do_cubes(void);
void triangulate(unsigned char lbl, unsigned int *e);
void normalize(Vector *a);
float interpolate(float c1, float v1,
                  float c2, float v2);
Surface *add_surface(char type, Surfdot  *dot, Triangle *tri,
   float contour, float *val, int npts, int ntri, char color, char *name );
void delete_surface(void);
void del_actualsurf(Surface *surf);
void approx_volumes(void);
void project_macu(int axis);
void free_projection(void);
void write_projection(void);
void getPropertyExtrema(Surface *sp);
void exact_volumes(Surface *sp);
double compute_volume(Surface *sp);
void defaultSurfaces(char *fileName);
void fast_surf(void);
void read_gcube(char *filename, int key);
void read_t41(char *fname, int key, char *density);
/*----------------------------------------*/
extern Mol *project_mol;
extern _2D_grid projection;
extern float cubemin, cubemax, cutoff;
extern float rectcol[4];
extern MacuHeader cubehead;

/*----------------------------------------*/
