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


/*
 *   Triangulation of a set of points describing a surface
 *   (Connolly-surface)
 */
#include <GLUT/glut.h>
/*#include <GLUI/glui.h>*/
#include <GL/glui.h>
#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "general.h"
#include "connect.h"
#include "macu.h"
#include "drawing.h"
#include "maininterf.h"
#include "utils.h"
#include "texture.h"
#include "chooseinterf.h"

#ifdef WIN32
#define srand48 srand
#define drand48 rand
#endif

/* no real improvement */
#define DIST(a,b) sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+\
                           (a[2]-b[2])*(a[2]-b[2]))
#define DOT_PRODUCT(a,b)   (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define DETERMINANT(a,b,c) (a[0]*b[1]*c[2]+b[0]*c[1]*a[2]+c[0]*a[1]*b[2]\
                           -a[2]*b[1]*c[0]-b[2]*c[1]*a[0]-c[2]*a[1]*b[0])
/* some improvement */
#define VSUB(a,b,c)  (c[0]=a[0]-b[0],c[1]=a[1]-b[1],c[2]=a[2]-b[2])
#define VADD(a,b,c)  (c[0]=a[0]+b[0],c[1]=a[1]+b[1],c[2]=a[2]+b[2])
#define VMEAN(a,b,c) (c[0]=(a[0]+b[0])*0.5,c[1]=(a[1]+b[1])*0.5,\
                      c[2]=(a[2]+b[2])*0.5)
#define TOOFAR(a,b,c) ((a[0]-b[0])>c)||((a[1]-b[1])>c)||((a[2]-b[2])>c)||\
                      ((b[0]-a[0])>c)||((b[1]-a[1])>c)||((b[2]-a[2])>c)
#define ABS(a) ((a)<0?-(a):(a))

#define MEMBLOCK 1024

typedef struct { Surfdot *a, *b, *opposite;
                 unsigned complete :  1;
               } Edge;
typedef struct { Surfdot **p;
                 Edge    **e;
                 short np, ne;
               } Subspace;

void read_sld2(char *, int, int);
Surface *readnewfile(FILE *, int, int);
Surface *readoldfile(FILE *fsld);
void do_subspace(void);
void tri_loop(void);
void seed_edge(void);
void make_edge(Surfdot *a, Surfdot *b, Surfdot *c);
int checking(Edge *actual_edge, int j, Surfdot *candidate);
void update_edgelist(Edge *actual_edge, int j, Surfdot *chose_one);
float distance(float *a, float *b);
void normalise(float *a);
float angle(float *a, float *b);
float dihedral(float *a, float *b, float *c);

float dotcubesize = 0.25, angle_degrees, dihedral_degrees;

static Surface *surfptr;
static Edge *iedg         = NULL;
static Surfdot *dot, *chosen_one;
static Triangle *tri;
static float *dotval;
static long  npts, npol = 0, nedg = 0, nempty = 0, nnextcubes;
static float area, distance_to_middle, inv;
static float cutoff_angle, cutoff_dihedral;
static int dot_mem;
static float edge_middle[3], potential_edge1[3], potential_edge2[3];
static float middle_normal[3];
static float xmin, xmax, ymin, ymax, zmin, zmax;
static int nx, ny, nz, ix, iy, iz, firstrun, nonormals;
static Subspace ***sspace;
static char dot_name[30];
static int isInside(float p[3], float a[3], float b[3], float c[3]);

void read_dots(char *fname)
{
   FILE *fp;
   char line[100], *s;
   register Surfdot *pt;
   float ar;
   int nitems, i;

   if(!actualmol){
      logprint("load a molecule first");
      update_logs();
      return;
   }

   if((fp = fopen(fname, "r")) == NULL){
      sprintf(line, "Can't open %s\n", fname);
      showinfobox(line);
      return;
   }

   if(!bit.addsurf){
      if(actualmol->firstsurf) delete_surface();
   }
   dot = NULL;
   tri = NULL;
   if((dot = (Surfdot *)malloc(MEMBLOCK*sizeof(Surfdot))) == NULL){
      showinfobox("Can't allocate enough memory for the dots");
      fclose(fp);
      return;
   }
   dot_mem = MEMBLOCK;
   npts = 0;
   pt = dot;
   xmin = ymin = zmin =  1.0e38;
   xmax = ymax = zmax = -1.0e38;
   area = 0;

/*  Reading of surface-point coordinates and normals */
   while(fgets(line, 99, fp)){
      if(line[40] == 'A') continue;
      nitems = sscanf(line+13, "%f%f%f%*s%f%f%f%f",
         pt->v, pt->v+1, pt->v+2, &ar, pt->n, pt->n+1, pt->n+2);
      if(ar > area) area = ar;
      npts++;
      pt++;
      if(npts == dot_mem){
         dot_mem += MEMBLOCK;
         if((dot = (Surfdot *)realloc(dot, dot_mem*sizeof(Surfdot))) == NULL){
            showinfobox("can't reallocate enough memory for the dots");
            fclose(fp);
            return;
         }
         pt = dot + npts;
      }
   }

   /*
    *  shake the points a little to remove perfect alignments of points, 
    *  which is a real bore for triangulation!
    */
    
    srand48(npts); /* same seed for same surface, for reproducability */
    for(pt = dot, i=0; i<npts; i++,  pt++) {

#ifdef WIN32
/*
	pt->v[0] += 0.0002 * (1/(rand()+1) - 0.5);
	pt->v[1] += 0.0002 * (1/(rand()+1) - 0.5);
	pt->v[2] += 0.0002 * (1/(rand()+1) - 0.5);
*/
	pt->v[0] += 0.0002 * ((rand()/RAND_MAX) - 0.5);
	pt->v[1] += 0.0002 * ((rand()/RAND_MAX) - 0.5);
	pt->v[2] += 0.0002 * ((rand()/RAND_MAX) - 0.5);
#else
	pt->v[0] += 0.0002 * (drand48() - 0.5);
	pt->v[1] += 0.0002 * (drand48() - 0.5);
	pt->v[2] += 0.0002 * (drand48() - 0.5);
#endif

      if(pt->v[0] < xmin) xmin = pt->v[0];
      if(pt->v[0] > xmax) xmax = pt->v[0];
      if(pt->v[1] < ymin) ymin = pt->v[1];
      if(pt->v[1] > ymax) ymax = pt->v[1];
      if(pt->v[2] < zmin) zmin = pt->v[2];
      if(pt->v[2] > zmax) zmax = pt->v[2];
   }

   if((dot = (Surfdot *)realloc(dot, npts*sizeof(Surfdot))) == NULL){
      showinfobox("Can't reallocate enough memory for the dots");
      fclose(fp);
      return;
   }

   if(nitems < 7){
      register int i;

      nonormals = 1;
      for(i=0, pt=dot; i<npts; i++, pt++) pt->n[0] = pt->n[1] = pt->n[2] = 0;
   }
   else nonormals = 0;

   if(area<0.1 || area>(zmax-zmin)/2) { /* make a crude estimation */
      area = 2.0*((xmax-xmin)*(ymax-ymin) + (xmax-xmin)*(zmax-zmin) +
                  (ymax-ymin)*(zmax-zmin))/npts;
   }
   dotcubesize = 0.8*sqrt(area);
   angle_degrees = 10;
   dihedral_degrees = 40;

//   drawfields();

#ifndef WIN32
   if((s = strrchr(fname, '/')) == NULL) s = fname;
#else
   if((s = strrchr(fname, '\\')) == NULL) s = fname;
#endif
   else s++;
   strcpy(dot_name, s);
   logprint(s);
   sprintf(line, "   %d points", npts);
   logprint(line);
   update_logs();

   surfptr = add_surface(DOT_SURFACE, dot, NULL, 0, NULL,
                               npts, 0, 0, dot_name);

   glutPostRedisplay();

   got_dots = 1;
   fclose(fp);
}

void read_dots_with_val(char *fname)
{
/* reads dots with an associated value to form a surface
 * 
 * expected format: x y z coordinated "associated area" "associated value"
 * %f %f %f %f %f
 * if no area values are present set the values in the column to 0.0.
 * lines at the top of the file beginning with # are treated as comment
 */

   FILE *fp;
   char line[256], *s;
   register Surfdot *pt;
   float ar, *valp;
   int nitems, i;

   if(!actualmol){
      logprint("load a molecule first");
      update_logs();
      return;
   }

   if((fp = fopen(fname, "r")) == NULL){
      sprintf(line, "Can't open %s\n", fname);
      showinfobox(line);
      return;
   }

   if(!bit.addsurf){
      if(actualmol->firstsurf) delete_surface();
   }

   dot = NULL;
   tri = NULL;
   npts = 0;

   while(fgets(line, 255, fp)){
      if(line[0] == '#') continue;
      npts++;
   }

   if((dot = (Surfdot *)malloc(npts*sizeof(Surfdot))) == NULL){
      showinfobox("Can't allocate enough memory for the dots");
      fclose(fp);
      return;
   }
   if((dotval = (float*) malloc(npts*sizeof(float))) == NULL) {
      showinfobox("can't allocate the value array");
      return;
   }

   pt = dot;
   valp = dotval;
   xmin = ymin = zmin =  1.0e38;
   xmax = ymax = zmax = -1.0e38;
   area = 0;

   rewind(fp);
/*  Reading of surface-point coordinates area and value*/
   while(fgets(line, 255, fp)){
      if(line[0] == '#') continue;
      nitems = sscanf(line, "%f%f%f%f%f", pt->v, pt->v+1, pt->v+2, &ar, valp);
      if(ar > area) area = ar;
      pt++;
      valp++;
   }

   /*
    *  shake the points a little to remove perfect alignments of points, 
    *  which is a real bore for triangulation!
    */
    
    srand48(npts); /* same seed for same surface, for reproducability */
    for(pt = dot, i=0; i<npts; i++,  pt++) {

#ifdef WIN32
/*
	pt->v[0] += 0.0002 * (1/(rand()+1) - 0.5);
	pt->v[1] += 0.0002 * (1/(rand()+1) - 0.5);
	pt->v[2] += 0.0002 * (1/(rand()+1) - 0.5);
*/
	pt->v[0] += 0.0002 * ((rand()/RAND_MAX) - 0.5);
	pt->v[1] += 0.0002 * ((rand()/RAND_MAX) - 0.5);
	pt->v[2] += 0.0002 * ((rand()/RAND_MAX) - 0.5);
#else
	pt->v[0] += 0.0002 * (drand48() - 0.5);
	pt->v[1] += 0.0002 * (drand48() - 0.5);
	pt->v[2] += 0.0002 * (drand48() - 0.5);
#endif

      if(pt->v[0] < xmin) xmin = pt->v[0];
      if(pt->v[0] > xmax) xmax = pt->v[0];
      if(pt->v[1] < ymin) ymin = pt->v[1];
      if(pt->v[1] > ymax) ymax = pt->v[1];
      if(pt->v[2] < zmin) zmin = pt->v[2];
      if(pt->v[2] > zmax) zmax = pt->v[2];
   }

   nonormals = 1;
   for(i=0, pt=dot; i<npts; i++, pt++) pt->n[0] = pt->n[1] = pt->n[2] = 0;

   if(area<0.1 || area>(zmax-zmin)/2) { /* make a crude estimation */
      area = 2.0*((xmax-xmin)*(ymax-ymin) + (xmax-xmin)*(zmax-zmin) +
                  (ymax-ymin)*(zmax-zmin))/npts;
   }
   dotcubesize = 0.8*sqrt(area);
   angle_degrees = 10;
   dihedral_degrees = 40;

#ifndef WIN32
   if((s = strrchr(fname, '/')) == NULL) s = fname;
#else
   if((s = strrchr(fname, '\\')) == NULL) s = fname;
#endif
   else s++;
   strcpy(dot_name, s);
   logprint(s);
   sprintf(line, "   %d points", npts);
   logprint(line);
   update_logs();

   surfptr = add_surface(DOT_SURFACE, dot, NULL, 0, dotval,
                               npts, 0, 0, dot_name);

   glutPostRedisplay();

   got_dots = 1;
   fclose(fp);

}

void set_triangulation_vars(Surface *sp)
/* for surface-dots which have been loaded from quality-files.
   see load_qual in molpro.h
*/
{
   int i;
   Surfdot *pt;

   surfptr = sp;
   dot = sp->dot;
   tri = sp->tri;
   npts = sp->npts;

   xmin = ymin = zmin =  1.0e38;
   xmax = ymax = zmax = -1.0e38;
   area = 0;

   for(pt=dot, i=0; i<npts; pt++, i++){
      if(pt->v[0] < xmin) xmin = pt->v[0];
      if(pt->v[0] > xmax) xmax = pt->v[0];
      if(pt->v[1] < ymin) ymin = pt->v[1];
      if(pt->v[1] > ymax) ymax = pt->v[1];
      if(pt->v[2] < zmin) zmin = pt->v[2];
      if(pt->v[2] > zmax) zmax = pt->v[2];
      pt->n[0] = pt->n[1] = pt->n[2] = 0;
   }

   /* make a crude estimation */
   area = 2.0*((xmax-xmin)*(ymax-ymin) + (xmax-xmin)*(zmax-zmin) +
               (ymax-ymin)*(zmax-zmin))/npts;

   dotcubesize = 0.8*sqrt(area);
   angle_degrees = 10;
   dihedral_degrees = 40;

//   drawfields();
   nonormals = 1;
   got_dots = 1;
}



void tridot(void)
{
   register int i;
   Subspace *sp;
   Edge *ep;
   char str[60];

   inv = 1./dotcubesize;
   nx = (int)fceil((xmax - xmin)*inv);
   ny = (int)fceil((ymax - ymin)*inv);
   nz = (int)fceil((zmax - zmin)*inv);

   if(!got_dots){
      logprint("dot triangulation :");
      logprint(" load a dot-surface first");
      update_logs();
      return;
   }
   if(tri) free(tri);
   if(surfptr->trinorm) free(surfptr->trinorm);
   surfptr->trinorm = NULL;
   tri = NULL;
   npol = nedg = 0;

/* why 3.5* and 2.5* malloc Peter? */
   if((iedg = (Edge *)malloc((int)(3.5*npts*sizeof(Edge)))) == NULL){
      showinfobox("Can't allocate enough memory for the edges");
      return;
   }
   if((tri  = (Triangle *)malloc((int)(2.5*npts*sizeof(Triangle)))) == NULL){
      showinfobox("Can't allocate enough memory for the triangles");
      return;
   }
   if((sspace = (Subspace ***)alloc_3D(nx, ny, nz, sizeof(Subspace))) == NULL){
      showinfobox("Can't allocate enough memory for the subspace");
      return;
   }

   sprintf(str, "nx = %d ny = %d nz = %d", nx, ny, nz);
   logprint(str);

   do_subspace();

/* display interactively in frontbuffer */
   glDrawBuffer(GL_FRONT);
   glPushMatrix();
   global_move();
   glPushMatrix();
   individual_move(actualmol);

   glDisable(GL_LIGHTING);

   set_color(0xff);
   seed_edge();

   nnextcubes = 1;
   firstrun = 1;
   cutoff_angle    = angle_degrees*M_PI/180.;
   cutoff_dihedral = dihedral_degrees*M_PI/180.;
   set_color(0xffff);
   tri_loop();

   if(nempty){
      nnextcubes = 2;
      firstrun = 0;
      cutoff_dihedral += 5.*M_PI/180.;
      set_color(0xff00);
      tri_loop();
   }

   nempty = 0;
   for(ep = iedg; (ep - iedg) < nedg; ep++) if(!ep->complete) nempty++;
   printf("%d free edges\n", nempty);

   free(iedg);
   iedg = NULL;
   for(i=0, sp=sspace[0][0]; i<nx*ny*nz; i++, sp++){
      if(sp->np) free(sp->p);
      if(sp->ne) free(sp->e);
   }
   free(sspace[0][0]);
   free(sspace[0]);
   free(sspace);

   glPopMatrix();
   glPopMatrix();
   glEnable(GL_LIGHTING);
   glDrawBuffer(GL_BACK);

   surfptr->tri = tri;
   surfptr->ntri = npol;

   if(!firsttexture) init_texture();

   surfptr->texture = firsttexture;
   surfptr->texenv = GL_MODULATE;
   surfptr->textype = TEX_MAP;
   set_vmin_vmax(surfptr);

   sprintf(str, "%d triangles %d edges\n", npol, nedg - 1);
   logprint(str); printf(str);
   logprint("Euler characteristics");
   printf("Euler characteristics");
   sprintf(str, " (Pts+Polygons-Edges): %d\n",
          npts + npol - nedg + 1);
   logprint(str); printf(str);
   if(nempty == 0) {
      sprintf(str, "enclosed volume: %f\n", compute_volume(surfptr));
      logprint(str); printf(str);
   }
   update_logs();
}




void do_subspace(void)
{
   register Surfdot *pt;
   register int i;
   int ix, iy, iz, j;
   register Subspace *sp;

   sp = sspace[0][0];
   for(i=0; i<nx*ny*nz; sp++, i++) sp->np = sp->ne = 0;

   for(pt=dot, i=0; i<npts; pt++, i++){
      ix = (int)ffloor((pt->v[0] - xmin)*inv);
      iy = (int)ffloor((pt->v[1] - ymin)*inv);
      iz = (int)ffloor((pt->v[2] - zmin)*inv);
      sp = sspace[iz][iy] + ix;

      if(!sp->np){
         if((sp->p = (Surfdot **)malloc(sizeof(Surfdot *))) == NULL){
            fprintf(stderr, "no place for sp->p\n");
            return;
         }
      }
      else {
         if((sp->p = (Surfdot **)realloc(sp->p, (sp->np+1)*sizeof(Surfdot *))) == NULL){
            fprintf(stderr, "no realloc place for sp->p\n");
            return;
         }
      }

      j = sp->np++;
      sp->p[j] = pt;
   }
}



void tri_loop(void)
{
   float edge_vector[3], middle_to_dot[3], tri_middle[3], new_normal[3];
   register int i, j;
   int jx, jy, jz;
   int jx1, jx2, jy1, jy2, jz1, jz2;
   register Edge *actual_edge;
   register Surfdot *candidate;
   Subspace *sp;
   float quality, alpha, beta, previous_quality = 0;

   for(actual_edge=iedg, j=0; j<nedg; actual_edge++, j++){
      if(actual_edge->complete) continue;
/*
      while (qtest()){
         if((qread(&key) == LEFTMOUSE) && key){
            strcpy(buffer, "Y");
            if(question(&xselection, &yselection,
               "Do you want to abort the triangulation ?", buffer)
               && (buffer[0] == 'Y')) {
               nempty = 0;
               winset(render_window);
               return;
            }
            winset(render_window);
         }
      }
*/

      VSUB(actual_edge->b->v, actual_edge->a->v, edge_vector);
      VMEAN(actual_edge->a->v, actual_edge->b->v, edge_middle);
      VSUB(edge_middle, actual_edge->opposite->v, tri_middle);
      cross_product(edge_vector, tri_middle, middle_normal);
      normalise(middle_normal);
      if(DOT_PRODUCT(middle_normal, middle_normal) == 0.0) {
	  fprintf(stderr, "zero length middle normal");
      }

      ix = (int)ffloor((edge_middle[0] - xmin)*inv);
      iy = (int)ffloor((edge_middle[1] - ymin)*inv);
      iz = (int)ffloor((edge_middle[2] - zmin)*inv);

      jx1 = MAX(ix - nnextcubes, 0);
      jx2 = MIN(ix + nnextcubes + 1, nx);
      jy1 = MAX(iy - nnextcubes, 0);
      jy2 = MIN(iy + nnextcubes + 1, ny);
      jz1 = MAX(iz - nnextcubes, 0);
      jz2 = MIN(iz + nnextcubes + 1, nz);

      chosen_one = NULL;
      for(jz = jz1; jz < jz2; jz++){
       for(jy = jy1; jy < jy2; jy++){
        sp = sspace[jz][jy]+jx1;
        for(jx = jx1; jx < jx2; jx++, sp++){
         for(i=0; i<sp->np; i++){
            candidate = sp->p[i];
            if(candidate == actual_edge->a || 
               candidate == actual_edge->b ||
               candidate == actual_edge->opposite)   continue;
            VSUB(candidate->v, edge_middle, middle_to_dot);
            VSUB(candidate->v, actual_edge->a->v, potential_edge1);
            VSUB(candidate->v, actual_edge->b->v, potential_edge2);

            if((beta = ABS(dihedral(middle_normal, edge_vector, middle_to_dot)))
                                                > cutoff_dihedral) continue;

            if(firstrun &&
                ((alpha = angle(potential_edge1, potential_edge2))
                                                < cutoff_angle)) continue;

            if(!(checking(actual_edge, j, candidate))) continue;
/*
	    quality = alpha;
            quality = distance_to_middle - (alpha + beta)*0.2;
	    quality = aspectRatio(actual_edge->a->v, actual_edge->b->v,
				  candidate->v);
*/
            distance_to_middle = length(middle_to_dot);
            quality = distance_to_middle;

            if((!chosen_one) || (quality < previous_quality)){
               chosen_one = candidate;
               previous_quality = quality;
            }
         }
        }
       }
      }
      if(!chosen_one) {
         nempty++;
         continue;
      }
/*      fprintf(stdout, "%7d %7d %7d\n",
             actual_edge->a - dot, actual_edge->b - dot, chosen_one - dot); */

      tri[npol].p1 = actual_edge->a - dot;
      tri[npol].p2 = actual_edge->b - dot;
      tri[npol].p3 = chosen_one - dot;

      if(nonormals) {
         VSUB(chosen_one->v, edge_middle, middle_to_dot);
         cross_product(edge_vector, middle_to_dot, new_normal);
         normalise(new_normal);
         chosen_one->n[0] += new_normal[0];
         chosen_one->n[1] += new_normal[1];
         chosen_one->n[2] += new_normal[2];
         actual_edge->a->n[0] += new_normal[0];
         actual_edge->a->n[1] += new_normal[1];
         actual_edge->a->n[2] += new_normal[2];
         actual_edge->b->n[0] += new_normal[0];
         actual_edge->b->n[1] += new_normal[1];
         actual_edge->b->n[2] += new_normal[2];
      }

      glBegin(GL_LINE_LOOP);
      glVertex3fv(actual_edge->a->v);
      glVertex3fv(actual_edge->b->v);
      glVertex3fv(chosen_one->v);
      glEnd();
      npol++;

      update_edgelist(actual_edge, j, chosen_one);

   }

   if(nonormals){
      for(i=0, candidate = dot; i<npts; i++, candidate++)
         normalise(candidate->n);
   }
}




void seed_edge(void)
{
   register int i, jx;
   int jy, jz;
   int jx1, jx2, jy1, jy2, jz1, jz2;
   float first, second, crit;
   float vector1[3], vector2[3], middle[3];
   static float normal[3] = {0, 0, 1};
   Subspace *sp;
   Surfdot *candidate, *point_one, *point_two;
   static Surfdot phantom1, phantom2;

/* the first point : the highest one in z */
   first = -1.e38;
   for(jx = 0, candidate = dot; jx < npts; jx++, candidate++){
      if(candidate->v[2] > first){
         point_one = candidate;
         first = candidate->v[2];
      }
   }

/* the second point : its nearest neighbour */
   second = 1.e38;
   ix = (int)ffloor((point_one->v[0] - xmin)*inv);
   iy = (int)ffloor((point_one->v[1] - ymin)*inv);
   iz = (int)ffloor((point_one->v[2] - zmin)*inv);
   jx1 = MAX(ix - 2, 0);
   jx2 = MIN(ix + 3, nx);
   jy1 = MAX(iy - 2, 0);
   jy2 = MIN(iy + 3, ny);
   jz1 = MAX(iz - 2, 0);
   jz2 = MIN(iz + 3, nz);

   for(jz = jz1; jz < jz2; jz++){
      for(jy = jy1; jy < jy2; jy++){
         sp = sspace[jz][jy]+jx1;
         for(jx = jx1; jx < jx2; jx++, sp++){
            for(i=0; i<sp->np; i++){
               if((candidate = sp->p[i]) == point_one) continue;
               if((crit = distance(point_one->v, candidate->v)) < second &&
                   crit > 0){
                  second = crit;
                  point_two = candidate;
               }
            }
         }
      }
   }

   VSUB(point_two->v, point_one->v, vector1);
   VMEAN(point_one->v, point_two->v, middle);



   cross_product(vector1, normal, vector2);
   VADD(middle, vector2, phantom1.v);
   VSUB(middle, vector2, phantom2.v);

   npol = 0;
   glBegin(GL_LINES);
   glVertex3fv(point_one->v);
   glVertex3fv(point_two->v);
   glEnd();

   nedg = nempty = 0;
   make_edge(point_one, point_two, &phantom2);
   make_edge(point_two, point_one, &phantom1);

}




void make_edge(Surfdot *a, Surfdot *b, Surfdot *c)
{
   float middle[3];
   int ix, iy, iz;
   Subspace *sp;

   iedg[nedg].a = b;
   iedg[nedg].b = a;
   iedg[nedg].opposite = c;
   iedg[nedg].complete = 0;

   VMEAN(b->v, a->v, middle);

   ix = (int)ffloor((middle[0] - xmin)*inv);
   iy = (int)ffloor((middle[1] - ymin)*inv);
   iz = (int)ffloor((middle[2] - zmin)*inv);

   sp = sspace[iz][iy]+ix;

   if(!sp->ne){
      if((sp->e = (Edge **)malloc(sizeof(Edge *))) == NULL){
         showinfobox("no place for sp->e");
         return;
      }
   }
   else {
      if((sp->e = (Edge **)realloc(sp->e, (sp->ne + 1)*sizeof(Edge *))) == NULL){
         showinfobox("can't reallocate sp->e");
         return;
      }
   }

   sp->e[sp->ne] = iedg + nedg;
   sp->ne++;
   nedg++;
}




int checking(Edge *actual_edge, int j, Surfdot *candidate)
{
   register int k;
   register Edge *other_edge;
   float test_edge[3], e1[3], e2[3], new_normal[3];
   float bcd1, bcd2, x1, x2, x3; 
   int jx, jy, jz;
   int jx1, jx2, jy1, jy2, jz1, jz2;
   Subspace *sp;

   jx1 = MAX(ix - nnextcubes, 0);
   jx2 = MIN(ix + nnextcubes + 1, nx);
   jy1 = MAX(iy - nnextcubes, 0);
   jy2 = MIN(iy + nnextcubes + 1, ny);
   jz1 = MAX(iz - nnextcubes, 0);
   jz2 = MIN(iz + nnextcubes + 1, nz);

/*** first check whether a, b and candidate are on a straight line ***/

    VSUB(actual_edge->b->v, actual_edge->a->v, e1);
    VSUB(candidate->v, actual_edge->a->v, e2);
    cross_product(e1, e2, new_normal);
    if(DOT_PRODUCT(new_normal, new_normal) == 0.0) {
	printf(" rejected %d : straight line\n", candidate - dot);
	return NULL;
    }


/*** then check if any other point is inside the potential triangle ***/

    for(jz = jz1; jz < jz2; jz++){
	for(jy = jy1; jy < jy2; jy++){
	    sp = sspace[jz][jy]+jx1;
	    for(jx = jx1; jx < jx2; jx++, sp++){

		Surfdot *point;

		for(k=0; k<sp->np; k++){
		    point = sp->p[k];
	    
		    if(point == actual_edge->a || 
			point == actual_edge->b ||
			point == candidate) continue;
		    
		    if(isInside(point->v, actual_edge->a->v, 
				actual_edge->b->v, candidate->v)) {
			return NULL;
		    }
		}
	    }
	}
    }


/*** then check if any existing edge would cross one of the potential edges ***/

   for(jz = jz1; jz < jz2; jz++){
    for(jy = jy1; jy < jy2; jy++){
     sp = sspace[jz][jy]+jx1;
     for(jx = jx1; jx < jx2; jx++, sp++){
      for(k=0; k<sp->ne; k++){
         other_edge = sp->e[k];
         if(actual_edge->a == other_edge->a &&
            candidate == other_edge->b) return NULL;
         if(actual_edge->b == other_edge->b &&
            candidate == other_edge->a) return NULL;

         if(other_edge->complete){
            if(actual_edge->a == other_edge->b &&
               candidate == other_edge->a) return NULL;
            if(actual_edge->b == other_edge->a &&
               candidate == other_edge->b) return NULL;
         }

         if(candidate == other_edge->a ||
            candidate == other_edge->b) continue;

         VSUB(other_edge->b->v, other_edge->a->v, test_edge);


	 if(actual_edge->a == other_edge->b ||
	    actual_edge->b == other_edge->a) continue;


/*
 * Determination of the intersection of the plane potential_edge,
 * middle_normal with the straight line through other_egde :
 *   plane         : Z + X1*potential_edge + X2*middle_normal
 *   straight line : other_edge->a + X3*test_edge
 *
 * then : X1*potential_edge + X2*middle_normal - X3*test_edge = E  
 *        with E = other_edge->a - Z
 */
         VSUB(other_edge->a->v, actual_edge->a->v, e1);
         VSUB(other_edge->a->v, actual_edge->b->v, e2);

         if(bcd1 = determinant(potential_edge1, middle_normal, test_edge)){
            x2 = determinant(potential_edge1, e1, test_edge)/bcd1;
            if(ABS(x2) <= dotcubesize*0.5) {
               x1 = determinant(e1, middle_normal, test_edge)/bcd1;
               x3 = determinant(potential_edge1, middle_normal, e1)/bcd1;
               if((x1 >  0.001) && (x1 <  .999) &&
                  (x3 > -.999) && (x3 <  -0.001)) return NULL;
            }
         }
         if(bcd2 = determinant(potential_edge2, middle_normal, test_edge)){
            x2 = determinant(potential_edge2, e2, test_edge)/bcd2;
            if(ABS(x2) <= dotcubesize*0.5) {
               x1 = determinant(e2, middle_normal, test_edge)/bcd2;
               x3 = determinant(potential_edge2, middle_normal, e2)/bcd2;
               if((x1 >  0.001) && (x1 <  .999) &&
                  (x3 > -.999) && (x3 <  -0.001)) return NULL;
            }
         }
      }
     }
    }
   }

   return 1;
}


static int isInside(float p[3], float a[3], float b[3], float c[3])
/*
 * returns 1 if p is inside (or close above or below) the triangle abc
 * otherwise 0
 */
{
    float ab[3], ac[3], ap[3], n[3], x0, x1, x2;
    float det;

    VSUB(b, a, ab);
    VSUB(c, a, ac);
    VSUB(p, a, ap);
    cross_product(ab, ac, n);
    normalise(n);
    
    det = determinant(ab, ac, n);
    if(det == 0.0) {
	fprintf(stderr, "warning: null determinant\n");
	return 1;
    }
    x2 = determinant(ab, ac, ap)/det;
    if(ABS(x2) > dotcubesize*0.5) return 0;
    
    x0 = determinant(ap, ac, n)/det;
    x1 = determinant(ab, ap, n)/det;
    
    if(x0 < 0.0 || x1 < 0.0) return 0;
    if(x0 + x1 > 1.0) return 0;

    return 1;
}



void update_edgelist(Edge *actual_edge, int j, Surfdot *chosen_one)
{
   register int i;
   register Edge *other_edge;
   int jx, jy, jz;
   int ix, iy, iz;
   int jx1, jx2, jy1, jy2, jz1, jz2;
   Subspace *sp;

   ix = (int)ffloor((edge_middle[0] - xmin)*inv);
   iy = (int)ffloor((edge_middle[1] - ymin)*inv);
   iz = (int)ffloor((edge_middle[2] - zmin)*inv);

   jx1 = MAX(ix - 2, 0);
   jx2 = MIN(ix + 3, nx);
   jy1 = MAX(iy - 2, 0);
   jy2 = MIN(iy + 3, ny);
   jz1 = MAX(iz - 2, 0);
   jz2 = MIN(iz + 3, nz);

   for(jz = jz1; jz < jz2; jz++){
      for(jy = jy1; jy < jy2; jy++){
         sp = sspace[jz][jy]+jx1;
         for(jx = jx1; jx < jx2; jx++, sp++){
            for(i=0; i<sp->ne; i++){
               other_edge = sp->e[i];
               if((actual_edge->b == other_edge->a) &&
                  (chosen_one == other_edge->b)){
                  other_edge->complete = 1;
                  goto X;
               }
            }
         }
      }
   }

   make_edge(actual_edge->b, chosen_one, actual_edge->a);

X :
   for(jz = jz1; jz < jz2; jz++){
      for(jy = jy1; jy < jy2; jy++){
         sp = sspace[jz][jy]+jx1;
         for(jx = jx1; jx < jx2; jx++, sp++){
            for(i=0; i<sp->ne; i++){
               other_edge = sp->e[i];
               if((actual_edge->a == other_edge->b) &&
                  (chosen_one == other_edge->a)){
                  other_edge->complete = 1;
                  goto Y;
               }
            }
         }
      }
   }

   make_edge(chosen_one, actual_edge->a, actual_edge->b);

Y :

   actual_edge->complete = 1;
}





/*
 *  vector-functions
 */


float distance(float *a, float *b)
{
   return sqrt((a[0]-b[0])*(a[0]-b[0]) +  (a[1]-b[1])*(a[1]-b[1]) + 
               (a[2]-b[2])*(a[2]-b[2]));
}



void normalise(float *a)
{
   float l, invl;

   if(!(l = length(a))) return;
   invl = 1.0/l;

   a[0] *= invl;
   a[1] *= invl;
   a[2] *= invl;
}




float angle(float *a, float *b)
{
   return facos(dot_product(a, b)/(length(a)*length(b)));
}



float dihedral(float *n1, float *b, float *c)
{
   float n2[3];

   cross_product(b, c, n2);
   if(DOT_PRODUCT(n2, n2) == 0.0) {
       fprintf(stderr, "zero length cross product\n");
       return 180.0;
   }

   return angle(n1, n2);
}




/*
 * aspectRatio returns twice the ratio of inscribed circle radius by
 * circumscribed circle ratio of a triangle as a triangle shape
 * quality descriptor. Equilateral triangles get 1.0
 * (the bigger, the better).
 *
 *		    pff january 97 
 */

float aspectRatio(float A[3], float B[3], float C[3])
{
    float a, b, c, s;
    float alpha;
    float umkreisRadius, inkreisRadius;
    
    a = DIST(B, C); b = DIST(A, C); c = DIST(A, B);

    alpha = acos((b*b + c*c - a*a)/(2.0*b*c));
    umkreisRadius = a / (2.0 * sin(alpha));

    s = 0.5*(a+b+c);
    inkreisRadius = sqrt((s-a)*(s-b)*(s-c)/s);
    
    return 2.0 * inkreisRadius / umkreisRadius;
}



/*****************************/


#define SLD_QUAL 0x00000001
#define SLD_NORM 0x00000001
#define SLD_ATOM 0x00000002
#define SLD_RESD 0x00000004
#define SLD_GRAD 0x00000008
#define QLEN 30
#define OLD_MAGIC -4711
#define NEW_MAGIC -4712


void read_sld(char *sldfile)
{
   read_sld2(sldfile, 0, 0);
}




void read_sld2(char *sldfile, int prop1, int prop2)
{
/* even it looks like you can read a second property, only one can be read
 * what OLD_MAGIC was is unknown
 * proper behaviour is tested for files generated with the OpenGL version of molekel
 * older files are not tested
*/

   FILE *fsld;
   int magic;
   Surface *sp;
   char str[150];

   if ((fsld = fopen(sldfile, "rb")) == NULL) {
      sprintf(str, "Can't open %s \n", sldfile);
      showinfobox(str);
      return;
   }

   fread(&magic, sizeof(int), 1, fsld);

   switch(magic){
      case OLD_MAGIC : sp = readoldfile(fsld);
                       break;
      case NEW_MAGIC : sp = readnewfile(fsld, prop1, prop2);
                       break;
      default        : showinfobox("wrong magic number");
   }

   exact_volumes(sp);

   glutPostRedisplay();
}


   
Surface *readnewfile(FILE *fsld, int prop1, int prop2)
{
   struct Header { int headlen; int trilen; int dotlen;
                   int ntri; int npts; int ncon; int natm; int nqual;
                   int t_flag; } *phead, header;
   char line[120];
   float *val, *val1, *valp, *valp1, fdum;
   int choice, choice1, nobjects, actual, second;
   short ato;
   int npospts, nnegpts, npostri, nnegtri, type, ntri;
   register int i, j;
   Triangle *trip, *tri2;
   Surfdot *dotp, *dot2;
   Surface *s1, *s2;

   if(!actualmol){
      logprint("Load a molecule first!");
      update_logs();
      return NULL;
   }

   if(!bit.addsurf){
      if(actualmol && actualmol->firstsurf) delete_surface();
   }
   dot = NULL;
   tri = NULL;

   /* The header */

   phead = &header;
   fread(phead, sizeof(struct Header), 1, fsld);

   sprintf(line, " %d dots and %d triangles \n", phead->npts, phead->ntri);
   logprint(line);
   update_logs();

/* allocation */

   if((tri = trip = (Triangle *)malloc(phead->ntri*sizeof(Triangle))) == NULL){
      showinfobox(" Can't allocate memory for the triangles !");
      return NULL;
   }
   if((dot = dotp = (Surfdot *)malloc(phead->npts*sizeof(Surfdot))) == NULL){
      showinfobox(" Can't allocate memory for the points !");
      return NULL;
   }

   if(phead->nqual > 0){
      if((val = valp = (float *)malloc(phead->npts*sizeof(float))) == NULL){
         showinfobox(" Can't allocate memory for the values !");
         return NULL;
      }
      val1 = NULL;
   }
   else val = val1 = NULL;

/*** read quality strings */

   choice = choice1 = -1;

   if(phead->nqual > 1){
      register int i;
      char **strptr;

      if(prop1) {
         choice = prop1 - 1;
         if(prop2) choice1 = prop2 - 1;
         for(i=0; i<phead->nqual; i++)
            fread(line, sizeof(char), QLEN, fsld);
         goto BED;
      }

      if((strptr = (char **)malloc((phead->nqual+1)*sizeof(char *))) == NULL){
         showinfobox("Can't allocate string-pointer");
         return NULL;
      }
      for(i=0; i<phead->nqual; i++){
         fread(line, sizeof(char), QLEN, fsld);
         line[QLEN-1] = '\0';
         if((strptr[i] = strdup(line)) == NULL){
            showinfobox("Can't allocate selector-string");
            return NULL;
         }
      }
      if((strptr[i] = strdup("None")) == NULL){
         showinfobox("Can't allocate selector-string");
         return NULL;
      }
/*
      choice = selection(&xselection, &yselection,
                         "Choose a first quality", phead->nqual+1, strptr);

      if((choice >= 0) && (choice < phead->nqual)) {
         choice1 = selection(&xselection, &yselection,
                         "Choose a second quality", phead->nqual+1, strptr);
      }
*/
      choice = -1; /* to be fixed */

      for(i=0; i<=phead->nqual; i++) free(strptr[i]);
      free(strptr);

BED:
      if(choice < 0) return NULL;
      if(choice >= phead->nqual) return NULL; /* "None" */

   }
   else if(phead->nqual == 1){
      fread(line, sizeof(char), QLEN, fsld);
      line[QLEN-1] = '\0';
      printf(" loading %s\n", line);
      choice = 0;
   }


   if((choice1 >= 0) && (choice1 < phead->nqual)) {
      if((val1 = valp1 = (float *)malloc(phead->npts*sizeof(float))) == NULL) {
         showinfobox(" Can't allocate memory for the second values !");
         return NULL;
      }
   }

   for (i=0; i<phead->ntri; i++)
      fread(trip++, sizeof(Triangle), 1, fsld);

   nobjects = second = actual = 0;
   for (i=0; i<phead->npts; i++, dotp++){
      fread (dotp->v, sizeof(float), 3, fsld);
      if(phead->t_flag & SLD_NORM) fread (dotp->n, sizeof(float), 3, fsld);
      if(phead->t_flag & SLD_ATOM){
         fread (&ato, sizeof(short), 1, fsld);
         if(ato != actual){
            nobjects++;
            actual = ato;
            if((!second) && (nobjects == 2)) second = i;
         }
      }
      for(j=0; j<phead->nqual; j++) {
         if(j == choice) {
            fread (valp, sizeof(float), 1, fsld);
            valp++;
         }
         else if(j == choice1) {
            fread (valp1, sizeof(float), 1, fsld);
            valp1++;
         }
         else fread (&fdum, sizeof(float), 1, fsld);
      }
   }

   if (nobjects == 2) {
      for(trip = tri, i = 0; i<phead->ntri ; i++, trip++) {
         if((trip->p1 == second) || (trip->p2 == second) ||
            (trip->p3 == second)) {
            npostri = i;
            nnegtri = phead->ntri - i;
            break;
         }
      }
      npospts = second;
      nnegpts = phead->npts - npospts;
      type = MACU_SURFACE;

      if((tri2 = (Triangle *)malloc(nnegtri*sizeof(Triangle))) == NULL) {
         showinfobox(" Can't allocate memory for the second set of triangles !");
         return NULL;
      }
      if((dot2 = (Surfdot *)malloc(nnegpts*sizeof(Surfdot))) == NULL) {
         showinfobox(" Can't allocate memory for the second set of points !");
         return NULL;
      }
      memcpy(tri2, tri + npostri, nnegtri*sizeof(Triangle));
      memcpy(dot2, dot + npospts, nnegpts*sizeof(Surfdot));

      if((tri = (Triangle *)realloc(tri, npostri*sizeof(Triangle))) == NULL) {
         showinfobox(" Can't reallocate memory for the first set of triangles !");
         return NULL;
      }
      if((dot = (Surfdot *)realloc(dot, npospts*sizeof(Surfdot))) == NULL) {
         showinfobox(" Can't reallocate memory for the first set of points !");
         return NULL;
      }

      for(trip=tri2, i=0; i<nnegtri; i++, trip++) {
         trip->p1 -= npospts;
         trip->p2 -= npospts;
         trip->p3 -= npospts;
      }

      s1 = add_surface(type, dot, tri, 0, val, npospts, npostri, 0, ".sld");
      s2 = add_surface(type, dot2, tri2, 0, val, nnegpts, nnegtri, 1, ".sld");
      s1->second = s2;
      s2->second = s1;
   }
   else {
      ntri = phead->ntri;
      npts = phead->npts;
      type = DOT_SURFACE;
      s1 = add_surface(type, dot, tri, 0, val, npts, ntri, 0, ".sld");
      s1->val1 = val1;
   }

   fclose(fsld);

   if(val) {
      s1->vmin = s1->vmax = *val;
      for(valp=val+1, i=0; i<npts; valp++, i++){
         if(*valp < s1->vmin) s1->vmin = *valp;
         if(*valp > s1->vmax) s1->vmax = *valp;
      }
         if(!firsttexture) init_texture();
         s1->texture = firsttexture;
         s1->texenv = GL_MODULATE;
         s1->textype = TEX_MAP;
   }
   if(val1) {
      s1->vmin1 = s1->vmax1 = *val1;
      for(valp1=val1+1, i=0; i<npts; valp1++, i++) {
         if(*valp1 < s1->vmin1) s1->vmin1 = *valp1;
         if(*valp1 > s1->vmax1) s1->vmax1 = *valp1;
      }
   }

   return s1;
}



Surface *readoldfile(FILE *fsld)
{
   struct Header { int headlen; int trilen; int dotlen;
                   int ntri; int npts; int ncon; int natm; int t_flag; };
   struct Header *phead, header;
   char line[120];
   float *val, *valp, fdum;
   int nobjects, actual, second;
   short ato;
   int npospts, nnegpts, npostri, nnegtri, type, ntri;
   register int i;
   Triangle *trip, *tri2;
   Surfdot *dotp, *dot2;
   Surface *s1, *s2;

   if(!actualmol){
      logprint("Load a molecule first!");
      update_logs();
      return NULL;
   }

   if(!bit.addsurf){
      if(actualmol->firstsurf) delete_surface();
   }
   dot = NULL;
   tri = NULL;

   /* The header */

   phead = &header;
   fread(phead, sizeof(struct Header), 1, fsld);

   sprintf(line, " %d dots and %d triangles \n", phead->npts, phead->ntri);
   logprint(line);
   update_logs();

/* allocation */

   if((tri = trip = (Triangle *)malloc(phead->ntri*sizeof(Triangle))) == NULL){
      showinfobox(" Can't allocate memory for the triangles !");
      return NULL;
   }
   if((dot = dotp = (Surfdot *)malloc(phead->npts*sizeof(Surfdot))) == NULL){
      showinfobox(" Can't allocate memory for the points !");
      return NULL;
   }

   if(phead->t_flag & SLD_QUAL){
      if((val = valp = (float *)malloc(phead->npts*sizeof(float))) == NULL){
         showinfobox(" Can't allocate memory for the values !");
         return NULL;
      }
   }
   else val = NULL;

   for (i=0; i<phead->ntri; i++)
      fread(trip++, sizeof(Triangle), 1, fsld);

   nobjects = second = actual = 0;
   for (i=0; i<phead->npts; i++, dotp++){
      fread (dotp->v, sizeof(float), 3, fsld);
      fread (dotp->n, sizeof(float), 3, fsld);
      if (phead->t_flag & SLD_QUAL)
         fread (valp++, sizeof(float), 1, fsld);
      if (phead->t_flag & SLD_ATOM){
         fread (&ato, sizeof(short), 1, fsld);
         if(ato != actual){
            nobjects++;
            actual = ato;
            if((!second) && (nobjects == 2)) second = i;
         }
      }
      if (phead->t_flag & SLD_RESD) fread (&ato, sizeof(short), 1, fsld);
      if (phead->t_flag & SLD_GRAD) fread (&fdum, sizeof(float), 1, fsld);
   }

   if (nobjects == 2){
      for(trip = tri, i = 0; i<phead->ntri ; i++, trip++){
         if((trip->p1 == second) || (trip->p2 == second) ||
            (trip->p3 == second)){
            npostri = i;
            nnegtri = phead->ntri - i;
            break;
         }
      }
      npospts = second;
      nnegpts = phead->npts - npospts;
      type = MACU_SURFACE;

      if((tri2 = (Triangle *)malloc(nnegtri*sizeof(Triangle))) == NULL){
         showinfobox(" Can't allocate memory for the second set of triangles !");
         return NULL;
      }
      if((dot2 = (Surfdot *)malloc(nnegpts*sizeof(Surfdot))) == NULL){
         showinfobox(" Can't allocate memory for the second set of points !");
         return NULL;
      }
      memcpy(tri2, tri + npostri, nnegtri*sizeof(Triangle));
      memcpy(dot2, dot + npospts, nnegpts*sizeof(Surfdot));

      if((tri = (Triangle *)realloc(tri, npostri*sizeof(Triangle))) == NULL){
         showinfobox(" Can't reallocate memory for the first set of triangles !");
         return NULL;
      }
      if((dot = (Surfdot *)realloc(dot, npospts*sizeof(Surfdot))) == NULL){
         showinfobox(" Can't allocate memory for the first set of points !");
         return NULL;
      }

      for(trip=tri2, i=0; i<nnegtri; i++, trip++){
         trip->p1 -= npospts;
         trip->p2 -= npospts;
         trip->p3 -= npospts;
      }

      s1 = add_surface(type, dot, tri, 0, val, npospts, npostri, 0, ".sld");
      s2 = add_surface(type, dot2, tri2, 0, val, nnegpts, nnegtri, 1, ".sld");
      s1->second = s2;
      s2->second = s1;
   }
   else {
      ntri = phead->ntri;
      npts = phead->npts;
      type = DOT_SURFACE;
      s1 = add_surface(type, dot, tri, 0, val, npts, ntri, 0, ".sld");
   }

   fclose(fsld);

   if(val) {
      s1->vmin = s1->vmax = *val;
      for(valp=val+1; --npts; valp++){
         if(*valp < s1->vmin) s1->vmin = *valp;
         if(*valp > s1->vmax) s1->vmax = *valp;
      }
   }

   return s1;
}

