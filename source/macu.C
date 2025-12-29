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
   functions concerning the marching cubes algorithm
*/
#ifndef WIN32
#include <sys/param.h>
#endif

#include "main.h"
#include "constant.h"
#include "molekel.h"
#include "general.h"
#include "glutwin.h"
#include "macu.h"
#include "patterns.h"
#include "utils.h"
#include "drawing.h"
#include "material.h"
#include "maininterf.h"
#include "surfaceinterf.h"
#include "texture.h"
#include "manip.h"
#include "connect.h"
#include "readadf.h"

/*----------------------------------------*/

#define MEMBLOCK 1024
#define ORBITAL_CUTOFF 0.05
#define DENSITY_CUTOFF 0.01

/*----------------------------------------*/
/* Globals */
float cutoff = 0, cubemin = 0, cubemax = 0;
Mol *project_mol;
/*----------------------------------------*/
static Surfdot *dot;
static Triangle *tri;
int nMacuPlanes;
char  macu_name[30];

MacuHeader    cubehead;
_2D_grid      projection = {NULL};

static float         xmin, xmax, ymin, ymax, zmin, zmax;
static float         za, zb, dx, dy, dz, invdx, invdy, invdz;
static int           nx, ny, nz, iz;
static int           dotmem, trimem;
static int           npts, ntri, surface_number = 0;
float                rectcol[] = {1.0, 1.0, 0.0, 1.0};
static tri_pattern   **pattern;
static unsigned char *npat;
unsigned char        **ccase;
static Plane         *A, *B, *C;

/* stuff for read_t41; not very elegant */
FILE *fpt;
static long previous_line = 0, preprevious = 0, preprepre = 0;
static char line[256];

char *find_string_t41(char *s);

/* _A : patterns around maximum (points outside cutoff are connected)
   _B : patterns around minimum (points inside cutoff are connected) */

/*----------------------------------------*/

void free_3D(Mol *mp)
{
   if(!mp) return;

   mp->got_macufile = 0;
   if(mp->cube_value){
      free(mp->cube_value[0][0]);
      free(mp->cube_value[0]);
      free(mp->cube_value);
      mp->cube_value = NULL;
   }
   if(ccase){
      free(ccase[0]);
      free(ccase);
      ccase = NULL;
   }

   mp->cubemin = mp->cubemax = mp->cutoff = 0;
   strcpy(macu_name, "");

   cubemin = cubemax = cutoff = 0;
   if(mp->plane) {
      if(mp->plane->plane_point) {
         free(mp->plane->plane_point[0]);
         free(mp->plane->plane_point);
      }
      if(mp->plane->tri) free(mp->plane->tri);
      free(mp->plane);
      mp->plane = NULL;
   }

   if(surfaceglui) {
      surfaceglui->sync_live();
   }
   if(planeglui) {
      planeglui->sync_live();
   }

}




Plane *plane_alloc(int nx, int ny)
{
   Point  *pointarray, *pp;
   Plane  *temp;
   int i;

   if((temp = (Plane*) malloc(sizeof(Plane))) == NULL) return NULL;

   if((pointarray = (Point*) malloc(nx*ny*sizeof(Point))) == NULL) return NULL;
   if((temp->point = (Point**) malloc(ny*sizeof(Point *))) == NULL) return NULL;

   for(i=0; i<ny; i++){
      temp->point[i]  = pointarray + i*nx;
   }

   for(i=0, pp = pointarray; i<nx*ny; i++, pp++){
      pp->weight = pp->index = 0;
      pp->ix = pp->iy = pp->iz = 0;
   }

   return temp;
}



/*
 * allocate dynamically a 3D-array  (or a 2D-array if nz == 0)
 *   which can be used through indices like an array that has been
 *   declared as array[nz][ny][nx] (or array[ny][nx])
 */
void *alloc_3D(int nx, int ny, int nz, size_t size)
{
   void ***temp;
   void **pointerarray;
   char *array;
   int  _2D;
   register int i;

   if(_2D = !nz) nz = 1;
   if((array = (char*) malloc(nx*ny*nz*size)) == NULL) return NULL;
      /* array will hold the data */
   if((pointerarray = (void**) malloc(ny*nz*sizeof(char *))) == NULL) return NULL;
      /* pointerarray will hold the pointers to the rows of data */
   for(i=0; i<ny*nz; i++) pointerarray[i] = array + i * nx * size;

   if(_2D) return pointerarray;

   if((temp = (void***) malloc(nz*sizeof(char **))) == NULL) return NULL;
      /* temp will hold the pointers to the planes */
   for(i=0; i<nz; i++)    temp[i] = pointerarray + i * ny;

   return temp;
}



void free_plane(Plane *p)
{
   free(p->point[0]);
   free(p->point);
   free(p);
}



void cubes(void)
{
   char str[50];
   int  npass;
   Surface *s1, *s2;
   Vector ***gradientCube;

   dot = NULL;
   tri = NULL;

   xmin = actualmol->box.x1;
   xmax = actualmol->box.x2;
   ymin = actualmol->box.y1;
   ymax = actualmol->box.y2;
   zmin = actualmol->box.z1;
   zmax = actualmol->box.z2;
   nx =  actualmol->box.nx;
   ny =  actualmol->box.ny;
   nz =  actualmol->box.nz;

   if(init_cubes()) return;

   npass = (bit.both_signs && cutoff) ? 2 : 1;


   dx = (xmax-xmin)/((float)nx-1.);
   dy = (ymax-ymin)/((float)ny-1.);
   dz = (zmax-zmin)/((float)nz-1.);
/*
   dx = (actualmol->box.x2-actualmol->box.x1)/((float)actualmol->box.nx-1);
   dy = (actualmol->box.y2-actualmol->box.y1)/((float)actualmol->box.ny-1);
   dz = (actualmol->box.z2-actualmol->box.z1)/((float)actualmol->box.nz-1);
*/
   invdx = 0.5/dx;
   invdy = 0.5/dy;
   invdz = 0.5/dz;

/* display interactively in frontbuffer */
   glDrawBuffer(GL_FRONT);
   glPushMatrix();
   global_move();
   glPushMatrix();
   individual_move(actualmol);


   if(!bit.addsurf){
      if(actualmol->firstsurf) delete_surface();
   }

/*** big loop ***/

/* to be fixed
   lmbind(MATERIAL, actualmol->nsurfs%21 + 100);
*/
   glCallList(surface_mat_base + actualmol->nsurfs%21);
   if(cutoff >= 0) {
      pattern = pat_B;
      npat    = npat_B;
   }
   else {
      pattern = pat_A;
      npat    = npat_A;
   }
   ntri = npts = 0;
   if(!alloc_tridot()) return;
   zb = zmin - dz;
   memset(ccase[0], 0, (nx+1)*(ny+1));
   for(iz = 0; iz < nz; iz++){
      shift_planes();
      do_points();
      do_cubes();
   }
   if(gradientCube = do_normals()) {
      surfaceNormals(gradientCube);
   }
   s1 = add_surface(MACU_SURFACE, dot, tri, cutoff, NULL,
                       npts, ntri, 0, macu_name );
   s2 = NULL;
   if(npass == 2 && -cutoff > actualmol->cubemin){
/* to be fixed
      lmbind(MATERIAL, actualmol->nsurfs%21 + 122);
*/
      glCallList(surface_mat_base + actualmol->nsurfs%21 + 22);
      cutoff  = -cutoff;
      if(cutoff >= 0) {
         pattern = pat_B;
         npat    = npat_B;
      }
      else {
         pattern = pat_A;
         npat    = npat_A;
      }
      ntri = npts = 0;
      if(!alloc_tridot()) return;
      zb = zmin - dz;
      memset(ccase[0], 0, (nx+1)*(ny+1));
      for(iz = 0; iz < nz; iz++){
         shift_planes();
         do_points();
         do_cubes();
      }
      if(gradientCube) {
         surfaceNormals(gradientCube);
      }
      
      s2 = add_surface(MACU_SURFACE, dot, tri, cutoff, NULL,
                       npts, ntri, 1, macu_name );
      s1->second = s2;
      s2->second = s1;
//      cutoff = -cutoff;
   }

/***********  end of big loop ***********/


   free_plane(A);
   free_plane(B);
   free_plane(C);
   
   if(gradientCube) {
      free(gradientCube[0][0]);
      free(gradientCube[0]);
      free(gradientCube);
   }

   glPopMatrix();
   glPopMatrix();
   glDrawBuffer(GL_BACK);

   if(s1){
      logprint("");
      sprintf(str, "surface at %f:", cutoff);
      logprint(str);
      sprintf(str, "%d points, %d triangles", s1->npts, s1->ntri);
      logprint(str);
   }
   if(s2){
      logprint("");
      sprintf(str, "surface at %f:", -cutoff);
      logprint(str);
      sprintf(str, "%d points, %d triangles", s2->npts, s2->ntri);
      logprint(str);
   }

   approx_volumes();
   exact_volumes(s1);
   
   glutSetWindow(mainwin);
   glutPostRedisplay();
}
/**** end of cubes ****/



int alloc_tridot(void)
{
   if((dot = (Surfdot*) malloc(MEMBLOCK * sizeof(Surfdot))) == NULL){
      showinfobox("can't allocate enough memory (dot)");
      return 0;
   }
   dotmem = MEMBLOCK;
   if((tri = (Triangle*) malloc(2 * MEMBLOCK * sizeof(Triangle))) == NULL){
      showinfobox("can't allocate enough memory (tri)");
      return 0;
   }
   trimem = 2 * MEMBLOCK;

   return 1;
}



void approx_volumes(void)
{
   register int i;
   int inside, outside, npts;
   register float *p;
   float dv, invol, outvol, inpercent, outpercent;
   char str[60];

   inside = outside = 0;
   npts = nx*ny*nz;
   dv   = dx*dy*dz;
   for(i=0, p=actualmol->cube_value[0][0]; i<npts; i++, p++){
      if(*p >= cutoff) inside++;
      if(*p < -cutoff) outside++;
   }
   invol  = dv*inside;
   inpercent = 100.0*inside/npts;
   logprint("");
   logprint("approx. volume above cutoff:");
   sprintf(str, "%.3f (= %.2f%%)", invol, inpercent);
   logprint(str);

   if(bit.both_signs){
      outvol  = dv*outside;
      outpercent = 100.0*outside/npts;
      logprint("approx. volume below");
      logprint("neg. cutoff:");
      sprintf(str, "%.3f (= %.2f%%)", outvol, outpercent);
      logprint(str);
   }
}


void exact_volumes(Surface *sp)
{
   double vol1, vol2;
   char str[120];

   if(!sp) return;

   vol1 = compute_volume(sp);

   if(!sp->second) {
      sprintf(str, "Exact volume = %f\n", vol1);
      logprint("");
      logprint(str);
      update_logs();
   }
   else {
      vol2 = compute_volume(sp->second);
      logprint("");
      logprint("Exact volume of");
      sprintf(str, "first surface = %f\n", vol1);
      logprint(str);
      logprint("Exact volume of");
      sprintf(str, "second surface = %f\n", vol2);
      logprint(str);
      update_logs();
   }
}




int init_cubes(void)
{
   static char *errmsg[] =  { "Marching Cubes :",
                              "load a .macu-file first",
                              "cutoff below minimal value!",
                              "cutoff above maximal value!",
                              "can't allocate enough memory for plane A\n",
                              "can't allocate enough memory for plane B\n",
                              "can't allocate enough memory for plane C\n" 
                            };
   char str[100];

   if(!actualmol->got_macufile){
      sprintf(str, "%s\n%s", errmsg[0], errmsg[1]); 
      showinfobox(str);
      return 1;
   }
   actualmol->cutoff = cutoff;
   if(cutoff < actualmol->cubemin){
      sprintf(str, "%s\n%s", errmsg[0], errmsg[2]); 
      showinfobox(str);
      return 2;
   }
   if(cutoff >\
         ((bit.both_signs && (-1*actualmol->cubemin) > actualmol->cubemax) ?\
          (-1*actualmol->cubemin) : actualmol->cubemax)){
      sprintf(str, "%s\n%s", errmsg[0], errmsg[3]); 
      showinfobox(str);
      return 3;
   }
   if((A = plane_alloc(nx, ny)) == NULL){
      sprintf(str, "%s\n%s", errmsg[0], errmsg[4]); 
      showinfobox(str);
      return 4;
   }
   if((B = plane_alloc(nx, ny)) == NULL){
      sprintf(str, "%s\n%s", errmsg[0], errmsg[5]); 
      showinfobox(str);
      return 5;
   }
   if((C = plane_alloc(nx, ny)) == NULL){
      sprintf(str, "%s\n%s", errmsg[0], errmsg[6]); 
      showinfobox(str);
      return 6;
   }

   return 0;
}


/*
 *   shift the three actual planes one plane upwards in the box
 */
void shift_planes(void)
{
   register short ix, iy;
   Plane *temp;
   Point *pt;

   temp = A;
   A = B;
   B = C;
   C = temp;

   for(iy=0; iy<ny+1; iy++){
      for(ix=0; ix<nx+1; ix++) ccase[iy][ix] >>= 4;
   }

   for(iy=0; iy<ny; iy++){
      for(ix=0, pt = C->point[iy]; ix<nx; ix++, pt++) {
         pt->weight = pt->index = pt->ix = pt->iy = pt->iz = 0;
      }
   }

   C->value = actualmol->cube_value[iz+1];
   B->value = actualmol->cube_value[iz];
   if(iz) A->value = actualmol->cube_value[iz-1];

   za = zb;
   zb += dz;
   rectangle();
}

void rectangle(void)
{
   float corner[4][3];
   register short i;

   corner[0][0] = xmin; corner[0][1] = ymin; corner[0][2] = zb;
   corner[1][0] = xmax; corner[1][1] = ymin; corner[1][2] = zb;
   corner[2][0] = xmax; corner[2][1] = ymax; corner[2][2] = zb;
   corner[3][0] = xmin; corner[3][1] = ymax; corner[3][2] = zb;

   glDisable(GL_LIGHTING);
   glColor4fv(rectcol);
   glBegin(GL_LINE_LOOP);
   for(i=0; i<4; i++) glVertex3fv(corner[i]);
   glEnd();
   glEnable(GL_LIGHTING);
}


/*
 *    calculate the normals
 *    as the normalized gradient vectors at each grid-point
 */
Vector ***do_normals(void)
{
   register short ix, iy, iz;
   Vector ***gc;
   float ***v;
   
   v = actualmol->cube_value;
   
   if((gc = (Vector***)alloc_3D(nx, ny, nz, sizeof(Vector))) == NULL) {
       fprintf(stderr, "can't allocate cube of gradients\n");
       return 0;
   }

   for(iz=0; iz<nz; iz++) {
     for(iy=0; iy<ny; iy++){
       for(ix=0; ix<nx; ix++){
         if(ix == 0)           gc[iz][iy][ix].x = (v[iz][iy][ix]   - v[iz][iy][ix+1]) * 2. * invdx;
         else if(ix == (nx-1)) gc[iz][iy][ix].x = (v[iz][iy][ix-1] - v[iz][iy][ix])   * 2. * invdx;
         else                  gc[iz][iy][ix].x = (v[iz][iy][ix-1] - v[iz][iy][ix+1]) * invdx;

         if(iy == 0)           gc[iz][iy][ix].y = (v[iz][iy][ix]   - v[iz][iy+1][ix]) * 2. * invdy;
         else if(iy == (ny-1)) gc[iz][iy][ix].y = (v[iz][iy-1][ix] - v[iz][iy][ix])   * 2. * invdy;
         else                  gc[iz][iy][ix].y = (v[iz][iy-1][ix] - v[iz][iy+1][ix]) * invdy;

         if(iz == 0)           gc[iz][iy][ix].z = (v[iz][iy][ix]   - v[iz+1][iy][ix])  * 2. * invdz;
         else if(iz == (nz-1)) gc[iz][iy][ix].z = (v[iz-1][iy][ix] - v[iz][iy][ix])    * 2. * invdz;
         else                  gc[iz][iy][ix].z = (v[iz-1][iy][ix] - v[iz+1][iy][ix])  * invdz;

         normalize(&gc[iz][iy][ix]);
       }
     }
   }
   return gc;
}



int addPoint(Point *p, float *coord,  int *np)
{
    if(!p->weight) {
	p->index = *np;
	*np += 1;
	dot[p->index].v[0] = coord[0];
	dot[p->index].v[1] = coord[1];
	dot[p->index].v[2] = coord[2];
	dot[p->index].n[0] = 0.0;
	dot[p->index].n[1] = 0.0;
	dot[p->index].n[2] = 1.0;
	p->weight = 1;
    }
    else {
	dot[p->index].v[0] = ((float)p->weight * dot[p->index].v[0] + coord[0])/(float)(p->weight + 1);
	dot[p->index].v[1] = ((float)p->weight * dot[p->index].v[1] + coord[1])/(float)(p->weight + 1);
	dot[p->index].v[2] = ((float)p->weight * dot[p->index].v[2] + coord[2])/(float)(p->weight + 1);
	p->weight += 1;
    }
    
    return p->index;
}


/*
 *  determine the surface-points
 *     as the intersection of the contour with the grid
 *     using linear interpolation for position and for the surface-normal
 */
void do_points(void)
{
   register short ix, iy;
   float x, y, value;
   float coord[3], difference1, difference2;
   int index;

   for(iy=0; iy<ny; iy++){
      y = ymin + iy*dy;
      for(ix=0; ix<nx; ix++){
         x = xmin + ix*dx;
	 value =  B->value[iy][ix];
         if(ix<(nx-1)){
            if((value >= cutoff) ^ (B->value[iy][ix+1] >= cutoff)) {
               coord[0] = interpolate(x, value, x+dx, B->value[iy][ix+1]);
               coord[1] = y;
               coord[2] = zb;
	       
	       difference1 = fabsf(value - cutoff);
	       difference2 = fabsf(B->value[iy][ix+1] - cutoff);
	       if(difference1 < difference2) {
		   index = addPoint(&B->point[iy][ix], coord, &npts);
	       }
	       else {
		   index = addPoint(&B->point[iy][ix+1], coord, &npts);
	       }
	       B->point[iy][ix].ix = index;
            }
         }
         if(iy<(ny-1)){
            if((value >= cutoff) ^ (B->value[iy+1][ix] >= cutoff)) {
               coord[0] = x;
               coord[1] = interpolate(y, value, y+dy, B->value[iy+1][ix]);
               coord[2] = zb;
	       
	       difference1 = fabsf(value - cutoff);
	       difference2 = fabsf(B->value[iy+1][ix] - cutoff);
	       if(difference1 < difference2) {
		   index = addPoint(&B->point[iy][ix], coord, &npts);
	       }
	       else {
		   index = addPoint(&B->point[iy+1][ix], coord, &npts);
	       }
	       B->point[iy][ix].iy = index;
            }
         }
         if (iz){
            if((value >= cutoff) ^ (A->value[iy][ix] >= cutoff)) {
               coord[0] = x;
               coord[1] = y;
               coord[2] = interpolate(zb, value, za, A->value[iy][ix]);
	       difference1 = fabsf(value - cutoff);
	       difference2 = fabsf(A->value[iy][ix] - cutoff);
	       if(difference1 < difference2) {
		   index = addPoint(&B->point[iy][ix], coord, &npts);
	       }
	       else {
		   index = addPoint(&A->point[iy][ix], coord, &npts);
	       }
	       B->point[iy][ix].iz = index;
            }
         }
	 
         if(6*(npts + 3) > dotmem){
            if((dot  = (Surfdot*) realloc(dot,
               (dotmem + MEMBLOCK) * sizeof(float))) == NULL){
               fprintf(stderr, " Can't reallocate enough memory ");
               return;
            }
            dotmem += MEMBLOCK;
         }
      }
   }
}



/*
 *  Determine the high-order bits of the cube-label,
 *  then triangulate the cubes
 */
void do_cubes(void)
{
   register short ix, iy;
   unsigned int e[12];
   
   for(iy=0; iy<ny; iy++){
      for(ix=0; ix<nx; ix++){
         if(B->value[iy][ix] >= cutoff){
            ccase[iy+1][ix+1] |=  16;
            ccase[iy+1][ix]   |=  32;
            ccase[iy][ix+1]   |= 128;
            ccase[iy][ix]     |=  64;
         }
      }
   }

   if(!iz) return;

   for(iy=1; iy<ny; iy++){
      for(ix=1; ix<nx; ix++){
         if((ccase[iy][ix]) && (ccase[iy][ix] != 255)){
            e[0]  = A->point[iy-1][ix-1].ix;
            e[1]  = A->point[iy-1][ix].iy;
            e[2]  = A->point[iy][ix-1].ix;
            e[3]  = A->point[iy-1][ix-1].iy;
            e[4]  = B->point[iy-1][ix-1].ix; 
            e[5]  = B->point[iy-1][ix].iy;
            e[6]  = B->point[iy][ix-1].ix;
            e[7]  = B->point[iy-1][ix-1].iy;
            e[8]  = B->point[iy-1][ix-1].iz;
            e[9]  = B->point[iy-1][ix].iz;
            e[10] = B->point[iy][ix].iz;
            e[11] = B->point[iy][ix-1].iz;

	    triangulate(ccase[iy][ix], e);
         }
      }
   }
}



/*
 *  linear inerpolation of the grid-points 1 and 2
 *  (coordinate c, vector a, value v) to get the normal at the surface-point
 *  returns the coordinate of the surface-point
 */
float interpolate(float c1, float v1,
                  float c2, float v2)
{
   float r;

   r = (cutoff-v1)/(v2-v1);
/*
   dot[npts].n[0] = a1.x + (a2.x - a1.x)*r;
   dot[npts].n[1] = a1.y + (a2.y - a1.y)*r;
   dot[npts].n[2] = a1.z + (a2.z - a1.z)*r;
*/
   return(c1 + (c2-c1)*r);
}


void normalize(Vector *a)
{
   float length, inv;

   length = sqrt(a->x*a->x + a->y*a->y + a->z*a->z);
   if(length == 0) return;
   inv = 1.0/length;
   if(cutoff < 0) inv *= -1.0;
   a->x *= inv;
   a->y *= inv;
   a->z *= inv;
}



void triangulate(unsigned char lbl, unsigned int *e)
{
   register int i;
   int pt[3];
   tri_pattern *tp;

   tp = pattern[lbl];

   for (i=0; i<npat[lbl]; i++, tp++){
      pt[0] = e[tp->edge0];
      pt[1] = e[tp->edge1];
      pt[2] = e[tp->edge2];
      
      if((pt[0] != pt[1]) && (pt[0] != pt[2]) && (pt[1] != pt[2])) {

         tri[ntri].p1 = pt[0];
         tri[ntri].p2 = pt[1];
         tri[ntri].p3 = pt[2];

         ntri++;
      }
   }

   if((ntri+6) > trimem){
      if((tri = (Triangle*) realloc(tri, (trimem + 2*MEMBLOCK)*sizeof(Triangle))) == NULL){
         fprintf(stderr, " Can't reallocate enough memory ");
         return;
      }
      trimem += 2 * MEMBLOCK;
   }
}



void surfaceNormals(Vector ***cu)
{
    float x, y, z, px, py, pz;
    int i, ix, iy, iz;
    Surfdot *pt;
    Vector norm;

    invdx = (nx-1.)/(xmax-xmin);
    invdy = (ny-1.)/(ymax-ymin);
    invdz = (nz-1.)/(zmax-zmin);

    for(i=0, pt=dot; i<npts; i++, pt++) {
	int jx, jy, jz, jx1, jy1, jz1;
    
	x = pt->v[0]; y = pt->v[1]; z = pt->v[2];

	ix = (int)floorf((x - xmin)*invdx);
	iy = (int)floorf((y - ymin)*invdy);
	iz = (int)floorf((z - zmin)*invdz);

	if(ix < 0 || ix >= nx) continue;
	if(iy < 0 || iy >= ny) continue;
	if(iz < 0 || iz >= nz) continue; /* outside the cube */

	px = (x-xmin)*invdx - ix; /* ratio */
	py = (y-ymin)*invdy - iy; /* ratio */
	pz = (z-zmin)*invdz - iz; /* ratio */

	jx = ix; jx1 = ix<nx-1 ? ix+1 : ix;
	jy = iy; jy1 = iy<ny-1 ? iy+1 : iy;
	jz = iz; jz1 = iz<nz-1 ? iz+1 : iz;

	cubeGradient(px, py, pz, &norm, 
	    cu[jz][jy]+jx, cu[jz][jy]+jx1, cu[jz][jy1]+jx1, cu[jz][jy1]+jx, 
	    cu[jz1][jy]+jx, cu[jz1][jy]+jx1, cu[jz1][jy1]+jx1, cu[jz1][jy1]+jx);

	pt->n[0] = norm.x;
	pt->n[1] = norm.y;
	pt->n[2] = norm.z;
    }
}

/****

    	     
	     P011________________P111
	       /|		/|
	      /	|	       / |
	  p001________________/p101
	     |  |	      |  |
	     |  |q01	      |  |q11
	r0-> | /________X_____|_/| <-r1
	  q00|/	|  	   q10./ |
	     |  | 	      |  |
	     |  |p010_________|__|p100
	     | /	      | /
	     |/_______________|/
	 p000		       p100

****/



/*
static float cubeInterpolGradient(float dx, float dy, float dz, Vector *grad, 
       float p000, float p100, float p110, float p010,
       float p001, float p101, float p111, float p011)
{
   float q00, q01, q10, q11, r0, r1, value, term1, term2;

   q00 = p000 + (p001-p000)*dz;
   q10 = p100 + (p101-p100)*dz;
   q01 = p010 + (p011-p010)*dz;
   q11 = p110 + (p111-p110)*dz;

   r0 = q00 + (q01-q00)*dy;
   r1 = q10 + (q11-q10)*dy;

   value = r0 + (r1-r0)*dx;
   
    grad->x = -(r1-r0)*invdx;
    grad->y = -((q01-q00)*(1.0-dx) + (q11-q10)*dx) * invdy;

    term1 =  (p001-p000)*(1.0-dy) + (p011-p010)*dy;
    term2 =  (p101-p100)*(1.0-dy) + (p111-p110)*dy;
    
    grad->z = -(term1*(1.0-dx) + term2*dx)*invdz;

    normalize(grad);

   return value;
}
*/





static void cubeGradient(float dx, float dy, float dz, Vector *grad, 
       Vector *p000, Vector *p100, Vector *p110, Vector *p010,
       Vector *p001, Vector *p101, Vector *p111, Vector *p011)
{
    Vector q00, q01, q10, q11, r0, r1;

    q00.x = p000->x + (p001->x-p000->x)*dz;
    q10.x = p100->x + (p101->x-p100->x)*dz;
    q01.x = p010->x + (p011->x-p010->x)*dz;
    q11.x = p110->x + (p111->x-p110->x)*dz;
    q00.y = p000->y + (p001->y-p000->y)*dz;
    q10.y = p100->y + (p101->y-p100->y)*dz;
    q01.y = p010->y + (p011->y-p010->y)*dz;
    q11.y = p110->y + (p111->y-p110->y)*dz;
    q00.z = p000->z + (p001->z-p000->z)*dz;
    q10.z = p100->z + (p101->z-p100->z)*dz;
    q01.z = p010->z + (p011->z-p010->z)*dz;
    q11.z = p110->z + (p111->z-p110->z)*dz;

    r0.x = q00.x + (q01.x - q00.x)*dy;
    r1.x = q10.x + (q11.x - q10.x)*dy;
    r0.y = q00.y + (q01.y - q00.y)*dy;
    r1.y = q10.y + (q11.y - q10.y)*dy;
    r0.z = q00.z + (q01.z - q00.z)*dy;
    r1.z = q10.z + (q11.z - q10.z)*dy;

    grad->x = r0.x + (r1.x - r0.x)*dx;
    grad->y = r0.y + (r1.y - r0.y)*dx;
    grad->z = r0.z + (r1.z - r0.z)*dx;
   
    normalize(grad);
}

int add_cutplane_tri(Mol *mp, unsigned int tri[3][2])
{
   if(!mp->plane->tri) {
      if((mp->plane->tri = (unsigned int(*)[3][2])malloc(sizeof(unsigned int[3][2]))) == NULL) {
         showinfobox("Can't allocate memory for plane");
         return 0;
      }
   } else
   {
      if((mp->plane->tri = 
         (unsigned int(*)[3][2])realloc(mp->plane->tri, (mp->plane->ntri+1)*sizeof(unsigned int[3][2]))) == NULL) {
         showinfobox("Can't allocate memory for plane");
         return 0;
      }
   }
   mp->plane->tri[mp->plane->ntri][0][0] = tri[0][0];
   mp->plane->tri[mp->plane->ntri][0][1] = tri[0][1];
   mp->plane->tri[mp->plane->ntri][1][0] = tri[1][0];
   mp->plane->tri[mp->plane->ntri][1][1] = tri[1][1];
   mp->plane->tri[mp->plane->ntri][2][0] = tri[2][0];
   mp->plane->tri[mp->plane->ntri][2][1] = tri[2][1];
//printf("ntri: %d\n", mp->plane->ntri);
//printf("%d %d %d %d %d %d\n", mp->plane->tri[mp->plane->ntri][0][0], mp->plane->tri[mp->plane->ntri][0][1], mp->plane->tri[mp->plane->ntri][1][0], mp->plane->tri[mp->plane->ntri][1][1], mp->plane->tri[mp->plane->ntri][2][0], mp->plane->tri[mp->plane->ntri][2][1]);
   mp->plane->ntri++;
   return 1;
}

void tri_cutplane(Mol *mp)
{
   int size, jj=0;
   int register i, j;
   unsigned int tri[3][2];

   size = mp->plane->npts;
   for(j=0; j<size-1; j++) {
//printf("j: %d\n", j);
       i = 0;
      jj = 0;
      while(i<size) {
//printf("i: %d\n", i);
         if(mp->plane->plane_point[j][i][6]) {
            tri[0][0]=j;
            tri[0][1]=i;
            i++;
            if(i<size && mp->plane->plane_point[j][i][6]) {
               tri[1][0]=j;
               tri[1][1]=i;
               while(jj < size) {
//printf("jj1: %d\n", jj);
                  if(mp->plane->plane_point[j+1][jj][6]) {
                     tri[2][0]=j+1;
                     tri[2][1]=jj;
// add trangle
//printf("position 1\n");
                     if(!add_cutplane_tri(mp, tri)) return;
                     break;
                  } else {
                     jj++;
                  }
               }
            }
            else {
               if(mp->plane->ntri == 0) i--;
            }
            while (jj<i) {
//printf("jj2: %d\n", jj);
               if(!mp->plane->plane_point[j+1][jj][6]) {
                  jj++;
                  continue;
               }
               tri[0][0]=j;
               tri[0][1]=i;
               tri[1][0]=j+1;
               tri[1][1]=jj;
               jj++;
               if(jj < size && mp->plane->plane_point[j+1][jj][6]) {
                  tri[2][0]=j+1;
                  tri[2][1]=jj;
//printf("position 2\n");
                  if(!add_cutplane_tri(mp, tri)) return;
               }
            }
            if(i+1<size && !mp->plane->plane_point[j][i+1][6]) {
               jj++;
               while (jj<size) {
//printf("jj3: %d\n", jj);
                  if(!mp->plane->plane_point[j+1][jj][6]) break;
                  tri[0][0]=j;
                  tri[0][1]=i;
                  tri[1][0]=j+1;
                  tri[1][1]=jj-1;
                  tri[2][0]=j+1;
                  tri[2][1]=jj;
//printf("position 3\n");
                  if(!add_cutplane_tri(mp, tri)) return;
                  jj++;
               }
            }
            if(mp->plane->ntri != 0) i--;
         }
         i++;
      }
   }
}

float interpolate_cutplane_value(Mol *mp, float p[3])
{
   int i, ix, iy, iz;
   float ddx, ddy, ddz, x, y, z;
   float p000, p100, p110, p010, p001, p101, p111, p011;
   float t[8];
   float v = 0;

   ix = (int)ffloor((p[0] - xmin) / dx);
   iy = (int)ffloor((p[1] - ymin) / dy);
   iz = (int)ffloor((p[2] - zmin) / dz);
   ddx = (p[0] - xmin) - dx*ix;
   ddy = (p[1] - ymin) - dy*iy;
   ddz = (p[2] - zmin) - dz*iz;
   ix = (ix<nx) ? ix : ix-1;
   iy = (iy<ny) ? iy : iy-1;
   iz = (iz<nz) ? iz : iz-1;
   x = ddx/dx;
   y = ddy/dy;
   z = ddz/dz;

   p000 = mp->cube_value[iz][iy][ix];
   p100 = mp->cube_value[iz][iy][ix+1];
   p010 = mp->cube_value[iz][iy+1][ix];
   p001 = mp->cube_value[iz+1][iy][ix];
   p110 = mp->cube_value[iz][iy+1][ix+1];
   p101 = mp->cube_value[iz+1][iy][ix+1];
   p011 = mp->cube_value[iz+1][iy+1][ix];
   p111 = mp->cube_value[iz+1][iy+1][ix+1];

   t[0] = p000*(1-x)*(1-y)*(1-z);
   t[1] = p100*x*(1-y)*(1-z);
   t[2] = p010*(1-x)*y*(1-z);
   t[3] = p001*(1-x)*(1-y)*z;
   t[4] = p101*x*(1-y)*z;
   t[5] = p011*(1-x)*y*z;
   t[6] = p110*x*y*(1-z);
   t[7] = p111*x*y*z;

   for(i=0; i<8; i++) {
      v += t[i];
   }

   return v;
}

/*
my idea did not work out!
float interpolate_cutplane_value(Mol *mp, float p[3])
{
   int ix, iy, iz;
   float ddx, ddy, ddz;
   float p000, p100, p110, p010, p001, p101, p111, p011;
   float q005, q050, q500, q105, q150, q015, q051, q505, q511, q151, q115, q550;
   float exyf, exyb, exzf, exzb, eyzf, eyzb;
   float v = 0;



   ix = (int)((p[0] - xmin) / dx);
   iy = (int)((p[1] - ymin) / dy);
   iz = (int)((p[2] - zmin) / dz);
   ddx = (p[0] - xmin) - dx*ix;
   ddy = (p[1] - ymin) - dy*iy;
   ddz = (p[2] - zmin) - dz*iz;
   ix = (ix<nx) ? ix : ix-1;
   iy = (iy<ny) ? iy : iy-1;
   iz = (iz<nz) ? iz : iz-1;

   p000 = mp->cube_value[iz][iy][ix];
   p100 = mp->cube_value[iz][iy][ix+1];
   p010 = mp->cube_value[iz][iy+1][ix];
   p001 = mp->cube_value[iz+1][iy][ix];
   p110 = mp->cube_value[iz][iy+1][ix+1];
   p101 = mp->cube_value[iz+1][iy][ix+1];
   p011 = mp->cube_value[iz+1][iy+1][ix];
   p111 = mp->cube_value[iz+1][iy+1][ix+1];
   
   q005 = p000 + (p001-p000)*ddz;
   q050 = p000 + (p010-p000)*ddy;
   q500 = p000 + (p100-p000)*ddx;
   q105 = p100 + (p101-p100)*ddz;
   q150 = p100 + (p110-p100)*ddy;
   q015 = p010 + (p011-p010)*ddz;
   q051 = p001 + (p011-p001)*ddy;
   q505 = p001 + (p101-p001)*ddx;
   q511 = p011 + (p111-p011)*ddx;
   q151 = p101 + (p111-p101)*ddy;
   q115 = p110 + (p111-p110)*ddz;
   q550 = p010 + (p110-p010)*ddx;

   exyf = (dx*(q505+(q511-q505)*ddy)+dy*(q051+(q151-q051)*ddx))/(dx+dy);
   exyb = (dx*(q500+(q550-q500)*ddy)+dy*(q050+(q150-q050)*ddx))/(dx+dy);
   exzf = (dz*(q005+(q105-q005)*ddx)+dx*(q500+(q505-q500)*ddz))/(dx+dz);
   exzb = (dz*(q015+(q115-q015)*ddx)+dx*(q550+(q511-q550)*ddz))/(dx+dz);
   eyzf = (dy*(q150+(q151-q150)*ddz)+dz*(q105+(q115-q105)*ddy))/(dy+dz);
   eyzb = (dy*(q050+(q051-q050)*ddz)+dz*(q005+(q051-q005)*ddy))/(dy+dz);

   v = ((dx+dy)*(exyb+(exyf-exyb)*ddz)+(dx+dz)*(exzf+(exzb-exzf)*ddy)+(dy+dz)*(eyzb+(eyzf-eyzb)*ddz))
         /(dx+dy+dz);

   return v;
}
*/

void set_cutplane(Mol *mp)
{
   float m1[3][3], m2[3][3], n1[3], n2[3], n3[3];
   float alpha, beta;
   int register i, j;

   if(mp->plane->ntri) {
      free(mp->plane->tri);
      mp->plane->tri = NULL;
      mp->plane->ntri = 0;
   }
   xmin = mp->box.x1;
   xmax = mp->box.x2;
   ymin = mp->box.y1;
   ymax = mp->box.y2;
   zmin = mp->box.z1;
   zmax = mp->box.z2;
   nx =  mp->box.nx;
   ny =  mp->box.ny;
   nz =  mp->box.nz;

   dx = (xmax-xmin)/((float)nx-1.);
   dy = (ymax-ymin)/((float)ny-1.);
   dz = (zmax-zmin)/((float)nz-1.);

   alpha = mp->plane->alpha/180*M_PI;  
   beta = mp->plane->beta/180*M_PI;
   n1[0] = n1[1] = 0;
   n1[2] = 1;
   m1[0][0] = 1;
   m1[1][0] = 0;
   m1[2][0] = 0;
   m1[0][1] = 0;
   m1[1][1] = cos(beta);
   m1[2][1] = sin(beta);
   m1[0][2] = 0;
   m1[1][2] = -sin(beta);
   m1[2][2] = cos(beta);
   m2[0][0] = cos(alpha);
   m2[1][0] = 0;
   m2[2][0] = sin(alpha);
   m2[0][1] = 0;
   m2[1][1] = 1;
   m2[2][1] = 0;
   m2[0][2] = -sin(alpha);
   m2[1][2] = 0;
   m2[2][2] = cos(alpha);

   for(j=0; j<mp->plane->npts; j++) {
      for(i=0; i<mp->plane->npts; i++) {
      n1[0] = mp->plane->plane_point[j][i][0];
      n1[1] = mp->plane->plane_point[j][i][1];
      n1[2] = 0;
      vec3mat(n1, m1, n2);
      vec3mat(n2, m2, n3);
      n3[0] = n3[0] + mp->plane->a[0];
      n3[1] = n3[1] + mp->plane->a[1];
      n3[2] = n3[2] + mp->plane->a[2];
      mp->plane->plane_point[j][i][6] = 0;
      if(n3[0] >= mp->box.x1 && n3[0] <= mp->box.x2) {
         if(n3[1] >= mp->box.y1 && n3[1] <= mp->box.y2) {
            if(n3[2] >= mp->box.z1 && n3[2] <= mp->box.z2) {
               mp->plane->plane_point[j][i][2] = n3[0];
               mp->plane->plane_point[j][i][3] = n3[1];
               mp->plane->plane_point[j][i][4] = n3[2];
               mp->plane->plane_point[j][i][5] = interpolate_cutplane_value(mp, n3);
               mp->plane->plane_point[j][i][6] = 1;
            }
         }
      }
      }
   }
   tri_cutplane(mp);
}

void init_cutplane(Mol *mp)
{
   int size, k, i, j;
   float zerox, zeroy, dd;

   if((mp->plane = (Cutplane *)malloc(sizeof(Cutplane))) ==NULL) {
      showinfobox("Can't allocate memory for plane");
      return;
   }
   dx = (xmax-xmin)/((float)nx-1.);
   dy = (ymax-ymin)/((float)ny-1.);
   dz = (zmax-zmin)/((float)nz-1.);
   mp->plane->a[0] = xmin + dx*(nx-1)*0.5;
   mp->plane->a[1] = ymin + dy*(ny-1)*0.5;
   mp->plane->a[2] = zmin + dz*(nz-1)*0.5;
   mp->plane->alpha = 0.0;
   mp->plane->beta = 0.0;
   dd = 0.5*sqrt(dx*dx + dy*dy + dz*dz);
   mp->plane->npts = 0;
   mp->plane->ntri = 0;
   mp->plane->tri = NULL;
   size = (int)(sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)+(zmax-zmin)*(zmax-zmin))/dd + 1);
   mp->plane->npts = size;
   if((mp->plane->plane_point = (float(**)[7])alloc_3D(size, size, 0, sizeof(float[7]))) ==NULL) {
      showinfobox("Can't allocate memory for plane");
      return;
   }
   zerox = mp->plane->a[0] - size/2 * dd;
   zeroy = mp->plane->a[1] - size/2 * dd;
   for(j=0; j<size; j++) {  // fill row per row
      for(i=0; i<size; i++) {
         mp->plane->plane_point[j][i][0] = zerox + i*dd;
         mp->plane->plane_point[j][i][1] = zeroy + j*dd;
      }
   }
   set_cutplane(mp);
}

void readmacu(char *fname, int key)
{
   FILE *fp;
   Box *bp;
   int len, dum;
   float *p, *q;
   register int h, i;
   static char str[50], *s;
   float *dummy_plane, sum, dv = 0;

   if(!actualmol){
      logprint("readmacu : no molecule loaded");
      update_logs();
      return;
   }

   if((fp = fopen(fname, "rb")) == NULL){
       sprintf(str, "can't open %s \n", fname);
       showinfobox(str);
       return;
   }

   fread(&dum,  sizeof(int),   1, fp);
   switch(dum){
      case 36      : //sprintf(str, ".macu-file : %s \n",fname);
                     //logprint(str);
                     break;
      default      : fclose(fp);
                     showinfobox("wrong magic number!");
                     return;
   }

   fread(&cubehead, sizeof(MacuHeader), 1, fp);
   bp = &actualmol->box;
   if(key == ADD_MACU || key == SUBST_MACU){
      if((nx != cubehead.nx) || (ny != cubehead.ny) || (nz != cubehead.nz)){
         showinfobox("The two macu-files don't correspond!");
         fclose(fp);
         return;
      }
   }
   xmin = bp->x1 = cubehead.xmin;
   xmax = bp->x2 = cubehead.xmax;
   ymin = bp->y1 = cubehead.ymin;
   ymax = bp->y2 = cubehead.ymax;
   zmin = bp->z1 = cubehead.zmin;
   zmax = bp->z2 = cubehead.zmax;
   nx   = bp->nx = cubehead.nx;
   ny   = bp->ny = cubehead.ny;
   nz   = bp->nz = cubehead.nz;
   nMacuPlanes = nz;

   logprint("");
   if((s = strrchr(fname, '/')) == NULL) s = fname;
   else s++;
   logprint(s);
   sprintf(str, "xmin =  ymin =  zmin =");
   logprint(str);
   sprintf(str, "%6.2f  %6.2f  %6.2f", xmin, ymin, zmin);
   logprint(str);
   sprintf(str, "xmax =  ymax =  zmax =");
   logprint(str);
   sprintf(str, "%6.2f  %6.2f  %6.2f", xmax, ymax, zmax);
   logprint(str);
   sprintf(str, "nx = %3d  ny = %3d  nz = %3d", nx, ny, nz); 
   logprint(str);

   strncpy(macu_name, s, 29);

   switch(key){
      case LOAD_MACU :
         free_3D(actualmol);
         if((actualmol->cube_value = (float***) alloc_3D(nx, ny, nz, sizeof(float))) == NULL){
            showinfobox(" can't allocate enough memory for the cube-values\n");
            return;
         }
         if((ccase = (unsigned char**) alloc_3D(nx+1, ny+1, 0, sizeof(char))) == NULL){
            showinfobox(" can't allocate enough memory for the cube-labels\n");
            return;
         }
            /* the cube-label array is (nx+1) * (ny+1),    *
             * so there is no need to test for the borders */
         for(i=0; i<nz; i++){
            fread(&len, sizeof(int), 1, fp);              /* FORTRAN - file */
            fread(actualmol->cube_value[i][0], sizeof(float), nx*ny, fp);
            fread(&len, sizeof(int), 1, fp);              /* FORTRAN - file */
         }
         break;
      case ADD_MACU :
         if((dummy_plane = (float*) malloc(nx*ny*sizeof(float))) == NULL){
            showinfobox(" can't allocate enough memory for the dummy plane\n");
            return;
         }
         for(i=0; i<nz; i++){
            fread(&len, sizeof(int), 1, fp);              /* FORTRAN - file */
            fread(dummy_plane, sizeof(float), nx*ny, fp);
            fread(&len, sizeof(int), 1, fp);              /* FORTRAN - file */
            for(p=dummy_plane, q=actualmol->cube_value[i][0], h=0; h<nx*ny; p++,q++,h++){
               *q += *p;
            }
         }
         free(dummy_plane);
         break;
      case SUBST_MACU :
         if((dummy_plane = (float*) malloc(nx*ny*sizeof(float))) == NULL){
            showinfobox(" can't allocate enough memory for the dummy plane\n");
            return;
         }
         for(i=0; i<nz; i++){
            fread(&len, sizeof(int), 1, fp);              /* FORTRAN - file */
            fread(dummy_plane, sizeof(float), nx*ny, fp);
            fread(&len, sizeof(int), 1, fp);              /* FORTRAN - file */
            for(p=dummy_plane,q=actualmol->cube_value[i][0],h=0; h<nx*ny; p++,q++,h++){
               *q -= *p;
            }
         }
         free(dummy_plane);
         break;
   }
   fclose(fp);
   p = actualmol->cube_value[0][0];
   cubemin = cubemax = *p;
   sum = 0;
   for(i=1; i < (nx*ny*nz); i++, p++){
      if(*p > cubemax) cubemax = *p;
      if(*p < cubemin) cubemin = *p;
      sum += *p;
   }
   if(cubemin < 0.000001 && cubemin > -0.000001) cubemin = 0.0;
   if(cubemax < 0.000001 && cubemax > -0.000001) cubemax = 0.0;
   actualmol->cubemin = cubemin;
   actualmol->cubemax = cubemax;
   actualmol->got_macufile = 1;
   actualmol->cutoff = cutoff = 0.00;
   actualmol->plane_iz = cubehead.nz/2;
   actualmol->plane_iy = cubehead.ny/2;
   actualmol->plane_ix = cubehead.nx/2;
   init_cutplane(actualmol);
   dv = (xmax-xmin) * _1_BOHR / (nx-1.) *
        (ymax-ymin) * _1_BOHR / (ny-1.) *
        (zmax-zmin) * _1_BOHR / (nz-1.);
   sum *= dv;
   logprint("sum of voxel contributions");
   sprintf(str, " (a.u.) = %f", sum);
   logprint(str);
   update_logs();

   if(surfaceglui) {
      surfaceglui->sync_live();
   }
   if(planeglui) {
      planeglui->sync_live();
   }

}




Surface *add_surface(char type, Surfdot  *dot, Triangle *tri,
   float contour, float *val, int npts, int ntri, char color, char *name )
{
   Surface    *temp, *last;

   if((temp = (Surface*) malloc(sizeof(Surface))) == NULL){
      showinfobox("Can't allocate enough memory for the surface\n");
      exit(-1);
   }

   if(actualmol->firstsurf) {
      last = actualmol->firstsurf;
      while(last->next) last = last->next;
      last->next = temp;
   }
   else {
      actualmol->firstsurf = temp;
   }
   actualsurf = temp;

    {
/*
	Triangle *tp;
    
	for(i=0, tp = tri; i<ntri; i++, tp++) {
	    if(tp->p1 < 0 || tp->p1 >= npts) {
		fprintf(stderr, "triangle %d, p1 = %d\n", i, tp->p1);
	    }
	    if(tp->p2 < 0 || tp->p2 >= npts) {
		fprintf(stderr, "triangle %d, p2 = %d\n", i, tp->p2);
	    }
	    if(tp->p3 < 0 || tp->p3 >= npts) {
		fprintf(stderr, "triangle %d, p3 = %d\n", i, tp->p3);
	    }
	}
*/
    }

   temp->dot = dot;
   temp->tri = tri;
   temp->contour = contour;
   temp->val = val;
   temp->val1 = NULL;
   temp->vmin1 = temp->vmax1 = 0;
   temp->npts = npts;
   temp->ntri = ntri;
   temp->type = type;
   temp->texture = NULL;
   temp->texenv = 0;
   temp->textype = 0;
   temp->alpha = 0;
   temp->clipplane_tvec[0] = temp->clipplane_tvec[1] = temp->clipplane_tvec[2] = 0.0;
   temp->clipplane_rvec[0] = temp->clipplane_rvec[1] = temp->clipplane_rvec[2] = 0.0;
   temp->clipplane_rvec[3] = 1.0;
   temp->surf_transp    = bit.surf_transp;
   temp->surf_clip      = 0;
   temp->chickenwire    = bit.chickenwire;
   temp->flatshade      = bit.flatshade;
   temp->dot_surf       = bit.dot_surf;
   if((temp->name = strdup(name)) == NULL){
      fprintf(stderr, "Can't allocate enough memory for the surface_name\n");
      exit(-1);
   }
   temp->next = NULL;
   temp->second = NULL;
   temp->trinorm = NULL;
   if(color) temp->matindex = actualmol->nsurfs%21 + surface_mat_base + 22;
   else      temp->matindex = actualmol->nsurfs%21 + surface_mat_base;
   surface_number++;
   actualmol->nsurfs++;

   if(val){
      getPropertyExtrema(temp);
   }

   update_interface_flags();
   return temp;
}



void getPropertyExtrema(Surface *sp)
{
   register int i;
   float *val, vmin, vmax;

   val = sp->val;
   vmax = vmin = *val++;
      
   for(i=1; i<sp->npts; i++, val++){
      if(*val > vmax) vmax = *val;
      if(*val < vmin) vmin = *val;
   }

   sp->vmax = vmax;
   sp->vmin = vmin;

   if(sp->val1) {
      val = sp->val1;
      vmax = vmin = *val++;
      
      for(i=1; i<sp->npts; i++, val++){
         if(*val > vmax) vmax = *val;
         if(*val < vmin) vmin = *val;
      }

      sp->vmax1 = vmax;
      sp->vmin1 = vmin;
   }

   sp->texture  = actual_texture;
   sp->textype = TEX_MAP;
}

void removeSurfFromList(Surface *x)
{
   Surface *i;

   if(!actualmol->firstsurf) return;

   if(x == actualmol->firstsurf) actualmol->firstsurf = x->next;
   else {
      for(i=actualmol->firstsurf; i; i=i->next) {
         if(x == i->next) {
            i->next = x->next;
            break;
         }
      }
   }
   x->next = NULL;
}

void delete_surface(void)
{
   Surface *sp, *surf;

   if(!actualmol || !actualmol->firstsurf) return;

   sp = surf = actualmol->firstsurf;
   if(sp->next) {
      while(sp->next->next) sp = sp->next;
      surf = sp->next;
      sp->next = NULL;
   }

   if(surf->dot) free(surf->dot);
   if(surf->tri) free(surf->tri);
   if(surf->val) free(surf->val);
   if(surf->name) free(surf->name);
   if(surf->trinorm) free(surf->trinorm);
   if(surf->second){
      surf->second->second = NULL;
      delete_surface();
   }
   free(surf);

   if(surf == actualmol->firstsurf) actualmol->firstsurf = NULL;
   actualsurf = actualmol->firstsurf;

   surface_number--;
   actualmol->nsurfs--;
}

void del_actualsurf(Surface *surf)
{
   Surface *second;
   if(!surf) return;

   if(surf->dot) free(surf->dot);
   if(surf->tri) free(surf->tri);
   if(surf->val) free(surf->val);
   if(surf->name) free(surf->name);
   if(surf->trinorm) free(surf->trinorm);
   if(surf->second){
      second = surf->second;
      surf->second->second = NULL;
      del_actualsurf(second);
   }
/*
   if(surf == actualmol->firstsurf) {
      actualmol->firstsurf = NULL;
   }
*/
   removeSurfFromList(surf);
   actualsurf = actualmol->firstsurf;
   free(surf);
   
   surface_number--;
   actualmol->nsurfs--;
   
}



/* sum all the values of the cube along one axis ( for spin-density-maps) */
void project_macu(int axis)
{
   int dim0, dim1, dim2;
   register short i, j, k;
   float *p;

   if(!actualmol->cube_value){
      logprint("Laod a .macu-file first");
      update_logs();
      return;
   }
   

   switch(axis){
      case 0 : dim0 = actualmol->box.nx; dim1 = actualmol->box.ny; dim2 = actualmol->box.nz;
         projection.p2[0]=actualmol->box.x1;
         projection.p2[1]=actualmol->box.y2; 
         projection.p2[2]=actualmol->box.z1;
         projection.p3[0]=actualmol->box.x1; 
         projection.p3[1]=actualmol->box.y1; 
         projection.p3[2]=actualmol->box.z2;
         break;
      case 1 : dim0 = actualmol->box.ny; dim1 = actualmol->box.nx; dim2 = actualmol->box.nz;
         projection.p2[0]=actualmol->box.x2; 
         projection.p2[1]=actualmol->box.y1; 
         projection.p2[2]=actualmol->box.z1;
         projection.p3[0]=actualmol->box.x1; 
         projection.p3[1]=actualmol->box.y1; 
         projection.p3[2]=actualmol->box.z2;
         break;
      case 2 : dim0 = actualmol->box.nz; dim1 = actualmol->box.nx; dim2 = actualmol->box.ny;
         projection.p2[0]=actualmol->box.x2; 
         projection.p2[1]=actualmol->box.y1; 
         projection.p2[2]=actualmol->box.z1;
         projection.p3[0]=actualmol->box.x1; 
         projection.p3[1]=actualmol->box.y2; 
         projection.p3[2]=actualmol->box.z1;
         break;
   }
   projection.p1[0]=actualmol->box.x1; 
   projection.p1[1]=actualmol->box.y1; 
   projection.p1[2]=actualmol->box.z1;
   projection.d1 = dim1;  projection.d2 = dim2;

   if(projection.value){
      free(projection.value[0]);
      free(projection.value);
   }
   if((projection.value = (float**) alloc_3D(dim1, dim2, 0, sizeof(float))) == NULL){
      showinfobox(" can't allocate plane for the projection");
      return;
   }

   for(i=0; i<dim2; i++){
      for(j=0; j<dim1; j++){
         projection.value[i][j] = 0;
         for(k=0; k<dim0; k++){
            switch(axis){
               case 0: projection.value[i][j] += actualmol->cube_value[i][j][k]; break;
               case 1: projection.value[i][j] += actualmol->cube_value[i][k][j]; break;
               case 2: projection.value[i][j] += actualmol->cube_value[k][i][j]; break;
            }
         }
      }
   }

   projection.vmin = projection.vmax = projection.value[0][0];
   for(i=0, p=projection.value[0]; i<dim1*dim2; i++, p++){
      if(*p < projection.vmin) projection.vmin = *p;
      if(*p > projection.vmax) projection.vmax = *p;
   }

   project_mol = actualmol;
//   write_projection();
}


void free_projection(void)
{
   if(projection.value){
      free(projection.value[0]);
      free(projection.value);
      projection.value = NULL;
      project_mol = NULL;
   }
}



void write_patterns(void)
{
   FILE *fp;
   int i, j, k;

   if((fp = fopen("pattern.a", "w")) == NULL){
      fprintf(stderr, "can't open pattern.a\n");
      return;
   }

   for(i=0, k=0; i<256; i++){
      for(j=0; j<npat_A[i]; j++, k++){
         fprintf(fp, "{%2d,%2d,%2d}, ", __A[k].edge0, __A[k].edge2, __A[k].edge1);
      }
      fprintf(fp, "\n");
   }
}




void defaultSurfaces(char *fileName)
{
    readmacu(fileName, LOAD_MACU);

    if(actualmol->cubemin >= 0.0) {
	cutoff = actualmol->cutoff = DENSITY_CUTOFF;
	bit.both_signs = 0;
    }
    else {
	cutoff = actualmol->cutoff = ORBITAL_CUTOFF;
	bit.both_signs = 1;
    }
    
    cubes();
}

double compute_volume(Surface *sp)
{
   int i;
   Triangle *tp;
   float *p1, *p2, *p3, u[2], v[2];
   double dv, area, volume;

/*
   screenprint("Warning: the volume computing method works only for\
 CLOSED and CORRECTLY TRIANGULATED surfaces!");
*/

   volume = 0;

   for(i=0, tp=sp->tri; i<sp->ntri; i++, tp++) {
      p1 = sp->dot[tp->p1].v;
      p2 = sp->dot[tp->p2].v;
      p3 = sp->dot[tp->p3].v;
      u[0] = p2[0] - p1[0]; u[1] = p2[1] - p1[1];
      v[0] = p3[0] - p1[0]; v[1] = p3[1] - p1[1];
      area = (u[0]*v[1] - u[1]*v[0]) * 0.5;

/* area: surface of the projection on the xy-plane, negative for downwards
         pointing triangles */

      dv = area * (p1[2] + p2[2] + p3[2]) * 0.333333333333333;

      volume += dv;
   }

   return volume;
}


/********************/

#define DIST(a,b) sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+\
                       (a[2]-b[2])*(a[2]-b[2]))
#define POW6(a) ((a)*(a)*(a)*(a)*(a)*(a))

void fast_surf(void)
{
   int i, j, k;
   int i1, i2, j1, j2, k1, k2, dk;
   float xyz[3], dx, dy, dz;
   float *p, tmp;
   AtoM *ap;

   if(!actualmol) return;
   if(!dotcubesize) return;

   free_3D(actualmol);

   ap = actualmol->firstatom;
   xmin = xmax = ap->coord[0];
   ymin = ymax = ap->coord[1];
   zmin = zmax = ap->coord[2];
   for(ap=ap->next; ap; ap=ap->next){
      if(ap->coord[0] < xmin) xmin = ap->coord[0];
      if(ap->coord[0] > xmax) xmax = ap->coord[0];
      if(ap->coord[1] < ymin) ymin = ap->coord[1];
      if(ap->coord[1] > ymax) ymax = ap->coord[1];
      if(ap->coord[2] < zmin) zmin = ap->coord[2];
      if(ap->coord[2] > zmax) zmax = ap->coord[2];
   }
   xmin -= (3.0 + dotcubesize);
   xmax += (3.0 + dotcubesize);
   ymin -= (3.0 + dotcubesize);
   ymax += (3.0 + dotcubesize);
   zmin -= (3.0 + dotcubesize);
   zmax += (3.0 + dotcubesize);

   nx = (int)fceil((xmax - xmin)/dotcubesize) + 1;
   ny = (int)fceil((ymax - ymin)/dotcubesize) + 1;
   nz = (int)fceil((zmax - zmin)/dotcubesize) + 1;

   actualmol->box.x1 = xmin;
   actualmol->box.x2 = xmax;
   actualmol->box.y1 = ymin;
   actualmol->box.y2 = ymax;
   actualmol->box.z1 = zmin;
   actualmol->box.z2 = zmax;
   actualmol->box.nx = nx;
   actualmol->box.ny = ny;
   actualmol->box.nz = nz;

/*
   printf("   xmin =%6.2f  ymin =%6.2f  zmin =%6.2f\n", xmin, ymin, zmin);
   printf("   xmax =%6.2f  ymax =%6.2f  zmax =%6.2f\n", xmax, ymax, zmax);
   printf("   nx   =%3d     ny   =%3d     nz   =%3d\n", nx, ny, nz); 
*/

   if((actualmol->cube_value = (float***)  alloc_3D(nx, ny, nz, sizeof(float))) == NULL){
       showinfobox(" can't allocate enough memory for the cube-values\n");
       return;
   }
   if((ccase = (unsigned char**) alloc_3D(nx+1, ny+1, 0, sizeof(char))) == NULL){
      showinfobox(" can't allocate enough memory for the cube-labels\n");
      return;
   }

   actualmol->got_macufile = 1;

   dx = (xmax - xmin)/(nx-1.);
   dy = (ymax - ymin)/(ny-1.);
   dz = (zmax - zmin)/(nz-1.);

   p = actualmol->cube_value[0][0];
   for(i=1; i < (nx*ny*nz); i++, p++) *p = 0.0;

   for(ap=actualmol->firstatom; ap; ap=ap->next){
      k = (int)((ap->coord[2]-zmin)/dz);
      dk = (int)(element[ap->ord].rvdw / dotcubesize * 1.7 + 1);
      k1 = MAX(k-dk, 0); k2 = MIN(k+dk+1, nz);
      for(k=k1, xyz[2]=zmin+k1*dz; k<k2; k++, xyz[2]+=dz) {
         j = (int)((ap->coord[1]-ymin)/dy);
         j1 = MAX(j-dk, 0); j2 = MIN(j+dk+1, ny);
         for(j=j1, xyz[1]=ymin+j1*dy; j<j2; j++, xyz[1]+=dy) {
            i = (int)((ap->coord[0]-xmin)/dx);
            i1 = MAX(i-dk, 0); i2 = MIN(i+dk+1, nx);
            for(i=i1, xyz[0]=xmin+i1*dx; i<i2; i++, xyz[0]+=dx) {
               tmp = element[ap->ord].rvdw/DIST(xyz, ap->coord);
               actualmol->cube_value[k][j][i] += POW6(tmp);
            }
         }
      }
   }

   p = actualmol->cube_value[0][0];
   cubemin = cubemax = *p;
   for(i=1; i < (nx*ny*nz); i++, p++){
      if(*p > cubemax) cubemax = *p;
      if(*p < cubemin) cubemin = *p;
   }
   if(cubemax > 1000) cubemax = 1000;
   actualmol->cubemin = cubemin;
   actualmol->cubemax = cubemax;
   actualmol->got_macufile = 1;
   actualmol->cutoff = cutoff = 1;
   actualmol->plane_iz = nz/2;
   actualmol->plane_iy = ny/2;
   actualmol->plane_ix = nx/2;

   bit.both_signs = 0;
   cubes();

   free_3D(actualmol);
}



void read_gcube(char *fname, int key)
{
   FILE *fp;
   Box *bp;
   register float *p;
   register int h, i, j, k;
   static char str[50], *s;
   float sum, dv;

   int isorb = 0, natoms, dimx, dimy, dimz;
   int norb, ijnk1;
   float xi, yi, zi, xf, yf, zf, xinc, yinc, zinc;
   float rjnk1, rjnk2, rjnk3;
   char buffer[132];
   char *bgnPtr, *endPtr;
   double value;

   if(!actualmol){
      logprint("readgcube: no molecule loaded");
      update_logs();
      return;
   }

   if((fp = fopen(fname, "r")) == NULL){
       sprintf(str, "can't open %s \n", fname);
       showinfobox(str);
       return;
   }

   fgets(buffer, 132, fp);
   fgets(buffer, 132, fp);
   fgets(buffer, 132, fp);
   sscanf(buffer, "%d%f%f%f", &natoms, &xi, &yi, &zi);
   if( natoms < 0 ) { natoms=-natoms;isorb=1;}
   fgets(buffer, 132, fp);
   sscanf(buffer, "%d%f%f%f", &dimx, &xinc, &rjnk2, &rjnk3);
   xf=(dimx-1)*xinc+xi;
   fgets(buffer, 132, fp);
   sscanf(buffer, "%d%f%f%f", &dimy, &rjnk1, &yinc, &rjnk3);
   yf=(dimy-1)*yinc+yi;
   fgets(buffer, 132, fp);
   sscanf(buffer, "%d%f%f%f", &dimz, &rjnk1, &rjnk2, &zinc);
   zf=(dimz-1)*zinc+zi;

   bp = &actualmol->box;
   if(key == ADD_MACU || key == SUBST_MACU){
      if((nx != dimx) || (ny != dimy) || (nz != dimz)){
         showinfobox("The two gcube-files don't correspond!");
         fclose(fp);
         return;
      }
   }
   xmin = bp->x1 = xi * BOHR;
   xmax = bp->x2 = xf * BOHR;
   ymin = bp->y1 = yi * BOHR;
   ymax = bp->y2 = yf * BOHR;
   zmin = bp->z1 = zi * BOHR;
   zmax = bp->z2 = zf * BOHR;
   nx   = bp->nx = dimx;
   ny   = bp->ny = dimy;
   nz   = bp->nz = dimz;
   nMacuPlanes = nz;

   i=0;
   while (i<natoms) {
     i=i+1;
     fgets(buffer, 132, fp);
   }

   if( isorb == 1) {
      fgets(buffer, 132, fp);
      sscanf(buffer, "%d%d", &ijnk1, &norb);
      if( ijnk1 > 1 ) { 
         showinfobox("Can do only one orbital at a time");
         fclose(fp);
         return;
      }
      sprintf(str," this is orbital number %d\n",norb);
      logprint(str);
   }


   logprint("");
   if((s = strrchr(fname, '/')) == NULL) s = fname;
   else s++;
   logprint(s);
   sprintf(str, "xmin =  ymin =  zmin =");
   logprint(str);
   sprintf(str, "%6.2f  %6.2f  %6.2f", xmin, ymin, zmin);
   logprint(str);
   sprintf(str, "xmax =  ymax =  zmax =");
   logprint(str);
   sprintf(str, "%6.2f  %6.2f  %6.2f", xmax, ymax, zmax);
   logprint(str);
   sprintf(str, "nx = %3d  ny = %3d  nz = %3d", nx, ny, nz); 
   logprint(str);

   strncpy(macu_name, s, 29);

   switch(key){
      case LOAD_MACU :
         free_3D(actualmol);
         if((actualmol->cube_value = (float***) alloc_3D(nx, ny, nz, sizeof(float))) == NULL){
            showinfobox(" can't allocate enough memory for the cube-values\n");
            return;
         }
         if((ccase = (unsigned char**) alloc_3D(nx+1, ny+1, 0, sizeof(char))) == NULL){
            showinfobox(" can't allocate enough memory for the cube-labels\n");
            return;
         }
         i = j = k = 0;
         while (fgets(buffer, 132, fp) != NULL) {
             bgnPtr = buffer;
             for (h = 0; h < 100; h++) {
                 value = strtod(bgnPtr, &endPtr);
                 if ((value == 0.0) && (endPtr == bgnPtr)) {
                     continue;
                 }
                 actualmol->cube_value[k][j][i] = (float)value;
                 bgnPtr = endPtr;
                 if(k<nz-1) k++;
                 else if(j<ny-1) {
                    j++; k = 0;
                 }
                 else if(i<nx-1) {
                    i++; j = 0; k = 0;
                 }
             }
         }
         break;
      case ADD_MACU :
         i = j = k = 0;
         while (fgets(buffer, 132, fp) != NULL) {
             bgnPtr = buffer;
             for (h = 0; h < 100; h++) {
                 value = strtod(bgnPtr, &endPtr);
                 if ((value == 0.0) && (endPtr == bgnPtr)) {
                     continue;
                 }
                 actualmol->cube_value[k][j][i] += (float)value;
                 bgnPtr = endPtr;
                 if(k<nz-1) k++;
                 else if(j<ny-1) {
                    j++; k = 0;
                 }
                 else if(i<nx-1) {
                    i++; j = 0; k = 0;
                 }
             }
         }
         break;
      case SUBST_MACU :
         i = j = k = 0;
         while (fgets(buffer, 132, fp) != NULL) {
             bgnPtr = buffer;
             for (h = 0; h < 100; h++) {
                 value = strtod(bgnPtr, &endPtr);
                 if ((value == 0.0) && (endPtr == bgnPtr)) {
                     continue;
                 }
                 actualmol->cube_value[k][j][i] -= (float)value;
                 bgnPtr = endPtr;
                 if(k<nz-1) k++;
                 else if(j<ny-1) {
                    j++; k = 0;
                 }
                 else if(i<nx-1) {
                    i++; j = 0; k = 0;
                 }
             }
         }
         break;
   }
   fclose(fp);

   p = actualmol->cube_value[0][0];
   cubemin = cubemax = *p;
   sum = 0;
   for(i=1; i < (nx*ny*nz); i++, p++){
//printf("%d %f\n", i, *p);
      if(*p > cubemax) cubemax = *p;
      if(*p < cubemin) cubemin = *p;
      sum += *p;
   }
   if(cubemin < 0.000001 && cubemin > -0.000001) cubemin = 0.0;
   if(cubemax < 0.000001 && cubemax > -0.000001) cubemax = 0.0;
   actualmol->cubemin = cubemin;
   actualmol->cubemax = cubemax;
   actualmol->got_macufile = 1;
   actualmol->cutoff = cutoff = 0.00;
   actualmol->plane_iz = cubehead.nz/2;
   actualmol->plane_iy = cubehead.ny/2;
   actualmol->plane_ix = cubehead.nx/2;
   init_cutplane(actualmol);

   dv = (xmax-xmin) * _1_BOHR / (nx-1.) *
        (ymax-ymin) * _1_BOHR / (ny-1.) *
        (zmax-zmin) * _1_BOHR / (nz-1.);
   sum *= dv;
   logprint("sum of voxel contributions");
   sprintf(str, " (a.u.) = %f", sum);
   logprint(str);
   update_logs();

   if(surfaceglui) {
      surfaceglui->sync_live();
   }
   if(planeglui) {
      planeglui->sync_live();
   }

}

static int getIntValue(char *str, int *ival)
{
    if(!find_string_t41(str)) {
	fprintf(stderr, "can't find %s\n", str);
	return 0;
    }
    fgets(line, 255, fpt);
    fgets(line, 255, fpt);
    sscanf(line, "%d", ival);
    return 1;
}

void read_t41(char *fname, int key, char *density)
{
   long fpos, fpos2;
   Box *bp;
   char str[50], *s, *spin = NULL;
   float dim[6], x, y, z, value = 0, value2 = 0, dv = 0, sum = 0;
   register float *p;
   int n[3], i, j, k;

   if(!actualmol){
      logprint("read_t41: no molecule loaded");
      update_logs();
      return;
   }

   if((fpt = fopen(fname, "r")) == NULL){
       sprintf(str, "can't open %s \n", fname);
       showinfobox(str);
       return;
   }

   if(!find_string_t41("Start_point")) {
	showinfobox("can't find start-point\n");
	return;
   }
   fgets(line, 255, fpt);
   fgets(line, 255, fpt);
   sscanf(line, "%f%f%f", dim, dim+2, dim+4);

   if(!getIntValue("nr of points x", n)) return;
   if(!getIntValue("nr of points y", n+1)) return;
   if(!getIntValue("nr of points z", n+2)) return;

   if(!find_string_t41("x-vector")) {
	showinfobox("can't find x-vector\n");
	return;
   }
   fgets(line, 255, fpt);
   fgets(line, 255, fpt);
   sscanf(line, "%f", &x);

   if(!find_string_t41("y-vector")) {
	showinfobox("can't find y-vector\n");
	return;
   }
   fgets(line, 255, fpt);
   fgets(line, 255, fpt);
   sscanf(line, "%*f%f", &y);

   if(!find_string_t41("z-vector")) {
	showinfobox("can't find z-vector\n");
	return;
   }
   fgets(line, 255, fpt);
   fgets(line, 255, fpt);
   sscanf(line, "%*f%*f%f", &z);

   dim[1] = dim[0] + (n[0]-1)*x;
   dim[3] = dim[2] + (n[1]-1)*y;
   dim[5] = dim[4] + (n[2]-1)*z;
   for(i=0; i<6; i++) dim[i] *= BOHR;

   bp = &actualmol->box;
   if(key == ADD_MACU || key == SUBST_MACU){
      if((nx != n[0]) || (ny != n[1]) || (nz != n[2])){
         showinfobox("The two t41-files don't correspond!");
         fclose(fpt);
         return;
      }
   }
   xmin = bp->x1 = dim[0];
   xmax = bp->x2 = dim[1];
   ymin = bp->y1 = dim[2];
   ymax = bp->y2 = dim[3];
   zmin = bp->z1 = dim[4];
   zmax = bp->z2 = dim[5];
   nx   = bp->nx = n[0];
   ny   = bp->ny = n[1];
   nz   = bp->nz = n[2];
   nMacuPlanes = nz;

   logprint("");
   if((s = strrchr(fname, '/')) == NULL) s = fname;
   else s++;
   logprint(s);
   sprintf(str, "xmin =  ymin =  zmin =");
   logprint(str);
   sprintf(str, "%6.2f  %6.2f  %6.2f", xmin, ymin, zmin);
   logprint(str);
   sprintf(str, "xmax =  ymax =  zmax =");
   logprint(str);
   sprintf(str, "%6.2f  %6.2f  %6.2f", xmax, ymax, zmax);
   logprint(str);
   sprintf(str, "nx = %3d  ny = %3d  nz = %3d", nx, ny, nz); 
   logprint(str);

   strncpy(macu_name, s, 29);

   logprint("");
   logprint("reading");
   logprint(density);
   if(strstr(density, "Spin")) spin = strcpy(density, "Density");
//   if(!find_string_t41(density)) {
   if(!find_string_t41(strcat(density," "))) {
	showinfobox("can't find Density\n");
	return;
   }

   if(key == LOAD_MACU){
      free_3D(actualmol);
      if((actualmol->cube_value = (float***) alloc_3D(nx, ny, nz, sizeof(float))) == NULL){
         showinfobox(" can't allocate enough memory for the cube-values\n");
         return;
      }
      if((ccase = (unsigned char**) alloc_3D(nx+1, ny+1, 0, sizeof(char))) == NULL){
         showinfobox(" can't allocate enough memory for the cube-labels\n");
         return;
      }
   }
   if(strstr(density, "Density")) {
       if(strstr(line, "Density_A")) { /* in unrestricted calc. density is */
                                       /* printed for spin A and spin b seperately */
                                       /* add values for density, substract for */
                                       /* spin density */
                                       /* read A and B at the same time, jump with */
                                       /* two pointers. Maybe not very efficient but */
                                       /* straight forward to implement */

         fgets(line, 255, fpt);
         fpos = ftell(fpt);
         if(!find_string_t41("Density_B")) {
             showinfobox("can't find Density\n");
             return;
         }
         fgets(line, 255, fpt);
         fpos2 = ftell(fpt);

         for(i=0; i<n[2]; i++) {
            for(j=0; j<n[1]; j++) {
               for(k=0; k<n[0]; k++) {
                 fseek(fpt, fpos, SEEK_SET);
                 if(!fscanf(fpt, "%f", &value)) {
                     sprintf(str, "can't read value %d in plane %d\n", j, i);
                     showinfobox(str);
                     return;
                 }
                 fpos = ftell(fpt);
                 fseek(fpt, fpos2, SEEK_SET);
                 if(!fscanf(fpt, "%f", &value2)) {
                     sprintf(str, "can't read value %d in plane %d\n", j, i);
                     showinfobox(str);
                     return;
                 }
                 fpos2 = ftell(fpt);
                 if(spin) { /* subtract -> spin */
                    value -= value2;
                 }
                 else { /* add -> density */
                    value += value2;
                 }
                 switch(key){
                    case LOAD_MACU:
                       actualmol->cube_value[i][j][k] = value;
                    break;
                    case ADD_MACU:
                       actualmol->cube_value[i][j][k] += value;
                    break;
                    case SUBST_MACU:
                       actualmol->cube_value[i][j][k] -= value;
                    break;
                 }
               }
             }
         }

         fclose(fpt);
         p = actualmol->cube_value[0][0];
         cubemin = cubemax = *p;
         sum = 0;
         for(i=1; i < (nx*ny*nz); i++, p++){
            if(*p > cubemax) cubemax = *p;
            if(*p < cubemin) cubemin = *p;
            sum += *p;
         }
         if(cubemin < 0.000001 && cubemin > -0.000001) cubemin = 0.0;
         if(cubemax < 0.000001 && cubemax > -0.000001) cubemax = 0.0;
         actualmol->cubemin = cubemin;
         actualmol->cubemax = cubemax;
         actualmol->got_macufile = 1;
         actualmol->cutoff = cutoff = 0.00;
         actualmol->plane_iz = cubehead.nz/2;
         actualmol->plane_iy = cubehead.ny/2;
         actualmol->plane_ix = cubehead.nx/2;
         init_cutplane(actualmol);

         dv = (xmax-xmin) * _1_BOHR / (nx-1.) *
              (ymax-ymin) * _1_BOHR / (ny-1.) *
              (zmax-zmin) * _1_BOHR / (nz-1.);
         sum *= dv;
         logprint("sum of voxel contributions");
         sprintf(str, " (a.u.) = %f", sum);
         logprint(str);
         update_logs();

         if(surfaceglui) {
            surfaceglui->sync_live();
         }
         if(planeglui) {
            planeglui->sync_live();
         }
         return;
       }
       fgets(line, 255, fpt);
   }
   else {
       do {
           fpos = ftell(fpt);
       } while(find_string_t41(density));
       fseek(fpt, fpos, SEEK_SET);
       fgets(line, 255, fpt);
       fgets(line, 255, fpt);
   }

   for(i=0; i<n[2]; i++) {
      for(j=0; j<n[1]; j++) {
         for(k=0; k<n[0]; k++) {
	    if(!fscanf(fpt, "%f", &value)) {
                sprintf(str, "can't read value %d in plane %d\n", j, i);
		showinfobox(str);
		return;
	    }
            switch(key){
               case LOAD_MACU:
                  actualmol->cube_value[i][j][k] = value;
               break;
               case ADD_MACU:
                  actualmol->cube_value[i][j][k] += value;
               break;
               case SUBST_MACU:
                  actualmol->cube_value[i][j][k] -= value;
               break;
            }
         }
      }
   }

   fclose(fpt);

   p = actualmol->cube_value[0][0];
   cubemin = cubemax = *p;
   sum = 0;
   for(i=1; i < (nx*ny*nz); i++, p++){
      if(*p > cubemax) cubemax = *p;
      if(*p < cubemin) cubemin = *p;
      sum += *p;
   }
   if(cubemin < 0.000001 && cubemin > -0.000001) cubemin = 0.0;
   if(cubemax < 0.000001 && cubemax > -0.000001) cubemax = 0.0;
   actualmol->cubemin = cubemin;
   actualmol->cubemax = cubemax;
   actualmol->got_macufile = 1;
   actualmol->cutoff = cutoff = 0.00;
   actualmol->plane_iz = cubehead.nz/2;
   actualmol->plane_iy = cubehead.ny/2;
   actualmol->plane_ix = cubehead.nx/2;
   init_cutplane(actualmol);

   dv = (xmax-xmin) * _1_BOHR / (nx-1.) *
        (ymax-ymin) * _1_BOHR / (ny-1.) *
        (zmax-zmin) * _1_BOHR / (nz-1.);
   sum *= dv;
   logprint("sum of voxel contributions");
   sprintf(str, " (a.u.) = %f", sum);
   logprint(str);
   update_logs();

   if(surfaceglui) {
      surfaceglui->sync_live();
   }
   if(planeglui) {
      planeglui->sync_live();
   }
}


char *find_string_t41(char *s)
{
   previous_line = ftell(fpt);
   do {
      if(!fgets(line, 255, fpt)) return NULL;
      if(strstr(line, s)) return line;
      preprepre = preprevious;
      preprevious = previous_line;
      previous_line = ftell(fpt);
   } while (1);
}
