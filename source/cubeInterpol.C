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


#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "general.h"
#include "cubeInterpol.h"
#include "macu.h"
#include "drawing.h"
#include "maininterf.h"
#include "chooseinterf.h"
#include "texture.h"

/*----------------------------------------*/
static void gridValuesOnSurface(Surface *sp, float ***cube, Box *bp);
static float getValueOnSurface(Surfdot *pt);
static float cubeInterpol(float dx, float dy, float dz,
       float p000, float p100, float p110, float p010,
       float p001, float p101, float p111, float p011);

static float xmin, xmax, ymin, ymax, zmin, zmax;
static float invdx, invdy, invdz;
static int   nx, ny, nz;
static float ***cu;

/*----------------------------------------*/

void gridValues(Surface *picked) /* callback of the GRID VALUE button */
{
   if(!actualmol){
      showinfobox("no molecule");
      return;
   }
   if(!actualmol->firstsurf){
      showinfobox("no surface");
      return;
   }
   if(!actualmol->cube_value){
      showinfobox("no 3D lattice");
      return;
   }

   gridValuesOnSurface(picked, actualmol->cube_value, &actualmol->box);
   getPropertyExtrema(picked);

   if(!firsttexture) init_texture();

   picked->texture = firsttexture;
   picked->texenv = GL_MODULATE;
   picked->textype = TEX_MAP;
   set_vmin_vmax(picked);

}

void gridValues_no_pick(void) /* callback of the GRID VALUE button */
{
   Surface *picked; 
   int count;

   if(!actualmol){
      showinfobox("no molecule");
      return;
   }
   if(!actualmol->firstsurf){
      showinfobox("no surface");
      return;
   }
   if(!actualmol->cube_value){
      showinfobox("no 3D lattice");
      return;
   }

/*
   if((picked = select_one_surface()) == NULL) return;

   gridValuesOnSurface(picked, actualmol->cube_value, &actualmol->box);

   getPropertyExtrema(picked);
*/
 
   picked = actualmol->firstsurf;
   for(count=1; count < surface_live_var; count++) picked = picked->next;

   gridValuesOnSurface(picked, actualmol->cube_value, &actualmol->box);
   getPropertyExtrema(picked);

//   drawit();
}



static void gridValuesOnSurface(Surface *sp, float ***cube, Box *bp)
{
   register int i;
   Surfdot *j;
   float *val;

   if((val = (float*) malloc(sp->npts*sizeof(float))) == NULL) {
      showinfobox("can't allocate the value array");
      return;
   }

   cu = cube;

   xmin = bp->x1; ymin = bp->y1; zmin = bp->z1;
   xmax = bp->x2; ymax = bp->y2; zmax = bp->z2;
   nx = bp->nx; ny = bp->ny; nz = bp->nz;
   invdx = (nx-1.)/(xmax-xmin);
   invdy = (ny-1.)/(ymax-ymin);
   invdz = (nz-1.)/(zmax-zmin);

   for(i=0, j=sp->dot; i<sp->npts; i++, j++){
      val[i] = getValueOnSurface(j);
   }

   if(sp->val) free(sp->val);
   sp->val = val;
}



static float getValueOnSurface(Surfdot *pt)
{
   float x, y, z, px, py, pz;
   int ix, iy, iz;

   x = pt->v[0]; y = pt->v[1]; z = pt->v[2];

   if(x < xmin || x > xmax) return 0;
   if(y < ymin || y > ymax) return 0;
   if(z < zmin || z > zmax) return 0; /* outside the cube */

   ix = (int)floorf((x - xmin)*invdx);
   iy = (int)floorf((y - ymin)*invdy);
   iz = (int)floorf((z - zmin)*invdz);

   px = (x-xmin)*invdx - ix; /* ratio */
   py = (y-ymin)*invdy - iy; /* ratio */
   pz = (z-zmin)*invdz - iz; /* ratio */

   return cubeInterpol(px, py, pz,
      cu[iz][iy][ix], cu[iz][iy][ix+1], cu[iz][iy+1][ix+1], cu[iz][iy+1][ix],
      cu[iz+1][iy][ix], cu[iz+1][iy][ix+1], cu[iz+1][iy+1][ix+1], cu[iz+1][iy+1][ix]);
}



static float cubeInterpol(float dx, float dy, float dz,
       float p000, float p100, float p110, float p010,
       float p001, float p101, float p111, float p011)
{
   float q00, q01, q10, q11, r0, r1;

   q00 = p000 + (p001-p000)*dz;
   q10 = p100 + (p101-p100)*dz;
   q01 = p010 + (p011-p010)*dz;
   q11 = p110 + (p111-p110)*dz;

   r0 = q00 + (q01-q00)*dy;
   r1 = q10 + (q11-q10)*dy;

   return r0 + (r1-r0)*dx;
}


