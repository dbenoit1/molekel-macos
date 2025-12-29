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
#include "constant.h"
#include "molekel.h"
#include "general.h"
#include "box.h"
#include "drawing.h"

extern unsigned long my_pick(short key);

/*----------------------------------------*/

Box box;
float orbcubesize = 0.25;

/*----------------------------------------*/

void create_box(void)
{
   register AtoM *atom;

   if(!actualmol) return;

   box.x1 = box.x2 = actualmol->firstatom->coord[0];
   box.y1 = box.y2 = actualmol->firstatom->coord[1];
   box.z1 = box.z2 = actualmol->firstatom->coord[2];

   for(atom = actualmol->firstatom->next; atom; atom = atom->next){
      if(atom->coord[0] < box.x1) box.x1 = atom->coord[0];
      if(atom->coord[0] > box.x2) box.x2 = atom->coord[0];
      if(atom->coord[1] < box.y1) box.y1 = atom->coord[1];
      if(atom->coord[1] > box.y2) box.y2 = atom->coord[1];
      if(atom->coord[2] < box.z1) box.z1 = atom->coord[2];
      if(atom->coord[2] > box.z2) box.z2 = atom->coord[2];
   }

   box.x1 -= 1.5; box.x2 += 1.5;
   box.y1 -= 1.5; box.y2 += 1.5;
   box.z1 -= 1.5; box.z2 += 1.5;
   box.cubesize = orbcubesize;

}


void reset_box(void)
{
   create_box();
}


void define_box(float *dim, int *ncub)
{
   float x, y, z;

   if(!actualmol) return;

   dim[0] = box.x1;   dim[1] = box.x2;
   dim[2] = box.y1;   dim[3] = box.y2;
   dim[4] = box.z1;   dim[5] = box.z2;

   x = dim[1] - dim[0]; y = dim[3] - dim[2]; z = dim[5] - dim[4];

   ncub[0] = (int)fceil(x/box.cubesize) + 1;
   ncub[1] = (int)fceil(y/box.cubesize) + 1;
   ncub[2] = (int)fceil(z/box.cubesize) + 1;
}


int correct_box(void)
{
   if(!actualmol) return 0;

   return (
      order(&box.x1, &box.x2) ||
      order(&box.y1, &box.y2) ||
      order(&box.z1, &box.z2) );
}




int order(float *small, float *big)
{
   float temp;

   if(*small <= *big)return 0;

   temp = *small;
   *small = *big;
   *big = temp;
   return 1;
}

