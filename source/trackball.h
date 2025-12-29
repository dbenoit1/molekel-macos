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
 * trackball.h
 * A virtual trackball implementation
 * Written by Gavin Bell for Silicon Graphics, November 1988.
 */

/*
 * NOTE: If, for some reason, gl.h shouldn't be included, one could
 * just define the Matrix type, like so:
 * typedef float Matrix[4][4];
 */
/* #include <gl/gl.h> */
/* typedef float Matrix[4][4]; in molekel.h */

/*
 *	Pass the x and y coordinates of the last and current positions of
 * the mouse, scaled so they are from (-1.0 ... 1.0).
 *
 * if ox,oy is the window's center and sizex,sizey is its size, then
 * the proper transformation from screen coordinates (sc) to world
 * coordinates (wc) is:
 * wcx = (2.0 * (scx-ox)) / (float)sizex - 1.0
 * wcy = (2.0 * (scy-oy)) / (float)sizey - 1.0
 *
 * For a really easy interface to this see 'ui.h'.
 */
void trackball(float *e, float p1x, float p1y, float p2x, float p2y);

/*
 *	Given two sets of Euler paramaters, add them together to get an
 * equivalent third set.  When incrementally adding them, the first
 * argument here should be the new rotation, the secon and third the
 * total rotation (which will be over-written with the resulting new
 * total rotation).
 */
void add_eulers(float *e1, float *e2, float *dest);

/*
 *	A useful function, builds a rotation matrix in Matrix based on
 * given Euler paramaters.
 */
void build_rotmatrix(Matrix m, float *e);

/*
 * This function computes the Euler paramaters given an xyz axis (the
 * first argument, 3 floats) and angle (expressed in radians, the
 * second argument).  The result is put into the third argument, which
 * must be an array of 4 floats.
 */
void axis_to_euler(float *a, float phi, float *e);

void invert_euler(float *a, float *e);

/** originally from vect.h **/

/*
 *	Definition for vector math.  Vectors are just arrays of 3 floats.
 */

#ifndef VECTDEF
#define VECTDEF

#ifndef _POLY9
#include <math.h>
#endif

 float *vnew(void);
 float *vclone(float *v);
 void  vcopy(float *v1, float *v2);
 void  vprint(float *v);
 void  vset(float *v, float x, float y, float z);
 void  vzero(float *v);
 void  vnormal(float *v);
 float vlength(float *v);
 void  vscale(float *v, float div);
 void  vmult(float *src1, float *src2, float *dst);
 void  vadd(float *src1, float *src2, float *dst);
 void  vsub(float *src1, float *src2, float *dst);
 void  vhalf(float *v1, float *v2, float *half);
 float vdot(float *v1, float *v2);
 void  vcross(float *v1, float *v2, float *cross);
 void  vdirection(float *v1, float *dir);
 void  vreflect(float *in, float *mirror, float *out);
 void  vmultmatrix(float m1[][4], float m2[][4], float prod[][4]);
 void  vtransform(float *v, float mat[][4], float *vt);

#define VNEW()           ((float *)malloc(3*sizeof(float)))
#define VCOPY(a,b)       (b[0]=a[0],b[1]=a[1],b[2]=a[2])
#define VSET(v,x,y,z)    (v[0]=(x),v[1]=(y),v[2]=(z))
#define VZERO(v)         (v[0]=v[1]=v[2]=0.0)
#define VLENGTH(v)	 (fsqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]))
#define VSCALE(v,d)      (v[0]*=(d),v[1]*=(d),v[2]*=(d))
#define VMULT(a,b,c)     (c[0]=a[0]*b[0],c[1]=a[1]*b[1],c[2]=a[2]*b[2])
#define VADD(a,b,c)      (c[0]=a[0]+b[0],c[1]=a[1]+b[1],c[2]=a[2]+b[2])
#define VSUB(a,b,c)      (c[0]=a[0]-b[0],c[1]=a[1]-b[1],c[2]=a[2]-b[2])
#define VDOT(a,b)        (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
#define VCROSS(a,b,c)    (c[0]=a[1]*b[2]-a[2]*b[1],\
                          c[1]=a[2]*b[0]-a[0]*b[2],\
                          c[2]=a[0]*b[1]-a[1]*b[0])


#endif   /* VECTDEF */

/*
$Log: trackball.h,v $
Revision 1.2  2000/03/31 13:40:35  portmann
macu, surfaces, clipplane implemented

Revision 1.1  1999/12/21 13:31:15  portmann
Add start files for glutwin!

Revision 1.1  1999/12/08 12:46:01  portmann
Conversion XGLmolekel to glutmolekel.
Functionality maintained!

Revision 1.1  1999/10/14 13:13:13  portmann
basic mouse and keybord manipulation done!

 * Revision 1.2  1995/09/13  09:20:00  flukiger
 * *** empty log message ***
 *
 * Revision 1.2  1995/09/13  09:20:00  flukiger
 * *** empty log message ***
 *
 * Revision 1.1  94/06/16  12:07:11  flukiger
 * Initial revision
 * 
*/
