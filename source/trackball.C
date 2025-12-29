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
 *	Implementation of a virtual trackball.  See trackball.h for the
 * interface to these routines.
 *	Implemented by Gavin Bell, lots of ideas from Thant Tessman and
 * the August '88 issue of Siggraph's "Computer Graphics," pp. 121-129.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include "molekel.h"
#include "trackball.h"


/* Size of virtual trackball, in relation to window size */
float trackballsize = 0.4;


/* Function prototypes for local functions */
float tb_project_to_sphere(float r, float x, float y);
void normalize_euler(float *e);

/*
 * Ok, simulate a track-ball.  Project the points onto the virtual
 * trackball, then figure out the axis of rotation, which is the cross
 * product of P1 P2 and O P1 (O is the center of the ball, 0,0,0)
 * Note:  This is a deformed trackball-- is a trackball in the center,
 * but is deformed into a hyperbolic solid of rotation away from the
 * center.
 * 
 * It is assumed that the arguments to this routine are in the range
 * (-1.0 ... 1.0)
 */
void trackball(float *e, float p1x, float p1y, float p2x, float p2y)
{
	float p1[3];
	float p2[3];
	float d[3];
	float phi;	/* how much to rotate about axis */
	float a[3];	/* Axis of rotation */

	VZERO(a);

	if (p1x == p2x && p1y == p2y)
	{
		VZERO(e);	/* Zero rotation */
		e[3] = 1.0;
		return ;
	}

/*
 * First, figure out z-coordinates for projection of P1 and P2 to
 * deformed sphere
 */
	VSET(p1, p1x, p1y, tb_project_to_sphere(trackballsize, p1x, p1y));
	VSET(p2, p2x, p2y, tb_project_to_sphere(trackballsize, p2x, p2y));

/*
 *	Now, we want the cross product of P1 and P2
 *  (Or cross product of p2 and p1... this was determined by trial and
 * error, so there may be a compensating mathematical boo-boo
 * somewhere else).
 */
	VCROSS(p2, p1, a);

/*
 *	Figure out how much to rotate around that axis.
 */
	VSUB(p1, p2, d);
	phi = fsin(VLENGTH(d) / (2.0 * trackballsize));

	axis_to_euler(a, phi, e);
}

/*
 *	Given an axis and angle, compute euler paramaters
 */
void axis_to_euler(float *a, float phi, float *e)
{
	vnormal(a);	/* Normalize axis */
	VCOPY(a, e);
	VSCALE(e, fsin(phi/2.0));
	e[3] = fcos(phi/2.0);
}

void invert_euler(float *a, float *e)
{
   e[0] = a[0] * -1;
   e[1] = a[1] * -1;
   e[2] = a[2] * -1;
   e[3] = a[3];
}

/*
 * Project an x,y pair onto a sphere of radius r OR a hyperbolic sheet
 * if we are away from the center of the sphere.
 */
float tb_project_to_sphere(float r, float x, float y)
{
	float d;

/*
	float t, z;

	d = sqrt(x*x + y*y);

	if (d < r*M_SQRT1_2)	* Inside sphere *
	{
		z = sqrt(r*r - d*d);
	}
	else	* On hyperbola *
	{
		t = r / M_SQRT2;
		z = t*t / d;
	}

   return z;

*/

/* Gauss-curve instead of sphere and hyperbola */

   d = (x*x + y*y)/r/r;
   return exp(-d);

}

/*
 *	Given two rotations, e1 and e2, expressed as Euler paramaters,
 * figure out the equivalent single rotation and stuff it into dest.
 * 
 * This routine also normalizes the result every COUNT times it is
 * called, to keep error from creeping in.
 */
#define COUNT 100
void add_eulers(float *e1, float *e2, float *dest)
{
	static int count=0;
	register int i;
	float t1[3], t2[3], t3[3];
	float tf[4];

	VCOPY(e1, t1); VSCALE(t1, e2[3]);
	VCOPY(e2, t2); VSCALE(t2, e1[3]);
	VCROSS(e2, e1, t3);
	VADD(t1, t2, tf);
	VADD(t3, tf, tf);
	tf[3] = e1[3] * e2[3] - VDOT(e1, e2);

	for (i = 0 ; i < 4 ;i++)
	{
		dest[i] = tf[i];
	}

	if (++count > COUNT)
	{
		count = 0;
		normalize_euler(dest);
	}
}

/*
 * Euler paramaters always obey:  a^2 + b^2 + c^2 + d^2 = 1.0
 * We'll normalize based on this formula.  Also, normalize greatest
 * component, to avoid problems that occur when the component we're
 * normalizing gets close to zero (and the other components may add up
 * to more than 1.0 because of rounding error).
 */
void normalize_euler(float *e)
{	/* Normalize result */
	int which, i;
	float gr;

	which = 0;
	gr = e[which];
	for (i = 1 ; i < 4 ; i++)
	{
		if (fabs(e[i]) > fabs(gr))
		{
			gr = e[i];
			which = i;
		}
	}

	e[which] = 0.0;

	e[which] = fsqrt(1.0 - (e[0]*e[0] + e[1]*e[1] +
		e[2]*e[2] + e[3]*e[3]));

	/* Check to see if we need negative square root */
	if (gr < 0.0)
		e[which] = -e[which];
}

/*
 * Build a rotation matrix, given Euler paramaters.
 */
void build_rotmatrix(Matrix m, float *e)
{
	m[0][0] = 1 - 2.0 * (e[1] * e[1] + e[2] * e[2]);
	m[0][1] = 2.0 * (e[0] * e[1] - e[2] * e[3]);
	m[0][2] = 2.0 * (e[2] * e[0] + e[1] * e[3]);
	m[0][3] = 0.0;

	m[1][0] = 2.0 * (e[0] * e[1] + e[2] * e[3]);
	m[1][1] = 1 - 2.0 * (e[2] * e[2] + e[0] * e[0]);
	m[1][2] = 2.0 * (e[1] * e[2] - e[0] * e[3]);
	m[1][3] = 0.0;

	m[2][0] = 2.0 * (e[2] * e[0] - e[1] * e[3]);
	m[2][1] = 2.0 * (e[1] * e[2] + e[0] * e[3]);
	m[2][2] = 1 - 2.0 * (e[1] * e[1] + e[0] * e[0]);
	m[2][3] = 0.0;

	m[3][0] = 0.0;
	m[3][1] = 0.0;
	m[3][2] = 0.0;
	m[3][3] = 1.0;
}

/** originally from the file vect.c **/

/*
 * vect -
 *	Various functions to support operations on vectors.
 *
 * David M. Ciemiewicz, Mark Grossman, Henry Moreton, and Paul Haeberli
 *
 * Modified for my own nefarious purposes-- Gavin Bell
 *                     prototyping : pff
 */

float *vnew(void)
{
	register float *v;

	v = (float *) malloc(sizeof(float)*3);
	return v;
}

float *vclone(float *v)
{
	register float *c;

	c = VNEW();
	VCOPY(v, c);
	return c;
}

void vcopy(float *v1, float *v2)
{
	register int i;
	for (i = 0 ; i < 3 ; i++)
		v2[i] = v1[i];
}

void vprint(float *v)
{
	printf("x: %f y: %f z: %f\n",v[0],v[1],v[2]);
}

void vset(float *v, float x, float y, float z)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
}

void vzero(float *v)
{
	v[0] = 0.0;
	v[1] = 0.0;
	v[2] = 0.0;
}

void vnormal(float *v)
{
	vscale(v,1.0/VLENGTH(v));
}

float vlength(float *v)
{
	return fsqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void vscale(float *v, float div)
{
	v[0] *= div;
	v[1] *= div;
	v[2] *= div;
}

void vmult(float *src1, float *src2, float *dst)
{
	dst[0] = src1[0] * src2[0];
	dst[1] = src1[1] * src2[1];
	dst[2] = src1[2] * src2[2];
}

void vadd(float *src1, float *src2, float *dst)
{
	dst[0] = src1[0] + src2[0];
	dst[1] = src1[1] + src2[1];
	dst[2] = src1[2] + src2[2];
}

void vsub(float *src1, float *src2, float *dst)
{
	dst[0] = src1[0] - src2[0];
	dst[1] = src1[1] - src2[1];
	dst[2] = src1[2] - src2[2];
}

void vhalf(float *v1, float *v2, float *half)
{
	float len;

	VADD(v2,v1,half);
	len = VLENGTH(half);
	if(len>0.0001)
		VSCALE(half,1.0/len);
	else
		*half = *v1;
}

float vdot(float *v1, float *v2)
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

void vcross(float *v1, float *v2, float *cross)
{
	float temp[3];

	temp[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
	temp[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
	temp[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
	VCOPY(temp, cross);
}

void vdirection(float *v1, float *dir)
{
	*dir = *v1;
	vnormal(dir);
}

void vreflect(float *in, float *mirror, float *out)
{
	float temp[3];

	VCOPY(mirror, temp);
	VSCALE(temp,vdot(mirror,in));
	VSUB(temp,in,out);
	VADD(temp,out,out);
}

void vmultmatrix(float m1[][4], float m2[][4], float prod[][4])
{
	register int row, col;
	float temp[4][4];

	for(row=0 ; row<4 ; row++) 
		for(col=0 ; col<4 ; col++)
			temp[row][col] = m1[row][0] * m2[0][col]
						   + m1[row][1] * m2[1][col]
						   + m1[row][2] * m2[2][col]
						   + m1[row][3] * m2[3][col];
	for(row=0 ; row<4 ; row++) 
		for(col=0 ; col<4 ; col++)
		prod[row][col] = temp[row][col];
}

void vtransform(float *v, float mat[][4], float *vt)
{
	float t[3];

	t[0] = v[0]*mat[0][0] + v[1]*mat[1][0] + v[2]*mat[2][0] + mat[3][0];
	t[1] = v[0]*mat[0][1] + v[1]*mat[1][1] + v[2]*mat[2][1] + mat[3][1];
	t[2] = v[0]*mat[0][2] + v[1]*mat[1][2] + v[2]*mat[2][2] + mat[3][2];
	VCOPY(t, vt);
}


