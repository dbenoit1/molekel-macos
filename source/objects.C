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


/* graphical objects sick and sphere */

#include "main.h"
#include "molekel.h"
#include "general.h"
#include "glutwin.h"
#include "objects.h"

/*----------------------------------------*/

GLuint stick_object, stick_object1, stick_object2;
GLuint rough_stick_object, rough_stick_object1, rough_stick_object2;
GLuint stick_object_thin, rough_stick_object_thin;
GLuint dot_stick_object;
GLuint cone_object, rough_cone_object;
GLuint ugly_sphere, rough_sphere, medium_sphere, fine_sphere;

float bond_factor = 0.15, vdW_factor = 0.2;

/*----------------------------------------*/

void generate_objects(void)
{
        stick_object = glGenLists(1);
        stick_object1 = glGenLists(1);
        stick_object2 = glGenLists(1);
        stick_object_thin = glGenLists(1);
        dot_stick_object = glGenLists(1);
        cone_object = glGenLists(1);
        rough_stick_object = glGenLists(1);
        rough_stick_object1 = glGenLists(1);
        rough_stick_object2 = glGenLists(1);
        rough_stick_object_thin = glGenLists(1);
        rough_cone_object = glGenLists(1);
        generate_stick_object();

        generate_sphere_objects();
}


/*----------------------------------------*/
/* sticks and spheres using gluCylinder and gluSphere */

void generate_sphere_objects(void)
{
        GLUquadricObj *qobj;
        unsigned long col = 0xff007fff;

        ugly_sphere   = glGenLists(1);
        if(ugly_sphere != 0) {
           glNewList(ugly_sphere, GL_COMPILE);
           qobj = gluNewQuadric();
           glColor4ub(col&0xff, col>>8&0xff, col>>16&0xff, col>>24&0xff);
           gluDisk( qobj, 0.,  1.0, 32, 1);
           gluDeleteQuadric(qobj);
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        rough_sphere   = glGenLists(1);
        if(rough_sphere != 0) {
           glNewList(rough_sphere, GL_COMPILE);
           qobj = gluNewQuadric();
           gluSphere( qobj, 1.0,  8, 4);
           gluDeleteQuadric(qobj);
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        medium_sphere   = glGenLists(1);
        if(medium_sphere != 0) {
           glNewList(medium_sphere, GL_COMPILE);
           qobj = gluNewQuadric();
/*           gluSphere( qobj, 1.0,  15, 15); */
           gluSphere( qobj, 1.0,  12, 6);
           gluDeleteQuadric(qobj);
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        fine_sphere   = glGenLists(1);
        if(fine_sphere != 0) {
           glNewList(fine_sphere, GL_COMPILE);
           qobj = gluNewQuadric();
/*           gluSphere( qobj, 1.0,  30, 30); */
           gluSphere( qobj, 1.0,  30, 15);
           gluDeleteQuadric(qobj);
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

}


void generate_stick_object(void)
{
        GLUquadricObj *qobj;
        int nfacets, i;



        nfacets = 12;
        if(stick_object != 0) {
           glNewList(stick_object, GL_COMPILE);
           qobj = gluNewQuadric();
           gluCylinder(qobj, 0.55 * bond_factor, 0.55 * bond_factor, 1.0, nfacets, 1);
           gluDeleteQuadric(qobj);
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        if(stick_object1 != 0) {
           glNewList(stick_object1, GL_COMPILE);
           qobj = gluNewQuadric();
           gluCylinder(qobj, 0.5 * bond_factor, 0.5 * bond_factor, 0.5, nfacets, 1);
           gluDeleteQuadric(qobj);
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        if(stick_object2 != 0) {
           glNewList(stick_object2, GL_COMPILE);
           glPushMatrix();
           glTranslatef(0.0, 0.0, 0.5);
           glCallList(stick_object1);
           glPopMatrix();
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        if(stick_object_thin != 0) {
           glNewList(stick_object_thin, GL_COMPILE);
           qobj = gluNewQuadric();
           gluCylinder(qobj, 0.3 * bond_factor, 0.3 * bond_factor, 1.0, nfacets, 1);
           gluDeleteQuadric(qobj);
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        if(dot_stick_object != 0) {
           glNewList(dot_stick_object, GL_COMPILE);
           glPushMatrix();
           glTranslatef(0.0, 0.0, 0.1);
           for(i=0; i<5; i++) {
              qobj = gluNewQuadric();
              gluSphere( qobj, 0.4*bond_factor,  12, 6);
              gluDeleteQuadric(qobj);
              glTranslatef(0.0, 0.0, 0.2);
           }
           glPopMatrix();
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        if(cone_object != 0) {
           glNewList(cone_object, GL_COMPILE);
           qobj = gluNewQuadric();
           gluCylinder(qobj, 0.6 * bond_factor, 0.0 * bond_factor, 1.0, nfacets, 1);
           gluDeleteQuadric(qobj);
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        nfacets =  6;
        if(rough_stick_object != 0) {
           glNewList(rough_stick_object, GL_COMPILE);
           qobj = gluNewQuadric();
           gluCylinder(qobj, 0.55 * bond_factor, 0.55 * bond_factor, 1.0, nfacets, 1);
           gluDeleteQuadric(qobj);
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        if(rough_stick_object1 != 0) {
           glNewList(rough_stick_object1, GL_COMPILE);
           qobj = gluNewQuadric();
           gluCylinder(qobj, 0.5 * bond_factor, 0.5 * bond_factor, 0.5, nfacets, 1);
           gluDeleteQuadric(qobj);
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        if(rough_stick_object2 != 0) {
           glNewList(rough_stick_object2, GL_COMPILE);
           glPushMatrix();
           glTranslatef(0.0, 0.0, 0.5);
           glCallList(rough_stick_object1);
           glPopMatrix();
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        if(rough_stick_object_thin != 0) {
           glNewList(rough_stick_object_thin, GL_COMPILE);
           qobj = gluNewQuadric();
           gluCylinder(qobj, 0.3 * bond_factor, 0.3 * bond_factor, 1.0, nfacets, 1);
           gluDeleteQuadric(qobj);
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        if(rough_cone_object!= 0) {
           glNewList(rough_cone_object, GL_COMPILE);
           qobj = gluNewQuadric();
           gluCylinder(qobj, 0.6 * bond_factor, 0.0 * bond_factor, 1.0, nfacets, 1);
           gluDeleteQuadric(qobj);
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }
}


/*------------------------------------------------------------*/
void generate_stick_object_orig(void)
/* original stick objects, currently not in use */
{
        register short i;
        int nfacets;
        float konst, vect[3], norm[3];

        if(bit.high_qual) nfacets = 30;
        else nfacets =  5;

        konst = 2.0*M_PI/nfacets;

        if(stick_object != 0) {
           glNewList(stick_object, GL_COMPILE);

           glBegin(GL_TRIANGLE_STRIP);
           for(i=0; i<=nfacets; i++){
                norm[0] = sin(konst*i);
                norm[1] = cos(konst*i);
                norm[2] = 0;
                glNormal3fv(norm);
                vect[0] = 0.5 * bond_factor * norm[0];
                vect[1] = 0.5 * bond_factor * norm[1];
                vect[2] = 0;
                glVertex3fv(vect);
                vect[2] = 1;
                glVertex3fv(vect);
           }
           glEnd();
           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        if(stick_object1 != 0) {
           glNewList(stick_object1, GL_COMPILE);

           glBegin(GL_TRIANGLE_STRIP);
           for(i=0; i<=nfacets; i++){
                norm[0] = sin(konst*i);
                norm[1] = cos(konst*i);
                norm[2] = 0;
                glNormal3fv(norm);
                vect[0] = 0.5 * bond_factor * norm[0];
                vect[1] = 0.5 * bond_factor * norm[1];
                vect[2] = 0;
                glVertex3fv(vect);
                vect[2] = 0.5;
                glVertex3fv(vect);
           }
           glEnd();

           glEndList();
        }
        else {
           printf("can't get new listindex");
        }

        if(stick_object2 != 0) {
           glNewList(stick_object2, GL_COMPILE);

           glBegin(GL_TRIANGLE_STRIP);
           for(i=0; i<=nfacets; i++){
                norm[0] = sin(konst*i);
                norm[1] = cos(konst*i);
                norm[2] = 0;
                glNormal3fv(norm);
                vect[0] = 0.5 * bond_factor * norm[0];
                vect[1] = 0.5 * bond_factor * norm[1];
                vect[2] = 0.5;
                glVertex3fv(vect);
                vect[2] = 1;
                glVertex3fv(vect);
           }
           glEnd();

           glEndList();
        }
        else {
           printf("can't get new listindex");
        }
}

/*------------------------------------------------------------*/

