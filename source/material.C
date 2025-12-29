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
#include "glutwin.h"
#include "material.h"
#include "objects.h"
#include "drawing.h"

GLfloat diffuse[][4] = {{ 0.7, 0.7, 0.7, 1.0 },    /* white */
                       { 0.7, 0.0, 0.0, 1.0 },     /* red */
                       { 0.0, 0.4, 0.0, 1.0 },     /* green */
                       { 0.0, 0.0, 0.6, 1.0 },     /* blue */
                       { 0.8, 0.8, 0.0, 1.0 },     /* yellow */
                       { 0.0, 0.7, 0.7, 1.0 },     /* cyan */
                       { 0.6, 0.0, 0.6, 1.0 },     /* magenta */
                       { 0.25, 0.25, 0.25, 1.0 },  /* grey */ 
                       { 0.3, 0.2, 0.1, 1.0 },     /* brown */
                       { 0.45, 0.40, 0.25, 1.0 },  /* gold */
                       { 0.7, 0.2, 0.1, 1.0 }};    /* red_hot */

GLfloat ambient[] = { 0.0, 0.0, 0.0, 1.0 };
GLfloat specular[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat specular_gold[] = { 1.0, 1.0, 0.6, 1.0 };
GLfloat specular_red_hot[] = { 0.8, 0.7, 0.7, 1.0 };
GLfloat emission[] = { 0.0, 0.0, 0.0, 1.0 };
GLfloat emission_red_hot[] = { 0.8, 0.3, 0.1, 1.0 };
GLfloat shininess = 100.0;
GLfloat shininess_gold = 30.0;
GLfloat shininess_red_hot = 15.0;


GLfloat diffuse_l0[] = {1.0, 1.0, 1.0, 1.0};
GLfloat specular_l0[] = {1.0, 1.0, 1.0, 1.0};
GLfloat ambient_l0[] = {0.1, 0.1, 0.1, 1.0};
GLfloat position_l0[] = {0.4, 0.4, 1.0, 0.0};

GLfloat diffuse_l1[] = {0.5, 0.5, 0.5, 1.0};
GLfloat specular_l1[] = {0.5, 0.5, 0.5, 1.0};
GLfloat ambient_l1[] = {0.1, 0.1, 0.1, 1.0};
GLfloat position_l1[] = {-1.0, 1.0, 1.0, 0.0};

GLfloat diffuse_l2[] = {0.7, 0.7, 0.7, 1.0};
GLfloat specular_l2[] = {0.7, 0.7, 0.7, 1.0};
GLfloat ambient_l2[] = {0.1, 0.1, 0.1, 1.0};
GLfloat position_l2[] = {0.0, 0.0, -1.0, 0.0};

GLfloat ambient_model[] = {0.2, 0.2, 0.2, 1.0};

GLfloat diffuse_surface[][4] = {{ 0.7, 0.7, 0.7, 1.0 },
                                { 0.8, 0.6, 0.0, 1.0 },
                                { 0.8, 0.45, 0.0, 1.0 },
                                { 0.75, 0.35, 0.0, 1.0 },
                                { 0.7, 0.25, 0.0, 1.0 },
                                { 0.65, 0.15, 0.0, 1.0 },
                                { 0.6, 0.0, 0.0, 1.0 },
                                { 0.6, 0.0, 0.15, 1.0 },
                                { 0.6, 0.0, 0.3, 1.0 },
                                { 0.6, 0.0, 0.45, 1.0 },
                                { 0.6, 0.0, 0.6, 1.0 },
                                { 0.0, 0.0, 0.6, 1.0 },
                                { 0.0, 0.2, 0.6, 1.0 },
                                { 0.0, 0.4, 0.6, 1.0 },
                                { 0.0, 0.5, 0.6, 1.0 },
                                { 0.0, 0.6, 0.6, 1.0 },
                                { 0.0, 0.6, 0.45, 1.0 },
                                { 0.0, 0.6, 0.3, 1.0 },
                                { 0.0, 0.6, 0.15, 1.0 },
                                { 0.0, 0.6, 0.0, 1.0 },
                                { 0.2, 0.6, 0.0, 1.0 },
                                { 0.4, 0.6, 0.0, 1.0 }};

GLfloat shininess_surface = 50.0;

GLubyte halftone1[] = {
0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55};

GLubyte halftone2[] = {
0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA};

GLubyte halftone3[] = {
0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC,
0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33,
0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC,
0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33,
0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC,
0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33,
0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC,
0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33,
0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC,
0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33,
0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33};

GLubyte halftone4[] = {
0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA,
0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA,
0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA, 0xAA,
0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55};

GLubyte halftone5[] = {
0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33,
0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC,
0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC,
0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33,
0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33,
0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC,
0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC,
0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33,
0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0x33, 0x33, 0x33, 0x33,
0x33, 0x33, 0x33, 0x33, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC, 0xCC,
0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33};
                      
GLuint def_mat_base;
GLuint element_mat_base;
GLuint surface_mat_base;
GLuint bond_mat;
GLuint arrow_mat;
GLuint pick_mat;
GLuint box_mat;
GLuint line_style_base;
GLubyte *halftone[4] = { halftone2, halftone3, halftone4, halftone5 };

void gen_def_materials(void)
{
   register int i;

   for(i=0; i<9; i++) {
      glNewList(def_mat_base + i, GL_COMPILE);
        glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
        glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse[i]);
        glMaterialfv(GL_FRONT, GL_EMISSION, emission);
        glMaterialf(GL_FRONT, GL_SHININESS, shininess);
      glEndList();
   }
      glNewList(def_mat_base + 9, GL_COMPILE);
        glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
        glMaterialfv(GL_FRONT, GL_SPECULAR, specular_gold);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse[9]);
        glMaterialfv(GL_FRONT, GL_EMISSION, emission);
        glMaterialf(GL_FRONT, GL_SHININESS, shininess_gold);
      glEndList();
      glNewList(def_mat_base + 10, GL_COMPILE);
        glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
        glMaterialfv(GL_FRONT, GL_SPECULAR, specular_red_hot);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse[10]);
        glMaterialfv(GL_FRONT, GL_EMISSION, emission_red_hot);
        glMaterialf(GL_FRONT, GL_SHININESS, shininess_red_hot);
      glEndList();
}

void gen_materials(void)
{
   register int i;

   def_mat_base = glGenLists(11);
   gen_def_materials();

   element_mat_base = glGenLists(100);
   for(i=0; i<99; i++){
      if(element[i].symbol != NULL) {
         glNewList(element_mat_base + i, GL_COMPILE);
            if(element[i].col < 11) 
               glCallList(def_mat_base + element[i].col);
            else glCallList(def_mat_base + 7);
         glEndList();
      }
   }

   surface_mat_base = glGenLists(44);
   for(i=0; i<22; i++) {
      glNewList(surface_mat_base + i, GL_COMPILE);
        glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
        glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_surface[i]);
        glMaterialfv(GL_FRONT, GL_EMISSION, emission);
        glMaterialf(GL_FRONT, GL_SHININESS, shininess_surface);
      glEndList();
      glNewList(surface_mat_base + 22 + (i+11)%21, GL_COMPILE);
        glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
        glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_surface[i]);
        glMaterialfv(GL_FRONT, GL_EMISSION, emission);
        glMaterialf(GL_FRONT, GL_SHININESS, shininess_surface);
      glEndList();
   }

   bond_mat = glGenLists(1);
   glNewList(bond_mat, GL_COMPILE);
      glCallList(def_mat_base + 7);
   glEndList();

   arrow_mat = glGenLists(1);
   glNewList(arrow_mat, GL_COMPILE);
      glCallList(def_mat_base + 1);
   glEndList();

   pick_mat = glGenLists(1);
   glNewList(pick_mat, GL_COMPILE);
      glCallList(def_mat_base + 10);
   glEndList();

   box_mat = glGenLists(1);
   glNewList(box_mat, GL_COMPILE);
      glCallList(def_mat_base +3);
   glEndList();

   line_style_base = glGenLists(6);
   glNewList(line_style_base, GL_COMPILE);
     glLineWidth(1.0);
     glDisable(GL_LINE_STIPPLE);
   glEndList();
   glNewList(line_style_base + 1, GL_COMPILE);
     glLineWidth(3.0);
     glLineStipple(1, 0x0f0f);
     glEnable(GL_LINE_STIPPLE);
   glEndList();
   glNewList(line_style_base + 2, GL_COMPILE);
     glLineWidth(4.0);
     glDisable(GL_LINE_STIPPLE);
   glEndList();
   glNewList(line_style_base + 3, GL_COMPILE);
     glLineWidth(1.0);
     glDisable(GL_LINE_STIPPLE);
   glEndList();
   glNewList(line_style_base + 4, GL_COMPILE);
     glLineWidth(2.0);
     glDisable(GL_LINE_STIPPLE);
   glEndList();
   glNewList(line_style_base + 5, GL_COMPILE);
     glLineWidth(1.0);
     glDisable(GL_LINE_STIPPLE);
   glEndList();
}

void scale_line(int w, int h)
{
   int win_fac;
   win_fac = w > h ? w/736 : h/736;
   glNewList(line_style_base, GL_COMPILE);
     glLineWidth(1.0*sfac*10*win_fac);
     glDisable(GL_LINE_STIPPLE);
   glEndList();
   glNewList(line_style_base + 1, GL_COMPILE);
     glLineWidth(3.0*sfac*10*win_fac);
     glLineStipple(1, 0x0f0f);
     glEnable(GL_LINE_STIPPLE);
   glEndList();
    glNewList(line_style_base + 2, GL_COMPILE);
     glLineWidth(4.0*sfac*10*win_fac);
     glDisable(GL_LINE_STIPPLE);
   glEndList();
   glNewList(line_style_base + 3, GL_COMPILE);
     glLineWidth(1.0*sfac*20*win_fac);
     glDisable(GL_LINE_STIPPLE);
   glEndList();
   glNewList(line_style_base + 4, GL_COMPILE);
     glLineWidth(2.0*sfac*10*win_fac*bond_factor/0.15);
     glDisable(GL_LINE_STIPPLE);
   glEndList();
   glNewList(line_style_base + 5, GL_COMPILE);
     glLineWidth(1.0*sfac*10*win_fac*bond_factor/0.15);
     glDisable(GL_LINE_STIPPLE);
   glEndList();
}

void set_light(void)
{
   glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse_l0);
   glLightfv(GL_LIGHT0, GL_SPECULAR,  specular_l0);
   glLightfv(GL_LIGHT0, GL_AMBIENT, ambient_l0);
   glLightfv(GL_LIGHT0, GL_POSITION, position_l0);

   glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse_l1);
   glLightfv(GL_LIGHT1, GL_SPECULAR,  specular_l1);
   glLightfv(GL_LIGHT1, GL_AMBIENT, ambient_l1);
   glLightfv(GL_LIGHT1, GL_POSITION, position_l1);

   glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuse_l2);
   glLightfv(GL_LIGHT2, GL_SPECULAR,  specular_l2);
   glLightfv(GL_LIGHT2, GL_AMBIENT, ambient_l2);
   glLightfv(GL_LIGHT2, GL_POSITION, position_l2);
  
   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_model);

   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_LIGHT1);
   glEnable(GL_LIGHT2);
}
