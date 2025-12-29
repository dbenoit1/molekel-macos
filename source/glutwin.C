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
 * generate GLUT window
*/

#include "main.h"
#include "molekel.h"
#include "manip.h"
#include "glutwin.h"
#include "menu.h"
#include "general.h"
#include "objects.h"
#include "material.h"
#include "drawing.h"
#include "browser.h"
#include "render.h"
#include "chooseinterf.h"

float x_offset = 0, y_offset = 0;
int xorig = 0, yorig = 0, xsize, ysize;
int mainwin, filewin, orbwin, freqwin;
GLdouble clipping_plane[] = {0, 0, -1, 0};
GLfloat clearcolor[4] = {0.0, 0.0, 0.0, 0.0};
GLfloat fogstart =  -0.5, fogend = 1.51;

void glutwin(int *pargc, char **argv)
{
  glutInit(pargc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
//  xsize = ysize = glutGet(GLUT_SCREEN_HEIGHT) - 8; 
//  xsize = ysize = 768 - 8;
  xsize = ysize = 736;
  glutInitWindowSize(xsize, ysize);
  glutInitWindowPosition(xorig, yorig);
  mainwin = glutCreateWindow("molekel main window");
  glutDisplayFunc(drawit);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(spec_keyboard);
  glutIdleFunc(NULL);
  make_main_menu();
  init_GL();
}

void init_GL(void)
{
    GLfloat fogcolor[4] = {0.0, 0.0, 0.0, 1.0};
/* setup OpenGL state */
    glClearDepth(1.0);
    glShadeModel(GL_SMOOTH);
    glFogi(GL_FOG_MODE, GL_LINEAR);
    glFogfv(GL_FOG_COLOR, fogcolor);
    glHint(GL_FOG_HINT, GL_DONT_CARE);
    glFogf(GL_FOG_START, fogstart);
    glFogf(GL_FOG_END, fogend);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);

//    init_persp(); seems to be handled by GLUT!

    set_light();
    gen_materials();

    generate_objects();
}


void init_persp(void)
{
   float ratio;

   ratio = (float)xsize/ysize;

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   if(bit.persp){
      gluPerspective(30.0,  ratio,  0.23, 20);
      gluLookAt(0.0 + x_offset, 0.0 + y_offset, 3.73,
             0.0 + x_offset, 0.0 + y_offset, 0.0, 0.0, 1.0, 0.0);
   }
   else {
      glOrtho(-ratio + x_offset, ratio + x_offset,
            -1.0 + y_offset, 1.0 + y_offset, -3.5, 17);

   }
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}


void clip_persp(void)
{
/* not used!
 * just used in case of no clipplanes are available, not the case in OGL
*/
   float ratio;

   ratio = (float)xsize/ysize;
   glPushMatrix();
   glMatrixMode(GL_PROJECTION);
   if(bit.persp){
      gluPerspective(30.0,  ratio,  3.73, 5.23);

      gluLookAt(0.0 + x_offset, 0.0 + y_offset, 3.73,
             0.0 + x_offset, 0.0 + y_offset, 0.0, 0.0, 1.0, 0.0);
   }
   else {
       glLoadIdentity ();
       glOrtho(-ratio + x_offset, ratio + x_offset,
            -1.0 + y_offset, 1.0 + y_offset, 0, 1.5);
   }
   glMatrixMode(GL_MODELVIEW);
   glPopMatrix();
}


void reshape(int w, int h)
{
   xsize = w;
   ysize = h;
   output_width = (xsize / 4) * 4;
   output_height = (ysize / 4) * 4;
   scale_line(xsize, ysize);
   glViewport(0, 0, (GLint) xsize, (GLint) ysize);
   init_persp();
   if(renderglui) renderglui->sync_live();
}


void dummy_gluiMouseFunc(int button, int state, int x, int y)
{
/* dummy function for GLUI mouse funtion */
/* a simple NULL argument causes a crash if the texture file is selected with the middle mouse button */
}

void dummy_gluiKeyboardFunc(unsigned char key, int x, int y)
{
}

void dummy_gluiSpecialFunc(int key, int x, int y)
{
}
