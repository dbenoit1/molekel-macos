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


#include <GLUT/glut.h>
/*#include <GLUI/glui.h>*/
#include <GL/glui.h>

#include "molekel.h"
#include "constant.h"
#include "general.h"
#include "glutwin.h"
#include "colorinterf.h"
#include "browser.h"
#include "material.h"
#include "objects.h"
#include "drawing.h"
#include "pick.h"
#include "maininterf.h"

void  vcopy(float *v1, float *v2); /* from trackball.h */

int colorinterfwin = 0;
static float ratio;

GLUI *colorinterfglui_bottom = NULL, *colorinterfglui_right;

GLfloat diffuse_live[4], specular_live[4], ambient_live[4], emission_live[4];
GLfloat shininess_live;
GLint   def_live = 0;

void colorinterfglui_cb(int key)
{
   unsigned short r, g, b;

   switch (key) {
      case DEFCOLOR:
         glutSetWindow(mainwin);
         glCallList(def_mat_base + def_live);
         glGetMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_live);
	 glGetMaterialfv(GL_FRONT, GL_SPECULAR, specular_live);
         glGetMaterialfv(GL_FRONT, GL_AMBIENT, ambient_live);
         glGetMaterialfv(GL_FRONT, GL_EMISSION, emission_live);
         glGetMaterialfv(GL_FRONT, GL_SHININESS, &shininess_live);
         if(colorinterfglui_bottom) colorinterfglui_bottom->sync_live();
         glutSetWindow(colorinterfwin);
      break;
      case ATOMCOLOR:
      case GET_ATOMCOLOR:
         if(!firstmol) {
            logprint("No molecule loaded");
            update_logs();
            return;
         }
         glutSetWindow(mainwin);
         enter_pick(key, "assign color to atom type");
      break;
      case BONDCOLOR:
         glutSetWindow(mainwin);
         glNewList(bond_mat, GL_COMPILE);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_live);
            glMaterialfv(GL_FRONT, GL_SPECULAR, specular_live);
            glMaterialfv(GL_FRONT, GL_AMBIENT, ambient_live);
            glMaterialfv(GL_FRONT, GL_EMISSION, emission_live);
            glMaterialf(GL_FRONT, GL_SHININESS, shininess_live);
         glEndList();
         glutPostRedisplay();
      break;
      case SURFCOLOR:
      case GET_SURFCOLOR:
         if(!firstmol) {
            logprint("No molecule loaded");
            update_logs();
            return;
         }
         glutSetWindow(mainwin);
         enter_pick(key, "assign color to surface");
      break;
      case ARROWCOLOR:
         glutSetWindow(mainwin);
         glNewList(arrow_mat, GL_COMPILE);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_live);
            glMaterialfv(GL_FRONT, GL_SPECULAR, specular_live);
            glMaterialfv(GL_FRONT, GL_AMBIENT, ambient_live);
            glMaterialfv(GL_FRONT, GL_EMISSION, emission_live);
            glMaterialf(GL_FRONT, GL_SHININESS, shininess_live);
         glEndList();
         glutPostRedisplay();
      break;
      case BOXCOLOR:
         glutSetWindow(mainwin);
         glNewList(box_mat, GL_COMPILE);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_live);
            glMaterialfv(GL_FRONT, GL_SPECULAR, specular_live);
            glMaterialfv(GL_FRONT, GL_AMBIENT, ambient_live);
            glMaterialfv(GL_FRONT, GL_EMISSION, emission_live);
            glMaterialf(GL_FRONT, GL_SHININESS, shininess_live);
         glEndList();
         glutPostRedisplay();
      break;
      case LABELCOLOR:
         vcopy(diffuse_live, labelcolor);
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case DISTCOLOR:
         vcopy(diffuse_live, distcolor);
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case ANGELCOLOR:
         vcopy(diffuse_live, anglecolor);
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case DIHEDCOLOR:
         vcopy(diffuse_live, dihedcolor);
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case BACKCOLOR1:
         r = 255 * diffuse_live[0];
         g = 255 * diffuse_live[1];
         b = 255 * diffuse_live[2];
         backcolor1 = 0xff000000 | b<<16 | g<<8 | r;
         glutSetWindow(colorinterfwin);
         glutPostRedisplay();
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case BACKCOLOR2:
         r = 255 * diffuse_live[0];
         g = 255 * diffuse_live[1];
         b = 255 * diffuse_live[2];
         backcolor2 = 0xff000000 | b<<16 | g<<8 | r;
         glutSetWindow(colorinterfwin);
         glutPostRedisplay();
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case ALPHABLENDING:
         glutSetWindow(colorinterfwin);
         glutPostRedisplay();
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case GET_BONDCOLOR:
         glutSetWindow(mainwin);
         glCallList(bond_mat);
         glGetMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_live);
	 glGetMaterialfv(GL_FRONT, GL_SPECULAR, specular_live);
         glGetMaterialfv(GL_FRONT, GL_AMBIENT, ambient_live);
         glGetMaterialfv(GL_FRONT, GL_EMISSION, emission_live);
         glGetMaterialfv(GL_FRONT, GL_SHININESS, &shininess_live);
         if(colorinterfglui_bottom) colorinterfglui_bottom->sync_live();
         glutSetWindow(colorinterfwin);
         glutPostRedisplay();
      break;
      case GET_ARROWCOLOR:
         glutSetWindow(mainwin);
         glCallList(arrow_mat);
         glGetMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_live);
	 glGetMaterialfv(GL_FRONT, GL_SPECULAR, specular_live);
         glGetMaterialfv(GL_FRONT, GL_AMBIENT, ambient_live);
         glGetMaterialfv(GL_FRONT, GL_EMISSION, emission_live);
         glGetMaterialfv(GL_FRONT, GL_SHININESS, &shininess_live);
         if(colorinterfglui_bottom) colorinterfglui_bottom->sync_live();
         glutSetWindow(colorinterfwin);
         glutPostRedisplay();
      break;
      case GET_BOXCOLOR:
         glutSetWindow(mainwin);
         glCallList(box_mat);
         glGetMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_live);
	 glGetMaterialfv(GL_FRONT, GL_SPECULAR, specular_live);
         glGetMaterialfv(GL_FRONT, GL_AMBIENT, ambient_live);
         glGetMaterialfv(GL_FRONT, GL_EMISSION, emission_live);
         glGetMaterialfv(GL_FRONT, GL_SHININESS, &shininess_live);
         if(colorinterfglui_bottom) colorinterfglui_bottom->sync_live();
         glutSetWindow(colorinterfwin);
         glutPostRedisplay();
      break;
      case CANCEL:
         if(colorinterfwin) {
            glutSetWindow(colorinterfwin);
            glutHideWindow();
         }
      break;
   }
}

void reshape_colorinterf(int w, int h)
{
   int tx, ty, tw, th;

   glutReshapeWindow(360, 660);
   GLUI_Master.get_viewport_area(&tx, &ty, &tw, &th);
   glViewport(tx, ty, tw, th);
   ratio = (float)tw/(float)th;
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glOrtho(-ratio, ratio, -1.0, 1.0, -3.5, 17);
   glMatrixMode(GL_MODELVIEW);
/* flexible window size
   glViewport(tx, ty, tw, th);
   ratio = (float)tw/(float)th;
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   if(ratio < 1)
      glOrtho(-1.0, 1.0, -1/ratio, 1/ratio, -3.5, 17);
   else glOrtho(-ratio, ratio, -1.0, 1.0, -3.5, 17);
   glMatrixMode(GL_MODELVIEW);
*/
}


void draw_colorinterf(void)
{
   GLUquadricObj *qobj;

   glClearColor(0.0, 0.0, 0.0, 1.0);
   glClear(GL_DEPTH_BUFFER_BIT);
   glDisable(GL_LIGHTING);
   if(backcolor1 == backcolor2) {
      set_clear_color(backcolor1);
      glClear(GL_COLOR_BUFFER_BIT);
   }
   else {
      glClear(GL_COLOR_BUFFER_BIT);
      glDepthMask(GL_FALSE);
      glBegin(GL_QUADS);
      if(ratio < 1) {
         set_color(backcolor1);
         glVertex3f(-1.0, -1.0/ratio, -3.5);
         glVertex3f( 1.0, -1.0/ratio, -3.5);
         set_color(backcolor2);
         glVertex3f( 1.0,  1.0/ratio, -3.5);
         glVertex3f(-1.0,  1.0/ratio, -3.5);
      }
      else {
         set_color(backcolor1);
         glVertex3f(-ratio, -1.0, -3.5);
         glVertex3f( ratio, -1.0, -3.5);
         set_color(backcolor2);
         glVertex3f( ratio,  1.0, -3.5);
         glVertex3f(-ratio,  1.0, -3.5);
      }
      glEnd();
      glDepthMask(GL_TRUE);
      glEnable(GL_LIGHTING);
   }
   glColor4fv(diffuse_live);
   glRectf(0.70, -0.60, 0.95, -0.85);
   glEnable(GL_LIGHTING);
   glMaterialfv(GL_FRONT, GL_SPECULAR, specular_live);
   glMaterialfv(GL_FRONT, GL_AMBIENT, ambient_live);
   glMaterialfv(GL_FRONT, GL_EMISSION, emission_live);
   glMaterialf(GL_FRONT, GL_SHININESS, shininess_live);
   if(alphablending) {
      glEnable(GL_BLEND);
      glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse[10]);
      qobj = gluNewQuadric();
      gluSphere( qobj, 0.3,  24, 24);
      gluDeleteQuadric(qobj); 
   }
   glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_live);
   qobj = gluNewQuadric();
   gluSphere( qobj, 0.8,  24, 24);
   gluDeleteQuadric(qobj); 
   if(alphablending) glDisable(GL_BLEND);
   glutSwapBuffers();
}


void colorinterf(void)
{
   if(colorinterfwin) {
      glutSetWindow(colorinterfwin);
      glutPopWindow();
      glutShowWindow();
   }
   else {
      glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
      glutInitWindowSize(360, 660);
      glutInitWindowPosition(winposition[0], winposition[1]);
      colorinterfwin = glutCreateWindow("color interface");
      glutDisplayFunc(draw_colorinterf);
      GLUI_Master.set_glutReshapeFunc(reshape_colorinterf);
      GLUI_Master.set_glutMouseFunc(dummy_gluiMouseFunc);
      GLUI_Master.set_glutKeyboardFunc(dummy_gluiKeyboardFunc);
      GLUI_Master.set_glutSpecialFunc(dummy_gluiSpecialFunc);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glShadeModel(GL_SMOOTH);
      glEnable(GL_NORMALIZE);
      glEnable(GL_DEPTH_TEST);
      set_light();
      colorinterfglui_cb(DEFCOLOR);

#ifdef WIN32
/* on WIN32 order of drawing the subwindows needs to be changed */
      colorinterfglui_right = GLUI_Master.create_glui_subwindow(colorinterfwin, GLUI_SUBWINDOW_RIGHT);
      GLUI_Panel *assign_mat_panel = colorinterfglui_right->add_panel("assign material to");
      colorinterfglui_right->add_statictext_to_panel(assign_mat_panel, "");
      GLUI_Button *atomtype = colorinterfglui_right->add_button_to_panel
                 (assign_mat_panel, "atom type", ATOMCOLOR, colorinterfglui_cb);
      GLUI_Button *bond = colorinterfglui_right->add_button_to_panel\
                 (assign_mat_panel, "bond", BONDCOLOR, colorinterfglui_cb);
      GLUI_Button *surface = colorinterfglui_right->add_button_to_panel\
                 (assign_mat_panel, "surface", SURFCOLOR, colorinterfglui_cb);
      GLUI_Button *arrow = colorinterfglui_right->add_button_to_panel\
                 (assign_mat_panel, "arrow", ARROWCOLOR, colorinterfglui_cb);
      GLUI_Button *box = colorinterfglui_right->add_button_to_panel\
                 (assign_mat_panel, "box", BOXCOLOR, colorinterfglui_cb);
      colorinterfglui_right->add_statictext("");
      GLUI_Panel *assign_col_panel = colorinterfglui_right->add_panel("assign diffuse color to");
      colorinterfglui_right->add_statictext_to_panel(assign_col_panel, "");
      colorinterfglui_right->add_statictext_to_panel(assign_col_panel, "background");
      GLUI_Button *back_top = colorinterfglui_right->add_button_to_panel\
                 (assign_col_panel, "top", BACKCOLOR2, colorinterfglui_cb);
      GLUI_Button *back_bot = colorinterfglui_right->add_button_to_panel\
                 (assign_col_panel, "bottom", BACKCOLOR1, colorinterfglui_cb);
      colorinterfglui_right->add_statictext_to_panel(assign_col_panel, "");
      colorinterfglui_right->add_statictext_to_panel(assign_col_panel, "labels");
      GLUI_Button *labelc = colorinterfglui_right->add_button_to_panel\
                 (assign_col_panel, "label", LABELCOLOR, colorinterfglui_cb);
      GLUI_Button *distc = colorinterfglui_right->add_button_to_panel\
                 (assign_col_panel, "distance", DISTCOLOR, colorinterfglui_cb);
      GLUI_Button *angelc = colorinterfglui_right->add_button_to_panel\
                 (assign_col_panel, "angel", ANGELCOLOR, colorinterfglui_cb);
      GLUI_Button *dihedc = colorinterfglui_right->add_button_to_panel\
                 (assign_col_panel, "dihedral", DIHEDCOLOR, colorinterfglui_cb);
      colorinterfglui_right->add_statictext("");
      GLUI_Checkbox *ablending = colorinterfglui_right->add_checkbox\
             ("enable alpha blending", &alphablending, ALPHABLENDING, colorinterfglui_cb);
      colorinterfglui_right->add_statictext("      (surface only)");
      colorinterfglui_right->add_statictext("");
      GLUI_Panel *get_mat_panel = colorinterfglui_right->add_panel("get material from");
      colorinterfglui_right->add_statictext_to_panel(get_mat_panel, "");
      GLUI_Button *g_atomtype = colorinterfglui_right->add_button_to_panel
                 (get_mat_panel, "atom type", GET_ATOMCOLOR, colorinterfglui_cb);
      GLUI_Button *g_bond = colorinterfglui_right->add_button_to_panel\
                 (get_mat_panel, "bond", GET_BONDCOLOR, colorinterfglui_cb);
      GLUI_Button *g_surface = colorinterfglui_right->add_button_to_panel\
                 (get_mat_panel, "surface", GET_SURFCOLOR, colorinterfglui_cb);
      GLUI_Button *g_arrow = colorinterfglui_right->add_button_to_panel\
                 (get_mat_panel, "arrow", GET_ARROWCOLOR, colorinterfglui_cb);
      GLUI_Button *g_box = colorinterfglui_right->add_button_to_panel\
                 (get_mat_panel, "box", GET_BOXCOLOR, colorinterfglui_cb);
#endif

      colorinterfglui_bottom = GLUI_Master.create_glui_subwindow(colorinterfwin, GLUI_SUBWINDOW_BOTTOM);
      colorinterfglui_bottom->set_main_gfx_window(colorinterfwin);
      GLUI_Panel *material_panel = colorinterfglui_bottom->add_panel("materials");
      GLUI_Spinner *defcolor_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (material_panel, "default materials", GLUI_SPINNER_INT, &def_live, DEFCOLOR, colorinterfglui_cb);
      defcolor_spinner->set_int_limits(0, 10);
      defcolor_spinner->set_alignment(GLUI_ALIGN_CENTER);
      GLUI_Panel *mat_panel = colorinterfglui_bottom->add_panel_to_panel(material_panel, "");
      colorinterfglui_bottom->add_statictext_to_panel(mat_panel, "diffuse color");
      GLUI_Spinner *diff_r_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (mat_panel, "r", GLUI_SPINNER_FLOAT, &(diffuse_live[0]));
      diff_r_spinner->set_float_limits(0, 1.0);
      GLUI_Spinner *diff_g_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (mat_panel, "g", GLUI_SPINNER_FLOAT, &(diffuse_live[1]));
      diff_g_spinner->set_float_limits(0, 1.0);
      GLUI_Spinner *diff_b_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (mat_panel, "b", GLUI_SPINNER_FLOAT, &(diffuse_live[2]));
      diff_b_spinner->set_float_limits(0, 1.0);
      colorinterfglui_bottom->add_statictext_to_panel(mat_panel, "specular color");
      GLUI_Spinner *specular_r_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (mat_panel, "r", GLUI_SPINNER_FLOAT, &(specular_live[0]));
      specular_r_spinner->set_float_limits(0, 1.0);
      GLUI_Spinner *specular_g_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (mat_panel, "g", GLUI_SPINNER_FLOAT, &(specular_live[1]));
      specular_g_spinner->set_float_limits(0, 1.0);
      GLUI_Spinner *specular_b_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (mat_panel, "b", GLUI_SPINNER_FLOAT, &(specular_live[2]));
      specular_b_spinner->set_float_limits(0, 1.0);
      colorinterfglui_bottom->add_statictext_to_panel(mat_panel, "ambient color");
      GLUI_Spinner *ambient_r_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (mat_panel, "r", GLUI_SPINNER_FLOAT, &(ambient_live[0]));
      ambient_r_spinner->set_float_limits(0, 1.0);
      GLUI_Spinner *ambient_g_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (mat_panel, "g", GLUI_SPINNER_FLOAT, &(ambient_live[1]));
      ambient_g_spinner->set_float_limits(0, 1.0);
      GLUI_Spinner *ambient_b_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (mat_panel, "b", GLUI_SPINNER_FLOAT, &(ambient_live[2]));
      ambient_b_spinner->set_float_limits(0, 1.0);
      colorinterfglui_bottom->add_statictext_to_panel(mat_panel, "emission color");
      GLUI_Spinner *emission_r_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (mat_panel, "r", GLUI_SPINNER_FLOAT, &(emission_live[0]));
      emission_r_spinner->set_float_limits(0, 1.0);
      GLUI_Spinner *emission_g_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (mat_panel, "g", GLUI_SPINNER_FLOAT, &(emission_live[1]));
      emission_g_spinner->set_float_limits(0, 1.0);
      GLUI_Spinner *emission_b_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (mat_panel, "b", GLUI_SPINNER_FLOAT, &(emission_live[2]));
      emission_b_spinner->set_float_limits(0, 1.0);
      colorinterfglui_bottom->add_statictext_to_panel(mat_panel, "");
      GLUI_Spinner *shininess_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (mat_panel, "shininess", GLUI_SPINNER_FLOAT, &shininess_live);
      shininess_spinner->set_float_limits(0, 100.0);
      GLUI_Spinner *alpha_spinner = colorinterfglui_bottom->add_spinner_to_panel\
          (mat_panel, "alpha value", GLUI_SPINNER_FLOAT, &(diffuse_live[3]));
      alpha_spinner->set_float_limits(0, 1.0);
      
      colorinterfglui_bottom->add_statictext("");
      GLUI_Button *cancel = colorinterfglui_bottom->add_button("cancel", CANCEL, colorinterfglui_cb);
      cancel->set_alignment(GLUI_ALIGN_CENTER);

#ifndef WIN32
      colorinterfglui_right = GLUI_Master.create_glui_subwindow(colorinterfwin, GLUI_SUBWINDOW_RIGHT);
      GLUI_Panel *assign_mat_panel = colorinterfglui_right->add_panel("assign material to");
      colorinterfglui_right->add_statictext_to_panel(assign_mat_panel, "");
      GLUI_Button *atomtype = colorinterfglui_right->add_button_to_panel
                 (assign_mat_panel, "atom type", ATOMCOLOR, colorinterfglui_cb);
      GLUI_Button *bond = colorinterfglui_right->add_button_to_panel\
                 (assign_mat_panel, "bond", BONDCOLOR, colorinterfglui_cb);
      GLUI_Button *surface = colorinterfglui_right->add_button_to_panel\
                 (assign_mat_panel, "surface", SURFCOLOR, colorinterfglui_cb);
      GLUI_Button *arrow = colorinterfglui_right->add_button_to_panel\
                 (assign_mat_panel, "arrow", ARROWCOLOR, colorinterfglui_cb);
      GLUI_Button *box = colorinterfglui_right->add_button_to_panel\
                 (assign_mat_panel, "box", BOXCOLOR, colorinterfglui_cb);
      colorinterfglui_right->add_statictext("");
      GLUI_Panel *assign_col_panel = colorinterfglui_right->add_panel("assign diffuse color to");
      colorinterfglui_right->add_statictext_to_panel(assign_col_panel, "");
      colorinterfglui_right->add_statictext_to_panel(assign_col_panel, "background");
      GLUI_Button *back_top = colorinterfglui_right->add_button_to_panel\
                 (assign_col_panel, "top", BACKCOLOR2, colorinterfglui_cb);
      GLUI_Button *back_bot = colorinterfglui_right->add_button_to_panel\
                 (assign_col_panel, "bottom", BACKCOLOR1, colorinterfglui_cb);
      colorinterfglui_right->add_statictext_to_panel(assign_col_panel, "");
      colorinterfglui_right->add_statictext_to_panel(assign_col_panel, "labels");
      GLUI_Button *labelc = colorinterfglui_right->add_button_to_panel\
                 (assign_col_panel, "label", LABELCOLOR, colorinterfglui_cb);
      GLUI_Button *distc = colorinterfglui_right->add_button_to_panel\
                 (assign_col_panel, "distance", DISTCOLOR, colorinterfglui_cb);
      GLUI_Button *angelc = colorinterfglui_right->add_button_to_panel\
                 (assign_col_panel, "angel", ANGELCOLOR, colorinterfglui_cb);
      GLUI_Button *dihedc = colorinterfglui_right->add_button_to_panel\
                 (assign_col_panel, "dihedral", DIHEDCOLOR, colorinterfglui_cb);
      colorinterfglui_right->add_statictext("");
      GLUI_Checkbox *ablending = colorinterfglui_right->add_checkbox\
             ("enable alpha blending", &alphablending, ALPHABLENDING, colorinterfglui_cb);
      colorinterfglui_right->add_statictext("      (surface only)");
      colorinterfglui_right->add_statictext("");
      GLUI_Panel *get_mat_panel = colorinterfglui_right->add_panel("get material from");
      colorinterfglui_right->add_statictext_to_panel(get_mat_panel, "");
      GLUI_Button *g_atomtype = colorinterfglui_right->add_button_to_panel
                 (get_mat_panel, "atom type", GET_ATOMCOLOR, colorinterfglui_cb);
      GLUI_Button *g_bond = colorinterfglui_right->add_button_to_panel\
                 (get_mat_panel, "bond", GET_BONDCOLOR, colorinterfglui_cb);
      GLUI_Button *g_surface = colorinterfglui_right->add_button_to_panel\
                 (get_mat_panel, "surface", GET_SURFCOLOR, colorinterfglui_cb);
      GLUI_Button *g_arrow = colorinterfglui_right->add_button_to_panel\
                 (get_mat_panel, "arrow", GET_ARROWCOLOR, colorinterfglui_cb);
      GLUI_Button *g_box = colorinterfglui_right->add_button_to_panel\
                 (get_mat_panel, "box", GET_BOXCOLOR, colorinterfglui_cb);
#endif
   }
}

void assign_color_to_atom(AtoM *ap)
{
   glutSetWindow(mainwin);
   glNewList(element_mat_base + ap->ord, GL_COMPILE);
      glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_live);
      glMaterialfv(GL_FRONT, GL_SPECULAR, specular_live);
      glMaterialfv(GL_FRONT, GL_AMBIENT, ambient_live);
      glMaterialfv(GL_FRONT, GL_EMISSION, emission_live);
      glMaterialf(GL_FRONT, GL_SHININESS, shininess_live);
   glEndList();
}

void assign_color_to_actualsurf(void)
{
   glutSetWindow(mainwin);
   glNewList(actualsurf->matindex, GL_COMPILE);
      glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_live);
      glMaterialfv(GL_FRONT, GL_SPECULAR, specular_live);
      glMaterialfv(GL_FRONT, GL_AMBIENT, ambient_live);
      glMaterialfv(GL_FRONT, GL_EMISSION, emission_live);
      glMaterialf(GL_FRONT, GL_SHININESS, shininess_live);
   glEndList();
   if(diffuse_live[3] < 0.999) actualsurf->alpha = diffuse_live[3]*254 + 1;
   else actualsurf->alpha = 0;
}

void get_color_from_atom(AtoM *ap)
{
   glutSetWindow(mainwin);
   glutPostRedisplay();
   glCallList(element_mat_base + ap->ord);
   glGetMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_live);
   glGetMaterialfv(GL_FRONT, GL_SPECULAR, specular_live);
   glGetMaterialfv(GL_FRONT, GL_AMBIENT, ambient_live);
   glGetMaterialfv(GL_FRONT, GL_EMISSION, emission_live);
   glGetMaterialfv(GL_FRONT, GL_SHININESS, &shininess_live);
   if(colorinterfglui_bottom) colorinterfglui_bottom->sync_live();
   glutSetWindow(colorinterfwin);
   glutPostRedisplay();
   glutSetWindow(mainwin);
}

void get_color_from_actualsurf(void)
{
   glutSetWindow(mainwin);
   glutPostRedisplay();
   glCallList(actualsurf->matindex);
   glGetMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_live);
   glGetMaterialfv(GL_FRONT, GL_SPECULAR, specular_live);
   glGetMaterialfv(GL_FRONT, GL_AMBIENT, ambient_live);
   glGetMaterialfv(GL_FRONT, GL_EMISSION, emission_live);
   glGetMaterialfv(GL_FRONT, GL_SHININESS, &shininess_live);
   if(colorinterfglui_bottom) colorinterfglui_bottom->sync_live();
   glutSetWindow(colorinterfwin);
   glutPostRedisplay();
   glutSetWindow(mainwin);
}

