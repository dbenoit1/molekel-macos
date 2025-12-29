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
#include "textureinterf.h"
#include "browser.h"
#include "texture.h"
#include "maininterf.h"
#include "pick.h"

int textureinterfwin = 0;
static float ratio;

int texprop_live[4];
int texenv_live = GL_MODULATE;

GLUI *textureinterfglui_bottom = NULL, *textureinterfglui_right;

void textureinterfglui_cb(int key)
{
   switch(key) {
      case TEXTURE_UP:
         texture_up();
         glutSetWindow(textureinterfwin);
         glutPostRedisplay();
      break;
      case TEXTURE_DOWN:
         texture_down();
         glutSetWindow(textureinterfwin);
         glutPostRedisplay();
      break;
      case LOAD_TEXTURE:
         file_select(key, rgb_ext);
      break;
      case DUPLEX_TEXTURE:
         glutSetWindow(mainwin);
         actual_texture = add_texture(actual_texture->image, actual_texture->idepth,
                           actual_texture->iwidth, actual_texture->iheight, actual_texture->texprop);
         glutSetWindow(textureinterfwin);
         glutPostRedisplay();
      break;
      case SURF_TEXTURE:
         if(!firstmol) {
            logprint("No molecule loaded");
            update_logs();
            return;
         }
         glutSetWindow(mainwin);
         enter_pick(key, "assign texture to surface");
      break;
      case SURF_REFLECT:
         if(!firstmol) {
            logprint("No molecule loaded");
            update_logs();
            return;
         }
         glutSetWindow(mainwin);
         enter_pick(key, "assign texture to surface");
      break;
      case SURF_PHONG:
         if(!firstmol) {
            logprint("No molecule loaded");
            update_logs();
            return;
         }
         glutSetWindow(mainwin);
         enter_pick(key, "assign texture to surface");
      break;
      case BOND_TEXTURE:
          if(assign_no_tex) bond_texture = NULL;
          else {
             glutSetWindow(mainwin);
             mod_texture((GLenum *)texprop_live);
             bond_texture = actual_texture;
             bond_texenv = (GLenum)texenv_live;
             bond_textype = TEX_MAP;
             glutSetWindow(textureinterfwin);
             glutPostRedisplay();
          }
          glutSetWindow(mainwin);
          glutPostRedisplay();
      break;
      case BOND_REFLECT:
          if(assign_no_tex) bond_texture = NULL;
          else {
             glutSetWindow(mainwin);
             mod_texture((GLenum *)texprop_live);
             bond_texture = actual_texture;
             bond_texenv = (GLenum)texenv_live;
             bond_textype = TEX_REFLECT;
             glutSetWindow(textureinterfwin);
             glutPostRedisplay();
          }
          glutSetWindow(mainwin);
          glutPostRedisplay();
      break;
      case ATOM_TEXTURE:
         if(!firstmol) {
            logprint("No molecule loaded");
            update_logs();
            return;
         }
         glutSetWindow(mainwin);
         enter_pick(key, "assign texture to atom type");
      break;
      case ATOM_REFLECT:
         if(!firstmol) {
            logprint("No molecule loaded");
            update_logs();
            return;
         }
         glutSetWindow(mainwin);
         enter_pick(key, "assign texture to atom type");
      break;
      case BACK_TEXTURE:
          if(assign_no_tex) back_texture = NULL;
          else {
             glutSetWindow(mainwin);
             mod_texture((GLenum *)texprop_live);
             back_texture = actual_texture;
             back_texenv = (GLenum)texenv_live;
             back_textype = TEX_REFLECT;
             glutSetWindow(textureinterfwin);
             glutPostRedisplay();
          }
          glutSetWindow(mainwin);
          glutPostRedisplay();
      break;
      case CANCEL:
         if(textureinterfwin) {
            glutSetWindow(textureinterfwin);
            glutHideWindow();
         }
      break;
   }
}

void reshape_textureinterf(int w, int h)
{
   int tx, ty, tw, th;

   glutReshapeWindow(400, 454);
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


void draw_textureinterf(void)
{
   glClearColor((float)200/255, (float)200/255, (float)200/255, 1.0);
   glClear(GL_COLOR_BUFFER_BIT);
   glColor4f(1.0, 1.0, 1.0, 1.0);
   glEnable(GL_TEXTURE_2D);
/* read_sgi_add_alpha() returns just idepth 2 and 4 */
   gluBuild2DMipmaps(GL_TEXTURE_2D,  actual_texture->idepth,
                      actual_texture->iwidth, actual_texture->iheight,
		      actual_texture->idepth == 4 ? GL_RGBA : GL_LUMINANCE_ALPHA,
                      GL_UNSIGNED_BYTE, actual_texture->image);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, actual_texture->texprop[0]);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, actual_texture->texprop[1]);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, actual_texture->texprop[2]);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, actual_texture->texprop[3]);
   glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
   texprop_live[0] = actual_texture->texprop[0];
   texprop_live[1] = actual_texture->texprop[1];
   texprop_live[2] = actual_texture->texprop[2];
   texprop_live[3] = actual_texture->texprop[3];
//   texenv_live = GL_MODULATE;
   glBegin(GL_QUADS);
   glTexCoord2f(0.0, 0.0);
   glVertex3f(-1.0, -1.0, -3.5);
   glTexCoord2f(1.0, 0.0);
   glVertex3f( 1.0, -1.0, -3.5);
   glTexCoord2f(1.0, 1.0);
   glVertex3f( 1.0,  1.0, -3.5);
   glTexCoord2f(0.0, 1.0);
   glVertex3f(-1.0,  1.0, -3.5);
   glEnd();
   glutSwapBuffers();
   textureinterfglui_bottom->sync_live();
}

void dummy_mouse(int button, int state, int x, int y)
{
/* dummy function for GLUI mouse funtion */
/* a simple NULL argument causes a crash if the texture file is selected with the middle mouse button */
}

void textureinterf(void)
{
   if (textureinterfwin) {
      glutSetWindow(textureinterfwin);
      glutPopWindow();
      glutShowWindow();
   }
   else {
      glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
      glutInitWindowSize(400, 454);
      glutInitWindowPosition(winposition[0], winposition[1]);
      textureinterfwin = glutCreateWindow("texture interface");
      glutDisplayFunc(draw_textureinterf);
      GLUI_Master.set_glutReshapeFunc(reshape_textureinterf);
      GLUI_Master.set_glutMouseFunc(dummy_gluiMouseFunc);
      GLUI_Master.set_glutKeyboardFunc(dummy_gluiKeyboardFunc);
      GLUI_Master.set_glutSpecialFunc(dummy_gluiSpecialFunc);
      if(!firsttexture) {
         glutSetWindow(mainwin);
         init_texture();
         glutSetWindow(textureinterfwin);
      }
      glutSetWindow(mainwin);
      add_contour_textures();
      glutSetWindow(textureinterfwin);
      texprop_live[0] = actual_texture->texprop[0];
      texprop_live[1] = actual_texture->texprop[1];
      texprop_live[2] = actual_texture->texprop[2];
      texprop_live[3] = actual_texture->texprop[3];
//      texenv_live = GL_MODULATE;
#ifdef WIN32
/* on WIN32 order of drawing the subwindows needs to be changed */
      textureinterfglui_right = GLUI_Master.create_glui_subwindow(textureinterfwin, GLUI_SUBWINDOW_RIGHT);
      GLUI_Button *up = textureinterfglui_right->add_button("tex list up", TEXTURE_UP, textureinterfglui_cb);
      GLUI_Button *down = textureinterfglui_right->add_button("tex list down", TEXTURE_DOWN, textureinterfglui_cb);
      textureinterfglui_right->add_statictext("");
      GLUI_Button *loadrgb = textureinterfglui_right->add_button("load rgb", LOAD_TEXTURE, textureinterfglui_cb);
      GLUI_Button *dupl = textureinterfglui_right->add_button("duplex texture", DUPLEX_TEXTURE, textureinterfglui_cb);
      textureinterfglui_right->add_statictext("");
      GLUI_Panel *add_tex_panel = textureinterfglui_right->add_panel("add texture to");
      GLUI_Button *add_surf = textureinterfglui_right->add_button_to_panel\
                  (add_tex_panel, "surface", SURF_TEXTURE, textureinterfglui_cb);
      GLUI_Button *add_bond = textureinterfglui_right->add_button_to_panel\
                  (add_tex_panel, "bond", BOND_TEXTURE, textureinterfglui_cb);
      GLUI_Button *add_atom = textureinterfglui_right->add_button_to_panel\
                  (add_tex_panel, "atom type", ATOM_TEXTURE, textureinterfglui_cb);
      GLUI_Button *bg = textureinterfglui_right->add_button_to_panel\
                  (add_tex_panel, "background", BACK_TEXTURE, textureinterfglui_cb);
      textureinterfglui_right->add_statictext("");
      GLUI_Panel *ref_tex_panel = textureinterfglui_right->add_panel("reflect texture from");
      GLUI_Button *ref_surf = textureinterfglui_right->add_button_to_panel\
                  (ref_tex_panel, "surface", SURF_REFLECT, textureinterfglui_cb);
      GLUI_Button *ref_bond = textureinterfglui_right->add_button_to_panel\
                  (ref_tex_panel, "bond", BOND_REFLECT, textureinterfglui_cb);
      GLUI_Button *ref_atom = textureinterfglui_right->add_button_to_panel\
                  (ref_tex_panel, "atom type", ATOM_REFLECT, textureinterfglui_cb);
      GLUI_Button *phong = textureinterfglui_right->add_button_to_panel\
                  (ref_tex_panel, "surface phong", SURF_PHONG, textureinterfglui_cb);
      textureinterfglui_right->add_statictext("");
      GLUI_Spinner *start_spinner = textureinterfglui_right->add_spinner\
                  ("start texture", GLUI_SPINNER_FLOAT, &start_texture);
      start_spinner->set_float_limits(0, 1.0);
      GLUI_Spinner *end_spinner = textureinterfglui_right->add_spinner\
                  ("end texture", GLUI_SPINNER_FLOAT, &end_texture);
      end_spinner->set_float_limits(0, 1.0);
#endif

      textureinterfglui_bottom = GLUI_Master.create_glui_subwindow(textureinterfwin, GLUI_SUBWINDOW_BOTTOM);
      GLUI_Checkbox *no_tex = textureinterfglui_bottom->add_checkbox("assign no texture", &assign_no_tex);
      textureinterfglui_bottom->add_statictext("");
      GLUI_Panel *texprop_panel = textureinterfglui_bottom->add_panel("texture properties");
      GLUI_Listbox *mag = textureinterfglui_bottom->add_listbox_to_panel\
                  (texprop_panel, "MAG_FILTER", &(texprop_live[2]));
      mag->add_item(GL_NEAREST, "GL_NEAREST");
      mag->add_item(GL_LINEAR, "GL_LINEAR");
      mag->set_int_val(texprop_live[2]);
      GLUI_Listbox *min = textureinterfglui_bottom->add_listbox_to_panel\
                  (texprop_panel, "MIN_FILTER", &(texprop_live[3]));
      min->add_item(GL_NEAREST, "GL_NEAREST");
      min->add_item(GL_LINEAR, "GL_LINEAR");
      min->add_item(GL_NEAREST_MIPMAP_NEAREST, "GL_N_MIPMAP_N");
      min->add_item(GL_NEAREST_MIPMAP_LINEAR, "GL_N_MIPMAP_L");
      min->add_item(GL_LINEAR_MIPMAP_NEAREST, "GL_L_MIPMAP_N");
      min->add_item(GL_LINEAR_MIPMAP_LINEAR, "GL_L_MIPMAP_L");
      min->set_int_val(texprop_live[3]);
      GLUI_Listbox *wraps = textureinterfglui_bottom->add_listbox_to_panel
                  (texprop_panel, "WRAP_S", &(texprop_live[0]));
      wraps->add_item(GL_CLAMP, "GL_CLAMP");
      wraps->add_item(GL_REPEAT, "GL_REPEAT");
      wraps->set_int_val(texprop_live[0]);
      GLUI_Listbox *wrapt = textureinterfglui_bottom->add_listbox_to_panel
                  (texprop_panel, "WRAP_T", &(texprop_live[1]));
      wrapt->add_item(GL_CLAMP, "GL_CLAMP");
      wrapt->add_item(GL_REPEAT, "GL_REPEAT");
      wrapt->set_int_val(texprop_live[1]);
      GLUI_Listbox *env = textureinterfglui_bottom->add_listbox_to_panel
                  (texprop_panel, "ENV_MOD", &texenv_live);
      env->add_item(GL_DECAL, "GL_DECAL");
      env->add_item(GL_REPLACE, "GL_REPLACE");
      env->add_item(GL_MODULATE, "GL_MODULATE");
      env->add_item(GL_BLEND, "GL_BLEND");
      env->set_int_val(texenv_live);
      textureinterfglui_bottom->add_statictext("");
      GLUI_Button *cancel = textureinterfglui_bottom->add_button("cancel", CANCEL, textureinterfglui_cb);
      cancel->set_alignment(GLUI_ALIGN_CENTER);

#ifndef WIN32
      textureinterfglui_right = GLUI_Master.create_glui_subwindow(textureinterfwin, GLUI_SUBWINDOW_RIGHT);
      GLUI_Button *up = textureinterfglui_right->add_button("tex list up", TEXTURE_UP, textureinterfglui_cb);
      GLUI_Button *down = textureinterfglui_right->add_button("tex list down", TEXTURE_DOWN, textureinterfglui_cb);
      textureinterfglui_right->add_statictext("");
      GLUI_Button *loadrgb = textureinterfglui_right->add_button("load rgb", LOAD_TEXTURE, textureinterfglui_cb);
      GLUI_Button *dupl = textureinterfglui_right->add_button("duplex texture", DUPLEX_TEXTURE, textureinterfglui_cb);
      textureinterfglui_right->add_statictext("");
      GLUI_Panel *add_tex_panel = textureinterfglui_right->add_panel("add texture to");
      GLUI_Button *add_surf = textureinterfglui_right->add_button_to_panel\
                  (add_tex_panel, "surface", SURF_TEXTURE, textureinterfglui_cb);
      GLUI_Button *add_bond = textureinterfglui_right->add_button_to_panel\
                  (add_tex_panel, "bond", BOND_TEXTURE, textureinterfglui_cb);
      GLUI_Button *add_atom = textureinterfglui_right->add_button_to_panel\
                  (add_tex_panel, "atom type", ATOM_TEXTURE, textureinterfglui_cb);
      GLUI_Button *bg = textureinterfglui_right->add_button_to_panel\
                  (add_tex_panel, "background", BACK_TEXTURE, textureinterfglui_cb);
      textureinterfglui_right->add_statictext("");
      GLUI_Panel *ref_tex_panel = textureinterfglui_right->add_panel("reflect texture from");
      GLUI_Button *ref_surf = textureinterfglui_right->add_button_to_panel\
                  (ref_tex_panel, "surface", SURF_REFLECT, textureinterfglui_cb);
      GLUI_Button *ref_bond = textureinterfglui_right->add_button_to_panel\
                  (ref_tex_panel, "bond", BOND_REFLECT, textureinterfglui_cb);
      GLUI_Button *ref_atom = textureinterfglui_right->add_button_to_panel\
                  (ref_tex_panel, "atom type", ATOM_REFLECT, textureinterfglui_cb);
      GLUI_Button *phong = textureinterfglui_right->add_button_to_panel\
                  (ref_tex_panel, "surface phong", SURF_PHONG, textureinterfglui_cb);
      textureinterfglui_right->add_statictext("");
      GLUI_Spinner *start_spinner = textureinterfglui_right->add_spinner\
                  ("start texture", GLUI_SPINNER_FLOAT, &start_texture);
      start_spinner->set_float_limits(0, 1.0);
      GLUI_Spinner *end_spinner = textureinterfglui_right->add_spinner\
                  ("end texture", GLUI_SPINNER_FLOAT, &end_texture);
      end_spinner->set_float_limits(0, 1.0);
#endif
   }
}

