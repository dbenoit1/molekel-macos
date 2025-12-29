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

#include "version.h"
#include "molekel.h"
#include "constant.h"
#include "general.h"
#include "glutwin.h"
#include "maininterf.h"
#include "browser.h"
#include "action.h"
#include "manip.h"
#include "macu.h"
#include "pick.h"
#include "chooseinterf.h"

GLUI *mainglui_top, *mainglui_bottom, *infoglui;
GLUI_StaticText *infotext[5];
int maininterfwin;
char active_mol_nbr[3];
char active_mol_name[30] = "";
char logstack[50][30];
int stacknr = 0, logscroll = 0, logsize = 8, oy;
int topbar[] = {520, 540}, bottombar[] = {625, 645};

int molprop_live_var = 1;

void mainglui_cb(int key)
/* Callback of maininterface */
{
   char infostring[500];

   switch (key) {
      case MOLMOD:
         switch (molprop_live_var) {
            case 0:
               Set_wire();
            break;
            case 1:
               Set_stick();
            break;
            case 2:
               Set_ball_and_stick();
            break;
            case 3:
               Set_spacefill();
            break;
         }
      break;
      case TOGGLE_MULT:
         if(actualmol) actualmol->multiple_bonds = bit.multiple_bonds;
      break;
      case BOND_ATTR:
          bond_interf();
      break;
      case DEL_MOL:
         if(actualmol) {
            logprint("Deleting molecule");
            free_mol(actualmol);
            update_interface_flags();
         }
         else logprint("No molecule loaded");
         update_logs();
      break;
      case MANIP_ALL:
          if(!bit.manip_all) {
             reset_globals();
          }
      break;
      case SWITCH_MOL:
         switch_mol();
      break;
      case PICK_MOL:
         if(actualmol && firstmol->next) {
            enter_pick(key, "Pick molecule");
         }
         else logprint("No or one molecule loaded");
         update_logs();
      break;
      case CENTER:
         if(actualmol) {
            logprint("Reseting rotation center");
            center();
         }
         else logprint("No molecule loaded");
         update_logs();
      break;
      case RESET:
         reset();
         if(actualmol) {
            logprint("Reseting orientation");
            update_logs();
         }
      break;
      case SPIN:
         Tgl_spinning();
      break;
      case FULLSCREEN:
         Tgl_fullscreen();
      break;
      case OPTION:
         select_option();
      break;
      case HELP:
         help();
      break;
      case INFO:
         sprintf(infostring, "MOLEKEL, Version %s, Date: %s, ",\
            MOLEKEL_VERSION, MOLEKEL_VERSION_DATE);
         strcat(infostring, "by Stefan Portmann, Copyright (C) 2002 CSCS/ETHZ\n");
         strcat(infostring, "(orig. IRIX GL implementation, concept and data structure by Peter F. Fluekiger, CSCS/UNI Geneva)\n");
         strcat(infostring, "MOLEKEL comes with ABSOLUTELY NO WARRANTY; for details see menu \"Warranty\". ");
         strcat(infostring, "The binary code is available\nfree of charge, but is not in the public domain. ");
         strcat(infostring, "See menu \"License\" for details on conditions and restrictions.\n");
         strcat(infostring, "Info: http://www.cscs.ch/molekel/");
/*
         strcat(infostring, "This is free software, and you are welcome to redistribute it ");
         strcat(infostring, "under certain conditions; menu \"License\" for details.\n");
         strcat(infostring, "Info: http://igc.ethz.ch/molekel/");
*/
/*
         sprintf(infostring, "MOLEKEL, Version %s, Date: %s\n",\
            MOLEKEL_VERSION, MOLEKEL_VERSION_DATE);
         strcat(infostring, "\n \nOpenGL/GLUT implementation by Stefan Portmann\n");
         strcat(infostring, "Original IRIX GL implementation by Peter F. Fluekiger\n");
         strcat(infostring, "Info: http://igc.ethz.ch/molekel/");
*/
         showinfobox(infostring);
      break;
   }
}

void update_logs(void)
/* Update logs in maininteface */
{
   int win = glutGetWindow();
   
   glutSetWindow(maininterfwin);
   glutPostRedisplay();
   glutSetWindow(win);
}
   
void logstring(int x, int y, char *string)
/* Print string in a GLUT window at position x, y */
{
   int len, i;
   void *font = GLUT_BITMAP_HELVETICA_12;

   glRasterPos2f(x, y);
   len = (int) strlen(string);
   for (i = 0; i < len; i++) {
      glutBitmapCharacter(font, string[i]);
  }
}

void logprint(char *string)
/* Set stackentry in logstack and step forward position by 1
 * Window update needs to be called seperately (update_logs)
*/
{
   strncpy(logstack[stacknr], string, 29);
   stacknr = (((stacknr + 1) % 50) == 0) ? 0 : stacknr + 1;
   logscroll = 0;
}

void draw_maininterf(void)
/* Draw GLUT maininterface */
{
   int nr = 1;
   int width = glutGet(GLUT_WINDOW_WIDTH);
   int i, tmp, lognr;
   float indent;
   Mol *mol;
   char *cpt;
   char tmpchar[30] = "";

/* Set actual molecule info */
   if(!firstmol) {
      nr = 0;
      strcpy(active_mol_name, "empty");
   }
   else {
      for (mol = firstmol; mol != actualmol; mol = mol->next) {
         nr++;
      }   
      if(actualmol->filename) {
#ifndef WIN32
         if(!(cpt = strrchr(actualmol->filename, '/')))
#else
         if(!(cpt = strrchr(actualmol->filename, '\\')))
#endif
            cpt = (char *)actualmol->filename;
         else cpt++;
         if(strlen(cpt) > 29) {
            strncpy(tmpchar, cpt, 26);
            sprintf(active_mol_name, "%s...", tmpchar);
         }
         else strncpy(active_mol_name, cpt, 29);
      }
      else strcpy(active_mol_name, "file name unknown!");
   }
   sprintf(active_mol_nbr, "%2d", nr);

/* Draw GLUT window */
   glClear(GL_COLOR_BUFFER_BIT);

   glColor4f(0.0, 0.0, 0.6, 1.0);
   logstring(10, 455, "actual molecule"); 
   glColor4f(0.0, 0.0, 0.0, 1.0);
   logstring(15, 475, "No."); 
   logstring(15, 490, active_mol_nbr); 
   logstring(60, 475, "file name");
   logstring(60, 490, active_mol_name);

   if(picking_mode) {
   glColor4f(0.0, 0.3, 0.0, 1.0);
   logstring(150, 455, "picking mode"); 
   glColor4f(0.0, 0.0, 0.0, 1.0);
   }

   if(bit.tgl_idle) {
   glColor4f(0.0, 0.3, 0.0, 1.0);
   logstring(170, 475, "GLUI Idle"); 
   glColor4f(0.0, 0.0, 0.0, 1.0);
   }

/* Draw separator ala GLUI */
   indent = width * .05;
   glLineWidth( 1.0 );
   glBegin( GL_LINES );
   glColor3f( .5, .5, .5 );
   glVertex2f( indent, 504);    
   glVertex2f( width-indent, 504);    

   glColor3f( 1., 1., 1. );
   glVertex2f( indent, 505);    
   glVertex2f( width-indent, 505);
   glEnd();

/* Draw scrollable log info */
   glColor4f(0.6, 0.0, 0.0, 1.0);
   logstring(10, 520, "log info");
   glColor4f(0.0, 0.0, 0.0, 1.0);

   tmp = stacknr - (logsize -1) + logscroll;
   lognr = (tmp <= 0) ? (tmp + 49) : tmp - 1;
   for (i = 0; i < logsize; i++) {
      logstring(15, 535 + 15*i, logstack[lognr]);
      lognr = (((lognr + 1) % 50) == 0) ? 0 : lognr + 1;
   }

   glutSwapBuffers();
}

void reshape_maininterf(int w, int h)
/* glutReshapeFunc for maininterface */
{
  glutReshapeWindow(1024 - 56 - xsize - xorig, ysize);
//  glViewport(0, 0, w, h);
  glViewport(0, 0, 1024 - 56 - xsize - xorig, ysize);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
//  gluOrtho2D(0, w, h, 0);
  gluOrtho2D(0, 1024 - 56 - xsize - xorig, ysize, 0);
  glMatrixMode(GL_MODELVIEW);
}

void mouse_maininterf(int button, int state, int x, int y)
/* glutMouseFunc for maininterface */
{
   if (state == GLUT_DOWN && y > bottombar[0] && y < bottombar[1]) {
      if (logscroll > -(50 - logsize)) logscroll--;
      glutPostRedisplay();
   }
   if (state == GLUT_DOWN && y > topbar[0] && y < topbar[1]) {
      if (logscroll < 0) logscroll++;
      glutPostRedisplay();
   }
   if (state == GLUT_DOWN && y > topbar[1] && y < bottombar[0] ) {
      oy = y;
   }
}

void motion_maininterf(int x, int y)
/* glutMotionFunc for maininterface */
{
   if (y > topbar[1] && y < bottombar[0] && (oy < y)) {
      if (logscroll > -(50 - logsize)) logscroll--;
      glutPostRedisplay();
   }
   if (y > topbar[1] && y < bottombar[0] && (oy > y)) {
      if (logscroll < 0) logscroll++;
      glutPostRedisplay();
   }

   if (y > topbar[0] && y < topbar[1]) glutSetCursor(GLUT_CURSOR_TOP_SIDE);
   else if (y > bottombar[0] && y < bottombar[1]) glutSetCursor(GLUT_CURSOR_BOTTOM_SIDE);
   else if (y > topbar[1] && y < bottombar[0]) glutSetCursor(GLUT_CURSOR_UP_DOWN);
   else glutSetCursor(GLUT_CURSOR_LEFT_ARROW);

   oy = y;
}

void passivemotion_maininterf(int x, int y)
/* glutPassiveMotionFunc for maininterface */
{
   if (y > topbar[0] && y < topbar[1]) glutSetCursor(GLUT_CURSOR_TOP_SIDE);
   else if (y > bottombar[0] && y < bottombar[1]) glutSetCursor(GLUT_CURSOR_BOTTOM_SIDE);
   else if (y > topbar[1] && y < bottombar[0]) glutSetCursor(GLUT_CURSOR_UP_DOWN);
   else glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}



void maininterf(void)
/* maininterface
 * Consists of one GLUT window and a top and bottom GLUI window
*/
{
   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
//   glutInitWindowSize(glutGet(GLUT_SCREEN_WIDTH) - xsize - 4*xorig, ysize);
//   glutInitWindowSize(1024 - xsize - 4*xorig, ysize);
//   glutInitWindowSize(1024 - 24 - xsize - 4*xorig, ysize);
   glutInitWindowSize(1024 - 56 - xsize - xorig, ysize);
//   glutInitWindowPosition(3*xorig + xsize , yorig);
   glutInitWindowPosition(24 + xorig + xsize , yorig);
   maininterfwin = glutCreateWindow("main interface");
   glClearColor((float)200/255, (float)200/255, (float)200/255, 1.0);
   glutDisplayFunc(draw_maininterf);
   GLUI_Master.set_glutReshapeFunc(reshape_maininterf);
   GLUI_Master.set_glutMouseFunc(mouse_maininterf);
   glutMotionFunc(motion_maininterf);
   glutPassiveMotionFunc(passivemotion_maininterf);
   GLUI_Master.set_glutKeyboardFunc(dummy_gluiKeyboardFunc);
   GLUI_Master.set_glutSpecialFunc(dummy_gluiSpecialFunc);
   glutSetCursor(GLUT_CURSOR_LEFT_ARROW);

/* top GLUI */
   mainglui_top = GLUI_Master.create_glui_subwindow(maininterfwin, GLUI_SUBWINDOW_TOP);
   mainglui_top->set_main_gfx_window(mainwin);

   GLUI_Panel *molprop_panel = mainglui_top->add_panel("molecular models");
   GLUI_RadioGroup *molprop_rbgrp = mainglui_top->add_radiogroup_to_panel\
             (molprop_panel, &molprop_live_var, MOLMOD, mainglui_cb);
   GLUI_RadioButton *wireframe = mainglui_top->add_radiobutton_to_group(molprop_rbgrp, "wire model");
   GLUI_RadioButton *sticks = mainglui_top->add_radiobutton_to_group(molprop_rbgrp, "sticks");
   GLUI_RadioButton *b_and_s = mainglui_top->add_radiobutton_to_group(molprop_rbgrp, "ball&sticks");
   GLUI_RadioButton *spacefill = mainglui_top->add_radiobutton_to_group(molprop_rbgrp, "spacefill");
   mainglui_top->add_column_to_panel(molprop_panel, true);
   GLUI_Checkbox *multiplebonds = mainglui_top->add_checkbox_to_panel\
             (molprop_panel,"multiplebonds", &(bit.multiple_bonds), TOGGLE_MULT, mainglui_cb);
   mainglui_top->add_statictext_to_panel(molprop_panel, "" );
   GLUI_Button *bondattrib = mainglui_top->add_button_to_panel\
             (molprop_panel, "bond attributes", BOND_ATTR, mainglui_cb);
   mainglui_top->add_statictext("");

   GLUI_Checkbox *manipall = mainglui_top->add_checkbox("manipulate all", &(bit.manip_all), MANIP_ALL, mainglui_cb);
   manipall->set_alignment(GLUI_ALIGN_CENTER);
   mainglui_top->add_statictext("");

   GLUI_Panel *switch_del_panel = mainglui_top->add_panel("", true);
   GLUI_Button *switchmol = mainglui_top->add_button_to_panel\
             (switch_del_panel, "switch mol", SWITCH_MOL, mainglui_cb);
   switchmol->set_w(92);
   mainglui_top->add_column_to_panel(switch_del_panel, false);
   GLUI_Button *delmol = mainglui_top->add_button_to_panel
             (switch_del_panel, "delete mol", DEL_MOL, mainglui_cb);
   delmol->set_w(92);
   mainglui_top->add_statictext("");

   GLUI_Button *delsurf = mainglui_top->add_button( "pick mol", PICK_MOL, mainglui_cb);
   mainglui_top->add_statictext("");

   GLUI_Panel *quality_panel = mainglui_top->add_panel("quality");
   GLUI_RadioGroup *quality_rbgrp = mainglui_top->add_radiogroup_to_panel\
             (quality_panel, &quality);
   GLUI_RadioButton *high = mainglui_top->add_radiobutton_to_group(quality_rbgrp, "high");
   GLUI_RadioButton *medium = mainglui_top->add_radiobutton_to_group(quality_rbgrp, "medium");
   GLUI_RadioButton *low = mainglui_top->add_radiobutton_to_group(quality_rbgrp, "low");
   mainglui_top->add_column_to_panel(quality_panel, true);
   GLUI_Checkbox *depthcue = mainglui_top->add_checkbox_to_panel\
             (quality_panel,"depthcue", &(bit.depthcue));
   GLUI_Checkbox *texture = mainglui_top->add_checkbox_to_panel\
             (quality_panel,"texture", &(bit.texture));
   GLUI_Checkbox *smooth = mainglui_top->add_checkbox_to_panel\
             (quality_panel,"smooth line", &(bit.smoothline));
   mainglui_top->add_statictext("");


   GLUI_Panel *view_panel = mainglui_top->add_panel("", true);
   GLUI_Button *center = mainglui_top->add_button_to_panel\
             (view_panel, "center", CENTER, mainglui_cb);
   center->set_w(92);
   GLUI_Button *reset = mainglui_top->add_button_to_panel\
             (view_panel, "reset", RESET, mainglui_cb);
   reset->set_w(92);
   mainglui_top->add_column_to_panel(view_panel, false);
   GLUI_Button *spin = mainglui_top->add_button_to_panel\
             (view_panel, "spin", SPIN, mainglui_cb);
   spin->set_w(92);
   GLUI_Button *fullscreen = mainglui_top->add_button_to_panel\
             (view_panel, "fullscreen", FULLSCREEN, mainglui_cb);
   fullscreen->set_w(92);
   mainglui_top->add_statictext("");

   GLUI_Button *option = mainglui_top->add_button\
             ("option", OPTION, mainglui_cb);

/* bottom GLUI */
   mainglui_bottom = GLUI_Master.create_glui_subwindow(maininterfwin, GLUI_SUBWINDOW_BOTTOM);
   mainglui_bottom->set_main_gfx_window(mainwin);
   GLUI_Panel *help_panel = mainglui_bottom->add_panel("", true);
   GLUI_Button *help = mainglui_bottom->add_button_to_panel\
             (help_panel, "help", HELP, mainglui_cb);
   help->set_w(92);
   mainglui_bottom->add_column_to_panel(help_panel, false);
   GLUI_Button *info = mainglui_bottom->add_button_to_panel\
             (help_panel, "info", INFO, mainglui_cb);
   info->set_w(92);
   mainglui_bottom->add_statictext("");

   GLUI_Button *quit = mainglui_bottom->add_button( "quit", 0,(GLUI_Update_CB)Quit );
}


void hideinfobox(int key)
{
   infoglui->close();
   infoglui = NULL;
   glutSetWindow(mainwin);
}

void showinfobox(char *s)
/* Show info box
 * string can contain max. 5 lines separated by \n 
 * \n\n is only one \n
 * \n \n gives a blank line
*/
{
   char *str;
   register int i;
   char *string;

   if(infoglui) {
      infoglui->close();
      infoglui = NULL;
   }
   makeinfowin();

   string = strdup(s);
   str = strtok(string, "\n");   
   infotext[0]->set_name(str); 
   for(i = 1; i < 5; i++) {
      if((str = strtok(NULL, "\n")) == NULL) infotext[i]->set_name("");
      else infotext[i]->set_name(str);
   }
}

void makeinfowin(void)
/* Create infobox
 * Contains 5 lines to write on
*/
{
   register int i;

//   infoglui = GLUI_Master.create_glui("info", 0, xsize/2, ysize/2);
   glutSetWindow(mainwin);
   infoglui = GLUI_Master.create_glui("info", 0, \
        winposition[0], winposition[1] + 300);
//   infoglui->set_main_gfx_window(mainwin); dont need redraw, no event to mainwin!
   GLUI_Panel *infotitle_panel = infoglui->add_panel("", true);
   infotitle_panel->set_alignment(GLUI_ALIGN_LEFT);
   infoglui->add_statictext_to_panel(infotitle_panel,"       I N F O !");
   GLUI_Panel *info_panel = infoglui->add_panel("", true);
   info_panel->set_alignment(GLUI_ALIGN_LEFT);
   infoglui->add_statictext_to_panel(info_panel, "");
   for (i = 0; i < 5; i++) {
      infotext[i] = infoglui->add_statictext_to_panel(info_panel, "information");
   }
   infoglui->add_statictext_to_panel(info_panel, "");
   GLUI_Button *ok = infoglui->add_button("ok", 0, hideinfobox);
}


void green_note(char *str)
{
   if(maininterfwin) glutSetWindow(maininterfwin);
   glDrawBuffer(GL_FRONT);
   glColor4f(0.0, 0.3, 0.0, 1.0);
   logstring(140, 455, str); 
   glColor4f(0.0, 0.0, 0.0, 1.0);
   glDrawBuffer(GL_BACK);
   if(mainwin) glutSetWindow(mainwin);
}

