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
#include <math.h>

#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "general.h"
#include "glutwin.h"
#include "maininterf.h"
#include "chooseinterf.h"
#include "menu.h"
#include "browser.h"
#include "calcdens.h"
#include "cubeInterpol.h"
#include "pick.h"
#include "objects.h"
#include "manip.h"
#include "drawing.h"
#include "box.h"
#include "connolly.h"
#include "connect.h"
#include "snap.h"
#include "colorinterf.h"
#include "textureinterf.h"
#include "action.h"
#include "render.h"
#include "material.h"
#include "utils.h"

#ifdef LINUX
#define NOSOCK
#endif

#ifdef WIN32
#define NOSOCK
#endif

GLUI *coef_or_matglui = NULL, *overwrite_filegui = NULL, *job_prioritygui = NULL;
GLUI *quitgui = NULL, *grid_or_dotglui = NULL, *select_surfgui = NULL;
GLUI *vmin_vmaxglui = NULL, *optionglui = NULL, *bondglui = NULL, *freqglui = NULL;
GLUI *playglui = NULL, *boxglui = NULL, *connollyglui = NULL, *triangglui = NULL, *labelglui = NULL;
GLUI *xyzformatglui = NULL, *t41contentglui = NULL, *renderglui = NULL, *dipoleglui = NULL;
GLUI_Spinner *xmin_spinner, *ymin_spinner, *zmin_spinner, *xmax_spinner, *ymax_spinner, *zmax_spinner;
GLUI_Spinner *zoom_spinner, *framenbr_spinner = NULL, *jpegqual_spinner;
GLUI_StaticText *freqtext;
GLUI_StaticText *dipoletext;
GLUI_StaticText *play_no;
GLUI_Button *freq_choose;
GLUI_Button *coord_cancel;
int datasource = USE_COEFFS;
int coef_or_mat_live_var = 0, grid_or_dot_live_var = 0;
int job_priority_live_var = 19, surface_live_var = 1, speed_live_var = 1;
int vibrating = 0, playing = 0;
int img_ff_live_var = 0, render_action_live_var = 0 , keep_ratio = 1;
int w_width = 736, w_height = 736;
float zoombox = 0;

void ok_cancel(GLUI *glui, void (*func) (int key))
{
   GLUI_Panel *ok_cancel_panel = glui->add_panel("");
   GLUI_Button *ok = glui->add_button_to_panel(ok_cancel_panel, "ok", OK, func);
   glui->add_column_to_panel(ok_cancel_panel, false);
   GLUI_Button *cancel =  glui->add_button_to_panel(ok_cancel_panel, "cancel", CANCEL, func);
}


void yes_no(GLUI *glui, void (*func) (int key))
{
   GLUI_Panel *yes_no_panel = glui->add_panel("");
   GLUI_Button *yes = glui->add_button_to_panel(yes_no_panel, "yes", YES, func);
   glui->add_column_to_panel(yes_no_panel, false);
   GLUI_Button *no =  glui->add_button_to_panel(yes_no_panel, "no", NO, func);
}

/* select_coeff_or_matrix interface */

void coef_or_mat_cb(int key)
{
    switch (key) {
      case OK:
         if (action_key == EL_DENS) {
            file_select(EL_DENS, macu_ext);
         }
         else if (action_key == SPIN_DENS) {
            file_select(SPIN_DENS, macu_ext);
         }
         else if (action_key == ATOM_SPIN) {
            spin_on_atoms();
            action_key = 0;
            mainglui_top->enable();
         }
         if(coef_or_matglui) {
            coef_or_matglui->close();
            coef_or_matglui = NULL;
         }
      break;
      case CANCEL:
         action_key = 0;
         coef_or_matglui->close();
         coef_or_matglui = NULL;
         mainglui_top->enable();
      break;
      case 0:
         switch (coef_or_mat_live_var) {
            case 0:
               datasource = USE_COEFFS;
            break;
            case 1:
               datasource = USE_MATRICES;
            break;
            default:
               datasource = USE_COEFFS;
         }
      break;
   }
}


void select_coeff_or_matrix(void)
{
   if(!actualmol->alphaDensity) {
      datasource = USE_COEFFS;
      coef_or_mat_cb(OK);
      return;
   }

   if(!coef_or_matglui) {
      coef_or_matglui = GLUI_Master.create_glui\
           ("coeff_or_matrix", 0, winposition[0] + 120, winposition[1] + 370);
//   coef_or_matglui ->set_main_gfx_window(mainwin);
      GLUI_Panel *title_panel = coef_or_matglui->add_panel("", true);
      coef_or_matglui->add_statictext_to_panel(title_panel,"Choose the data for the density");
      GLUI_Panel *choose_panel = coef_or_matglui->add_panel("", true);
      GLUI_RadioGroup *choose_rbgrp = coef_or_matglui->add_radiogroup_to_panel\
              (choose_panel, &coef_or_mat_live_var, 0, coef_or_mat_cb);
      GLUI_RadioButton *mo_coef = coef_or_matglui->add_radiobutton_to_group\
              (choose_rbgrp, "Use the MO-coefficients");
      GLUI_RadioButton *dens= coef_or_matglui->add_radiobutton_to_group\
              (choose_rbgrp, "Use the density-matrices");
      ok_cancel(coef_or_matglui, coef_or_mat_cb);
      mainglui_top->disable();
   }
   else {
      coef_or_matglui->show();
      glutSetWindow(coef_or_matglui->get_glut_window_id());
      glutPopWindow();
   }
}

/* select overwrite_file interface */

void overwrite_file_cb(int key)
{
    switch (key) {
      case YES:
         overwrite_filegui->hide();
         overwrite_filegui->close();
         overwrite_filegui = NULL;
         switch (action_key) {
            case EL_DENS:
            case CALC_ORB:
            case SPIN_DENS:
            case MEP:
               select_job_priority();
            break; 
            case  WRITE_XYZ_ORIG:
               logprint("");
               logprint("Saving XYZ (orig. orient.)");
               update_logs();
               glutTimerFunc(1, timer_cb, action_key);
               mainglui_top->enable();
               action_key = 0;
            break;
            case  WRITE_XYZ_CURRENT:
               logprint("");
               logprint("Saving XYZ (current orient.)");
               update_logs();
               glutTimerFunc(1, timer_cb, action_key);
               mainglui_top->enable();
               action_key = 0;
            break;
            case  WRITE_PDB_ORIG:
               logprint("");
               logprint("Saving PDB (orig. orient.)");
               update_logs();
               glutTimerFunc(1, timer_cb, action_key);
               mainglui_top->enable();
               action_key = 0;
            break;
            case  WRITE_PDB_CURRENT:
               logprint("");
               logprint("Saving PDB (current orient.)");
               update_logs();
               glutTimerFunc(1, timer_cb, action_key);
               mainglui_top->enable();
               action_key = 0;
            break;
            case  WRITE_MS:
               logprint("");
               logprint("Saving MS");
               update_logs();
               glutTimerFunc(1, timer_cb, action_key);
               mainglui_top->enable();
               action_key = 0;
            break;
            case WRITE_SLD:
               logprint("");
               logprint("Saving SLD");
               update_logs();
               glutTimerFunc(1, timer_cb, action_key);
               mainglui_top->enable();
               action_key = 0;
            break;
            case  WRITE_DOTVAL:
               logprint("");
               logprint("Saving dots with value");
               update_logs();
               glutTimerFunc(1, timer_cb, action_key);
               mainglui_top->enable();
               action_key = 0;
            break;
            case WRITE_RGB:
               logprint("");
               logprint("Saving RGB");
               update_logs();
               glutTimerFunc(1, timer_cb, action_key);
               mainglui_top->enable();
               action_key = 0;
            break;
            case WRITE_TIFF:
               logprint("");
               logprint("Saving TIFF");
               update_logs();
               glutTimerFunc(1, timer_cb, action_key);
               mainglui_top->enable();
               action_key = 0;
            break;
            case WRITE_TIFF_LZW:
               logprint("");
               logprint("Saving TIFF LZW");
               update_logs();
               glutTimerFunc(1, timer_cb, action_key);
               mainglui_top->enable();
               action_key = 0;
            break;
            case WRITE_TIFF_PACK:
               logprint("");
               logprint("Saving TIFF packbits");
               update_logs();
               glutTimerFunc(1, timer_cb, action_key);
               mainglui_top->enable();
               action_key = 0;
            break;
            case WRITE_JPEG:
               logprint("");
               logprint("Saving JPEG");
               update_logs();
               glutTimerFunc(1, timer_cb, action_key);
               mainglui_top->enable();
               action_key = 0;
            break;
            case TIF_CURR_IMG:
            case TIFLZW_CURR_IMG:
            case TIFPACK_CURR_IMG:
            case PPM_CURR_IMG:
            case RGB_CURR_IMG:
            case JPEG_CURR_IMG:
               logprint("");
               logprint("Rendering image");
               update_logs();
               glutTimerFunc(500, timer_cb, action_key); // 500ms required!
               mainglui_top->enable();
               action_key = 0;
            break;
         }
      break;
      case NO:
         overwrite_filegui->close();
         overwrite_filegui = NULL;
         mainglui_top->enable();
         if(pixels) {
            free(pixels);
            pixels = NULL;
         }
         action_key = 0;
      break;
    }
}

void quit_cb(int key)
{
    switch (key) {
      case OK:
/* trick to close windows, when a socket process is still active
 * a simple exit(0) leaves the windows on the screen, even tough they
 * are not active anymore
*/
#ifdef WIN32
         unlink("C:\\TEMP\\connolly.ms");
         unlink("C:\\TEMP\\input");
#endif         
         glutDestroyWindow(mainwin);
         glutDestroyWindow(maininterfwin);
         if(filewin) glutDestroyWindow(filewin);
         if(freqwin) glutDestroyWindow(freqwin);
         if(orbwin) glutDestroyWindow(orbwin);
         if(colorinterfwin) glutDestroyWindow(colorinterfwin);
         if(textureinterfwin) glutDestroyWindow(textureinterfwin);
         quitgui->close();
         quitgui = NULL;
         glutTimerFunc(1, exit, 0);
      break;
      case CANCEL:
         quitgui->close();
         quitgui = NULL;
         return;
    }
}


void select_overwrite_file(char *fname)
{
   if(!overwrite_filegui) {
      overwrite_filegui = GLUI_Master.create_glui\
           ("overwrite file", 0, winposition[0] + 15, winposition[1] + 385);
      GLUI_Panel *title_panel = overwrite_filegui->add_panel("", true);
      overwrite_filegui->add_statictext_to_panel(title_panel, fname);
      overwrite_filegui->add_statictext_to_panel(title_panel, "");
      overwrite_filegui->add_statictext_to_panel(title_panel, "File exists! Overwrite it?");
      yes_no(overwrite_filegui, overwrite_file_cb);
      mainglui_top->disable();
   }
   else {
      overwrite_filegui->show();
      glutSetWindow(overwrite_filegui->get_glut_window_id());
      glutPopWindow();
   }
}


void setup_job_cb(int key)
{
   job_prioritygui->close();
   job_prioritygui = NULL;
   mainglui_top->enable();
   calcdens(action_key);
   action_key = 0;
}

void select_job_priority(void)
{
#ifndef NOSOCK
   if(!job_prioritygui) {
      job_prioritygui = GLUI_Master.create_glui\
        ("job priority", 0, winposition[0] + 100, winposition[1] + 370);
      GLUI_Panel *title_panel = job_prioritygui->add_panel("", true);
      job_prioritygui->add_statictext_to_panel(title_panel, "Set job priority");
      job_prioritygui->add_statictext("");
      GLUI_Spinner *priority_spinner = job_prioritygui->add_spinner\
          ("job priority", GLUI_SPINNER_INT, &job_priority_live_var);
      priority_spinner->set_int_limits(0, 19);
      job_prioritygui->add_statictext("");
      GLUI_Button *ok = job_prioritygui->add_button("ok", OK, setup_job_cb);
   }
   else {
      job_prioritygui->show();
      glutSetWindow(job_prioritygui->get_glut_window_id());
      glutPopWindow();
   }
#else
   job_priority_live_var = 19;
   mainglui_top->enable();
   calcdens(action_key);
   action_key = 0;
#endif
}


void select_quit(void)
{
   if(!quitgui) {
      quitgui =  GLUI_Master.create_glui\
          ("quit", 0, winposition[0] + 240, winposition[1] + 300);
      GLUI_Panel *title_panel = quitgui->add_panel("", true);
      quitgui->add_statictext_to_panel(title_panel, "Are you sure you want to quit?");
      ok_cancel(quitgui, quit_cb);
   }
   else {
      quitgui->show();
      glutSetWindow(quitgui->get_glut_window_id());
      glutPopWindow();
   }
}


void select_surface_cb(int key)
{
    switch (key) {
       case MEP:
        if(select_surfgui) select_surfgui->close();
        select_surfgui = NULL;
        action_key = 0;
        mainglui_top->enable();
        glutSetWindow(mainwin);
        mep_dot_surface_no_pick(); 
       break;
       case GRID_VALUES:
        if(select_surfgui) select_surfgui->close();
        select_surfgui = NULL;
        action_key = 0;
        mainglui_top->enable();
        glutSetWindow(mainwin);
        gridValues_no_pick();
       break;
    }
}


void select_surface(int key)
{
/* not used anymore */
   Surface *surf;
   int limit;

   if(!(surf = actualmol->firstsurf)){
      action_key = 0;
      mainglui_top->enable();
      showinfobox("Load or generate a molecular surface first!");
      return;
   }
   for(surface_live_var = 0; surf; surf=surf->next) surface_live_var++;
   if(surface_live_var == 1) {
      select_surface_cb(key);
   }
   else {
      limit = surface_live_var;
      if (!select_surfgui) {
         select_surfgui = GLUI_Master.create_glui\
          ("select surface", 0, winposition[0], winposition[1]);
         GLUI_Panel *title_panel = select_surfgui->add_panel("", true);
         select_surfgui->add_statictext_to_panel(title_panel, "Select surface");
         select_surfgui->add_statictext("");
         GLUI_Spinner *priority_spinner = select_surfgui->add_spinner\
             ("surface number", GLUI_SPINNER_INT, &surface_live_var);
         priority_spinner->set_int_limits(1, limit);
         select_surfgui->add_statictext("");
         GLUI_Button *ok = select_surfgui->add_button("ok", key, select_surface_cb);
      }
      else {
         select_surfgui->show();
         glutSetWindow(select_surfgui->get_glut_window_id());
         glutPopWindow();
      }
   }
}

void grid_or_dot_cb(int key)
{
    switch (key) {
      case OK:
         if (action_key == MEP) {
            switch (grid_or_dot_live_var) {
               case 0:
                  pick_surface(MEP);
               break;
               case 1:
                  file_select(MEP, macu_ext);
               break;
            }
         }
         grid_or_dotglui->close();
         grid_or_dotglui = NULL;
      break;
      case CANCEL:
         action_key = 0;
         grid_or_dotglui->close();
         grid_or_dotglui = NULL;
         mainglui_top->enable();
      break;
   }
}

void select_grid_or_dot(void)
{
   if(!grid_or_dotglui) {
      grid_or_dotglui = GLUI_Master.create_glui\
          ("3D-grid or map to surface", 0, winposition[0], winposition[1]);
      GLUI_Panel *title_panel = grid_or_dotglui->add_panel("", true);
      grid_or_dotglui->add_statictext_to_panel(title_panel,"Calculate the MEP");
      GLUI_Panel *choose_panel = grid_or_dotglui->add_panel("", true);
      GLUI_RadioGroup *choose_rbgrp = grid_or_dotglui->add_radiogroup_to_panel\
              (choose_panel, &grid_or_dot_live_var, 0, grid_or_dot_cb);
      GLUI_RadioButton *mo_coef = grid_or_dotglui->add_radiobutton_to_group\
              (choose_rbgrp, "map to surface");
      GLUI_RadioButton *dens= grid_or_dotglui->add_radiobutton_to_group\
              (choose_rbgrp, "3D-grid");
      ok_cancel(grid_or_dotglui, grid_or_dot_cb);
      mainglui_top->disable();
   }
   else {
      grid_or_dotglui->show();
      glutSetWindow(grid_or_dotglui->get_glut_window_id());
      glutPopWindow();
   }
}

void vmin_vmax_cb(int key)
{
    char str[100];

    vmin_vmaxglui->close();
    vmin_vmaxglui = NULL;

    glutSetWindow(mainwin);
    glutPostRedisplay();

    sprintf(str, "vmin = %f", actualsurf->vmin);
    logprint("");
    logprint("Mapping grid values:");
    logprint(str);
    sprintf(str, "vmax = %f", actualsurf->vmax);
    logprint(str);
    update_logs();

}

 
void set_vmin_vmax(Surface *surf)
{
   if(!vmin_vmaxglui) {
      vmin_vmaxglui = GLUI_Master.create_glui\
         ("set vmin vmax", 0, winposition[0], winposition[1]);
         GLUI_Panel *title_panel = vmin_vmaxglui->add_panel("", true);
         vmin_vmaxglui->add_statictext_to_panel(title_panel, "set vmin vmax");
         vmin_vmaxglui->add_statictext("");
         GLUI_Spinner *vmax_spinner = vmin_vmaxglui->add_spinner\
             ("vmax", GLUI_SPINNER_FLOAT, &(surf->vmax));
// 28.10.02: no limits: make in some cases no sense, let the user fully decide
// and think
//         vmax_spinner->set_float_limits(surf->vmin, fabsf(10 * (surf->vmax)));
         GLUI_Spinner *vmin_spinner = vmin_vmaxglui->add_spinner\
             ("vmin", GLUI_SPINNER_FLOAT, &(surf->vmin));
//         vmin_spinner->set_float_limits(10 * surf->vmin, surf->vmax);
         vmin_vmaxglui->add_statictext("");
         GLUI_Button *ok = vmin_vmaxglui->add_button("ok", 0, vmin_vmax_cb);
 
   }
   else {
      vmin_vmaxglui->show();
      glutSetWindow(vmin_vmaxglui->get_glut_window_id());
      glutPopWindow();
   }
}


void option_cb(int key)
{
    switch (key) {
      case PERSP:
         init_persp();
      break;
      case DEPTHCUE:
         glFogf(GL_FOG_START, fogstart);
         glFogf(GL_FOG_END, fogend);
      break;
      case SET_WSIZE:
         glutReshapeWindow(w_width, w_height);
      break;
      case TRANSPARENT:
         if(actualmol) actualmol->transparent = bit.transparent;
      break;
      case NO_H:
         if(actualmol) actualmol->suppress_H = bit.suppress_H;
      break;
      case RESIDUES:
         if(actualmol) actualmol->residues = bit.residues;
      break;
      case TGL_IDLE:
         if(!bit.tgl_idle) {
            Set_glut_idle();
         }
         else {
            Set_glui_idle();
         }
         update_logs();
      break;
      case CANCEL_OPT:
         optionglui->close();
         optionglui = NULL;
      break;
    }
}


void select_option(void)
{
   if(!optionglui) {
      optionglui = GLUI_Master.create_glui\
         ("option", 0, winposition[0], winposition[1]);
      optionglui->set_main_gfx_window(mainwin);
      GLUI_Panel *proj_panel = optionglui->add_panel("projection", true);
      GLUI_RadioGroup *option_rbgrp = optionglui->add_radiogroup_to_panel\
              (proj_panel, &bit.persp, PERSP, option_cb);
      GLUI_RadioButton *ortho = optionglui->add_radiobutton_to_group\
              (option_rbgrp, "orthogonal projection");
      GLUI_RadioButton *persp = optionglui->add_radiobutton_to_group\
              (option_rbgrp, "perspective projection");
      optionglui->add_statictext("");
      GLUI_Checkbox *coord_sys = optionglui->add_checkbox("show x-, y-, z-axis", &(bit.coord_sys));
      coord_sys->set_alignment(GLUI_ALIGN_CENTER);
      optionglui->add_statictext("");
      GLUI_Panel *wsize_panel = optionglui->add_panel("window size", true);
      GLUI_Spinner *size_x_spinner = optionglui->add_spinner_to_panel\
             (wsize_panel, "width ", GLUI_SPINNER_INT, &w_width);
      size_x_spinner->set_int_limits(16, glutGet(GLUT_SCREEN_WIDTH) ? glutGet(GLUT_SCREEN_WIDTH) : 736);
      GLUI_Spinner *size_y_spinner = optionglui->add_spinner_to_panel\
             (wsize_panel, "height", GLUI_SPINNER_INT, &w_height);
      size_y_spinner->set_int_limits(16, glutGet(GLUT_SCREEN_HEIGHT) ? glutGet(GLUT_SCREEN_HEIGHT) : 736);
      optionglui->add_statictext_to_panel(wsize_panel, "");
      GLUI_Button *set= optionglui->add_button_to_panel(wsize_panel, "set", SET_WSIZE, option_cb);
      optionglui->add_statictext("");
      GLUI_Panel *depth_panel = optionglui->add_panel("depthcue range", true);
      GLUI_Spinner *start_spinner = optionglui->add_spinner_to_panel\
             (depth_panel, "start", GLUI_SPINNER_FLOAT, &(fogstart), DEPTHCUE, option_cb);
      GLUI_Spinner *end_spinner = optionglui->add_spinner_to_panel\
             (depth_panel, "end", GLUI_SPINNER_FLOAT, &(fogend), DEPTHCUE, option_cb);
      optionglui->add_statictext("");
      GLUI_Panel *molprop_panel = optionglui->add_panel("molecular properties", true);
      GLUI_Checkbox *transp_cpk = optionglui->add_checkbox_to_panel\
             (molprop_panel,"transparent spacefill", &(bit.transparent), TRANSPARENT, option_cb);
      GLUI_Checkbox *hydrogen = optionglui->add_checkbox_to_panel\
             (molprop_panel,"hide hydrogens", &(bit.suppress_H), NO_H, option_cb);
      GLUI_Checkbox *residue = optionglui->add_checkbox_to_panel\
             (molprop_panel,"backbone only", &(bit.residues), RESIDUES, option_cb);
      GLUI_StaticText *text = optionglui->add_statictext_to_panel(molprop_panel,"(pdb files only)");
      text->set_alignment(GLUI_ALIGN_RIGHT);
      optionglui->add_statictext("");
      GLUI_Panel *idle_panel = optionglui->add_panel("idle", true);
      GLUI_Checkbox *glui_idle = optionglui->add_checkbox_to_panel\
              (idle_panel, "enable GLUI idle function", &(bit.tgl_idle), TGL_IDLE, option_cb);
      GLUI_StaticText *idle_text1 = optionglui->add_statictext_to_panel\
              (idle_panel, "makes spinners really spinning");
      idle_text1->set_alignment(GLUI_ALIGN_RIGHT);
      GLUI_StaticText *idle_text2 = optionglui->add_statictext_to_panel\
              (idle_panel, "(uses resources even MOLEKEL is idle)");
      idle_text2->set_alignment(GLUI_ALIGN_RIGHT);
      
      optionglui->add_statictext("");
      GLUI_Button *cancel = optionglui->add_button("cancel", CANCEL_OPT, option_cb);
   }
   else {
      optionglui->show();
      glutSetWindow(optionglui->get_glut_window_id());
      glutPopWindow();
   }
} 

void bond_interf_cb(int key)
{
    switch (key) {
      case ADD_BOND:
         glutSetWindow(mainwin);
         enter_pick(key, "Add bond: pick two atoms");
      break;
      case DEL_BOND:
         glutSetWindow(mainwin);
         enter_pick(key, "Del bond: pick a bond");
      break;
      case SINGLEBOND:
         glutSetWindow(mainwin);
         enter_pick(key, "change to single: pick a bond");
      break;
      case DOUBLEBOND:
         glutSetWindow(mainwin);
         enter_pick(key, "change to double: pick a bond");
      break;
      case TRIPLEBOND:
         glutSetWindow(mainwin);
         enter_pick(key, "change to triple: pick a bond");
      break;
      case H_BOND:
         glutSetWindow(mainwin);
         enter_pick(key, "change to triple: pick a bond");
      break;
      case BOND_COL:
         if(actualmol) actualmol->bond_col = bit.bond_col;
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case SHOW_H_BOND:
         if(actualmol) actualmol->h_bond = bit.h_bond;
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case UPDATE_H_BOND:
         if(actualmol) create_Hbonds();
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case VDW_SCALE:
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case BOND_THICK:
         glutSetWindow(mainwin);
         generate_stick_object();
         scale_line(xsize, ysize);
         glutPostRedisplay();
      break;
      case CANCEL:
         bondglui->close();
         bondglui = NULL;
      break;
    }
}


void bond_interf(void)
{
   if(!bondglui) {
      bondglui = GLUI_Master.create_glui\
         ("bond attributes", 0, winposition[0], winposition[1]);
      GLUI_Panel *modbond_panel = bondglui->add_panel("", true);
      GLUI_Button *del_bond = bondglui->add_button_to_panel(modbond_panel, "del bond", DEL_BOND, bond_interf_cb);
      bondglui->add_column_to_panel(modbond_panel, false);
      GLUI_Button *add_bond = bondglui->add_button_to_panel(modbond_panel, "add bond", ADD_BOND, bond_interf_cb);

      GLUI_Panel *h_bond_panel = bondglui->add_panel("H-bonds", true);
      bondglui->add_statictext_to_panel(h_bond_panel, "");
      GLUI_Checkbox *show_h_bond= bondglui->add_checkbox_to_panel\
             (h_bond_panel, "show H-bonds", &(bit.h_bond), SHOW_H_BOND, bond_interf_cb);
      bondglui->add_statictext_to_panel(h_bond_panel, "");
      GLUI_Button *h_bond = bondglui->add_button_to_panel(h_bond_panel, "mod to H-bond", H_BOND, bond_interf_cb);
      bondglui->add_column_to_panel(h_bond_panel, true);
      GLUI_Spinner *h_b_dist_spinner = bondglui->add_spinner_to_panel\
             (h_bond_panel, "max. distance", GLUI_SPINNER_FLOAT, &max_h_bond);
      h_b_dist_spinner->set_float_limits(0.1, 10.0);
      GLUI_Spinner *h_b_angle_spinner = bondglui->add_spinner_to_panel\
             (h_bond_panel, "min. angle", GLUI_SPINNER_FLOAT, &min_h_angle);
      h_b_angle_spinner->set_float_limits(0.1, 180.0);
      bondglui->add_statictext_to_panel(h_bond_panel, "");
      GLUI_Button *h_b_recalc = bondglui->add_button_to_panel(h_bond_panel, "update H-bonds", UPDATE_H_BOND, bond_interf_cb);
      
      bondglui->add_statictext("");
      GLUI_Panel *ballstick_panel = bondglui->add_panel("ball and stick attributes", true);
      GLUI_Panel *ballstick_panel1 = bondglui->add_panel_to_panel(ballstick_panel, "", false);
      GLUI_Button *single_bond = bondglui->add_button_to_panel(ballstick_panel1, "mod to single", SINGLEBOND, bond_interf_cb);
      bondglui->add_column_to_panel(ballstick_panel1, false);
      GLUI_Button *double_bond = bondglui->add_button_to_panel(ballstick_panel1, "mod to double", DOUBLEBOND, bond_interf_cb);
      GLUI_Button *triple_bond = bondglui->add_button_to_panel(ballstick_panel, "mod to tripple", TRIPLEBOND, bond_interf_cb);
      bondglui->add_statictext_to_panel(ballstick_panel, "");
      GLUI_Checkbox *bondcolor= bondglui->add_checkbox_to_panel\
             (ballstick_panel, "own bond color", &(bit.bond_col), BOND_COL, bond_interf_cb);
      bondglui->add_statictext_to_panel(ballstick_panel, "");
      GLUI_Spinner *vdw_spinner = bondglui->add_spinner_to_panel\
             (ballstick_panel, "vdW scale-factor", GLUI_SPINNER_FLOAT, &(vdW_factor), VDW_SCALE, bond_interf_cb);
      vdw_spinner->set_float_limits(0, 1.0);
      GLUI_Spinner *bondthick_spinner = bondglui->add_spinner_to_panel\
             (ballstick_panel, "bond-thickness (A)", GLUI_SPINNER_FLOAT, &(bond_factor), BOND_THICK, bond_interf_cb);
      bondthick_spinner->set_float_limits(0, 1.0);
      bondglui->add_statictext("");
      GLUI_Button *cancel = bondglui->add_button("cancel", CANCEL, bond_interf_cb);
   }
   else {
      bondglui->show();
      glutSetWindow(bondglui->get_glut_window_id());
      glutPopWindow();
   }
}

void freq_interf_cb(int key)
{

    switch (key) {
      case VIBRATION:
         if(actualvib) {
            if(!save_coords(actualmol)) return;
            vibrating = 1;
            mainglui_top->disable();
            freq_choose->disable();
            glutSetWindow(mainwin);
            glutSetMenu(truncF_menu_id);
            glutAttachMenu(GLUT_RIGHT_BUTTON);
            if(bit.tgl_idle) GLUI_Master.set_glutIdleFunc(idle_freq);
            else glutIdleFunc(idle_freq);
         }
         else {
            logprint("No vibration selected");
            update_logs();
            return;
         }
      break;
      case STOP:
         if(!vibrating) return;
         vibrating = 0;
         freq_pos = 0;
         restore_coords(actualmol);
         mainglui_top->enable();
         freq_choose->enable();
         glutSetWindow(mainwin);
         glutSetMenu(main_menu_id);
         glutAttachMenu(GLUT_RIGHT_BUTTON);
         if(bit.tgl_idle) GLUI_Master.set_glutIdleFunc(NULL);
         else glutIdleFunc(NULL);
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case FREQ_ARROW:
         if(actualmol) actualmol->show_freq_arrow = bit.show_freq_arrow;
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case FREQ_SCALE:
         if(actualmol) actualmol->sc_freq_ar = 5 * sc_freq_ar;
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case FREQ:
         if(!actualmol || !actualmol->vibration || !actualmol->n_frequencies){
            logprint("No vibrations loaded");
            update_logs();
            return;
         }
         select_vib();
      break;
      case SOLID_ARROW:
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case CANCEL:
         freqglui->close();
         freqglui = NULL;
      break;
    }
}


void freq_interf(void)
{
   char str[100];

   if(!actualmol) {
      logprint("No molecule loaded");
      update_logs();
      return;
   }
   if(!freqglui) {
      freqglui = GLUI_Master.create_glui\
         ("frequencies", 0, winposition[0], winposition[1]);
      GLUI_Panel *frame_panel = freqglui->add_panel("", true);
      GLUI_Panel *animfreq_panel = freqglui->add_panel_to_panel(frame_panel, "", true);
      GLUI_Button *animate = freqglui->add_button_to_panel(animfreq_panel, "animate", VIBRATION, freq_interf_cb);
      freqglui->add_column_to_panel(animfreq_panel, false);
      GLUI_Button *stop = freqglui->add_button_to_panel(animfreq_panel, "stop", STOP, freq_interf_cb);
      GLUI_Panel *arrowfreq_panel = freqglui->add_panel_to_panel(frame_panel, "", true);
      GLUI_Checkbox *arrow = freqglui->add_checkbox_to_panel\
             (arrowfreq_panel, "show arrows", &(bit.show_freq_arrow), FREQ_ARROW, freq_interf_cb);
      GLUI_Checkbox *sld_arrow = freqglui->add_checkbox_to_panel\
             (arrowfreq_panel, "solid arrows", &(bit.solid_arrow), SOLID_ARROW, freq_interf_cb);
      freqglui->add_column_to_panel(arrowfreq_panel, false);
      GLUI_Spinner *freq_sc_spinner = freqglui->add_spinner_to_panel\
             (arrowfreq_panel, "arrow scale-factor", GLUI_SPINNER_FLOAT, &(sc_freq_ar), FREQ_SCALE, freq_interf_cb);
      freq_sc_spinner->set_float_limits(0, 1.0);
      GLUI_Spinner *speed_spinner = freqglui->add_spinner_to_panel\
             (arrowfreq_panel, "slow down", GLUI_SPINNER_INT, &(speed_live_var));
      speed_spinner->set_float_limits(1, 100);
      GLUI_Panel *choose_panel = freqglui->add_panel_to_panel(frame_panel, "", true);
      if(actualvib) sprintf(str, "%5.0f/cm %s", actualvib->frequency, actualvib->type); 
      else sprintf(str, "no frequency selected");
      freqtext = freqglui->add_statictext_to_panel(choose_panel, str);
      freqglui->add_column_to_panel(choose_panel, false);
      freq_choose = freqglui->add_button_to_panel(choose_panel, "choose", FREQ, freq_interf_cb);
      if(vibrating) freq_choose->disable();

      freqglui->add_statictext("");
      GLUI_Button *cancel = freqglui->add_button("cancel", CANCEL, freq_interf_cb);
   }
   else {
      freqglui->show();
      glutSetWindow(freqglui->get_glut_window_id());
      glutPopWindow();
   }
}

void play_interf_cb(int key)
{
    char str[50];
    if(key == CANCEL) {
         playglui->close();
         playglui = NULL;
         mainglui_top->enable();
         glutSetWindow(mainwin);
         glutSetMenu(main_menu_id);
         glutAttachMenu(GLUT_RIGHT_BUTTON);
         return;
    }
      
    if(!dynamics.trajectory){
       showinfobox("no structures to animate!"); return;
    }
    if(actualmol != dynamics.molecule){
       showinfobox("wrong molecule!"); return;
    }
    if(!dynamics.end) dynamics.end = dynamics.ntotalsteps - 1;

    switch (key) {
      case PLAY_F:
         playing = 1;
         dynamics.direction =  1;
         drawingtime = (drawingtime > 10) ? drawingtime : 10;
         dynamics.stepsize = dynamics.timestep = 100.0/drawingtime;
         coord_cancel->disable();
         glutSetWindow(mainwin);
         if(bit.tgl_idle) GLUI_Master.set_glutIdleFunc(idle_coord);
         else glutIdleFunc(idle_coord);
      break;
      case PLAY_B:
         playing = 1;
         dynamics.direction =  -1;
         drawingtime = (drawingtime > 10) ? drawingtime : 10;
         dynamics.stepsize = dynamics.timestep = 100.0/drawingtime;
         coord_cancel->disable();
         glutSetWindow(mainwin);
         if(bit.tgl_idle) GLUI_Master.set_glutIdleFunc(idle_coord);
         else glutIdleFunc(idle_coord);
      break;
      case PLAY_FIRST:
         dynamicsGoto(dynamics.current = dynamics.start);
         if(playglui){
            sprintf(str, "structure # = %d", dynamics.current + 1);
            play_no->set_text(str);
         }
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case PLAY_LAST:
         dynamicsGoto(dynamics.current = dynamics.end);
         if(playglui){
            sprintf(str, "structure # = %d", dynamics.current + 1);
            play_no->set_text(str);
         }
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case NEXT_F:
         if(dynamics.current < dynamics.end) 
               dynamicsGoto(++dynamics.current);
         if(playglui){
            sprintf(str, "structure # = %d", dynamics.current + 1);
            play_no->set_text(str);
         }
         glutSetWindow(mainwin);
         glutPostRedisplay();
         break;
         case NEXT_B:
         if(dynamics.current > dynamics.start) 
               dynamicsGoto(--dynamics.current);
         if(playglui){
            sprintf(str, "structure # = %d", dynamics.current + 1);
            play_no->set_text(str);
         }
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case STOP:
         playing = 0;
         if(playglui){
            if(dynamics.current < 0) sprintf(str, "structure # = %d", 1);
            else if (dynamics.current >= dynamics.ntotalsteps) 
                 sprintf(str, "structure # = %d", dynamics.ntotalsteps);
            else sprintf(str, "structure # = %d", dynamics.current + 1);
            play_no->set_text(str);
         }
         coord_cancel->enable();
         glutSetWindow(mainwin);
         if(bit.tgl_idle) GLUI_Master.set_glutIdleFunc(NULL);
         else glutIdleFunc(NULL);
      break;
    }
}

void play_interf(void)
{
   char str[50];

   if(!actualmol) {
      logprint("No molecule loaded");
      update_logs();
      return;
   }
   if(!dynamics.trajectory){
      showinfobox("no structures to animate!"); return;
   }
   if(actualmol != dynamics.molecule){
      showinfobox("wrong active molecule!"); return;
   }

   if(!playglui) {
      playglui = GLUI_Master.create_glui\
         ("play", 0, winposition[0], winposition[1]);
      GLUI_Panel *frame_panel = playglui->add_panel("");
      sprintf(str, "structure # = %d", dynamics.current + 1);
      play_no = playglui->add_statictext_to_panel(frame_panel, str);
      play_no->set_alignment(GLUI_ALIGN_CENTER);
      GLUI_Panel *play_panel = playglui->add_panel_to_panel(frame_panel, "", true);
      GLUI_Panel *play_panel2 = playglui->add_panel_to_panel(play_panel, "", false);
      GLUI_Button *play_f = playglui->add_button_to_panel(play_panel2, "play >", PLAY_F, play_interf_cb);
      GLUI_Button *play_first = playglui->add_button_to_panel(play_panel2, "first", PLAY_FIRST, play_interf_cb);
      GLUI_Button *next_f = playglui->add_button_to_panel(play_panel2, "next >|", NEXT_F, play_interf_cb);
      playglui->add_column_to_panel(play_panel2, false);
      GLUI_Button *play_b = playglui->add_button_to_panel(play_panel2, "play <", PLAY_B, play_interf_cb);
      GLUI_Button *play_last = playglui->add_button_to_panel(play_panel2, "last", PLAY_LAST, play_interf_cb);
      GLUI_Button *next_b = playglui->add_button_to_panel(play_panel2, "next |<", NEXT_B, play_interf_cb);
      GLUI_Button *stop = playglui->add_button_to_panel(play_panel, "stop", STOP, play_interf_cb);
      GLUI_RadioGroup *play_rbgrp = playglui->add_radiogroup_to_panel\
              (play_panel, &dynamics.runtype);
      play_rbgrp->set_int_val(0);
      GLUI_RadioButton *single = playglui->add_radiobutton_to_group\
              (play_rbgrp, "single");
      GLUI_RadioButton *loop = playglui->add_radiobutton_to_group\
              (play_rbgrp, "loop");
      GLUI_RadioButton *swing = playglui->add_radiobutton_to_group\
              (play_rbgrp, "swing");
      GLUI_Panel *attr_panel = playglui->add_panel_to_panel(frame_panel, "", true);
      GLUI_Checkbox *keep_b = playglui->add_checkbox_to_panel\
             (attr_panel, "keep bonds", &(bit.keepbonds));
      playglui->add_column_to_panel(attr_panel, false);
      GLUI_Checkbox *superimpose = playglui->add_checkbox_to_panel\
             (attr_panel, "superimpose", &(bit.superimpose));

      playglui->add_statictext("");
      coord_cancel = playglui->add_button("cancel", CANCEL, play_interf_cb);
      
   }
   else {
      playglui->show();
      glutSetWindow(playglui->get_glut_window_id());
      glutPopWindow();
   }

   mainglui_top->disable();
   glutSetWindow(mainwin);
   glutSetMenu(truncC_menu_id);
   glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void box_interf_cb(int key)
{
    static float zoom = 0;
    float zoomfac = 0, x, y, z, min;

    x = box.x2 - box.x1;
    y = box.y2 - box.y1;
    z = box.z2 - box.z1;
    min = (x < y) ? x : y;
    min = (min < z) ? min : z;
    switch (key) {
      case ZOOMBOX:
         zoomfac = (zoombox - zoom);
         if(zoombox != zoom) {
            if(zoomfac < 0 && min < 0.1) zoomfac = 0;
            box.x1 -= zoomfac;
            box.x2 += zoomfac;
            box.y1 -= zoomfac;
            box.y2 += zoomfac;
            box.z1 -= zoomfac;
            box.z2 += zoomfac;
            zoom = zoombox;
         }
      break;
      case RESETBOX:
         reset_box();
         zoombox = zoom = 0;
         boxglui->sync_live();
      break;
      case BOX:
         if(box.x2 - 0.1 <= box.x1) xmax_spinner->set_float_val(box.x1 + 0.1);
         if(box.y2 - 0.1 <= box.y1) ymax_spinner->set_float_val(box.y1 + 0.1);
         if(box.z2 - 0.1 <= box.z1) zmax_spinner->set_float_val(box.z1 + 0.1);
      break;
      case BOXVAL_UPDATE:
          boxglui->sync_live();
      break;
      case CANCEL:
         boxglui->close();
         boxglui = NULL;
         mainglui_top->enable();
         glutSetWindow(mainwin);
         glutSetMenu(main_menu_id);
         glutAttachMenu(GLUT_RIGHT_BUTTON);
      break;
    }
}


void box_interf(void)
{
   if(!actualmol) {
      logprint("No molecule loaded");
      update_logs();
      return;
   }
   if(!boxglui) {
      boxglui = GLUI_Master.create_glui\
         ("box", 0, winposition[0], winposition[1]);
      boxglui->set_main_gfx_window(mainwin);
      GLUI_Panel *frame_panel = boxglui->add_panel("");
      GLUI_Checkbox *showbox = boxglui->add_checkbox_to_panel\
             (frame_panel, "show box", &(actualmol->box_on));
      showbox->set_alignment(GLUI_ALIGN_CENTER);

      GLUI_Panel *box_panel = boxglui->add_panel_to_panel(frame_panel, "", false);

      GLUI_Panel *box_panel2 = boxglui->add_panel_to_panel(box_panel, "", false);
      GLUI_Spinner *cubesize_spinner = boxglui->add_spinner_to_panel\
             (box_panel2, "cubesize", GLUI_SPINNER_FLOAT, &(box.cubesize));
      cubesize_spinner->set_float_limits(0.01, 2.5);
      boxglui->add_column_to_panel(box_panel2, false);
      zoom_spinner = boxglui->add_spinner_to_panel\
             (box_panel2, "zoom box", GLUI_SPINNER_FLOAT, &zoombox, ZOOMBOX, box_interf_cb);

      GLUI_Button *boxval_update = boxglui->add_button_to_panel\
             (box_panel, "update values", BOXVAL_UPDATE, box_interf_cb);

      GLUI_Panel *box_panel3 = boxglui->add_panel_to_panel(box_panel, "", false);
      xmin_spinner = boxglui->add_spinner_to_panel\
             (box_panel3, "x-min", GLUI_SPINNER_FLOAT, &(box.x1), BOX, box_interf_cb);
      ymin_spinner = boxglui->add_spinner_to_panel\
             (box_panel3, "y-min", GLUI_SPINNER_FLOAT, &(box.y1), BOX, box_interf_cb);
      zmin_spinner = boxglui->add_spinner_to_panel\
             (box_panel3, "z-min", GLUI_SPINNER_FLOAT, &(box.z1), BOX, box_interf_cb);
      boxglui->add_column_to_panel(box_panel3, false);
      xmax_spinner = boxglui->add_spinner_to_panel\
             (box_panel3, "x-max", GLUI_SPINNER_FLOAT, &(box.x2), BOX, box_interf_cb);
      ymax_spinner = boxglui->add_spinner_to_panel\
             (box_panel3, "y-max", GLUI_SPINNER_FLOAT, &(box.y2), BOX, box_interf_cb);
      zmax_spinner = boxglui->add_spinner_to_panel\
             (box_panel3, "z-max", GLUI_SPINNER_FLOAT, &(box.z2), BOX, box_interf_cb);

      GLUI_Button *resetbox = boxglui->add_button_to_panel(frame_panel, "reset box", RESETBOX, box_interf_cb);
      boxglui->add_statictext("");
      GLUI_Button *cancel = boxglui->add_button("cancel", CANCEL, box_interf_cb);
   }
   else {
      boxglui->show();
      glutSetWindow(boxglui->get_glut_window_id());
      glutPopWindow();
   }
   mainglui_top->disable();
   glutSetWindow(mainwin);
   glutSetMenu(truncB_menu_id);
   glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void connolly_interf_cb(int key)
{
    switch(key)
    {
      case CONNOLLY:
         connolly();
      break;
      case CONNECT:
         triangglui_interf();
      break;
      case CANCEL:
         connollyglui->close();
         connollyglui = NULL;
      break;
    }
}

void connolly_interf(void)
{
   if(!actualmol) {
      logprint("No molecule loaded");
      update_logs();
      return;
   }
   if(!connollyglui) {
      connollyglui = GLUI_Master.create_glui\
         ("connolly", 0, winposition[0], winposition[1]);
      GLUI_Panel *frame_panel = connollyglui->add_panel("");
      GLUI_Spinner *dot_dens_spinner = connollyglui->add_spinner_to_panel\
          (frame_panel, "dot density", GLUI_SPINNER_FLOAT, &dot_density);
      GLUI_Spinner *probe_rad_spinner = connollyglui->add_spinner_to_panel\
          (frame_panel, "probe radius", GLUI_SPINNER_FLOAT, &probe_radius);
      GLUI_Spinner *r_incr_spinner = connollyglui->add_spinner_to_panel\
          (frame_panel, "R increment", GLUI_SPINNER_FLOAT, &radius_increment);
      connollyglui->add_statictext_to_panel(frame_panel, "");
      GLUI_Checkbox *burried = connollyglui->add_checkbox_to_panel\
          (frame_panel, "burried surf", &(conflag.bury));
      GLUI_Checkbox *longout = connollyglui->add_checkbox_to_panel\
          (frame_panel, "long output", &(conflag.lon));
      connollyglui->add_statictext_to_panel(frame_panel, "");
      GLUI_Button *go = connollyglui->add_button_to_panel\
          (frame_panel, "go", CONNOLLY, connolly_interf_cb);
      GLUI_Button *triangulate = connollyglui->add_button_to_panel\
          (frame_panel, "dot triangulate", CONNECT, connolly_interf_cb);
      connollyglui->add_statictext("");
      GLUI_Button *cancel = connollyglui->add_button("cancel", CANCEL, connolly_interf_cb);
   }
   else {
      connollyglui->show();
      glutSetWindow(connollyglui->get_glut_window_id());
      glutPopWindow();
   }
}


void triang_interf_cb(int key)
{
    switch(key)
    {
      case CONNECT:
         glutSetWindow(mainwin);
         tridot();
         glutPostRedisplay();
      break;
      case CANCEL:
         triangglui->close();
         triangglui = NULL;
      break;
    }
}


void triangglui_interf(void)
{
   if(!actualmol) {
      logprint("No molecule loaded");
      update_logs();
      return;
   }
   if(!triangglui) {
      triangglui = GLUI_Master.create_glui\
         ("triangulate", 0, winposition[0], winposition[1]);
      GLUI_Panel *frame_panel = triangglui->add_panel("");
      GLUI_Spinner *dotcubesize_spinner = triangglui->add_spinner_to_panel\
          (frame_panel, "dot cubesize", GLUI_SPINNER_FLOAT, &dotcubesize);
      dotcubesize_spinner->set_float_limits(0.5*dotcubesize, 5*dotcubesize);
      GLUI_Spinner *angle_spinner = triangglui->add_spinner_to_panel\
          (frame_panel, "angle", GLUI_SPINNER_FLOAT, &angle_degrees);
      angle_spinner->set_float_limits(1, 180);
      GLUI_Spinner *dihedral_spinner = triangglui->add_spinner_to_panel\
          (frame_panel, "dihedral", GLUI_SPINNER_FLOAT, &dihedral_degrees);
      dihedral_spinner->set_float_limits(1, 180);
      triangglui->add_statictext_to_panel(frame_panel, "");
      GLUI_Button *go = triangglui->add_button_to_panel(frame_panel, "go", CONNECT, triang_interf_cb);
      triangglui->add_statictext("");
      GLUI_Button *cancel = triangglui->add_button("cancel", CANCEL, triang_interf_cb);
   }
   else {
      triangglui->show();
      glutSetWindow(triangglui->get_glut_window_id());
      glutPopWindow();
   }
}


void label_interf_cb(int key)
{
    switch(key)
    {
       case LABELS:
          if(actualmol) actualmol->labels = bit.labels;
       break;
       case ATM_CHAR:
          if(actualmol) actualmol->atm_char = bit.atm_char;
       break;
       case ATM_SPIN:
          if(actualmol) actualmol->atm_spin = bit.atm_spin;
       break;
       case CANCEL:
          labelglui->close();
          labelglui = NULL;
       break;
    }
}

void label_interf(void)
{
   if(!actualmol) {
      logprint("No molecule loaded");
      update_logs();
      return;
   }
   if(!labelglui) {
      labelglui = GLUI_Master.create_glui\
         ("labels", 0, winposition[0], winposition[1]);
      labelglui->set_main_gfx_window(mainwin);
      GLUI_Panel *label_panel = labelglui->add_panel("", true);
      GLUI_Checkbox *labels = labelglui->add_checkbox_to_panel\
             (label_panel, "atom labels", &(bit.labels), LABELS, label_interf_cb);
      GLUI_Checkbox *charg = labelglui->add_checkbox_to_panel\
             (label_panel, "atom charges", &(bit.atm_char), ATM_CHAR, label_interf_cb);
      GLUI_Checkbox *spin = labelglui->add_checkbox_to_panel\
             (label_panel, "atom spin", &(bit.atm_spin), ATM_SPIN, label_interf_cb);
      labelglui->add_statictext_to_panel(label_panel, "");
      GLUI_Checkbox *top = labelglui->add_checkbox_to_panel\
             (label_panel, "labels on top", &(bit.lbl_on_top));
      GLUI_Checkbox *nomeasurlab = labelglui->add_checkbox_to_panel\
             (label_panel, "hide measure labels", &(bit.no_measure_lbl));
      labelglui->add_column_to_panel(label_panel, true);
      labelglui->add_statictext_to_panel(label_panel, "global label size");
      GLUI_RadioGroup *size_rbgrp = labelglui->add_radiogroup_to_panel\
              (label_panel, &label_size);
      GLUI_RadioButton *small = labelglui->add_radiobutton_to_group\
              (size_rbgrp, "S");
      GLUI_RadioButton *medium = labelglui->add_radiobutton_to_group\
              (size_rbgrp, "M");
      GLUI_RadioButton *large = labelglui->add_radiobutton_to_group\
              (size_rbgrp, "L");
      GLUI_RadioButton *xlarge = labelglui->add_radiobutton_to_group\
              (size_rbgrp, "XL");
      GLUI_RadioButton *scalable = labelglui->add_radiobutton_to_group\
              (size_rbgrp, "scalable");
      labelglui->add_statictext("");
      GLUI_Button *cancel = labelglui->add_button("cancel", CANCEL, label_interf_cb);
   }
   else {
      labelglui->show();
      glutSetWindow(labelglui->get_glut_window_id());
      glutPopWindow();
   }
}


void xyzformatglui_cb(int key)
{
    switch (key) {
       case OK:
          file_select(action_key, xyz_ext);
       case CANCEL:
          action_key = 0;
          xyzformatglui->close();
          xyzformatglui = NULL;
       break;
    }
}

void xyzformat_interf(void)
{
   if(!xyzformatglui) {
      xyzformatglui = GLUI_Master.create_glui\
         ("xyz format", 0, winposition[0], winposition[1]);
      GLUI_Panel *frame_panel = xyzformatglui->add_panel("");
      GLUI_RadioGroup *format_rbgrp = xyzformatglui->add_radiogroup_to_panel\
              (frame_panel, &xyzformat);
      GLUI_RadioButton *sfff = xyzformatglui->add_radiobutton_to_group\
              (format_rbgrp, "%s %f %f %f");
      GLUI_RadioButton *dfff = xyzformatglui->add_radiobutton_to_group\
              (format_rbgrp, "%d %f %f %f");
      GLUI_RadioButton *fffs = xyzformatglui->add_radiobutton_to_group\
              (format_rbgrp, "%f %f %f %s");
      GLUI_RadioButton *fffd = xyzformatglui->add_radiobutton_to_group\
              (format_rbgrp, "%f %f %f %d");
      xyzformatglui->add_column_to_panel(frame_panel, true);
      xyzformatglui->add_statictext_to_panel(frame_panel, "%s: atomic symbol");
      xyzformatglui->add_statictext_to_panel(frame_panel, "%d: atomic number");
      xyzformatglui->add_statictext_to_panel(frame_panel, "%f: xyz coordinates");
      ok_cancel(xyzformatglui, xyzformatglui_cb);
   }
   else {
      xyzformatglui->show();
      glutSetWindow(xyzformatglui->get_glut_window_id());
      glutPopWindow();
   }
}


void t41content_cb(int key)
{
    switch (key) {
       case OK:
          file_select(action_key, t41_ext);
       case CANCEL:
          action_key = 0;
          t41contentglui->close();
          t41contentglui = NULL;
       break;
    }
}

void t41content_interf(void)
{
   if(!t41contentglui) {
      t41contentglui = GLUI_Master.create_glui\
         ("t41 content", 0, winposition[0], winposition[1]);
      GLUI_Panel *frame_panel = t41contentglui->add_panel("");
      GLUI_EditText *t41content = t41contentglui->add_edittext_to_panel\
         (frame_panel, "tape41 content", GLUI_EDITTEXT_TEXT, t41cont);
      t41content->set_w(200);
      t41contentglui->add_column_to_panel(frame_panel, true);
      t41contentglui->add_statictext_to_panel(frame_panel, "enter:");
      t41contentglui->add_statictext_to_panel(frame_panel, "\"Density\",");
      t41contentglui->add_statictext_to_panel(frame_panel, "\"Spin\" if A and B density present");
      t41contentglui->add_statictext_to_panel(frame_panel, "or the orbital e.g. \"SCF_B2\"");
      t41contentglui->add_statictext_to_panel(frame_panel, "(check the tape41 file!)");
      ok_cancel(t41contentglui, t41content_cb);
   }
   else {
      t41contentglui->show();
      glutSetWindow(t41contentglui->get_glut_window_id());
      glutPopWindow();
   }
}


void render_cb(int key)
{
    float ratio;
    
    switch (key) {
       case ADJUST_IMG_W:
          ratio = (float)xsize/(float)ysize;
          output_width = (output_width / 4) * 4;
          if(keep_ratio) output_height = ((int)(output_width / ratio) / 4) * 4;
          renderglui->sync_live();
       break;
       case ADJUST_IMG_H:
          ratio = (float)xsize/(float)ysize;
          output_height = (output_height / 4) * 4;
          if(keep_ratio) output_width = ((int)(output_height * ratio) / 4) * 4;
          renderglui->sync_live();
       break;
       case RESET:
          output_width = xsize;
          output_height = ysize;
          renderglui->sync_live();
       break;
       case EN_DIS_ABLE:
          if (framenbr_spinner) {
             switch (render_action_live_var) {
                case 0:
                   framenbr_spinner->disable();
                break;
                case 1:
                case 2:
                   framenbr_spinner->enable();
                break;
             }
          }
          if (img_ff_live_var == 4) jpegqual_spinner->enable();
          else jpegqual_spinner->disable();
       break;
       case OK:
          switch (img_ff_live_var) {
             case 0:
                switch (render_action_live_var) {
                   case 0:
                      action_key = RGB_CURR_IMG;
                   break;
                   case 1:
                      action_key = RGB_FREQ_IMG;
                   break;
                   case 2:
                      action_key = RGB_STRC_IMG;
                   break;
                }
                file_select(action_key, rgb_ext);
             break;
             case 1:
                switch (render_action_live_var) {
                   case 0:
                      action_key = TIF_CURR_IMG;
                   break;
                   case 1:
                      action_key = TIF_FREQ_IMG;
                   break;
                   case 2:
                      action_key = TIF_STRC_IMG;
                   break;
                }
                file_select(action_key, tif_ext);
             break;
             case 2:
                switch (render_action_live_var) {
                   case 0:
                      action_key = TIFLZW_CURR_IMG;
                   break;
                   case 1:
                      action_key = TIFLZW_FREQ_IMG;
                   break;
                   case 2:
                      action_key = TIFLZW_STRC_IMG;
                   break;
                }
                file_select(action_key, tif_ext);
             break;
             case 3:
                switch (render_action_live_var) {
                   case 0:
                      action_key = TIFPACK_CURR_IMG;
                   break;
                   case 1:
                      action_key = TIFPACK_FREQ_IMG;
                   break;
                   case 2:
                      action_key = TIFPACK_STRC_IMG;
                   break;
                }
                file_select(action_key, tif_ext);
             break;
             case 4:
                switch (render_action_live_var) {
                   case 0:
                      action_key = JPEG_CURR_IMG;
                   break;
                   case 1:
                      action_key = JPEG_FREQ_IMG;
                   break;
                   case 2:
                      action_key = JPEG_STRC_IMG;
                   break;
                }
                file_select(action_key, jpeg_ext);
             break;
             case 5:
                switch (render_action_live_var) {
                   case 0:
                      action_key = PPM_CURR_IMG;
                   break;
                   case 1:
                      action_key = PPM_FREQ_IMG;
                   break;
                   case 2:
                      action_key = PPM_STRC_IMG;
                   break;
                }
                file_select(action_key, ppm_ext);
             break;

           }
          renderglui->close();
          renderglui = NULL;
          framenbr_spinner = NULL;
       break;
       case CANCEL:
          action_key = 0;
          renderglui->close();
          renderglui = NULL;
          framenbr_spinner = NULL;
       break;
    }
}

void render_interf(void)
{
   if(!renderglui) {
      renderglui = GLUI_Master.create_glui\
         ("render", 0, winposition[0], winposition[1]);

      GLUI_Panel *image_panel = renderglui->add_panel("", true);
      GLUI_Spinner *imgwidth_spinner = renderglui->add_spinner_to_panel\
          (image_panel, "image width ", GLUI_SPINNER_INT, &output_width, ADJUST_IMG_W, render_cb);
      imgwidth_spinner->set_int_limits(16, 32000);
      GLUI_Spinner *imgheight_spinner = renderglui->add_spinner_to_panel\
          (image_panel, "image height", GLUI_SPINNER_INT, &output_height, ADJUST_IMG_H, render_cb);
      imgheight_spinner->set_int_limits(16, 32000);
      renderglui->add_column_to_panel(image_panel, false);
      GLUI_Checkbox *keep = renderglui->add_checkbox_to_panel\
             (image_panel, "keep ratio", &keep_ratio);
      renderglui->add_statictext_to_panel(image_panel,"");
      GLUI_Button *reset = renderglui->add_button_to_panel(image_panel, "reset", RESET, render_cb);

      GLUI_Panel *render_panel = renderglui->add_panel("", true);
      renderglui->add_statictext_to_panel(render_panel,"file format");
      GLUI_RadioGroup *format_rbgrp = renderglui->add_radiogroup_to_panel\
              (render_panel, &img_ff_live_var, EN_DIS_ABLE, render_cb);
      GLUI_RadioButton *rgb = renderglui->add_radiobutton_to_group\
              (format_rbgrp, "rgb");
      GLUI_RadioButton *tiff = renderglui->add_radiobutton_to_group\
              (format_rbgrp, "tiff");
      GLUI_RadioButton *tifflzw = renderglui->add_radiobutton_to_group\
              (format_rbgrp, "tiff lzw");
      GLUI_RadioButton *tiffpack = renderglui->add_radiobutton_to_group\
              (format_rbgrp, "tiff packbits");
      GLUI_RadioButton *jpeg = renderglui->add_radiobutton_to_group\
              (format_rbgrp, "jpeg");
      GLUI_RadioButton *ppm = renderglui->add_radiobutton_to_group\
              (format_rbgrp, "ppm (low memory!)");
      renderglui->add_statictext_to_panel(render_panel,"");
      jpegqual_spinner = renderglui->add_spinner_to_panel\
          (render_panel, "jpeg quality ", GLUI_SPINNER_INT, &jpeg_qual);
      jpegqual_spinner->set_int_limits(10, 100);
      render_cb(EN_DIS_ABLE);

      renderglui->add_column_to_panel(render_panel, true);
      renderglui->add_statictext_to_panel(render_panel,"action");
      GLUI_RadioGroup *render_action_rbgrp = renderglui->add_radiogroup_to_panel\
              (render_panel, &render_action_live_var, EN_DIS_ABLE, render_cb);
      GLUI_RadioButton *current = renderglui->add_radiobutton_to_group\
              (render_action_rbgrp, "current image");
      GLUI_RadioButton *freq = renderglui->add_radiobutton_to_group\
              (render_action_rbgrp, "freq animation");
      GLUI_RadioButton *structure = renderglui->add_radiobutton_to_group\
              (render_action_rbgrp, "structure animation");
      renderglui->add_statictext_to_panel(render_panel,"");
      framenbr_spinner = renderglui->add_spinner_to_panel\
          (render_panel, "start frame numbering at ", GLUI_SPINNER_INT, &frame_nbr);
      framenbr_spinner->set_int_limits(1, 999);
      render_cb(EN_DIS_ABLE);
      ok_cancel(renderglui, render_cb);
   }
   else {
      renderglui->show();
      glutSetWindow(renderglui->get_glut_window_id());
      glutPopWindow();
   }
}


void dipole_interf_cb(int key)
{
   switch (key) {
      case SHOW_DIPOLE:
         if(actualmol) actualmol->show_dipole = bit.show_dipole;
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case DIPOLE_SCALE:
         if(actualmol) actualmol->sc_dipole_ar = 2 * sc_dipole_ar;
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case CANCEL:
         dipoleglui->close();
         dipoleglui = NULL;
      break;
   }
}

void dipole_interf(void)
{
   char str[30];

   if(!actualmol) {
      logprint("No molecule loaded");
      update_logs();
      return;
   }

   if(actualmol->dipole) sprintf(str, "dipole moment: %7.4f Debye", actualmol->dipole->absolute);
   else sprintf(str, "no dipole moment info");
   
   if(!dipoleglui) {
      dipoleglui = GLUI_Master.create_glui\
         ("dipole moment", 0, winposition[0], winposition[1]);
      GLUI_Panel *frame_panel = dipoleglui->add_panel("", true);
      dipoletext = dipoleglui->add_statictext_to_panel(frame_panel, str);
      dipoleglui->add_statictext_to_panel(frame_panel, "");
      GLUI_Checkbox *show = dipoleglui->add_checkbox_to_panel\
             (frame_panel, "show dipole moment", &(bit.show_dipole), SHOW_DIPOLE, dipole_interf_cb);
      dipoleglui->add_statictext_to_panel(frame_panel, "");
      GLUI_Spinner *dipole_sc_spinner = dipoleglui->add_spinner_to_panel\
             (frame_panel, "dipole scale-factor", GLUI_SPINNER_FLOAT, &(sc_dipole_ar), DIPOLE_SCALE, dipole_interf_cb);
      dipole_sc_spinner->set_float_limits(0, 1.0);

      dipoleglui->add_statictext("");
      GLUI_Button *cancel = dipoleglui->add_button("cancel", CANCEL, dipole_interf_cb);
   }
   else {
      dipoleglui->show();
      glutSetWindow(dipoleglui->get_glut_window_id());
      glutPopWindow();
   }
}
