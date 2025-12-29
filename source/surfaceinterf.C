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
#include "surfaceinterf.h"
#include "glutwin.h"
#include "browser.h"
#include "macu.h"
#include "chooseinterf.h"
#include "maininterf.h"
#include "manip.h"
#include "utils.h"
#include "pick.h"
#include "action.h"
#include "menu.h"

GLUI *surfaceglui = NULL, *planeglui = NULL, *surf_optglui = NULL, *dotsldglui = NULL;
GLUI *projectglui = NULL;
GLUI_Rotation *rot;
GLUI_Translation *trans;
GLUI_RadioGroup *format_rbgrp;
int clip_move_rot = 0;
int clipplane_moving = 0;
int projectaxis = 0;
float plane_trans = 0.0;

void surfaceglui_cb(int key)
{
   switch (key) {
      case CANCEL:
         if(planeglui) {
            planeglui->close();
            planeglui = NULL;
         }
         if(dotsldglui) {
            dotsldglui->close();
            dotsldglui = NULL;
         }
         if(projectglui) {
            projectglui->close();
            projectglui = NULL;
         }
         if(surf_optglui) {
            surf_optglui->close();
            surf_optglui = NULL;
         }
         clip_move_rot = 0;
         surfaceglui->close();
         surfaceglui = NULL;
      break;
      case CANCEL_PLANE:
         planeglui->close();
         planeglui = NULL;
      break;
      case CANCEL_DOTSLD:
         dotsldglui->close();
         dotsldglui = NULL;
      break;
      case CANCEL_SURF_OPT:
         surf_optglui->close();
         surf_optglui = NULL;
      break;
   }

   if(!actualmol) {
      logprint("No molecule loaded");
      update_logs();
      return;
   }
   switch (key) {
      case DUMMY:
         printf("Dummy of surfaceglui\n");
      break;
      case LOAD_MACU:
         bit.both_signs = 0;
      case ADD_MACU:
      case SUBST_MACU:
      switch (format_rbgrp->get_int_val()) {
         case 0:
             file_select(key, macu_ext);
         break;
         case 1:
             file_select(key, gcube_ext);
         break;
         case 2:
             action_key = key;
             t41content_interf();
         break;
      }
      break;
      case UNLOAD_MACU:
         free_3D(actualmol);
      break;
      case RESET:
         cubemin = actualmol->cubemin;
         cubemax = actualmol->cubemax;
         cutoff = 0.0;
         surfaceglui->sync_live();
      break;
      case PICK_SURF:
          enter_pick(key, "Pick surface");
      break;
      case SWITCH_SURF:
         switch_surf();
         update_interface_flags();
      break;
      case ADD_SURF:
         cubes();
         update_interface_flags();
      break;
      case DEL_SURF:
         del_actualsurf(actualsurf);
         update_interface_flags();
      break;
      case GRID_VALUES:
         pick_surface(key);
         update_interface_flags();
      break;
      case XYZ_PLANE:
         if(!actualmol->plane) {
            logprint("No grid data loaded");
            update_logs();
            return;
         }
         if(!planeglui)makeplaneglui();
         else {
            planeglui->show();
            glutSetWindow(planeglui->get_glut_window_id());
            glutPopWindow();
         }
      break;
      case DOTSLD:
         if(!dotsldglui)makedotsldglui();
         else {
            dotsldglui->show();
            glutSetWindow(dotsldglui->get_glut_window_id());
            glutPopWindow();
         }
      break;
      case SURF_OPT:
         if(!surf_optglui)makesurf_optglui();
         else {
            surf_optglui->show();
            glutSetWindow(surf_optglui->get_glut_window_id());
            glutPopWindow();
         }
      break;
      case SURF_CLIP:
         if(actualsurf) actualsurf->surf_clip = bit.surf_clip;
      break;
      case CONNECT:
         triangglui_interf();
      break;
      case LOAD_DOTS:
         file_select(key, ms_ext);
      break;
      case LOAD_SLD:
         file_select(key, sld_ext);
      break;
      case LOAD_DOTVAL:
         file_select(key, dotval_ext);
      break;
      case WRITE_MS:
         action_key = key;
         file_select(key, ms_ext);
      break;
      case WRITE_SLD:
         action_key = key;
         file_select(key, sld_ext);
      break;
      case WRITE_DOTVAL:
         action_key = key;
         file_select(key, dotval_ext);
      break;
      case PROJECT_MACU:
         if(!projectglui)makeprojectglui();
         else {
            projectglui->show();
            glutSetWindow(projectglui->get_glut_window_id());
            glutPopWindow();
         }
      break;
      case FREE_PROJECTION:
         free_projection();
      break;
   }
 
}


void surf_optglui_cb(int key)
{
   switch (key) {
      case TRANSPARENT:
         Set_surf_transp();
      break;
      case CHICKENWIRE:
         Set_chickenwire();
      break;
      case FLATSHADE:
         Set_flatshade();
      break;
      case DOT_SURFACE:
         Set_dot_surface();
      break;
   }
   surf_optglui->sync_live();
   glutSetWindow(mainwin);
   glutPostRedisplay();
}


void makesurfaceglui(void)
{
   surfaceglui = GLUI_Master.create_glui("surface", 0, winposition[0], winposition[1]);
   surfaceglui->set_main_gfx_window(mainwin);

   GLUI_Panel *manipmacu_panel = surfaceglui->add_panel(""); 
   GLUI_Panel *manipmacu_panel2 = surfaceglui->add_panel_to_panel\
             (manipmacu_panel, "", GLUI_PANEL_NONE); 
   GLUI_Button *loadmacu = surfaceglui->add_button_to_panel\
             (manipmacu_panel2, "load", LOAD_MACU , surfaceglui_cb);
   GLUI_Button *addmacu = surfaceglui->add_button_to_panel\
             (manipmacu_panel2, "add", ADD_MACU , surfaceglui_cb);
   surfaceglui->add_column_to_panel(manipmacu_panel2, false);
   GLUI_Button *unload = surfaceglui->add_button_to_panel
             (manipmacu_panel2, "unload", UNLOAD_MACU, surfaceglui_cb);
   GLUI_Button *subtmacu = surfaceglui->add_button_to_panel\
             (manipmacu_panel2, "subtract", SUBST_MACU , surfaceglui_cb);
   format_rbgrp = surfaceglui->add_radiogroup_to_panel(manipmacu_panel);
   GLUI_RadioButton *macu_button = \
             surfaceglui->add_radiobutton_to_group(format_rbgrp, "macu");
   GLUI_RadioButton *cube_button = \
             surfaceglui->add_radiobutton_to_group(format_rbgrp, "gaussian cube");
   GLUI_RadioButton *t41_button = \
             surfaceglui->add_radiobutton_to_group(format_rbgrp, "adf tape41");

   surfaceglui->add_statictext("");

   GLUI_Panel *surf_panel = surfaceglui->add_panel("");
   GLUI_Panel *surf_panel2 = surfaceglui->add_panel_to_panel\
             (surf_panel, "", GLUI_PANEL_NONE);
   GLUI_Checkbox *bothsign = surfaceglui->add_checkbox_to_panel\
             (surf_panel2, "both signs", &(bit.both_signs));
   GLUI_Checkbox *addsurf = surfaceglui->add_checkbox_to_panel\
             (surf_panel2, "add surface", &(bit.addsurf));
   GLUI_Panel *surf_panel5 = surfaceglui->add_panel_to_panel\
             (surf_panel2, "");
   GLUI_Spinner *Vminspinner =  surfaceglui->add_spinner_to_panel\
             (surf_panel5, "V-min", GLUI_SPINNER_FLOAT, &cubemin);
   GLUI_Spinner *Vmaxspinner =  surfaceglui->add_spinner_to_panel\
             (surf_panel5, "V-max", GLUI_SPINNER_FLOAT, &cubemax);
   GLUI_Spinner *cutoffspinner =  surfaceglui->add_spinner_to_panel\
             (surf_panel5, "cutoff", GLUI_SPINNER_FLOAT, &cutoff);
   GLUI_Button *reset = surfaceglui->add_button_to_panel
             (surf_panel5, "reset", RESET, surfaceglui_cb);
   surfaceglui->add_column_to_panel(surf_panel2, false);
   GLUI_Button *pick_surf = surfaceglui->add_button_to_panel
             (surf_panel2, "pick surface", PICK_SURF, surfaceglui_cb);
   pick_surf->set_alignment(GLUI_ALIGN_LEFT);
   surfaceglui->add_statictext_to_panel(surf_panel2, "");
   GLUI_Button *switch_surf = surfaceglui->add_button_to_panel
             (surf_panel2, "switch surface", SWITCH_SURF, surfaceglui_cb);
   switch_surf->set_alignment(GLUI_ALIGN_LEFT);
   surfaceglui->add_statictext_to_panel(surf_panel2, "");
   GLUI_Button *create_surf = surfaceglui->add_button_to_panel
             (surf_panel2, "create surface", ADD_SURF, surfaceglui_cb);
   create_surf->set_alignment(GLUI_ALIGN_LEFT);
   surfaceglui->add_statictext_to_panel(surf_panel2, "");
   GLUI_Button *del_surf = surfaceglui->add_button_to_panel\
             (surf_panel2, "delete surface", DEL_SURF , surfaceglui_cb);
   del_surf->set_alignment(GLUI_ALIGN_LEFT);

   GLUI_Panel *surf_panel3 = surfaceglui->add_panel_to_panel(surf_panel, "");
   GLUI_Button *grid = surfaceglui->add_button_to_panel
             (surf_panel3, "grid value", GRID_VALUES, surfaceglui_cb);
   GLUI_Button *proj = surfaceglui->add_button_to_panel
             (surf_panel3, "project", PROJECT_MACU, surfaceglui_cb);
   surfaceglui->add_column_to_panel(surf_panel3, false);
   GLUI_Button *planes = surfaceglui->add_button_to_panel\
             (surf_panel3, "planes", XYZ_PLANE , surfaceglui_cb);
   GLUI_Button *delproj = surfaceglui->add_button_to_panel
             (surf_panel3, "del projection", FREE_PROJECTION, surfaceglui_cb);

   surfaceglui->add_statictext_to_panel(surf_panel, "");
   GLUI_RadioGroup *legend_rbgrp = surfaceglui->add_radiogroup_to_panel(surf_panel, &maplegend);
   legend_rbgrp->set_alignment(GLUI_ALIGN_CENTER);
   GLUI_RadioButton *no_leg_button = \
             surfaceglui->add_radiobutton_to_group(legend_rbgrp, "map color legend off");
   GLUI_RadioButton *surf_leg_button = \
             surfaceglui->add_radiobutton_to_group(legend_rbgrp, "show surface map color legend");
   GLUI_RadioButton *plane_leg_button = \
             surfaceglui->add_radiobutton_to_group(legend_rbgrp, "show plane map color legend");

   surfaceglui->add_statictext("");

   GLUI_Panel *surf_panel4 = surfaceglui->add_panel("");
   GLUI_Button *option = surfaceglui->add_button_to_panel\
             (surf_panel4, "surface option", SURF_OPT, surfaceglui_cb);
   surfaceglui->add_column_to_panel(surf_panel4, false);
   GLUI_Button *dotsld = surfaceglui->add_button_to_panel\
             (surf_panel4, "dot/sld", DOTSLD, surfaceglui_cb);
   surfaceglui->add_statictext("");

   GLUI_Panel *clip_panel = surfaceglui->add_panel("");
   GLUI_Checkbox *on_off_clip = surfaceglui->add_checkbox_to_panel\
             (clip_panel, "surface clip", &bit.surf_clip, SURF_CLIP, surfaceglui_cb);
   on_off_clip->set_alignment(GLUI_ALIGN_CENTER);
   surfaceglui->add_column_to_panel(clip_panel, false);
   GLUI_Checkbox *move_rot_clip = surfaceglui->add_checkbox_to_panel\
             (clip_panel, "move rot clip", &clip_move_rot);
   move_rot_clip->set_alignment(GLUI_ALIGN_CENTER);
   surfaceglui->add_statictext("");

   GLUI_Button *cancel = surfaceglui->add_button("cancel", CANCEL, surfaceglui_cb);
}

void makedotsldglui(void)
{
   dotsldglui = GLUI_Master.create_glui\
        ("dot/sld", 0, winposition[0] + 350, winposition[1] + 450);

   GLUI_Panel *triang_panel = dotsldglui->add_panel("");
   GLUI_Panel *rw_panel = dotsldglui->add_panel_to_panel(triang_panel, "", GLUI_PANEL_NONE);
   GLUI_Button *load_ms = dotsldglui->add_button_to_panel
             (rw_panel, "load ms", LOAD_DOTS, surfaceglui_cb);
   GLUI_Button *load_sld = dotsldglui->add_button_to_panel
             (rw_panel, "load sld", LOAD_SLD, surfaceglui_cb);
   GLUI_Button *load_dotval = dotsldglui->add_button_to_panel
             (rw_panel, "load dots with value", LOAD_DOTVAL, surfaceglui_cb);
   dotsldglui->add_column_to_panel(rw_panel, false);
   GLUI_Button *save_ms = dotsldglui->add_button_to_panel
             (rw_panel, "save ms", WRITE_MS, surfaceglui_cb);
   GLUI_Button *save_sld = dotsldglui->add_button_to_panel
             (rw_panel, "save sld", WRITE_SLD, surfaceglui_cb);
   GLUI_Button *save_dotval = dotsldglui->add_button_to_panel
             (rw_panel, "save dots with value", WRITE_DOTVAL, surfaceglui_cb);
   GLUI_Button *dot_trang = dotsldglui->add_button_to_panel
             (triang_panel, "dot triangulate", CONNECT, surfaceglui_cb);
   dotsldglui->add_statictext("");
   GLUI_Button *cancel = dotsldglui->add_button("cancel", CANCEL_DOTSLD, surfaceglui_cb);
}

void makesurf_optglui(void)
{
   surf_optglui = GLUI_Master.create_glui\
        ("surface options", 0, winposition[0] + 350, winposition[1] + 450);

   GLUI_Panel *option_panel = surf_optglui->add_panel("");
   GLUI_Checkbox *transparent = surf_optglui->add_checkbox_to_panel\
             (option_panel, "transparent", &(bit.surf_transp), TRANSPARENT, surf_optglui_cb);
   GLUI_Checkbox *dots = surf_optglui->add_checkbox_to_panel\
             (option_panel, "dots", &(bit.dot_surf), DOT_SURFACE, surf_optglui_cb);
   surf_optglui->add_column_to_panel(option_panel, false);
   GLUI_Checkbox *chicken = surf_optglui->add_checkbox_to_panel\
             (option_panel, "chickenwire", &(bit.chickenwire), CHICKENWIRE, surf_optglui_cb);
   GLUI_Checkbox *flat = surf_optglui->add_checkbox_to_panel\
             (option_panel, "flatshade", &(bit.flatshade), FLATSHADE, surf_optglui_cb);
   surf_optglui->add_statictext("");
   GLUI_Button *cancel = surf_optglui->add_button("cancel", CANCEL_SURF_OPT, surfaceglui_cb);
 
}

void planeglui_cb(int key)
{
   static float old_plane_trans = 0.0;

   if(!actualmol) return;
   switch (key) {
      case CUBEPLANES:
         actualmol->cubeplanes = bit.cubeplanes;
      break;
      case RESET:
         actualmol->plane_ix = actualmol->box.nx/2;
         actualmol->plane_iy = actualmol->box.ny/2;
         actualmol->plane_iz = actualmol->box.nz/2;
         actualmol->plane->alpha = 0.0;
         actualmol->plane->beta = 0.0;
         actualmol->plane->a[0] = actualmol->box.x1 + 
           (actualmol->box.x2-actualmol->box.x1)/((float)actualmol->box.nx-1.)*
              (actualmol->box.nx-1)*0.5;
         actualmol->plane->a[1] = actualmol->box.y1 + 
           (actualmol->box.y2-actualmol->box.y1)/((float)actualmol->box.ny-1.)*
              (actualmol->box.ny-1)*0.5;
         actualmol->plane->a[2] = actualmol->box.z1 + 
           (actualmol->box.z2-actualmol->box.z1)/((float)actualmol->box.nz-1.)*
              (actualmol->box.nz-1)*0.5;
         set_cutplane(actualmol);
         planeglui->sync_live();
      break;
      case MOVEPLANE:
         switch (cubeplane_axis) {
            case 0:
               if(plane_trans > old_plane_trans) {
                  actualmol->plane_ix++;
               }
               else {
                  actualmol->plane_ix--;
               }
               if(actualmol->plane_ix < 0) actualmol->plane_ix = 0;
               else if(actualmol->plane_ix >= actualmol->box.nx) actualmol->plane_ix = actualmol->box.nx - 1;
            break;
            case 1:
               if(plane_trans > old_plane_trans) {
                  actualmol->plane_iy++;
               }
               else {
                  actualmol->plane_iy--;
               }
               if(actualmol->plane_iy < 0) actualmol->plane_iy = 0;
               else if(actualmol->plane_iy >= actualmol->box.ny) actualmol->plane_iy = actualmol->box.ny - 1;
            break;
            case 2:
               if(plane_trans > old_plane_trans) {
                  actualmol->plane_iz++;
               }
               else {
                  actualmol->plane_iz--;
               }
               if(actualmol->plane_iz < 0) actualmol->plane_iz = 0;
               else if(actualmol->plane_iz >= actualmol->box.nz) actualmol->plane_iz = actualmol->box.nz - 1;
            break;
         }
         old_plane_trans = plane_trans;
      break;
      case SET_CUTPLANE:
         set_cutplane(actualmol);
      break;
   }
}


void makeplaneglui(void)
{
   planeglui = GLUI_Master.create_glui\
        ("plane", 0, winposition[0] + 200, winposition[1] + 350);
   planeglui->set_main_gfx_window(mainwin);
   GLUI_Panel *frame_panel = planeglui->add_panel("");
   GLUI_Panel *plane_panel = planeglui->add_panel_to_panel(frame_panel,"", GLUI_PANEL_NONE); 
   GLUI_Panel *plane_panel2 = planeglui->add_panel_to_panel(plane_panel, ""); 
   GLUI_RadioGroup *axisgrp = planeglui->add_radiogroup_to_panel\
              (plane_panel2, &cubeplane_axis);
   GLUI_RadioButton *x = planeglui->add_radiobutton_to_group(axisgrp, "x-axis");
   GLUI_RadioButton *y = planeglui->add_radiobutton_to_group(axisgrp, "y-axis");
   GLUI_RadioButton *z = planeglui->add_radiobutton_to_group(axisgrp, "z-axis");
   GLUI_RadioButton *def = planeglui->add_radiobutton_to_group(axisgrp, "def. plane:");
   planeglui->add_statictext_to_panel(plane_panel2, "      plane center position:");
   GLUI_Spinner *ax_spinner =  planeglui->add_spinner_to_panel\
             (plane_panel2, "x", GLUI_SPINNER_FLOAT, &actualmol->plane->a[0], SET_CUTPLANE, planeglui_cb);
   GLUI_Spinner *ay_spinner =  planeglui->add_spinner_to_panel\
             (plane_panel2, "y", GLUI_SPINNER_FLOAT, &actualmol->plane->a[1], SET_CUTPLANE, planeglui_cb);
   GLUI_Spinner *az_spinner =  planeglui->add_spinner_to_panel\
             (plane_panel2, "z", GLUI_SPINNER_FLOAT, &actualmol->plane->a[2], SET_CUTPLANE, planeglui_cb);
   planeglui->add_column_to_panel(plane_panel2, false);
   GLUI_Checkbox *show = planeglui->add_checkbox_to_panel\
             (plane_panel2, "show plane", &(bit.cubeplanes), CUBEPLANES, planeglui_cb);
   planeglui->add_statictext_to_panel(plane_panel2, "");
   planeglui->add_statictext_to_panel(plane_panel2, "");
   planeglui->add_statictext_to_panel(plane_panel2, "");
   planeglui->add_statictext_to_panel(plane_panel2, "rotate around:");
   GLUI_Spinner *alpha_spinner =  planeglui->add_spinner_to_panel\
             (plane_panel2, "y-axis", GLUI_SPINNER_FLOAT, &actualmol->plane->alpha, SET_CUTPLANE, planeglui_cb);
   alpha_spinner->set_float_limits(0, 360, GLUI_LIMIT_WRAP);
   GLUI_Spinner *beta_spinner =  planeglui->add_spinner_to_panel\
             (plane_panel2, "z-axis", GLUI_SPINNER_FLOAT, &actualmol->plane->beta, SET_CUTPLANE, planeglui_cb);
   beta_spinner->set_float_limits(0, 360, GLUI_LIMIT_WRAP);
   planeglui->add_column_to_panel(plane_panel, false);
   GLUI_Translation *moveplane = planeglui->add_translation_to_panel
                  (plane_panel, "move axis plane", GLUI_TRANSLATION_Z, &plane_trans, MOVEPLANE, planeglui_cb);
   GLUI_Button *reset = planeglui->add_button_to_panel(frame_panel, "reset position", RESET, planeglui_cb);
   planeglui->add_statictext("");
   GLUI_Button *cancel = planeglui->add_button("cancel", CANCEL_PLANE, surfaceglui_cb);
}

void projectglui_cb(int key)
{
   switch (key) {
      case PROJECT_MACU:
         project_macu(projectaxis);
         glutSetWindow(mainwin);
         glutPostRedisplay();
      break;
      case CANCEL:
         projectglui->close();
         projectglui = NULL;
      break;
   }
}

void makeprojectglui(void)
{
   projectglui = GLUI_Master.create_glui\
        ("project", 0, winposition[0] + 350, winposition[1] + 450);
   GLUI_Panel *frame_panel = projectglui->add_panel("");
   GLUI_RadioGroup *axisgrp = projectglui->add_radiogroup_to_panel(frame_panel, &projectaxis);
   GLUI_RadioButton *x = projectglui->add_radiobutton_to_group(axisgrp, "x-axis");
   GLUI_RadioButton *y = projectglui->add_radiobutton_to_group(axisgrp, "y-axis");
   GLUI_RadioButton *z = projectglui->add_radiobutton_to_group(axisgrp, "z-axis");
   projectglui->add_statictext_to_panel(frame_panel, "");
   GLUI_Button *go = projectglui->add_button_to_panel(frame_panel, "go", PROJECT_MACU, projectglui_cb);
   projectglui->add_statictext("");
   GLUI_Button *cancel = projectglui->add_button("cancel", CANCEL, projectglui_cb);
}
