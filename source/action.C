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


/* event actions!
 * actions caused by the glutwin, menu or keyboard, start with a capital letter
*/

#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "general.h"
#include "browser.h"
#include "glutwin.h"
#include "manip.h"
#include "action.h"
#include "maininterf.h"
#include "chooseinterf.h"

unsigned int Toggle(int flag)
{
   if(flag) return 0;
   else return 1;   
}

void Quit(void)
{
   select_quit();
}

void Tgl_manip_all(void)
{
   if(bit.manip_all) {
      reset_globals();
      bit.manip_all  = 0;
   }
   else bit.manip_all  = 1;
}

void Tgl_depth(void)
{
   bit.depthcue = Toggle(bit.depthcue);
}

void Set_wire(void)
{
   bit.spacefill = bit.sticks = bit.ball_and_stick = 0;
   bit.wire_model = 1;
   if(actualmol){
      actualmol->spacefill = actualmol->sticks = actualmol->ball_and_stick = 0;
      actualmol->wire_model = 1;
   }
   molprop_live_var = 0;
}


void Set_ball_and_stick(void)
{
   bit.spacefill = bit.wire_model = bit.sticks = 0;
   bit.ball_and_stick = 1;
   if(actualmol){
      actualmol->spacefill = actualmol->wire_model = actualmol->sticks = 0;
      actualmol->ball_and_stick = 1;
   }
   molprop_live_var = 2;
}


void Set_stick(void)
{
   bit.spacefill = bit.wire_model = bit.ball_and_stick = 0;
   bit.sticks = 1;
   if(actualmol){
      actualmol->spacefill = actualmol->wire_model = actualmol->ball_and_stick = 0;
      actualmol->sticks = 1;
   }
   molprop_live_var = 1;
}


void Set_spacefill(void)
{
   bit.sticks = bit.wire_model = bit.ball_and_stick = 0;
   bit.spacefill = 1;
   if(actualmol){
      actualmol->sticks = actualmol->wire_model = actualmol->ball_and_stick = 0;
      actualmol->spacefill = 1;
   }
   molprop_live_var = 3;
}

void Tgl_fullscreen(void)
{
    static int fullscreen = 0;
    static int old_x;
    static int old_y;
    static int old_width;
    static int old_height;

    fullscreen = !fullscreen;
    if (fullscreen) {
       old_x = glutGet(GLUT_WINDOW_X);
       old_y = glutGet(GLUT_WINDOW_Y);
       old_width = glutGet(GLUT_WINDOW_WIDTH);
       old_height = glutGet(GLUT_WINDOW_HEIGHT);
       glutPopWindow();
       glutFullScreen();
    } else {
       glutReshapeWindow(old_width, old_height);
       glutPositionWindow(old_x, old_y);
    }
}

void Tgl_spinning(void)
{
   if(!actualmol) {
      logprint("No molecule loaded");
      update_logs();
      return;
   }
   if (!spinning) {
      logprint("Left mouse to rotate/spin");
      logprint("Click spin again to stop");
      update_logs();
      glutSetWindow(mainwin);
      glutMotionFunc(spinmotion);
      if(bit.tgl_idle) GLUI_Master.set_glutIdleFunc(idle_spin);
      else glutIdleFunc(idle_spin);
   }
   else {
      spinning = 0;
      glutSetWindow(mainwin);
      glutMotionFunc(motion);
      if(bit.tgl_idle) GLUI_Master.set_glutIdleFunc(NULL);
      else glutIdleFunc(NULL);
      logprint("Spinning stoped");
      update_logs();
   }
}


void Set_surf_transp(void)
{
   if(actualsurf) {
      actualsurf->surf_transp = bit.surf_transp;
      if(bit.surf_transp) {
         bit.chickenwire = bit.dot_surf = 0;
         actualsurf->chickenwire = actualsurf->dot_surf = 0;
      }
   }
   else {
      if(bit.surf_transp) {
         bit.chickenwire = bit.dot_surf = 0;
      }
   }
}

void Set_chickenwire(void)
{
   if(actualsurf) { 
      actualsurf->chickenwire =  bit.chickenwire;
      if(bit.chickenwire) {
         bit.surf_transp = bit.dot_surf = bit.flatshade = 0;
         actualsurf->surf_transp = actualsurf->dot_surf = actualsurf->flatshade = 0;
      }
   }
   else {
      if(bit.chickenwire){
         bit.surf_transp = bit.dot_surf = bit.flatshade = 0;
      }
   }
}

void Set_dot_surface(void)
{
   if(actualsurf) {
      actualsurf->dot_surf =  bit.dot_surf;
      if(bit.dot_surf) {
         bit.surf_transp = bit.chickenwire = bit.flatshade = 0;
         actualsurf->surf_transp = actualsurf->chickenwire = actualsurf->flatshade = 0;
      }
   }
   else {
      if(bit.dot_surf) {
         bit.chickenwire = bit.chickenwire = bit.flatshade = 0;
      }
   }
}

void Set_flatshade(void)
{
   if(actualsurf) {
      actualsurf->flatshade = bit.flatshade;
      if(bit.flatshade) {
         bit.chickenwire = bit.dot_surf = 0;
         actualsurf->chickenwire = actualsurf->dot_surf = 0;
      }
   }
   else {
      if(bit.flatshade) {
         bit.chickenwire = bit.dot_surf = 0;
      }
   }
}

void Tgl_idle(void)
{
   if(bit.tgl_idle) {
      bit.tgl_idle  = 0;
      Set_glut_idle();
   }
   else {
      bit.tgl_idle  = 1;
      Set_glui_idle();
   }
   update_logs();
}

void Set_glui_idle(void)
{
   if(vibrating) {
      GLUI_Master.set_glutIdleFunc(idle_freq);
   }
   else if(playing) {
      GLUI_Master.set_glutIdleFunc(idle_coord);
   }
   else if(spinning){
      GLUI_Master.set_glutIdleFunc(idle_spin);
   }
   else {
      GLUI_Master.set_glutIdleFunc(NULL);
   }
}

void Set_glut_idle(void)
{
   if(vibrating) {
      glutIdleFunc(idle_freq);
   }
   else if(playing) {
      glutIdleFunc(idle_coord);
   }
   else if(spinning){
      glutIdleFunc(idle_spin);
   }
   else {
      glutIdleFunc(NULL);
   }
}
