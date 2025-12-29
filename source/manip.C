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


/* glut window manipulations */
#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "general.h"
#include "glutwin.h"
#include "trackball.h"
#include "drawing.h"
#include "utils.h"
#include "macu.h"
#include "cubeInterpol.h"
#include "action.h"
#include "manip.h"
#include "maininterf.h"
#include "surfaceinterf.h"
#include "pick.h"
#include "chooseinterf.h"
#include "menu.h"
#include "material.h"

/*----------------------------------*/
static float trackballcenter[2];
int picking_mode, spinning = 0;
/*----------------------------------*/
/* test */
/*----------------------------------*/

static long omx, omy, nmx, nmy;
static float angle = 5*M_PI/180.;
static int current_button = -1, modifiers;

void mouse(int button, int state, int x, int y)
{
   char str[20];
   int dologprint = 0;

   if(!firstmol) return;

   if (current_button == -1) {
      if (state == GLUT_DOWN) {
         omx = x;
         omy = y;
         current_button = button;
         modifiers = glutGetModifiers();

         if(current_button == GLUT_LEFT_BUTTON && modifiers & GLUT_ACTIVE_SHIFT) {
            if(modifiers & GLUT_ACTIVE_ALT) { strcpy(str, "z-translation"); dologprint = 1; }
            else if(modifiers & GLUT_ACTIVE_CTRL) { strcpy(str, "scaling"); dologprint = 1; }
         }
         if(current_button == GLUT_MIDDLE_BUTTON) {
            if(modifiers & GLUT_ACTIVE_SHIFT) { strcpy(str, "z-translation"); dologprint = 1; }
            else if(modifiers & GLUT_ACTIVE_CTRL) { strcpy(str, "scaling"); dologprint = 1; }
         }
         if(dologprint) {
            green_note(str);
         }
      }
   } else if (current_button == button) {
      if (state == GLUT_UP) {
         current_button = -1;
         update_logs();
      }
   }

   if (state == GLUT_UP && rapidmove) {
      rapidmove = 0;
      glutPostRedisplay();
   }
}

void picking_mouse(int button, int state, int x, int y)
{
   if(!firstmol) return;

   modifiers = glutGetModifiers();
   if(modifiers & GLUT_ACTIVE_SHIFT) {
      if (current_button == -1) {
         if (state == GLUT_DOWN) {
            omx = x;
            omy = y;
            current_button = button;
         }
      } else if (current_button == button) {
         if (state == GLUT_UP) {
            current_button = -1;
         }
      }
      if (state == GLUT_UP && rapidmove) {
         rapidmove = 0;
         glutPostRedisplay();
      }
  }
   else {
      switch (button) {
          case GLUT_LEFT_BUTTON:
             if(state == GLUT_DOWN) {
                pick_action(pick_key, x, y);
             }
          break;
          case GLUT_MIDDLE_BUTTON:
             if(state == GLUT_DOWN) {
                pick_action(DONE, x, y);
             }
          break;
      }
   }


}


void motion(int x, int y)
{
   if(!firstmol) return;

   nmx = x;
   nmy = y;

   if(drawingtime > 600) rapidmove = 1;
   if(actualsurf && clip_move_rot && actualsurf->surf_clip) {
      switch (current_button) {
         case GLUT_LEFT_BUTTON:
            if(modifiers & GLUT_ACTIVE_CTRL) {
               mouse_rotation(&omx, &omy, &nmx, &nmy);
            }
            else if(modifiers & GLUT_ACTIVE_SHIFT) {
               t[0] = (float)(nmx-omx)/(float)xsize*2.0;
               t[1] = (float)(omy-nmy)/(float)ysize*2.0;
               t[2] = 0;
               clipplane_moving = 1;
               mouse_clip_translation();
            }
            else {
               mouse_clip_rotation(&omx, &omy, &nmx, &nmy);
            }
         break;
         case GLUT_MIDDLE_BUTTON:
            t[0] = (float)(nmx-omx)/(float)xsize*2.0;
            t[1] = (float)(omy-nmy)/(float)ysize*2.0;
            t[2] = 0;
            clipplane_moving = 1;
            mouse_clip_translation();
         break;
      }
   }
   else {
      switch (current_button) {
         case GLUT_LEFT_BUTTON:
            if(modifiers & GLUT_ACTIVE_SHIFT) {
               if (modifiers & GLUT_ACTIVE_ALT) {
                  t[0] = t[1] = 0;
                  t[2] = (float)(omy-nmy)/(float)ysize;
                  mouse_translation();
               }
               else if (modifiers & GLUT_ACTIVE_CTRL) {
                  sfac *= pow(2.0, 2.0*(nmy-omy)/(float)ysize);
                  scale_line(xsize, ysize);
               } 
               else {
                  t[0] = (float)(nmx-omx)/(float)xsize*2.0;
                  t[1] = (float)(omy-nmy)/(float)ysize*2.0;
                  t[2] = 0;
                  mouse_translation();
               }
            }
            else {
               mouse_rotation(&omx, &omy, &nmx, &nmy);
            }
         break;
         case GLUT_MIDDLE_BUTTON:
            if (modifiers & GLUT_ACTIVE_SHIFT) {
               t[0] = t[1] = 0;
               t[2] = (float)(omy-nmy)/(float)ysize;
               mouse_translation();
            }
            else if (modifiers & GLUT_ACTIVE_CTRL) {
               sfac *= pow(2.0, 2.0*(nmy-omy)/(float)ysize);
               scale_line(xsize, ysize);
            } 
            else {
               t[0] = (float)(nmx-omx)/(float)xsize*2.0;
               t[1] = (float)(omy-nmy)/(float)ysize*2.0;
               t[2] = 0;
               mouse_translation();
            }
         break;
      }
   }

   omx = nmx;
   omy = nmy;
   glutPostRedisplay();
   clipplane_moving = 0;
}


void picking_motion(int x, int y)
{
   if(!firstmol) return;

   nmx = x;
   nmy = y;

   if(modifiers & GLUT_ACTIVE_SHIFT) {
      if(drawingtime > 600) rapidmove = 1;
      switch (current_button) {
         case GLUT_LEFT_BUTTON:
               mouse_rotation(&omx, &omy, &nmx, &nmy);
               omx = nmx;
               omy = nmy;
               glutPostRedisplay();
         break;
      }
   }
}


void idle_spin(void)
{
   if(!firstmol) return;

   if(spinning) {
      if(bit.manip_all){
         add_eulers(spinrot, global_rot, global_rot);
      }
      else add_eulers(spinrot, actualmol->rvec, actualmol->rvec);
      glutSetWindow(mainwin);
      glutPostRedisplay();
   }
}


void idle_freq(void)
{
   static float step;
   float factor;
   int i;
   register AtoM *ap;
   Vector *vw;

   step = M_PI/(10 * speed_live_var);
   factor = actualmol->sc_freq_ar * sin(freq_pos) / 2.5;
   freq_pos += step;
   for(ap = actualmol->firstatom, vw=coord_backup, i=0; ap; ap=ap->next, vw++, i++){
      ap->coord[0] = vw->x + actualvib->coord[i].x * factor;
      ap->coord[1] = vw->y + actualvib->coord[i].y * factor;
      ap->coord[2] = vw->z + actualvib->coord[i].z * factor;
   }
   glutSetWindow(mainwin);
   glutPostRedisplay();
}

void idle_coord(void)
{
   static int increment;
   char str[50];

   dynamics.isrunning = 1;
   increment = dynamics.direction;
      if((dynamics.current > dynamics.end) && (increment > 0)) {
         switch(dynamics.runtype) {
            case 0 : /* single run */
               dynamics.isrunning = 0;
               dynamics.current = dynamics.end;
               increment = 0;
               break;
            case 1 : /* loop */
               dynamics.current = dynamics.start - increment;
               break;
            case 2 : /* swing */
               dynamics.current = dynamics.end + increment;
               dynamics.direction = -1;
               increment = -increment;
               break;
         }
      }
      else if((dynamics.current < dynamics.start) && (increment < 0)) {
         switch(dynamics.runtype) {
            case 0 : /* single run */
               dynamics.isrunning = 0;
               dynamics.current = dynamics.start;
               increment = 0;
               break;
            case 1 : /* loop */
               dynamics.current = dynamics.end - increment;
               break;
            case 2 : /* swing */
               dynamics.current = dynamics.start + increment;
               dynamics.direction = 1;
               increment = -increment;
               break;
         }
      }

      
      dynamics.timestep -= 1;
//      if(increment == 0 || dynamics.current < 0 || dynamics.current >= dynamics.ntotalsteps || dynamics.timestep < 0) {
      if(dynamics.timestep < 0) {
         dynamics.current += increment;
         dynamics.timestep = dynamics.stepsize;
      }
      dynamicsGoto(dynamics.current);

      if(!dynamics.isrunning) {
         glutSetWindow(mainwin);
         if(bit.tgl_idle) GLUI_Master.set_glutIdleFunc(NULL);
         else glutIdleFunc(NULL);
         coord_cancel->enable();
      }

      if(playglui){
         if(dynamics.current < 0) sprintf(str, "structure # = %d", 1);
         else if (dynamics.current >= dynamics.ntotalsteps) 
                 sprintf(str, "structure # = %d", dynamics.ntotalsteps);
         else sprintf(str, "structure # = %d", dynamics.current + 1);
         play_no->set_text(str);
         playglui->sync_live();
      }

      glutSetWindow(mainwin);
      glutPostRedisplay();
}


void spinmotion(int x, int y)
{
   if(!firstmol) return;

   nmx = x;
   nmy = y;


   if(drawingtime > 600) rapidmove = 1;
   switch (current_button) {
      case GLUT_LEFT_BUTTON:
         t[0] = t[1] = t[2] = 0;
         r[0] = r[1] = r[2] = 0;
         r[3] = 1.0;
         spinrot[0] = spinrot[1] = spinrot[2] = 0;
         spinrot[3] = 1.0;
         spinmouse_rotation(&omx, &omy, &nmx, &nmy);
         spinning = 1;
      break;
      case GLUT_MIDDLE_BUTTON:
         if (modifiers & GLUT_ACTIVE_SHIFT) {
            t[0] = t[1] = 0;
            t[2] = (float)(omy-nmy)/(float)ysize;
            mouse_translation();
         }
         else if (modifiers & GLUT_ACTIVE_CTRL) {
            sfac *= pow(2.0, 2.0*(nmy-omy)/(float)ysize);
            scale_line(xsize, ysize);
         } 
         else {
            t[0] = (float)(nmx-omx)/(float)xsize*2.0;
            t[1] = (float)(omy-nmy)/(float)ysize*2.0;
            t[2] = 0;
            mouse_translation();
         }
      break;
   }
   omx = nmx;
   omy = nmy;
   glutPostRedisplay();
}


void keyboard(unsigned char key, int x, int y)
{
   switch(key) {
      case 27: /* ESC */
         Tgl_fullscreen();
      break;
      case 'a':
         Tgl_manip_all();
      break;
      case 'd':
         Tgl_depth();
      break;
      case 'w':
         Set_wire();
      break;
      case 'b':
         Set_ball_and_stick();
      break;
      case 's':
         Set_stick();
      break;
      case 'f':
         Set_spacefill();
      break;
      case 't':
         switch_mol();
      break;
      case 'i':
         Tgl_idle();
      break;
   }
   mainglui_top->sync_live();
   if(optionglui) optionglui->sync_live();
}


void spec_keyboard(int key, int x, int y)
{
   int key_modifiers;
   key_modifiers = glutGetModifiers();

   if(!firstmol) return;

   switch(key) {
      case GLUT_KEY_UP:
         if(key_modifiers & GLUT_ACTIVE_CTRL) {
            sfac *= 1.04;
            scale_line(xsize, ysize);
         }
         else if(key_modifiers & GLUT_ACTIVE_SHIFT) {
            t[0] = 0;
            t[1] = 5/(float)ysize*2.0;
            t[2] = 0;
            mouse_translation();
         }
         else if(key_modifiers & GLUT_ACTIVE_ALT) {
            t[0] = 0;
            t[1] = 0;
            t[2] = -10/(float)ysize*2.0;
            mouse_translation();
         }
         else {
            flip('x', angle);
         }
      break;
      case GLUT_KEY_DOWN:
         if(key_modifiers & GLUT_ACTIVE_CTRL) {
            sfac *= 0.96;
            scale_line(xsize, ysize);
         }
         else if(key_modifiers & GLUT_ACTIVE_SHIFT) {
            t[0] = 0;
            t[1] = -5/(float)ysize*2.0;
            t[2] = 0;
            mouse_translation();
         }
         else if(key_modifiers & GLUT_ACTIVE_ALT) {
            t[0] = 0;
            t[1] = 0;
            t[2] = 10/(float)ysize*2.0;
            mouse_translation();
         }
         else {
            flip('x', -angle);
         }
      break;
      case GLUT_KEY_LEFT:
         if(key_modifiers & GLUT_ACTIVE_CTRL) {
            flip('z', -angle);
         }
         else if(key_modifiers & GLUT_ACTIVE_SHIFT) {
            t[0] = -5/(float)xsize*2.0;
            t[1] = 0;
            t[2] = 0;
            mouse_translation();
         }
         else {
            flip('y', angle);
         }
      break;
      case GLUT_KEY_RIGHT:
         if(key_modifiers & GLUT_ACTIVE_CTRL) {
            flip('z', angle);
         }
         else if(key_modifiers & GLUT_ACTIVE_SHIFT) {
            t[0] = 5/(float)xsize*2.0;
            t[1] = 0;
            t[2] = 0;
            mouse_translation();
         }
         else {
            flip('y', -angle);
         }
      break;
   }
   glutPostRedisplay();
}


void mouse_rotation(long *omx, long *omy, long *nmx, long *nmy)
{
   float p1x, p1y, p2x, p2y;

   if(*nmx == *omx && *nmy == *omy) return;

   new_trackballcenter(actualmol);
   
   u_to_world(*omx, *omy, &p1x, &p1y);
   u_to_world(*nmx, *nmy, &p2x, &p2y);
   trackball(r, p1x, p1y, p2x, p2y);

   if(bit.manip_all){
   add_eulers(r, global_rot, global_rot);
   }
   else add_eulers(r, actualmol->rvec, actualmol->rvec);
}

void mouse_clip_rotation(long *omx, long *omy, long *nmx, long *nmy)
{
   float p1x, p1y, p2x, p2y, r[4], t[4];
   Surface *sp;

   if(*nmx == *omx && *nmy == *omy) return;

   new_trackballcenter(actualmol);
   
   u_to_world(*omx, *omy, &p1x, &p1y);
   u_to_world(*nmx, *nmy, &p2x, &p2y);
   trackball(r, p1x, p1y, p2x, p2y);

   t[0] = t[1] = t[2] = t[3] = 0;
   if(bit.manip_all) {
      if(bit.independent) {
         for(sp=actualmol->firstsurf; sp; sp=sp->next) {
           set_clipplane_matrix(sp, actualmol, r, t);
         }
      }
   }
   else {
      set_clipplane_matrix(actualsurf, actualmol, r, t);
   }
}

void spinmouse_rotation(long *omx, long *omy, long *nmx, long *nmy)
{
   float p1x, p1y, p2x, p2y;

   if(*nmx == *omx && *nmy == *omy) return;

   new_trackballcenter(actualmol);
   
   u_to_world(*omx, *omy, &p1x, &p1y);
   u_to_world(*nmx, *nmy, &p2x, &p2y);
   trackball(r, p1x, p1y, p2x, p2y);

   spinrot[0] = r[0];
   spinrot[1] = r[1];
   spinrot[2] = r[2];
   spinrot[3] = r[3];

   if(bit.manip_all){
   add_eulers(r, global_rot, global_rot);
   }
   else add_eulers(r, actualmol->rvec, actualmol->rvec);
}


void mouse_translation(void)
{
   if(bit.manip_all) {
      vadd(t,  global_trans,  global_trans);
   }
   else {
      float inv = 1.0/sfac;
      t[0] *= inv; t[1] *= inv; t[2] *= inv;
      vadd(t, actualmol->tvec, actualmol->tvec);
   }

}

void mouse_clip_translation(void)
{
   Surface *sp;

   r[0] = r[1] = r[2] = 0;
   r[3] = 1.0;
   if(bit.manip_all) {
      for(sp=actualmol->firstsurf; sp; sp=sp->next) {
         set_clipplane_matrix(sp, actualmol, r, t);
      }
   }
   else set_clipplane_matrix(actualsurf, actualmol, r, t);
}

/* transformation from window-coordinates to world-coordinates */
void u_to_world(long sx, long sy, float *wx, float *wy)
{
      *wx = (2.0 * sx) / (float) xsize - 1.0 - trackballcenter[0];
      *wy = (-2.0 * sy) / (float) ysize + 1.0 - trackballcenter[1];
}

void new_trackballcenter(Mol *mp)
{
   Matrix m;
   float v[4], w[4];

   if(!mp || (bit.manip_all && !clip_move_rot)){
      trackballcenter[0] = trackballcenter[1] = 0;
      return;
   }

   glGetFloatv(GL_PROJECTION_MATRIX, &m[0][0]);
   v[0] = mp->tvec[0]*sfac;   v[1] = mp->tvec[1]*sfac;
   v[2] = mp->tvec[2]*sfac;   v[3] = 1.0;
   vecxmat(v, m, w);

   if(w[3] != 1.0){ /* normalize */
      w[0] /= w[3];   w[1] /= w[3];
      w[2] /= w[3];   w[3] = 1.0;
   }

   trackballcenter[0] = w[0];
   trackballcenter[1] = w[1];
}

void flip(char axis, float angle)
{
   Mol *mp;

   if(!actualmol){
      logprint("no molecule loaded!");
      update_logs();
      return;
   }

   r[0] = r[1] = r[2] = 0;

   switch(axis){
      case 'x' : r[0] = fsin(0.5*angle); break;
      case 'y' : r[1] = fsin(0.5*angle); break;
      case 'z' : r[2] = fsin(0.5*angle); break;
   }

   r[3] = fcos(0.5*angle);

   if(bit.manip_all){
      for(mp=firstmol; mp; mp=mp->next) add_eulers(r, mp->rvec, mp->rvec);
   }
   else add_eulers(r, actualmol->rvec, actualmol->rvec);
}


void switch_mol(void)
{
   register AtoM *ap;

   if(!firstmol || !firstmol->next) {
      logprint("No molecule to switch with");
      update_logs();
      return;
   }

   if(!actualmol->next) actualmol = firstmol;
   else actualmol = actualmol->next;

   actualsurf = actualmol->firstsurf;
   update_interface_flags();

   for(ap=actualmol->firstatom; ap; ap=ap->next) ap->picked = 1;
   picking_mode = 1;
   drawit();
   if(actualsurf) highLight((unsigned long)actualsurf, PICK_MOL, actualmol, NULL);
   picking_mode = 0;
   for(ap=actualmol->firstatom; ap; ap=ap->next) ap->picked = 0;

}


void switch_surf(void)
{
   if(!actualmol || !actualmol->firstsurf || !actualmol->firstsurf->next) {
      logprint("No surface to switch with");
      update_logs();
      return;
   }

   if(!actualsurf->next) actualsurf = actualmol->firstsurf;
   else actualsurf = actualsurf->next;
   update_interface_flags();
   highLight((unsigned long)actualsurf, PICK_SURF, actualmol, NULL);
}


void get_translation(Matrix m, float *t)
{
   t[0] = m[3][0];
   t[1] = m[3][1];
   t[2] = m[3][2];
}



void get_rotation(Matrix m, float *e)
{
   e[0] = 0.5*sqrt(fabs(1.+m[0][0]-m[1][1]-m[2][2]));
   e[1] = 0.5*sqrt(fabs(1.-m[0][0]+m[1][1]-m[2][2]));
   e[2] = 0.5*sqrt(fabs(1.-m[0][0]-m[1][1]+m[2][2]));
   e[3] = 0.5*sqrt(fabs(1.+m[0][0]+m[1][1]+m[2][2]));

   /* set the first sign to positive, the others as follows */
   if(m[0][1]+m[1][0] < 0) e[1] = -e[1];
   if(m[0][2]+m[2][0] < 0) e[2] = -e[2];
   if(m[2][1]-m[1][2] < 0) e[3] = -e[3];
}

void reset_globals(void)
{
/*
   set the global modelviev matrix to scaling only,
   and modify each molecule's matrix so there's no net change
*/
   Mol *mp;
   Matrix m;
   float invScale;

/*
   first construct the global modelview matrix without scaling
   (by premultiplying with the inverse scaling matrix)
                (Sg-1 * Sg) * Tg * Rg
*/
   glPushMatrix();
   glLoadIdentity();
   invScale = 1.0/sfac;
   glScalef(invScale, invScale, invScale);
   glTranslatef(global_trans[0], global_trans[1], global_trans[2]);
   glScalef(sfac, sfac, sfac);
   build_rotmatrix(m, global_rot);
   glMultMatrixf(&m[0][0]);

/*
   then, for each molecule multiply with its matrix, except the
   centering translation.
                Tg * Rg * Tm * Rm
   extract the rotation and translation components rvec and tvec
   (which build the new modelview matrix) from the
   resulting matrix
*/
   for(mp=firstmol; mp; mp=mp->next){
      glPushMatrix();
      glTranslatef(mp->tvec[0], mp->tvec[1], mp->tvec[2]);
      build_rotmatrix(m, mp->rvec);
      glMultMatrixf(&m[0][0]);
      glGetFloatv(GL_MODELVIEW_MATRIX, &m[0][0]);
      glPopMatrix();
      get_translation(m, mp->tvec);
      get_rotation(m, mp->rvec);
   }
   glPopMatrix();

/*
   finally reset the global translation and rotation
*/
   global_trans[0] = global_trans[1] = global_trans[2] = 0;
   global_rot[0] = global_rot[1] = global_rot[2] = 0;
   global_rot[3] = 1;
}


void reset(void)
/* reset the rotation matrix to its original value */
{
   Mol *mp;

   if(!actualmol){
      logprint("No molecule loaded");
      update_logs();
      return;
   }

   glLoadIdentity();
   if(bit.manip_all){
      for(mp=firstmol; mp; mp=mp->next){
         mp->tvec[0] = mp->tvec[1] = mp->tvec[2] = 0;
         mp->rvec[0] = mp->rvec[1] = mp->rvec[2] = 0;
         center_of_mass();
         mp->rvec[3] = 1.0;
      }
      global_trans[0] = global_trans[1] = global_trans[2] = 0;
      global_rot[0] = global_rot[1] = global_rot[2] = 0;
      global_rot[3] = 1.0;
   }
   else {
      actualmol->tvec[0] = actualmol->tvec[1] = actualmol->tvec[2] = 0;
      actualmol->rvec[0] = actualmol->rvec[1] = actualmol->rvec[2] = 0;
      center_of_mass();
      actualmol->rvec[3] = 1.0;
   }
}


void center_of_mass(void)
{
   register AtoM *a;
   float vect[3];

   vect[0] = vect[1] = vect[2] = 0;
   actualmol->mass = 0;
   for(a = actualmol->firstatom; a; a = a->next){
      vect[0] -= a->coord[0]*element[a->ord].mass;
      vect[1] -= a->coord[1]*element[a->ord].mass;
      vect[2] -= a->coord[2]*element[a->ord].mass;
      actualmol->mass += element[a->ord].mass;
   }
   if(!actualmol->mass){
      actualmol->centervec[0] = actualmol->centervec[1] =
      actualmol->centervec[2] = 0;
      return;
   }
   actualmol->centervec[0] = vect[0]/actualmol->mass;
   actualmol->centervec[1] = vect[1]/actualmol->mass;
   actualmol->centervec[2] = vect[2]/actualmol->mass;
}


void update_interface_flags(void)
{
   char str[100];
   got_dots = 0;
   if(actualmol) {
      bit.wire_model = actualmol->wire_model;
      bit.sticks = actualmol->sticks;
      bit.ball_and_stick = actualmol->ball_and_stick;
      bit.spacefill = actualmol->spacefill;
      bit.transparent = actualmol->transparent;
      bit.suppress_H = actualmol->suppress_H;
      bit.residues = actualmol->residues;
      bit.bond_col = actualmol->bond_col;
      bit.h_bond = actualmol->h_bond;
      bit.multiple_bonds = actualmol->multiple_bonds;
      bit.show_freq_arrow = actualmol->show_freq_arrow;
      bit.show_dipole = actualmol->show_dipole;
      bit.labels = actualmol->labels;
      bit.atm_char = actualmol->atm_char;
      bit.atm_spin = actualmol->atm_spin;
      bit.cubeplanes = actualmol->cubeplanes;
      sc_freq_ar = actualmol->sc_freq_ar / 5;
      sc_dipole_ar = actualmol->sc_dipole_ar / 2;
      if(actualsurf) {
         bit.surf_clip = actualsurf->surf_clip;
         bit.surf_transp = actualsurf->surf_transp;
         bit.dot_surf = actualsurf->dot_surf;
         bit.chickenwire = actualsurf->chickenwire;
         bit.flatshade = actualsurf->flatshade;
      }

      if(actualmol->wire_model) molprop_live_var = 0;
      else if (actualmol->sticks) molprop_live_var = 1;
      else if (actualmol->ball_and_stick) molprop_live_var = 2;
      else if (actualmol->spacefill) molprop_live_var =3;

      if(actualmol->cube_value) {
         cubemin = actualmol->cubemin;
         cubemax = actualmol->cubemax;
         cutoff = actualmol->cutoff;
      }
      else {
         cubemin = 0;
         cubemax = 0;
         cutoff  = 0;
      }
      actualvib = actualmol->freq_arrow;
   }

   mainglui_top->sync_live();
   if(surfaceglui) surfaceglui->sync_live();
   if(surf_optglui) surf_optglui->sync_live();
   if(optionglui) optionglui->sync_live();
   if(bondglui) bondglui->sync_live();
   if(labelglui) labelglui->sync_live();
   if(planeglui) {
      planeglui->close();
      planeglui = NULL;
   }
   if(dipoleglui) {
      if(actualmol && actualmol->dipole) 
         sprintf(str, "dipole moment: %8.4f Debye", actualmol->dipole->absolute);
      else sprintf(str, "no dipole moment info");
      dipoletext->set_text(str);
      dipoleglui->sync_live();
   }
   if(freqglui) {
      if(actualvib) sprintf(str, "%5.0f/cm %s", actualvib->frequency, actualvib->type); 
      else sprintf(str, "no frequency selected");
      freqtext->set_text(str);
      freqglui->sync_live();
   }
   update_logs();
}
