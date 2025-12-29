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


/* drawing routines */
#include "constant.h"
#include "main.h"
#include "molekel.h"
#include "glutwin.h"
#include "general.h"
#include "drawing.h"
#include "drawmol.h"
#include "drawprop.h"
#include "trackball.h"
#include "objects.h"
#include "drawprop.h"
#include "utils.h"
#include "texture.h"
#include "material.h"

/*----------------------------------------*/

int drawingtime;
float sfac = 0.1;
static Matrix m;
float t[3], r[4] = {0, 0, 0, 1}, s;
float spinrot[4] = {0, 0, 0, 1};
float global_trans[3]  = {0, 0, 0}, global_rot[4] = {0, 0, 0, 1};
float global_center[3] = {0, 0, 0}, global_scale  = 1.0;
float global_spin[4]   = {0, 0, 0, 1};
short rapidmove = 0;
float  scale_factor;

unsigned long backcolor1 = 0x00000000, backcolor2 = 0x00000000;
float distcolor[] = {1.0, 1.0, 0.0, 1.0};
float anglecolor[] = {0.0, 1.0, 1.0, 1.0};
float dihedcolor[] = {1.0, 0.0, 1.0, 1.0};
float labelcolor[] = {1.0, 1.0, 0.0, 1.0};
int halftone_pattern;

/*----------------------------------------*/
/* test */
GLuint testlist, testlist2;
/*----------------------------------------*/

/* draw the whole scene */
void drawit(void)
{
   int starttime, endtime;
   TRcontext *tr = NULL;

   starttime = glutGet(GLUT_ELAPSED_TIME);

   glClearColor(0.0, 0.0, 0.0, 0.0);
   glClear(GL_DEPTH_BUFFER_BIT);
   draw_back();

   if(bit.depthcue) glEnable(GL_FOG);
   else             glDisable(GL_FOG);

/* to be fixed!!!

   if(bit.rc_stereo){
      RGBwritemask(0xff, 0, 0);
      drawit2();
      RGBwritemask(0, 0xff, 0xff);
      flip('y', 5*M_PI/180);
      zclear();
      drawit2();
      RGBwritemask(0xff, 0xff, 0xff);
      flip('y', -5*M_PI/180);
   }

   else if(bit.acbuf  && (!rapidmove)){
      int i, j, ndiv;
      float delta;

** the square of ndiv is the oversampling rate **
      ndiv = 3;
      delta = 2.0/(float)ysize/(float)ndiv;


      readsource(SRC_AUTO);
      acbuf(AC_CLEAR, 0);

      for(i=0; i<ndiv; i++) {
         x_offset = i*delta;
         for(j=0; j<ndiv; j++) {
            y_offset = j*delta;
            init_persp();
            zclear(); draw_back();
            drawit2();
            acbuf(AC_ACCUMULATE, 1.0);
         }
      }

      acbuf(AC_RETURN, 1.0/(float)ndiv/(float)ndiv);
   }
   else {
       drawit2();
   }
*/

   drawit2();

/* to be fixed!!!
   if(getgdesc(GD_FOGVERTEX)) fogvertex(FG_OFF, fog_params);

   draw_text();

   if(bit.activ_logo) {
       draw_activ_mol();
       draw_logo();
   }
*/

   glutSwapBuffers();


   endtime = glutGet(GLUT_ELAPSED_TIME);
   drawingtime = endtime - starttime;
//printf("drawing time: %d milliseconds\n", drawingtime);
}


void drawit2(void)
{
   extern Int_dist *first_intermolecular_distance;
   Mol *mp;

   halftone_pattern = 0;
   glCallList(line_style_base);
   glPushMatrix();
   global_move();


   for(mp=firstmol; mp; mp=mp->next){
      glPushMatrix();
      individual_move(mp);
      draw_structure(mp);
      draw_properties(mp);
      if(mp == actualmol && bit.coord_sys) {
         draw_coord_system(mp);
      }
/* draw color map scale */
      if(maplegend) {
         draw_mapscale(mp);
      }

/* test: case of
      draw3();
      draw5();
      draw6();
      glCallList(testlist);
      glCallList(testlist2);
*/

      glPopMatrix();

      if((mp->labels || mp->atm_char || mp->atm_spin) && (!rapidmove)) draw_labels(m, mp);
   }

   if(first_intermolecular_distance) drawIntermolecularDistances();

   glPopMatrix();
/* to be fixed!!!
   subpixel(0);
*/
}


void set_color(unsigned long c)
{
   glColor4ub(c&0xff, c>>8&0xff, c>>16&0xff, c>>24&0xff);
}

void set_clear_color(unsigned long c)
{
   glClearColor((c&0xff) / 255.0, (c>>8&0xff) / 255.0, (c>>16&0xff) / 255.0, (c>>24&0xff) / 255.0);
}

void global_move(void)
{
   glTranslatef(global_trans[0], global_trans[1], global_trans[2]);
   glScalef(sfac, sfac, sfac);
   build_rotmatrix(m, global_rot);
   glMultMatrixf(&m[0][0]);
}


void individual_move(Mol *mp)
{
      glTranslatef(mp->tvec[0], mp->tvec[1], mp->tvec[2]);
      build_rotmatrix(m, mp->rvec);
      glMultMatrixf(&m[0][0]);

/* translation in object's coordinates */
      glTranslatef(mp->centervec[0], mp->centervec[1], mp->centervec[2]);
}


void draw_back(void)
{
   float ratio;
   unsigned long col1, col2;
   static float vect[][3] = {{ -1.5, -1.5, -1.49}, {  1.5, -1.5, -1.49},
                             {  1.5,  1.5, -1.49}, { -1.5,  1.5, -1.49} };
   static float t[][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};

   glDisable(GL_FOG);

   col1 = backcolor1; col2 = backcolor2;

/* to be fixed!!!
 * divide each component by 2 * 
   if(shadow_offset){
      col1 = (backcolor1>>1 & 0x7f7f7f) | 0xff000000;
      col2 = (backcolor2>>1 & 0x7f7f7f) | 0xff000000;
   }
   else {
      col1 = backcolor1; col2 = backcolor2;
   }
*/

   if(bit.texture && back_texture){
      enable_texture(back_texture, back_texenv, back_textype);
      glDisable(GL_TEXTURE_GEN_S);
      glDisable(GL_TEXTURE_GEN_T);
   }

   if((col1 == col2) && !(bit.texture && back_texture)) {
      set_clear_color(col1);
      glClear(GL_COLOR_BUFFER_BIT);
      return;
   }

   glClearColor(clearcolor[0], clearcolor[1], clearcolor[2], clearcolor[3]);
   glClear(GL_COLOR_BUFFER_BIT);

   ratio = 1.5*(float)xsize/ysize;
   vect[0][0] = vect[3][0] = -ratio;
   vect[1][0] = vect[2][0] =  ratio;

   glDisable(GL_LIGHTING);
   glDepthMask(GL_FALSE);
   glBegin(GL_TRIANGLE_STRIP);
      set_color(col1);
      glTexCoord2fv(vect[0]); glVertex3fv(vect[0]);
      glTexCoord2fv(vect[1]); glVertex3fv(vect[1]);
      set_color(col2);
      glTexCoord2fv(vect[3]); glVertex3fv(vect[3]);
      glTexCoord2fv(vect[2]); glVertex3fv(vect[2]);
   glEnd();
   glDepthMask(GL_TRUE);
   glEnable(GL_LIGHTING);

   if(bit.texture && back_texture) disable_texture();

   if(bit.depthcue) glEnable(GL_FOG);

}

void get_center_traj(int n, float *c, Vector *vp)
{
   int i;
   float x, y, z;

   if(!n) return;
   x = y = z = 0;

   for(i=0; i<n; i++, vp++){
      x += vp->x;
      y += vp->y;
      z += vp->z;
   }

   c[0] = x/n;   c[1] = y/n;   c[2] = z/n;

}

void get_center_mol(int n, float *c, Mol *mp)
{
   float x, y, z;
   AtoM *ap;

   if(!n) return;
   x = y = z = 0;

   for(ap=actualmol->firstatom; ap; ap=ap->next){
      x += ap->coord[0];
      y += ap->coord[1];
      z += ap->coord[2];
   }

   c[0] = x/n;   c[1] = y/n;   c[2] = z/n;

}

#define TOLERANCE   0.001
#define MAXIT     400
int superImposeTraj(long step, Matrix rt)
{
   int i, j, k, ict, ix, iy, iz, iflag;
   float xt, yt, zt, xx;
   float aa[3][3], da[3], db[3], cgm[3], cgf[3];
   float gam, gam2, sig, sig2, bb, cc, sg;
   AtoM *ap;
   Vector *vw;

   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) aa[i][j] = rt[i][j] = 0;
      rt[i][i] = 1.;
   }

   get_center_traj(actualmol->natoms, cgm, dynamics.trajectory[step]);
   get_center_mol(actualmol->natoms, cgf, actualmol);


   for(ap=actualmol->firstatom, vw=dynamics.trajectory[step]; ap; ap=ap->next, vw++){
      da[0] = vw->x - cgm[0];
      da[1] = vw->y - cgm[1];
      da[2] = vw->z - cgm[2];
      db[0] = ap->coord[0] - cgf[0];
      db[1] = ap->coord[1] - cgf[1];
      db[2] = ap->coord[2] - cgf[2];
      for(j=0; j<3; j++){
         xx = da[j];
         aa[j][0] +=  xx * db[0];
         aa[j][1] +=  xx * db[1];
         aa[j][2] +=  xx * db[2];
      }
   }

   ict = 0;
   do {
      iflag = 0;
      for(ix=0; ix<3; ix++){
         if(++ict > MAXIT) return 0;
         iy = (ix + 1)%3;
         iz = (ix + 2)%3;
         sig = aa[iz][iy] - aa[iy][iz];
         gam = aa[iy][iy] + aa[iz][iz];
         sig2 = sqrt(sig*sig);
         gam2 = sqrt(gam*gam);
         sg = sqrt(sig*sig + gam*gam);
         if((sg > 0.001) && sig2 > TOLERANCE*gam2){
            sg = 1./sg;
            for(k=0; k<3; k++){
               bb = gam*aa[iy][k] + sig*aa[iz][k];
               cc = gam*aa[iz][k] - sig*aa[iy][k];
               aa[iy][k] = bb * sg;
               aa[iz][k] = cc * sg;
               bb = gam*rt[iy][k] + sig*rt[iz][k];
               cc = gam*rt[iz][k] - sig*rt[iy][k];
               rt[iy][k] = bb * sg;
               rt[iz][k] = cc * sg;
            }
            iflag = 1;
         } /* endif */
      }
   } while(iflag);

   transpose(rt);
   for(i=0; i<3; i++){
      xt    = -rt[0][i] * cgm[0];
      yt    = -rt[1][i] * cgm[1];
      zt    = -rt[2][i] * cgm[2];
      rt[3][i] =  xt + yt + zt + cgf[i];
      rt[i][3]= 0;
   }
   rt[3][3] = 1;
   return 1;
}

void dynamicsGoto(long step)
{
   register AtoM *ap;
   register int i;
   Vector *vw;
   Matrix rt;
   float u[4], v[4];

   if(step < 0 || step >= dynamics.ntotalsteps) return;

   if(bit.superimpose) {
      if(!superImposeTraj(step, rt)) return;
      for(ap=actualmol->firstatom, vw=dynamics.trajectory[step]; ap; ap=ap->next, vw++){
            u[0] = vw->x;
            u[1] = vw->y;
            u[2] = vw->z;
            u[3] = 1.0;
            vecxmat(u, rt, v);
            ap->coord[0] = v[0];
            ap->coord[1] = v[1];
            ap->coord[2] = v[2];
      }
   }
   else {
      if(dynamics.freeat && step){
         for(i=0, vw=dynamics.trajectory[step]; i<dynamics.nfreat; i++, vw++){
            dynamics.freeat[i]->coord[0] = vw->x;
            dynamics.freeat[i]->coord[1] = vw->y;
            dynamics.freeat[i]->coord[2] = vw->z;
         }
      }
      else {
         for(ap=actualmol->firstatom, vw=dynamics.trajectory[step]; ap; ap=ap->next, vw++){
            ap->coord[0] = vw->x;
            ap->coord[1] = vw->y;
            ap->coord[2] = vw->z;
         }
      }
   }

   dynamics.current = step;
   
    if(!bit.keepbonds) {
	AtoM *ap;
	if(actualmol->firstbond) free_bonds(actualmol->firstbond);
	actualmol->nbonds = 0;
	actualmol->firstbond = NULL;
	for(ap=actualmol->firstatom; ap; ap=ap->next) {
	    ap->nbonds = 0;
	}
	createBonds(0);
        create_Hbonds();
    }

}


