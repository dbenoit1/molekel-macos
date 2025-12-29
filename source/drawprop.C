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
    property drawing routines
*/

#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "general.h"
#include "glutwin.h"
#include "drawprop.h"
#include "drawing.h"
#include "macu.h"
#include "utils.h"
#include "material.h"
#include "geometry.h"
#include "texture.h"
#include "maininterf.h"
#include "objects.h"
#include "drawmol.h"
#include "trackball.h"
#include "render.h"

static unsigned long property_col[] = 
   {0xff0000ff, 0xff000fff, 0xff001fff, 0xff002fff, 0xff003fff, 0xff004fff,
    0xff005fff, 0xff006fff, 0xff007fff, 0xff008fff, 0xff009fff, 0xff00afff,
    0xff00bfff, 0xff00cfff, 0xff00dfff, 0xff00efff, 0xff00ffff, 0xff00ffef,
    0xff00ffdf, 0xff00ffcf, 0xff00ffbf, 0xff00ffaf, 0xff00ff9f, 0xff00ff8f,
    0xff00ff7f, 0xff00ff6f, 0xff00ff5f, 0xff00ff4f, 0xff00ff3f, 0xff00ff2f,
    0xff00ff1f, 0xff00ff0f, 0xff00ff00, 0xff0fff00, 0xff1fff00, 0xff2fff00,
    0xff3fff00, 0xff4fff00, 0xff5fff00, 0xff6fff00, 0xff7fff00, 0xff8fff00,
    0xff9fff00, 0xffafff00, 0xffbfff00, 0xffcfff00, 0xffdfff00, 0xffefff00,
    0xffffff00, 0xffffef00, 0xffffdf00, 0xffffcf00, 0xffffbf00, 0xffffaf00,
    0xffff9f00, 0xffff8f00, 0xffff7f00, 0xffff6f00, 0xffff5f00, 0xffff4f00,
    0xffff3f00, 0xffff2f00, 0xffff1f00, 0xffff0f00, 0xffff0000 };
static unsigned long shadow_col[] = 
   {0xff00007f, 0xff00077f, 0xff000f7f, 0xff00177f, 0xff001f7f, 0xff00277f,
    0xff002f7f, 0xff00377f, 0xff003f7f, 0xff00477f, 0xff004f7f, 0xff00577f,
    0xff005f7f, 0xff00677f, 0xff006f7f, 0xff00777f, 0xff007f7f, 0xff007f77,
    0xff007f6f, 0xff007f57, 0xff007f5f, 0xff007f57, 0xff007f4f, 0xff007f47,
    0xff007f4f, 0xff007f37, 0xff007f2f, 0xff007f27, 0xff007f1f, 0xff007f17,
    0xff007f0f, 0xff007f07, 0xff007f00, 0xff077f00, 0xff0f7f00, 0xff177f00,
    0xff1f7f00, 0xff277f00, 0xff2f7f00, 0xff377f00, 0xff3f7f00, 0xff477f00,
    0xff4f7f00, 0xff577f00, 0xff5f7f00, 0xff677f00, 0xff6f7f00, 0xff777f00,
    0xff7f7f00, 0xff7f7700, 0xff7f6f00, 0xff7f6700, 0xff7f5f00, 0xff7f5700,
    0xff7f4f00, 0xff7f4700, 0xff7f3f00, 0xff7f2700, 0xff7f2f00, 0xff7f2700,
    0xff7f1f00, 0xff7f1700, 0xff7f0f00, 0xff7f0700, 0xff7f0000 };
unsigned long vector_color[44];


static Mol *molecule;
void draw_distances(Mol *mp);
void draw_angles(Mol *mp);
void draw_torsions(Mol *mp);
void get_norm(float *a, float *b, float *c, Vector *d);
int make_trinorm(Surface *sp);




void draw_properties(Mol *mp)
{
   molecule = mp;

/* to be fixed
   if(bit.depthcue && !getgdesc(GD_FOGVERTEX)) depthcue(1);

   if(position) draw_positions(mp);
*/

   if(bit.smoothline){
      glEnable(GL_POINT_SMOOTH);
      glEnable(GL_LINE_SMOOTH);
//      glBlendFunc(GL_SRC_ALPHA,  GL_ZERO);
      glBlendFunc(GL_SRC_ALPHA,  GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_BLEND);
   }

   if(mp->firstdist) draw_distances(mp);

   if(mp->firstang)  draw_angles(mp);

   if(mp->firsttor)  draw_torsions(mp);

   if(bit.smoothline){
      glDisable(GL_POINT_SMOOTH);
      glDisable(GL_LINE_SMOOTH);
      glBlendFunc(GL_ONE,  GL_ZERO);
      glDisable(GL_BLEND);
   }

   if(mp->firstsurf) draw_surfaces(mp);

   if(mp->cubeplanes && mp->cube_value){
      draw_cubeplanes(mp);
   }
   
   if(mp->show_freq_arrow && mp->freq_arrow && mp->vibration && mp->n_frequencies){
      if(bit.solid_arrow && !rapidmove) draw_solid_freq_arrow(mp);
      else draw_freq_arrow(mp);
   }

   if(mp->show_dipole && mp->dipole) {
      draw_dipole(mp);
   }

   if(mp==project_mol) draw_2Dgrid(mp);

}

void draw_string(float x, float y, float z, char *string)
/* Print string in a GLUT window at position x, y, z */
{
   int len, i;
   void *font;
   float mat[4][4];

   switch(label_size) {
      case 0:
         font = GLUT_BITMAP_HELVETICA_12;
      break;
      case 1:
         font = GLUT_BITMAP_9_BY_15;
      break;
      case 2:
         font = GLUT_BITMAP_HELVETICA_18;
      break;
      case 3:
         font = GLUT_BITMAP_TIMES_ROMAN_24;
      break;
      case 4:
         font = GLUT_STROKE_ROMAN;
      break;
   }
//   void *font = GLUT_BITMAP_HELVETICA_18;
//   void *font = GLUT_BITMAP_9_BY_15;
//   void *font = GLUT_BITMAP_TIMES_ROMAN_24;

   if(font == GLUT_STROKE_ROMAN) {
      glCallList(line_style_base + 3);
      glGetFloatv(GL_MODELVIEW_MATRIX, &mat[0][0]);
//for (i=0; i < 4; i++) {
//printf("%f %f %f %f\n", mat[0][i], mat[1][i], mat[2][i], mat[3][i]);
//}
      transpose(mat);
      for (i=0; i < 3; i++) {
         mat[i][3] = 0;
         mat[3][i] = 0;
      }
//for (i=0; i < 4; i++) {
//printf("%f %f %f %f\n", mat[0][i], mat[1][i], mat[2][i], mat[3][i]);
//}
      glPushMatrix();
      glTranslatef(x, y, z);
      glMultMatrixf(&mat[0][0]);
      glScalef(0.01, 0.01, 0.01);
      len = (int) strlen(string);
      for (i = 0; i < len; i++) {
         glutStrokeCharacter(font, string[i]);
      }
      glPopMatrix();
      glCallList(line_style_base + 1);
   }
   else {
      if(bit.render) trRasterPos3f(tr, x, y, z);
      else glRasterPos3f(x, y, z);
      len = (int) strlen(string);
      for (i = 0; i < len; i++) {
         glutBitmapCharacter(font, string[i]);
      }
   }
}


void draw_surfaces(Mol *mp)
{
   register Surface *sp;


/* to be fixed
   if(bit.smoothline) {
      if(mp->chickenwire || mp->wire_model) blendfunction(BF_SA, BF_ZERO);
   }
   else blendfunction(BF_ONE, BF_ZERO);


   for(sp=mp->firstsurf; sp; sp=sp->next){
      if(rapidmove)            draw_dots(sp);
      else if(mp->chickenwire) draw_chickenwire(sp);
      else if(mp->wire_model)   draw_dots(sp);
      else {
         if(mp->transparent) setpattern(2);
         if(mp->flatshade)   draw_flatshaded_surface(sp);
         else                draw_solid_surface(sp);
         setpattern(0);
      }
   }
*/

   for(sp=mp->firstsurf; sp; sp=sp->next){
      if(sp->surf_clip) clipplane_on(sp);
      if(sp->dot_surf || rapidmove) draw_dots(sp);
      else if(sp->chickenwire) draw_chickenwire(sp);
      else if(sp->surf_transp) {
         glEnable(GL_POLYGON_STIPPLE);
         glPolygonStipple(halftone[halftone_pattern]);
         if(halftone_pattern < 3) halftone_pattern++; else halftone_pattern = 0;
         if(sp->flatshade) {
            glShadeModel(GL_FLAT);
            draw_flatshaded_surface(sp);
            glShadeModel(GL_SMOOTH);
         }
         draw_solid_surface(sp);
         glDisable(GL_POLYGON_STIPPLE);
      }
      else {
         if(sp->flatshade) {
            glShadeModel(GL_FLAT);
            draw_flatshaded_surface(sp);
            glShadeModel(GL_SMOOTH);
         }
         draw_solid_surface(sp);
      }
      if(sp->surf_clip) clipplane_off();
   }

}



void draw_chickenwire(Surface *sp)
{
   register int i;
   Surfdot  *dot;
   Triangle *j;
   int col;
   float factor, val;
   float diffuse_color[4];

   if(j = sp->tri){
      if(bit.smoothline){
         glEnable(GL_POINT_SMOOTH);
         glEnable(GL_LINE_SMOOTH);
         glBlendFunc(GL_SRC_ALPHA,  GL_ZERO);
         glEnable(GL_BLEND);
      }
      glCallList(sp->matindex);
      glDisable(GL_LIGHTING);
      dot = sp->dot;
      if(sp->val){
         factor = 64.99 / (sp->vmax - sp->vmin);
         for(i=0; i<sp->ntri; i++, j++){
            val = (sp->val[j->p1] + sp->val[j->p2]
                   + sp->val[j->p3]) * .333333333333;
            if(val > sp->vmax) col = (int)64.99;
            else if(val< sp->vmin) col = (int)0.00;
            else col = (int)(0.01 + (val - sp->vmin) * factor);
            set_color(property_col[col]);

            glBegin(GL_LINE_LOOP);
            glVertex3fv(dot[j->p1].v);
            glVertex3fv(dot[j->p2].v);
            glVertex3fv(dot[j->p3].v);
            glEnd();

         }
      }
      else {
         glGetMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_color);
         glColor4fv(diffuse_color);
         for(i=0; i<sp->ntri; i++, j++){
            glBegin(GL_LINE_LOOP);
            glVertex3fv(dot[j->p1].v);
            glVertex3fv(dot[j->p2].v);
            glVertex3fv(dot[j->p3].v);
            glEnd();
         }
      }
      if(bit.smoothline){
         glDisable(GL_POINT_SMOOTH);
         glDisable(GL_LINE_SMOOTH);
         glBlendFunc(GL_ONE,  GL_ZERO);
         glDisable(GL_BLEND);
      }
      glEnable(GL_LIGHTING);
   }
   else draw_dots(sp);

}





void draw_dots(Surface *sp)
{
   register int i;
   register Surfdot *j;
   float factor, size;
   float diffuse_color[4];
   int col;

   glCallList(sp->matindex);
   glDisable(GL_LIGHTING);
   if(bit.smoothline){
      glEnable(GL_POINT_SMOOTH);
      glEnable(GL_LINE_SMOOTH);
      glBlendFunc(GL_SRC_ALPHA,  GL_ZERO);
      glEnable(GL_BLEND);
   }
   if(sp->val){
      factor = 64.99 / (sp->vmax - sp->vmin);
      glGetFloatv(GL_POINT_SIZE, &size);
      glPointSize(2.0);
      glBegin(GL_POINTS);
      for(i=0, j=sp->dot; i<sp->npts; i++, j++){
         if(sp->val[i] > sp->vmax) col = (int)64.99;
         else if(sp->val[i] < sp->vmin) col = (int)0.00;
         else col = (int)(0.01 + (sp->val[i] - sp->vmin) * factor);
         set_color(property_col[col]);
         glVertex3fv(j->v);
      }
      glEnd();
      glPointSize(size);
   }
   else {
      glGetMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_color);
      glColor4fv(diffuse_color);
      glGetFloatv(GL_POINT_SIZE, &size);
      glPointSize(2.0);
      glBegin(GL_POINTS);
      for(i=0, j = sp->dot; i<sp->npts; i++, j++) glVertex3fv(j->v); 
      glEnd();
      glPointSize(size);
   }
   if(bit.smoothline){
      glDisable(GL_POINT_SMOOTH);
      glDisable(GL_LINE_SMOOTH);
      glBlendFunc(GL_ONE,  GL_ZERO);
      glDisable(GL_BLEND);
   }
   glEnable(GL_LIGHTING);

}




void draw_solid_surface(Surface *sp)
{
   register int i;
   Surfdot *dot;
   Triangle *j;
   static short col;
   float factor, factor1, t[2];
   unsigned long *cp;

   if(j = sp->tri){

/* to be fixed
      if(sp->alpha && getgdesc(GD_BLEND))
         blendfunction(BF_SA, BF_MSA);
*/
      glEnable(GL_ALPHA_TEST);

      dot = sp->dot;
/* to be fixed
      lmbind(MATERIAL, sp->matindex + shadow_offset);
*/
      glCallList(sp->matindex);

      if(sp->val){
         if(bit.texture && sp->texture){
            factor = 1.0 / (sp->vmax - sp->vmin);
            factor1 = 1.0 / (sp->vmax1 - sp->vmin1);
            enable_texture(sp->texture, sp->texenv, sp->textype);
            glDisable(GL_TEXTURE_GEN_S);
            glDisable(GL_TEXTURE_GEN_T);
            t[1] = 0;
/* alpha blending */
               if(alphablending) {
                  glEnable(GL_BLEND);
                  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
               }
               for(i=0; i<sp->ntri; i++, j++){
                  glBegin(GL_TRIANGLES);
                  t[0] = (sp->val[j->p1] - sp->vmin) * factor;
                  if(sp->val1) t[1] = (sp->val1[j->p1] - sp->vmin1) * factor1;
                  glTexCoord2fv(t);   glNormal3fv(dot[j->p1].n);   glVertex3fv(dot[j->p1].v);

                  t[0] = (sp->val[j->p2] - sp->vmin) * factor;
                  if(sp->val1) t[1] = (sp->val1[j->p2] - sp->vmin1) * factor1;
                  glTexCoord2fv(t);   glNormal3fv(dot[j->p2].n);   glVertex3fv(dot[j->p2].v);

                  t[0] = (sp->val[j->p3] - sp->vmin) * factor;
                  if(sp->val1) t[1] = (sp->val1[j->p3] - sp->vmin1) * factor1;
                  glTexCoord2fv(t);   glNormal3fv(dot[j->p3].n);   glVertex3fv(dot[j->p3].v);
                  glEnd();
               }
               if(alphablending) glDisable(GL_BLEND);
            disable_texture();
         }
         else {
            glColorMaterial(GL_FRONT, GL_DIFFUSE);
            glEnable(GL_COLOR_MATERIAL);
            factor = 64.99 / (sp->vmax - sp->vmin);
/* to be fixed
            if(shadow_offset) cp = shadow_col;
            else              cp = property_col;
*/
            cp = property_col;
            if(sp->alpha){
               unsigned long test;
               test = 0x00ffffff | (sp->alpha <<24);
               for(i=0; i<65; i++) cp[i] &= test;
            }
/* alpha blending */
               if(alphablending) {
                  glEnable(GL_BLEND);
                  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
               }
               glBegin(GL_TRIANGLES);
               for(i=0; i<sp->ntri; i++, j++){
                  if(sp->val[j->p1] > sp->vmax) col = (int)64.99;
                  else if(sp->val[j->p1] < sp->vmin) col = (int)0.0;
                  else col = (int)((sp->val[j->p1] - sp->vmin) * factor);
                  set_color(cp[col]);
                  glNormal3fv(dot[j->p1].n);   glVertex3fv(dot[j->p1].v);
                  if(sp->val[j->p2] > sp->vmax) col = (int)64.99;
                  else if(sp->val[j->p2] < sp->vmin) col = (int)0.0;
                  else col = (int)((sp->val[j->p2] - sp->vmin) * factor);
                  set_color(cp[col]);
                  glNormal3fv(dot[j->p2].n);   glVertex3fv(dot[j->p2].v);
                  if(sp->val[j->p3] > sp->vmax) col = (int)64.99;
                  else if(sp->val[j->p3] < sp->vmin) col = (int)0.0;
                  else col = (int)((sp->val[j->p3] - sp->vmin) * factor);
                  set_color(cp[col]);
                  glNormal3fv(dot[j->p3].n);   glVertex3fv(dot[j->p3].v);
               }
               glEnd();
               if(alphablending) glDisable(GL_BLEND);
            glDisable(GL_COLOR_MATERIAL);
            if(sp->alpha) for(i=0; i<65; i++) cp[i] |= 0xff000000;
         }
      }
      else {
         if(bit.texture && sp->texture && sp->textype == TEX_PHONG){
            /* simulated phong */
            glDisable(GL_LIGHTING);
            glColor4f(1.0, 1.0, 1.0, 1.0);
            enable_texture(sp->texture, sp->texenv, sp->textype);
            if(alphablending){
               glEnable(GL_BLEND);
               glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            }
            get_pmat();
            //cpack(0xffffffff);
            glBegin(GL_TRIANGLES);
            for(i=0; i<sp->ntri; i++, j++){
               myt3f(dot[j->p1].n);   glVertex3fv(dot[j->p1].v);
               myt3f(dot[j->p2].n);   glVertex3fv(dot[j->p2].v);
               myt3f(dot[j->p3].n);   glVertex3fv(dot[j->p3].v);
            }
            glEnd();
            if(alphablending) glDisable(GL_BLEND);
            disable_texture();
            glEnable(GL_LIGHTING);
         }
         else {
            if(bit.texture && sp->texture) enable_texture(sp->texture, sp->texenv, sp->textype);
/* alpha blending */
            if(alphablending){
               glEnable(GL_BLEND);
               glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            }
            glBegin(GL_TRIANGLES);
            for(i=0; i<sp->ntri; i++, j++){
               glNormal3fv(dot[j->p1].n);   glVertex3fv(dot[j->p1].v);
               glNormal3fv(dot[j->p2].n);   glVertex3fv(dot[j->p2].v);
               glNormal3fv(dot[j->p3].n);   glVertex3fv(dot[j->p3].v);
            }
            glEnd();
            if(alphablending) glDisable(GL_BLEND);
            if(bit.texture && sp->texture) disable_texture();
         }
      }

/* to be fixed
      blendfunction(BF_ONE, BF_ZERO);
*/
      glDisable(GL_ALPHA_TEST);

   }

   else draw_dots(sp);
}






void draw_flatshaded_surface(Surface *sp)
{
   register int i;
   Surfdot *dot;
   Triangle *j;
   Vector *vp;
   static short col;
   float factor, t[2];
   unsigned long *cp;

   if(j = sp->tri){

      if(!sp->trinorm) {
         if(!make_trinorm(sp)) return;
      }

/* to be fixed
      if(sp->alpha && getgdesc(GD_BLEND))
         blendfunction(BF_SA, BF_MSA);
*/

      dot = sp->dot;
/* to be fixed
      lmbind(MATERIAL, sp->matindex + shadow_offset);
*/
      glCallList(sp->matindex);

      if(sp->val){
         if(bit.texture && sp->texture){
            factor = 1.0 / (sp->vmax - sp->vmin);
            enable_texture(sp->texture, sp->texenv, sp->textype);
            glDisable(GL_TEXTURE_GEN_S);
            glDisable(GL_TEXTURE_GEN_T);
            t[1] = 0;
/* alpha blending */
               if(alphablending) {
                  glEnable(GL_BLEND);
                  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
               }
               for(i=0, vp=sp->trinorm; i<sp->ntri; i++, j++, vp++){
                  glBegin(GL_TRIANGLES);
                  t[0] = (sp->val[j->p1] - sp->vmin) * factor;
                  glTexCoord2fv(t);  glNormal3fv((float *)vp);   glVertex3fv(dot[j->p1].v);

                  t[0] = (sp->val[j->p2] - sp->vmin) * factor;
                  glTexCoord2fv(t);  glVertex3fv(dot[j->p2].v);

                  t[0] = (sp->val[j->p3] - sp->vmin) * factor;
                  glTexCoord2fv(t);  glVertex3fv(dot[j->p3].v);
                  glEnd();
            }
            if(alphablending) glDisable(GL_BLEND);
            disable_texture();
         }
         else {
            glColorMaterial(GL_FRONT, GL_DIFFUSE);
            glEnable(GL_COLOR_MATERIAL);
            factor = 64.99 / (sp->vmax - sp->vmin);
/* to be fixed
            if(shadow_offset) cp = shadow_col;
            else              cp = property_col;
*/
            cp = property_col;
            if(sp->alpha){
               unsigned long test;
               test = 0x00ffffff | (sp->alpha <<24);
               for(i=0; i<65; i++) cp[i] &= test;
            }
               glBegin(GL_TRIANGLES);
               for(i=0, vp=sp->trinorm; i<sp->ntri; i++, j++, vp++){
                  col = (int)((sp->val[j->p1] - sp->vmin) * factor);
                  set_color(cp[col]);
                  glNormal3fv((float *)vp);   glVertex3fv(dot[j->p1].v);
                  col = (int)((sp->val[j->p2] - sp->vmin) * factor);
                  set_color(cp[col]);
                  glVertex3fv(dot[j->p2].v);
                  col = (int)((sp->val[j->p3] - sp->vmin) * factor);
                  set_color(cp[col]);
                  glVertex3fv(dot[j->p3].v);
               }
               glEnd();
            glDisable(GL_COLOR_MATERIAL);
            if(sp->alpha) for(i=0; i<65; i++) cp[i] |= 0xff000000;
         }
      }
      else {
         if(bit.texture && sp->texture && sp->textype == TEX_PHONG){
            /* simulated phong */
            glDisable(GL_LIGHTING);
            glColor4f(1.0, 1.0, 1.0, 1.0);
            enable_texture(sp->texture, sp->texenv, sp->textype);
            if(alphablending){
               glEnable(GL_BLEND);
               glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            }
            get_pmat();
            //cpack(0xffffffff);
            glBegin(GL_TRIANGLES);
            for(i=0; i<sp->ntri; i++, j++){
               myt3f(dot[j->p1].n);   glVertex3fv(dot[j->p1].v);
               myt3f(dot[j->p2].n);   glVertex3fv(dot[j->p2].v);
               myt3f(dot[j->p3].n);   glVertex3fv(dot[j->p3].v);
            }
            glEnd();
            if(alphablending) glDisable(GL_BLEND);
            disable_texture();
            glEnable(GL_LIGHTING);
         }
         else {
            if(bit.texture && sp->texture) enable_texture(sp->texture, sp->texenv, sp->textype);
/* alpha blending */
            if(alphablending){
               glEnable(GL_BLEND);
               glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            }
            glBegin(GL_TRIANGLES);
            for(i=0, vp=sp->trinorm; i<sp->ntri; i++, j++, vp++){
               glNormal3fv((float *)vp);   glVertex3fv(dot[j->p1].v);
               glVertex3fv(dot[j->p2].v);
               glVertex3fv(dot[j->p3].v);
            }
            glEnd();
            if(alphablending) glDisable(GL_BLEND);
            if(bit.texture && sp->texture) disable_texture();
         }
      }

/* to be fixed
      blendfunction(BF_ONE, BF_ZERO);
*/

   }

   else draw_dots(sp);
}

/* to be fixed

void draw_positions(Mol *mp)
{
   register int i, j;
   float fraction;
   register Vector *vp;
   extern Dynamics dynamics;;

   fraction = npositions/64.;
   for(j=0, i=0, vp=position; j<64; j++){
      set_color(property_col[j]);
      bgnline();
      for(; i<(j+1)*fraction; i++, vp++){
         v3f((float *)vp);
         if(i>= dynamics.current) { endline(); return; }
      }
      endline();
      i--; vp--;
   }
}

*/

void draw_cubeplanes(Mol *mp)
{
   short colnum;
   register short ix, iy, iz;
   float v[3], t[3], dx, dy, dz, factor;
   float tip[3];

/*
   if(mp->wire_model || rapidmove){
      draw_cubeplanegrid(mp);
      return;
   }
*/

   glDisable(GL_LIGHTING);
   glEnable(GL_POLYGON_STIPPLE);
   glPolygonStipple(halftone1);

   dx = (mp->box.x2 - mp->box.x1)/(mp->box.nx-1.);
   dy = (mp->box.y2 - mp->box.y1)/(mp->box.ny-1.);
   dz = (mp->box.z2 - mp->box.z1)/(mp->box.nz-1.);

   switch (cubeplane_axis) {
      case 0:
         v[0] = mp->plane_ix*dx + mp->box.x1;
         if(bit.texture && (!firsttexture)) init_texture();
         if(bit.texture && actual_texture){
            glColor4f(1.0, 1.0, 1.0, 1.0);
            enable_texture(actual_texture, GL_MODULATE, TEX_MAP);
            glDisable(GL_TEXTURE_GEN_S);
            glDisable(GL_TEXTURE_GEN_T);
            factor = 1.0 / (cubemax - cubemin);
            t[1] = 0;
            for(iz=0; iz<mp->box.nz-1; iz++){
               v[2] = mp->box.z1 + iz*dz;;
               glBegin(GL_TRIANGLE_STRIP);
               v[1] = mp->box.y1;
               for(iy=0; iy<mp->box.ny; iy++, v[1] += dy){
                  t[0] = (mp->cube_value[iz][iy][mp->plane_ix] - cubemin) * factor;
                  glTexCoord2fv(t);
                  glVertex3fv(v);
                  v[2] += dz;
                  t[0] = (mp->cube_value[iz+1][iy][mp->plane_ix] - cubemin) * factor;
                  glTexCoord2fv(t);
                  glVertex3fv(v);
                  v[2] -= dz;
               }
               glEnd();
            }
      disable_texture();
         }
         else {
            factor = 64.98 / (cubemax - cubemin);
            for(iz=0; iz<mp->box.nz-1; iz++){
               v[2] = mp->box.z1 + iz*dz;;
               glBegin(GL_TRIANGLE_STRIP);
               v[1] = mp->box.y1;
               for(iy=0; iy<mp->box.ny; iy++, v[1] += dy){
                  colnum = (int)(0.01 + 
                     (mp->cube_value[iz][iy][mp->plane_ix] - cubemin) * factor);
                  if(colnum<0) set_color(0x7f7fff);
                  else if(colnum >64) set_color(0xff7f7f);
                  else set_color(property_col[colnum]);
                  glVertex3fv(v);
                  v[2] += dz;
                  colnum = (int)(0.01 +
                     (mp->cube_value[iz+1][iy][mp->plane_ix] - cubemin) * factor);
                  if(colnum<0) set_color(0x7f7fff);
                  else if(colnum >64) set_color(0xff7f7f);
                  else set_color(property_col[colnum]);
                  glVertex3fv(v);
                  v[2] -= dz;
               }
               glEnd();
            }
         }
      break;
      case 1:
         v[1] = mp->plane_iy*dy + mp->box.y1;
         if(bit.texture && (!firsttexture)) init_texture();
         if(bit.texture && actual_texture){
            glColor4f(1.0, 1.0, 1.0, 1.0);
            enable_texture(actual_texture, GL_MODULATE, TEX_MAP);
            glDisable(GL_TEXTURE_GEN_S);
            glDisable(GL_TEXTURE_GEN_T);
            factor = 1.0 / (cubemax - cubemin);
            t[1] = 0;
            for(ix=0; ix<mp->box.nx-1; ix++){
               v[0] = mp->box.x1 + ix*dx;;
               glBegin(GL_TRIANGLE_STRIP);
               v[2] = mp->box.z1;
               for(iz=0; iz<mp->box.nz; iz++, v[2] += dz){
                  t[0] = (mp->cube_value[iz][mp->plane_iy][ix] - cubemin) * factor;
                  glTexCoord2fv(t);
                  glVertex3fv(v);
                  v[0] += dx;
                  t[0] = (mp->cube_value[iz][mp->plane_iy][ix+1] - cubemin) * factor;
                  glTexCoord2fv(t);
                  glVertex3fv(v);
                  v[0] -= dx;
               }
               glEnd();
            }
            disable_texture();
         }
         else {
            factor = 64.98 / (cubemax - cubemin);
            for(ix=0; ix<mp->box.nx-1; ix++){
               v[0] = mp->box.x1 + ix*dx;;
               glBegin(GL_TRIANGLE_STRIP);
               v[2] = mp->box.z1;
               for(iz=0; iz<mp->box.nz; iz++, v[2] += dz){
                  colnum = (int)(0.01 + 
                     (mp->cube_value[iz][mp->plane_iy][ix] - cubemin) * factor);
                  if(colnum<0) set_color(0x7f7fff);
                  else if(colnum >64) set_color(0xff7f7f);
                  else set_color(property_col[colnum]);
                  glVertex3fv(v);
                  v[0] += dx;
                  colnum = (int)(0.01 +
                     (mp->cube_value[iz][mp->plane_iy][ix+1] - cubemin) * factor);
                  if(colnum<0) set_color(0x7f7fff);
                  else if(colnum >64) set_color(0xff7f7f);
                  else set_color(property_col[colnum]);
                  glVertex3fv(v);
                  v[0] -= dx;
               }
               glEnd();
            }
         }
      break;
      case 2:
         v[2] = mp->plane_iz*dz + mp->box.z1;
         if(bit.texture && (!firsttexture)) init_texture();
         if(bit.texture && actual_texture){
            glColor4f(1.0, 1.0, 1.0, 1.0);
            enable_texture(actual_texture, GL_MODULATE, TEX_MAP);
            glDisable(GL_TEXTURE_GEN_S);
            glDisable(GL_TEXTURE_GEN_T);
            factor = 1.0 / (cubemax - cubemin);
            t[1] = 0;
            for(iy=0; iy<mp->box.ny-1; iy++){
               v[1] = mp->box.y1 + iy*dy;;
               glBegin(GL_TRIANGLE_STRIP);
               v[0] = mp->box.x1;
               for(ix=0; ix<mp->box.nx; ix++, v[0] += dx){
                  t[0] = (mp->cube_value[mp->plane_iz][iy][ix] - cubemin) * factor;
                  glTexCoord2fv(t);
                  glVertex3fv(v);
                  v[1] += dy;
                  t[0] = (mp->cube_value[mp->plane_iz][iy+1][ix] - cubemin) * factor;
                  glTexCoord2fv(t);
                  glVertex3fv(v);
                  v[1] -= dy;
               }
               glEnd();
            }
      disable_texture();
         }
         else {
            factor = 64.98 / (cubemax - cubemin);
            for(iy=0; iy<mp->box.ny-1; iy++){
               v[1] = mp->box.y1 + iy*dy;;
               glBegin(GL_TRIANGLE_STRIP);
               v[0] = mp->box.x1;
               for(ix=0; ix<mp->box.nx; ix++, v[0] += dx){
                  colnum = (int)(0.01 + 
                     (mp->cube_value[mp->plane_iz][iy][ix] - cubemin) * factor);
                  if(colnum<0) set_color(0x7f7fff);
                  else if(colnum >64) set_color(0xff7f7f);
                  else set_color(property_col[colnum]);
                  glVertex3fv(v);
                  v[1] += dy;
                  colnum = (int)(0.01 +
                     (mp->cube_value[mp->plane_iz][iy+1][ix] - cubemin) * factor);
                  if(colnum<0) set_color(0x7f7fff);
                  else if(colnum >64) set_color(0xff7f7f);
                  else set_color(property_col[colnum]);
                  glVertex3fv(v);
                  v[1] -= dy;
               }
               glEnd();
            }
         }
      break;
      case 3:
/*
         tip[0] = mp->plane->a[0] + mp->plane->n[0];
         tip[1] = mp->plane->a[1] + mp->plane->n[1];
         tip[2] = mp->plane->a[2] + mp->plane->n[2];
         glBegin(GL_LINES);
         glVertex3fv(mp->plane->a);
         glVertex3fv(tip);
         glEnd();
         glBegin(GL_POINTS);
         for(iy=0; iy<mp->plane->npts; iy++) {
            for(ix=0; ix<mp->plane->npts; ix++) {
            if(mp->plane->plane_point[iy][ix][6]) {
               v[0] = mp->plane->plane_point[iy][ix][2];
               v[1] = mp->plane->plane_point[iy][ix][3];
               v[2] = mp->plane->plane_point[iy][ix][4];
               glVertex3fv(v);
            }
            }
         }
         glEnd();
*/
         if(bit.texture && (!firsttexture)) init_texture();
         if(bit.texture && actual_texture){
            glColor4f(1.0, 1.0, 1.0, 1.0);
            enable_texture(actual_texture, GL_MODULATE, TEX_MAP);
            glDisable(GL_TEXTURE_GEN_S);
            glDisable(GL_TEXTURE_GEN_T);
            factor = 1.0 / (cubemax - cubemin);
            t[1] = 0;
            glBegin(GL_TRIANGLES);
//printf("ntri: %d\n", mp->plane->ntri);
            for(ix=0; ix<mp->plane->ntri; ix++) {
               for(iy=0; iy<3; iy++) {
//printf("%d %d ", mp->plane->tri[ix][iy][0], mp->plane->tri[ix][iy][1]);
               v[0] = mp->plane->plane_point[mp->plane->tri[ix][iy][0]][mp->plane->tri[ix][iy][1]][2];
               v[1] = mp->plane->plane_point[mp->plane->tri[ix][iy][0]][mp->plane->tri[ix][iy][1]][3];
               v[2] = mp->plane->plane_point[mp->plane->tri[ix][iy][0]][mp->plane->tri[ix][iy][1]][4];
               t[0] = (mp->plane->plane_point[mp->plane->tri[ix][iy][0]][mp->plane->tri[ix][iy][1]][5] - cubemin) * factor;
               glTexCoord2fv(t);
//printf("x y z %f %f %f", v[0], v[1], v[2]);
               glVertex3fv(v);
               }
//printf("\n");
            }
            glEnd();
            disable_texture();
         }
         else {
            factor = 64.98 / (cubemax - cubemin);
            glBegin(GL_TRIANGLES);
//printf("ntri: %d\n", mp->plane->ntri);
            for(ix=0; ix<mp->plane->ntri; ix++) {
               for(iy=0; iy<3; iy++) {
//printf("%d %d ", mp->plane->tri[ix][iy][0], mp->plane->tri[ix][iy][1]);
               v[0] = mp->plane->plane_point[mp->plane->tri[ix][iy][0]][mp->plane->tri[ix][iy][1]][2];
               v[1] = mp->plane->plane_point[mp->plane->tri[ix][iy][0]][mp->plane->tri[ix][iy][1]][3];
               v[2] = mp->plane->plane_point[mp->plane->tri[ix][iy][0]][mp->plane->tri[ix][iy][1]][4];
               colnum = (int)(0.01 + 
                  (mp->plane->plane_point[mp->plane->tri[ix][iy][0]][mp->plane->tri[ix][iy][1]][5] - cubemin) * factor);
               if(colnum<0) set_color(0x7f7fff);
               else if(colnum >64) set_color(0xff7f7f);
               else set_color(property_col[colnum]);
//printf("x y z %f %f %f", v[0], v[1], v[2]);
               glVertex3fv(v);
               }
//printf("\n");
            }
            glEnd();
         }
      break;
   }
   glDisable(GL_POLYGON_STIPPLE);
   glEnable(GL_LIGHTING);
}


/* to be fixed

void draw_cubeplanegrid(Mol *mp)
{
   short colnum;
   register short ix, iy;
   float v[3], dx, dy, dz;

   dx = (cubehead.xmax - cubehead.xmin)/(cubehead.nx-1.);
   dy = (cubehead.ymax - cubehead.ymin)/(cubehead.ny-1.);
   dz = (cubehead.zmax - cubehead.zmin)/(cubehead.nz-1.);

   v[2] = plane_iz*dz + cubehead.zmin;
   for(iy=0; iy<cubehead.ny; iy++){
      v[1] = cubehead.ymin + iy*dy;
      bgnline();
      v[0] = cubehead.xmin;
      for(ix=0; ix<cubehead.nx; ix++, v[0] += dx){
         colnum = 0.01 +
            (mp->cube_value[plane_iz][iy][ix] - mp->cubemin)/(mp->cubemax-mp->cubemin)*64.98;
         if(colnum<0) set_color(0x7f7faf);
         else if(colnum >64) set_color(0xaf7f7f);
         else set_color(property_col[colnum]);
         v3f(v);
      }
      endline();
   }
   for(ix=0; ix<cubehead.nx; ix++){
      v[0] = cubehead.xmin + ix*dx;
      bgnline();
      v[1] = cubehead.ymin;
      for(iy=0; iy<cubehead.ny; iy++, v[1] += dy){
         colnum = 0.01 +
            (mp->cube_value[plane_iz][iy][ix] - mp->cubemin)/(mp->cubemax-mp->cubemin)*64.98;
         if(colnum<0) set_color(0x7f7faf);
         else if(colnum >64) set_color(0xaf7f7f);
         else set_color(property_col[colnum]);
         v3f(v);
      }
      endline();
   }
}





extern float cutoff;
float splat_size;
float splatv[4][3];



void draw_splat(int colnum, float x, float y, float z)
{
   Matrix m;

   pushmatrix();
   translate(x, y, z);
   getmatrix(m);
   m[0][0] = m[1][1] = m[2][2] = scale_factor;
   m[0][1] = m[0][2] = m[1][2] = m[1][0] = m[2][0] = m[2][1] = 0.0;
   loadmatrix(m);

   cpack(property_col[colnum]);
   rectf(-splat_size, -splat_size, splat_size, splat_size);

   popmatrix();
}



void draw_volume(Mol *mp)
{
   register int i, j, k;
   int colnum;
   float x, y, z, dx, dy, dz;

   if(mp->cube_value == NULL) return;

   set_scale_factor();

   dx = (mp->box.x2-mp->box.x1)/(mp->box.nx-1.);
   dy = (mp->box.y2-mp->box.y1)/(mp->box.ny-1.);
   dz = (mp->box.z2-mp->box.z1)/(mp->box.nz-1.);

   splat_size = (dx+dy+dz)/3.;

   if(getgdesc(GD_BLEND)) {
      unsigned long test;

      blendfunction(BF_SA, BF_MSA);
      test = 0x17ffffff;
      for(i=0; i<65; i++) property_col[i] &= test;
   }
   zwritemask(0);

   for(k=0, z=mp->box.z1; k<mp->box.nz; k++, z+=dz){
      for(j=0, y=mp->box.y1; j<mp->box.ny; j++, y+=dy){
         for(i=0, x=mp->box.x1; i<mp->box.nx; i++, x+=dx){
            if(mp->cube_value[k][j][i] > cutoff){
               colnum = 0.01 + (mp->cube_value[k][j][i] - mp->cubemin)/
                         (mp->cubemax-mp->cubemin)*64.98;
               draw_splat(colnum, x, y, z);
            }
         }
      }
   }

   for(i=0; i<65; i++) property_col[i] |= 0xff000000;
   zwritemask(0xffffffff);
}



*/

void draw_2Dgrid(Mol *mp)
{
   short colnum;
   register short i, j, k;
   float v[3], r1[3], r2[3], factor, t[2];;

   for(i=0; i<3; i++){
      r1[i] = (projection.p2[i] - projection.p1[i])/(projection.d1-1.);
      r2[i] = (projection.p3[i] - projection.p1[i])/(projection.d2-1.);
   }

   glDisable(GL_LIGHTING);
   glEnable(GL_POLYGON_STIPPLE);
   glPolygonStipple(halftone1);

   if(bit.texture && (!firsttexture)) init_texture();
   if(bit.texture && actual_texture){
      factor = 1.0 / (projection.vmax-projection.vmin);
      glColor4f(1.0, 1.0, 1.0, 1.0);
      enable_texture(actual_texture, GL_MODULATE, TEX_MAP);
      glDisable(GL_TEXTURE_GEN_S);
      glDisable(GL_TEXTURE_GEN_T);
      t[1] = 0;
      for(j=0; j<projection.d2-1; j++){
         glBegin(GL_TRIANGLE_STRIP);
         for(i=0; i<projection.d1; i++){
            for(k=0; k<3; k++){
               v[k] = projection.p1[k] + i*r1[k] + j*r2[k];
            }
            t[0] = (projection.value[j][i] - projection.vmin) * factor;
            glTexCoord2fv(t);
            glVertex3fv(v);
            for(k=0; k<3; k++){
               v[k] += r2[k];
            }
            t[0] = (projection.value[j+1][i] - projection.vmin) * factor;
            glTexCoord2fv(t);
            glVertex3fv(v);
         }
         glEnd();
      }
      disable_texture();
   }
   else {
      for(j=0; j<projection.d2-1; j++){
         glBegin(GL_TRIANGLE_STRIP);
         for(i=0; i<projection.d1; i++){
            for(k=0; k<3; k++){
               v[k] = projection.p1[k] + i*r1[k] + j*r2[k];
            }
            colnum = (int)(0.01 + (projection.value[j][i] - projection.vmin)/
                            (projection.vmax-projection.vmin)*64.98);
            set_color(property_col[colnum]);
            glVertex3fv(v);
            for(k=0; k<3; k++){
               v[k] += r2[k];
            }
            colnum = (int)(0.01 + (projection.value[j+1][i] - projection.vmin)/
                            (projection.vmax-projection.vmin)*64.98);
            set_color(property_col[colnum]);
            glVertex3fv(v);
         }
         glEnd();
      }
   }
   glDisable(GL_POLYGON_STIPPLE);
   glEnable(GL_LIGHTING);
}


void draw_distances(Mol *mp)
{
   Mon_dist *dp;
   float m[3];
   char str[20];

   glColor4fv(distcolor);
   glDisable(GL_LIGHTING);
   glCallList(line_style_base + 1);


   for(dp=mp->firstdist; dp; dp=dp->next){

      glBegin(GL_LINES);
      glVertex3fv(dp->a->coord);
      glVertex3fv(dp->b->coord);
      glEnd();

      sprintf(str, " %.3f", dist(dp->a->coord, dp->b->coord));

      m[0] = (dp->a->coord[0] + dp->b->coord[0])*0.5;
      m[1] = (dp->a->coord[1] + dp->b->coord[1])*0.5;
      m[2] = (dp->a->coord[2] + dp->b->coord[2])*0.5;

      if(!bit.no_measure_lbl) {
         if(bit.lbl_on_top) {
             glDisable(GL_DEPTH_TEST);
             draw_string(m[0],  m[1],  m[2], str);
             glEnable(GL_DEPTH_TEST);
         }
         else draw_string(m[0],  m[1],  m[2], str);
      }

   }

   glCallList(line_style_base);
   glEnable(GL_LIGHTING);
}


void getMovedPoints(Int_dist *ip, float p[4], float q[4])
{
   Matrix ma, mb;
   float a[4], b[4];
   int i;

   for(i=0; i<3; i++) {
      a[i] = ip->a->coord[i];
      b[i] = ip->b->coord[i];
   }
   a[3] = b[3] = 1.0;

   glPushMatrix();
   glLoadIdentity();
   individual_move(ip->ma);
   glGetFloatv(GL_MODELVIEW_MATRIX, &ma[0][0]);
   glPopMatrix();
   
   vecxmat(a, ma, p);

   glPushMatrix();
   glLoadIdentity();
   individual_move(ip->mb);
   glGetFloatv(GL_MODELVIEW_MATRIX, &mb[0][0]);
   glPopMatrix();

   vecxmat(b, mb, q);
}



static void drawIntermolecularDistance(Int_dist *ip)
{
   float p[4], q[4], m[3];
   char str[50];

   getMovedPoints(ip, p, q);

   glColor4fv(distcolor);
   glDisable(GL_LIGHTING);
   glCallList(line_style_base + 1);

   glBegin(GL_LINES);
   glVertex3fv(p);
   glVertex3fv(q);
   glEnd();

   sprintf(str, "%.3f", dist(p, q));

   m[0] = (p[0] + q[0])*0.5;
   m[1] = (p[1] + q[1])*0.5;
   m[2] = (p[2] + q[2])*0.5;

   if(bit.lbl_on_top) {
      glDisable(GL_DEPTH_TEST);
      draw_string(m[0],  m[1],  m[2], str);
      glEnable(GL_DEPTH_TEST);
   }
   else draw_string(m[0],  m[1],  m[2], str);

   glCallList(line_style_base);
   glEnable(GL_LIGHTING);
}



void drawIntermolecularDistances(void)
{
   extern Int_dist *first_intermolecular_distance;
   Int_dist *i;

   for(i=first_intermolecular_distance; i; i=i->next) {
      drawIntermolecularDistance(i);
   }
}



static void getEndPoint(float *a, float *b, float l, float *x)
{
   float length, factor;

   length = dist(a, b);
   if(!length) { x[0] = x[1] = x[2] = 0; return; }
   factor = l/length;

   x[0] = a[0] + factor*(b[0] - a[0]);
   x[1] = a[1] + factor*(b[1] - a[1]);
   x[2] = a[2] + factor*(b[2] - a[2]);
}


static void getHalfAngle(float *a, float *b, float *c, float l, float *x)
{
   float sum[3];

   sum[0] = a[0] + c[0] - b[0];
   sum[1] = a[1] + c[1] - b[1];
   sum[2] = a[2] + c[2] - b[2];
   getEndPoint(b, sum, l, x);
}




static void drawAngle(float a[3], float b[3], float c[3], float radius, int sign)
{
   int i;
   float p[17][3];
   char str[20];
   float angle;

   angle = sign * valence_angle(a, b, c);
   sprintf(str, " %.1f", angle);

   getEndPoint(b, a, radius, p[0]);
   getEndPoint(b, c, radius, p[8]);
   getHalfAngle(p[0], b, p[8], radius, p[4]);
   getHalfAngle(p[0], b, p[4], radius, p[2]);
   getHalfAngle(p[0], b, p[2], radius, p[1]);
   getHalfAngle(p[2], b, p[4], radius, p[3]);
   getHalfAngle(p[4], b, p[8], radius, p[6]);
   getHalfAngle(p[4], b, p[6], radius, p[5]);
   getHalfAngle(p[6], b, p[8], radius, p[7]);

   glBegin(GL_LINE_LOOP);
   glVertex3fv(b);
   for(i=0; i<9; i++) glVertex3fv(p[i]);
   glEnd();

   if(!bit.no_measure_lbl) {
      if(bit.lbl_on_top) {
         glDisable(GL_DEPTH_TEST);
         draw_string(p[4][0], p[4][1], p[4][2], str);
         glEnable(GL_DEPTH_TEST);
      }
      else draw_string(p[4][0], p[4][1], p[4][2], str);
   }

}




void draw_angles(Mol *mp)
{
   static float angleRadius = 0.5;
   Mon_ang *ap;

   glColor4fv(anglecolor);
   glDisable(GL_LIGHTING);
   glCallList(line_style_base + 1);

   for(ap=mp->firstang; ap; ap=ap->next){
      drawAngle(ap->a->coord, ap->b->coord, ap->c->coord, angleRadius, 1);
   }

   glCallList(line_style_base);
   glEnable(GL_LIGHTING);
}


#define DOTPRODUCT(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

void draw_torsions(Mol *mp)
{
   static float torsionRadius = 0.5;
   Mon_tor *tp;
   float m[3];  
   double n[4]; 
   float aPrime[3], dPrime[3];  
   double mu, nu;
   float *a, *b, *c, *d;
   int sign;

   glColor4fv(dihedcolor);
   glDisable(GL_LIGHTING);
   glCallList(line_style_base + 1);

   for(tp=mp->firsttor; tp; tp=tp->next){
      a = tp->a->coord;
      b = tp->b->coord;
      c = tp->c->coord;
      d = tp->d->coord;
      sign = (dihedral_angle(a, b, c, d) < 0) ? -1 : 1;

      m[0] = (b[0] + c[0])*0.5;
      m[1] = (b[1] + c[1])*0.5;
      m[2] = (b[2] + c[2])*0.5;
      n[0] = (c[0] - b[0]);
      n[1] = (c[1] - b[1]);
      n[2] = (c[2] - b[2]);
      n[3] = -DOTPRODUCT(n, m);
      mu = -(DOTPRODUCT(a, n) + n[3])/DOTPRODUCT(n, n);
      aPrime[0] = a[0] + mu*n[0];
      aPrime[1] = a[1] + mu*n[1];
      aPrime[2] = a[2] + mu*n[2];
      nu = -(DOTPRODUCT(d, n) + n[3])/DOTPRODUCT(n, n);
      dPrime[0] = d[0] + nu*n[0];
      dPrime[1] = d[1] + nu*n[1];
      dPrime[2] = d[2] + nu*n[2];

      drawAngle(aPrime, m, dPrime, torsionRadius, sign);
   }

   glCallList(line_style_base);
   glEnable(GL_LIGHTING);
}





void get_norm(float *a, float *b, float *c, Vector *d)
{
   float ab[3], ac[3], n[3], len;

   ab[0] = b[0]-a[0]; ab[1] = b[1]-a[1]; ab[2] = b[2]-a[2];
   ac[0] = c[0]-a[0]; ac[1] = c[1]-a[1]; ac[2] = c[2]-a[2];

   n[0] = ab[1]*ac[2] - ab[2]*ac[1];
   n[1] = ab[2]*ac[0] - ab[0]*ac[2];
   n[2] = ab[0]*ac[1] - ab[1]*ac[0];

   len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
   if(!len) d->x = d->y = d->z = 0;

   d->x = n[0]/len;
   d->y = n[1]/len;
   d->z = n[2]/len;
}


int make_trinorm(Surface *sp)
{
   Surfdot *dp;
   Triangle *tp;
   Vector *vp;
   int i;

   if((sp->trinorm = (Vector *)malloc(sp->ntri * sizeof(Vector))) == NULL){
      showinfobox("can't allocate triangle normals");
      return 0;
   }

   dp = sp->dot;
   for(tp=sp->tri, i=0, vp=sp->trinorm; i<sp->ntri; tp++, i++, vp++){
      get_norm(dp[tp->p1].v, dp[tp->p2].v, dp[tp->p3].v, vp);
   }

   return 1;
}

void draw_freq_arrow(Mol *mp)
{
   AtoM *ap;
   float tip[3], arrowcolor[4];
   register int i = 0;

   if(bit.smoothline){
      glEnable(GL_POINT_SMOOTH);
      glEnable(GL_LINE_SMOOTH);
      glBlendFunc(GL_SRC_ALPHA,  GL_ZERO);
      glEnable(GL_BLEND);
   }
   glCallList(arrow_mat);
   glGetMaterialfv(GL_FRONT, GL_DIFFUSE, arrowcolor);
   glColor4fv(arrowcolor);
   glDisable(GL_LIGHTING);
   glCallList(line_style_base + 2);

   for(ap=mp->firstatom; ap; ap=ap->next){
      tip[0] = ap->coord[0] + mp->sc_freq_ar * mp->freq_arrow->coord[i].x;
      tip[1] = ap->coord[1] + mp->sc_freq_ar * mp->freq_arrow->coord[i].y;
      tip[2] = ap->coord[2] + mp->sc_freq_ar * mp->freq_arrow->coord[i].z;
      
      glBegin(GL_LINES);
      glVertex3fv(ap->coord);
      glVertex3fv(tip);
      glEnd();

      i++;
   }
   glCallList(line_style_base);
   glEnable(GL_LIGHTING);
   if(bit.smoothline){
      glDisable(GL_POINT_SMOOTH);
      glDisable(GL_LINE_SMOOTH);
      glBlendFunc(GL_ONE,  GL_ZERO);
      glDisable(GL_BLEND);
   }
}

void draw_solid_freq_arrow(Mol *mp)
{
   Matrix bondMatrix;
   AtoM *ap, *atom1, *atom2;
   Bond *bond;
   float v[3];
   register int i = 0;

   if((atom1 = (AtoM*) malloc(sizeof(AtoM))) == NULL){
      showinfobox("can't allocate enough memory\n");
      return;
   }
   if((atom2 = (AtoM*) malloc(sizeof(AtoM))) == NULL){
      showinfobox("can't allocate enough memory\n");
      return;
   }
   if((bond = (Bond*) malloc(sizeof(Bond))) == NULL){
      showinfobox("can't allocate enough memory\n");
      return;
   }

   glCallList(arrow_mat);

   for(ap=mp->firstatom; ap; ap=ap->next){
      atom1->coord[0] = ap->coord[0] + mp->sc_freq_ar * mp->freq_arrow->coord[i].x;
      atom1->coord[1] = ap->coord[1] + mp->sc_freq_ar * mp->freq_arrow->coord[i].y;
      atom1->coord[2] = ap->coord[2] + mp->sc_freq_ar * mp->freq_arrow->coord[i].z;
      
      bond->a = ap;
      bond->b = atom1;

      glPushMatrix();
      getBondMatrix(bond, bondMatrix);
      glMultMatrixf(&bondMatrix[0][0]);
      if(quality != 0) glCallList(rough_stick_object_thin);
      else glCallList(stick_object_thin);
      glPopMatrix();

      vsub(atom1->coord, ap->coord, v);
      vnormal(v);
      vscale(v, 0.3);
      vadd(atom1->coord, v, atom2->coord);

      bond->a = atom1;
      bond->b = atom2;

      glPushMatrix();
      getBondMatrix(bond, bondMatrix);
      glMultMatrixf(&bondMatrix[0][0]);
      if(quality != 0) glCallList(rough_cone_object);
      else glCallList(cone_object);
      glPopMatrix();

      i++;
   }
}

void draw_dipole(Mol *mp)
{
   Matrix bondMatrix;
   AtoM *ap, *atom1, *atom2;
   Bond *bond;
   float v[3];
   float vec[3], factor;
   register int i = 0;

   if((atom1 = (AtoM*) malloc(sizeof(AtoM))) == NULL){
      showinfobox("can't allocate enough memory\n");
      return;
   }
   if((atom2 = (AtoM*) malloc(sizeof(AtoM))) == NULL){
      showinfobox("can't allocate enough memory\n");
      return;
   }
   if((bond = (Bond*) malloc(sizeof(Bond))) == NULL){
      showinfobox("can't allocate enough memory\n");
      return;
   }

   glCallList(arrow_mat);

   vec[0] =  (mp->dipole->end[0] - mp->dipole->start[0])/2;
   vec[1] =  (mp->dipole->end[1] - mp->dipole->start[1])/2;
   vec[2] =  (mp->dipole->end[2] - mp->dipole->start[2])/2;
   factor = mp->sc_dipole_ar - 1.0;

   atom1->coord[0] = mp->dipole->end[0] + vec[0] * factor;
   atom1->coord[1] = mp->dipole->end[1] + vec[1] * factor;
   atom1->coord[2] = mp->dipole->end[2] + vec[2] * factor;
   atom2->coord[0] = mp->dipole->start[0] - vec[0] * factor;
   atom2->coord[1] = mp->dipole->start[1] - vec[1] * factor;
   atom2->coord[2] = mp->dipole->start[2] - vec[2] * factor;
      
   bond->a = atom2;
   bond->b = atom1;

   glPushMatrix();
   getBondMatrix(bond, bondMatrix);
   glMultMatrixf(&bondMatrix[0][0]);
   if(quality != 0) glCallList(rough_stick_object_thin);
   else glCallList(stick_object_thin);
   glPopMatrix();

   vsub(atom1->coord, atom2->coord, v);
   vnormal(v);
   vscale(v, 0.3);
   vadd(atom1->coord, v, atom2->coord);

   bond->a = atom1;
   bond->b = atom2;

   glPushMatrix();
   getBondMatrix(bond, bondMatrix);
   glMultMatrixf(&bondMatrix[0][0]);
   if(quality != 0) glCallList(rough_cone_object);
   else glCallList(cone_object);
   glPopMatrix();

   free(atom1);
   free(atom2);
   free(bond);
}

void draw_labels(Matrix m, Mol *mp)
{
   register AtoM *atom;
   char str[20];

/* to be fixed
   if(bit.depthcue && !getgdesc(GD_FOGVERTEX)) depthcue(1);
*/
   if(bit.smoothline){
      glEnable(GL_POINT_SMOOTH);
      glEnable(GL_LINE_SMOOTH);
      glBlendFunc(GL_SRC_ALPHA,  GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_BLEND);
   }

   glPushMatrix();
   glTranslatef(mp->tvec[0], mp->tvec[1], mp->tvec[2]);
   if(!mp->wire_model){
      if(mp->spacefill) glTranslatef(0, 0, 2.0);
      else              glTranslatef(0, 0, 0.5);
   }
   glMultMatrixf(&m[0][0]);
   glTranslatef(mp->centervec[0], mp->centervec[1], mp->centervec[2]);

/* no need anymore
   if(mp->wire_model){
      if(bit.depthcue && !getgdesc(GD_FOGVERTEX)) depthcue(1);
      set_color(labelcolor);
   }
   else cpack(labelcolor);
*/
   glColor4fv(labelcolor);
   glDisable(GL_LIGHTING);

   for(atom = mp->firstatom; atom; atom = atom->next){
      if(mp->suppress_H && atom->ord == 1) continue;
      if(mp->residues && !atom->main) continue;
/* to be fixed
      if(mp==dynamics.molecule && dynamics.freeat && !showFixed) {
	  if(atom->fixed) continue;
      }
*/

      strcpy(str, " ");
      if(mp->labels) sprintf(str, " %s %d", element[atom->ord].symbol, atom->name);
      if(mp->atm_char) sprintf(str, "%s %5.3f", str, atom->charge);
      if(mp->atm_spin) sprintf(str, "%s %5.3f", str, atom->spin);
/*      sprintf(str, "%s %d", element[atom->ord].symbol, atom->name); */

      if(bit.lbl_on_top) {
         glDisable(GL_DEPTH_TEST);
         draw_string(atom->coord[0], atom->coord[1], atom->coord[2], str);
         glEnable(GL_DEPTH_TEST);
      }
      else draw_string(atom->coord[0], atom->coord[1], atom->coord[2], str);

   }

   glPopMatrix();

   glCallList(line_style_base);
   glEnable(GL_LIGHTING);

   if(bit.smoothline){
      glDisable(GL_POINT_SMOOTH);
      glDisable(GL_LINE_SMOOTH);
      glBlendFunc(GL_ONE,  GL_ZERO);
      glDisable(GL_BLEND);
   }

/* to be fixed
   depthcue(0);
*/
}


void draw_highlight_surf(Surface *sp)
{
   Triangle *tp;
   float vec[3];
   int i;

   tp = sp->tri;
   glBegin(GL_POINTS);
   for(i=0; i<sp->ntri; i++, tp++){
      vec[0] = sp->dot[tp->p1].v[0] + 0.01 * sp->dot[tp->p1].n[0];
      vec[1] = sp->dot[tp->p1].v[1] + 0.01 * sp->dot[tp->p1].n[1];
      vec[2] = sp->dot[tp->p1].v[2] + 0.01 * sp->dot[tp->p1].n[2];
      glVertex3fv(vec);
      vec[0] = sp->dot[tp->p2].v[0] + 0.01 * sp->dot[tp->p2].n[0];
      vec[1] = sp->dot[tp->p2].v[1] + 0.01 * sp->dot[tp->p2].n[1];
      vec[2] = sp->dot[tp->p2].v[2] + 0.01 * sp->dot[tp->p2].n[2];
      glVertex3fv(vec);
      vec[0] = sp->dot[tp->p3].v[0] + 0.01 * sp->dot[tp->p3].n[0];
      vec[1] = sp->dot[tp->p3].v[1] + 0.01 * sp->dot[tp->p3].n[1];
      vec[2] = sp->dot[tp->p3].v[2] + 0.01 * sp->dot[tp->p3].n[2];
      glVertex3fv(vec);
      }
   glEnd();
}



void draw_coord_system(Mol *mp)
{
   Matrix bondMatrix;
   AtoM *atom1, *atom2;
   Bond *bond;
   float v[3];
   float vec[3], factor;
   char *xyz[3] = {"X", "Y", "Z"};
   register int i = 0;

   if((atom1 = (AtoM*) malloc(sizeof(AtoM))) == NULL){
      showinfobox("can't allocate enough memory\n");
      return;
   }
   if((atom2 = (AtoM*) malloc(sizeof(AtoM))) == NULL){
      showinfobox("can't allocate enough memory\n");
      return;
   }
   if((bond = (Bond*) malloc(sizeof(Bond))) == NULL){
      showinfobox("can't allocate enough memory\n");
      return;
   }

   for(i=0; i<3; i++) {
      glCallList(def_mat_base + i + 1);
      atom1->coord[0] = 0.0;
      atom1->coord[1] = 0.0;
      atom1->coord[2] = 0.0;
      atom1->coord[i] = 1.0;
      atom2->coord[0] = 0.0;
      atom2->coord[1] = 0.0;
      atom2->coord[2] = 0.0;
      bond->a = atom2;
      bond->b = atom1;

      glPushMatrix();
      getBondMatrix(bond, bondMatrix);
      glMultMatrixf(&bondMatrix[0][0]);
      if(quality != 0) glCallList(rough_stick_object_thin);
      else glCallList(stick_object_thin);
      glPopMatrix();

      vsub(atom1->coord, atom2->coord, v);
      vnormal(v);
      vscale(v, 0.3);
      vadd(atom1->coord, v, atom2->coord);

      bond->a = atom1;
      bond->b = atom2;

      glPushMatrix();
      getBondMatrix(bond, bondMatrix);
      glMultMatrixf(&bondMatrix[0][0]);
      if(quality != 0) glCallList(rough_cone_object);
      else glCallList(cone_object);
      glPopMatrix();

      glColor4fv(labelcolor);
      glDisable(GL_LIGHTING);
      if(bit.smoothline){
         glEnable(GL_POINT_SMOOTH);
         glEnable(GL_LINE_SMOOTH);
         glBlendFunc(GL_SRC_ALPHA,  GL_ONE_MINUS_SRC_ALPHA);
         glEnable(GL_BLEND);
      }
      if(bit.lbl_on_top) {
         glDisable(GL_DEPTH_TEST);
         draw_string(1.3*atom1->coord[0], 1.3*atom1->coord[1], 1.3*atom1->coord[2], xyz[i]);
         glEnable(GL_DEPTH_TEST);
      }
      else draw_string(1.3*atom1->coord[0], 1.3*atom1->coord[1], 1.3*atom1->coord[2], xyz[i]);
      if(bit.smoothline){
         glDisable(GL_POINT_SMOOTH);
         glDisable(GL_LINE_SMOOTH);
         glBlendFunc(GL_ONE,  GL_ZERO);
         glDisable(GL_BLEND);
      }
      glEnable(GL_LIGHTING);
   }

   free(atom1);
   free(atom2);
   free(bond);
}


void draw_mapscale(Mol *mp)
{
  GLint viewport[4];
  Surface *sp;
  char str[20];
  int cp_labelsize, i;
  float delta, min;
  Texture *texture = NULL;
  GLdouble left, right, top, bottom;

  
  switch (maplegend) {
     case 1:
        if(mp != actualmol) return;
        if(mp->firstsurf) {
           sp = mp->firstsurf;
           if(sp->val){
              if(bit.texture && sp->texture){
                texture = sp->texture; 
                delta = (sp->vmax - sp->vmin) / 20;
                min = sp->vmin;
              } else return;
           } else return;
        } else return;
     break;
     case 2:
        if(mp == project_mol) {
           texture = actual_texture;
           delta = (projection.vmax-projection.vmin) / 20;
           min = projection.vmin;
        } else if(mp == actualmol) {
           if(mp->cubeplanes && mp->cube_value){
              if(bit.texture && (!firsttexture)) init_texture();
              if(bit.texture && actual_texture){
                 texture = actual_texture;
                 delta = (cubemax - cubemin) / 20;
                 min = cubemin;
              } else return;
           } else return;
        } else return;
     break;
     default: return;
  }

  glGetIntegerv(GL_VIEWPORT, viewport);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
/* render does not work! for the mapscale */
  if(bit.render) {
     left = 736 * (GLdouble) output_width / (GLdouble) output_height
                  * (tr->CurrentColumn * tr->TileWidthNB - tr->TileBorder) / tr->ImageWidth;
     right = left + 736 * (GLdouble) output_width / (GLdouble) output_height
                  * tr->CurrentTileWidth / tr->ImageWidth;
     bottom = 736 * (tr->CurrentRow * tr->TileHeightNB - tr->TileBorder) / tr->ImageHeight;
     top = bottom + 736 * tr->CurrentTileHeight / tr->ImageHeight;
     gluOrtho2D(left, right, bottom, top);
  }
  else gluOrtho2D(0.0, 736 * (GLdouble) viewport[2] / (GLdouble) viewport[3], 0.0, 736);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glColor4f(0.7, 0.7, 0.7, 1.0);
  glEnable(GL_TEXTURE_2D);
#ifdef NOTEXBIND
  glCallList(texture->texname);
#else
  glBindTexture(GL_TEXTURE_2D, texture->texname);
#endif
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, texture->texprop[0]);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, texture->texprop[1]);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, texture->texprop[2]);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, texture->texprop[3]);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

  glDisable(GL_DEPTH_TEST);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0, 0.0);
  glVertex2f(10.0, 10.0);
  glTexCoord2f(0.0, 0.0);
  glVertex2f( 30.0, 10.0);
  glTexCoord2f(1.0, 0.0);
  glVertex2f( 30.0,  370.0);
  glTexCoord2f(1.0, 0.0);
  glVertex2f(10.0,  370.0);
  glEnd();
  glDisable(GL_TEXTURE_2D);
  cp_labelsize = label_size;
  if(label_size == 4) label_size = 3;
  glColor4fv(labelcolor);
  glDisable(GL_LIGHTING);
  for(i=0; i<21; i++) {
     sprintf(str, "%9.5f", min + i * delta);
     draw_string(35.0, 10.0 + i * 18, 0.0, str);
  }
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  label_size = cp_labelsize;
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
}
