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


#include "constant.h"
#include "main.h"
#include "molekel.h"
#include "general.h"
#include "glutwin.h"
#include "drawing.h"
#include "drawmol.h"
#include "utils.h"
#include "material.h"
#include "objects.h"
#include "manip.h"
#include "box.h"
#include "texture.h"

/*----------------------------------------*/

unsigned long drawcol[] =
  { 0xffffffff, 0xff0000ff, 0xff00ff00, 0xffff0000, 0xff00ffff,
    0xffffff00, 0xffff00ff, 0xff00ff00, 0xff00007f, 0xff000000};
//float atomcolor[4];

/*----------------------------------------*/

void getBondMatrix(register Bond *bond, Matrix matrix);

/*----------------------------------------*/

void draw_structure(Mol *mp)
{
   if(mp->wire_model || rapidmove) {
      draw_wiremodel(mp);
   }
   else draw_solidmodel(mp);
}


void draw_wiremodel(Mol *mp)
{
   register AtoM *atom;
   register Bond *bond;
   AtoM *atomA, *atomB;
   float middle[3];

   if(bit.smoothline){
      glEnable(GL_POINT_SMOOTH);
      glEnable(GL_LINE_SMOOTH);
      glBlendFunc(GL_SRC_ALPHA,  GL_ZERO);
      glEnable(GL_BLEND);
      glCallList(line_style_base + 4);
   }

   glDisable(GL_LIGHTING);
   glCallList(line_style_base + 5);

/*** draw non-bonded atoms ***/
   glBegin(GL_POINTS);
   for(atom = mp->firstatom; atom; atom = atom->next){
      if(atom->nbonds) continue;
      if(mp->suppress_H && (atom->ord == 1)) continue;
      if(mp==dynamics.molecule && dynamics.freeat && !bit.showfixed) {
          if(atom->fixed) continue;
      }
      if(mp->residues && !atom->main) continue;

//      glCallList(element_mat_base + atom->ord);
//      glGetMaterialfv(GL_FRONT, GL_DIFFUSE, atomcolor);
//      glColor4fv(atomcolor);
      set_color(drawcol[element[atom->ord].col]);
      glVertex3fv(atom->coord);
   }
   glEnd();


/*** draw the bonds ***/
   for(bond = mp->firstbond; bond; bond = bond->next){
        atomA = bond->a;
        atomB = bond->b;
      if(mp->suppress_H){
         if(atomA->ord == 1 || atomB->ord == 1) continue;
      }
      if(mp==dynamics.molecule && dynamics.freeat && !bit.showfixed) {
          if(atomA->fixed || atomB->fixed) continue;
      }
      if(mp->residues) {
         if(!atomA->main || !atomB->main) continue;
      }
      if(bond->type == H_BOND) {
         if(mp->h_bond) {
            glColor4fv(distcolor);
            glCallList(line_style_base + 1);
            glBegin(GL_LINES);
            glVertex3fv(atomA->coord);
            glVertex3fv(atomB->coord);
            glEnd();
            glCallList(line_style_base);
         }
      }
      else {
         middle[0] = (atomA->coord[0] + atomB->coord[0])*0.5;
         middle[1] = (atomA->coord[1] + atomB->coord[1])*0.5;
         middle[2] = (atomA->coord[2] + atomB->coord[2])*0.5;
//      glCallList(element_mat_base + atomA->ord);
//      glGetMaterialfv(GL_FRONT, GL_DIFFUSE, atomcolor);
//      glColor4fv(atomcolor);
         set_color(drawcol[element[atomA->ord].col]);
         glBegin(GL_LINES);
         glVertex3fv(atomA->coord);
         glVertex3fv(middle);
         glEnd();
//      glCallList(element_mat_base + atomB->ord);
//      glGetMaterialfv(GL_FRONT, GL_DIFFUSE, atomcolor);
//      glColor4fv(atomcolor);
         set_color(drawcol[element[atomB->ord].col]);
         glBegin(GL_LINES);
         glVertex3fv(middle);
         glVertex3fv(atomB->coord);
         glEnd();
      }
   }
   if(bit.smoothline){
      glDisable(GL_POINT_SMOOTH);
      glDisable(GL_LINE_SMOOTH);
      glBlendFunc(GL_ONE,  GL_ZERO);
      glDisable(GL_BLEND);
      glCallList(line_style_base);
   }

   if(mp->box_on) draw_wirebox();

/* to be fixed!!!
   depthcue(0);
*/
   glCallList(line_style_base);
   glEnable(GL_LIGHTING);

   if(picking_mode) {
      for(atom = mp->firstatom; atom; atom = atom->next) {
          if(atom->picked) draw_atom(mp, atom);
      }
   }
}

void draw_solidmodel(Mol *mp)
{
   register Bond *bond;

   set_scale_factor();

   drawAllAtoms(mp);

   if(mp->ball_and_stick){
      glCallList(bond_mat);
/* to be fixed
      lmbind(MATERIAL, 10 + shadow_offset);
*/
      if(bit.texture && bond_texture) enable_texture(bond_texture, bond_texenv, bond_textype);
      for(bond = mp->firstbond; bond; bond = bond->next){
         draw_bond(mp, bond);
      }
      if(bit.texture && bond_texture) disable_texture();
   }
   else if(mp->sticks || (mp->spacefill && mp->transparent)){
      if(bit.texture && bond_texture) enable_texture(bond_texture, bond_texenv, bond_textype);
      for(bond = mp->firstbond; bond; bond = bond->next){
         draw_color_stick(mp, bond);
      }
      if(bit.texture && bond_texture) disable_texture();
   }

/* to be fixed
   if(mp->box_on){
      lmbind(MATERIAL, 8 + shadow_offset);
      draw_solidbox();
   }
*/

   if(mp->box_on){
      glCallList(box_mat);
      draw_solidbox();
   }
}

void drawAllAtoms(Mol *mp)
{
   float factor;
   Matrix m;
   GLuint sphere_obj;
   AtomType *at;
   AtoM *ap;
   int i;

   for(at = mp->firstAtomType; at; at = at->next) {

      if(mp->suppress_H && at->ord == 1) continue;

      if(mp->spacefill) {
         factor = element[at->ord].rvdw;
         if(quality == 0) sphere_obj = fine_sphere;
         else if(quality == 1) sphere_obj = medium_sphere;
         else sphere_obj = rough_sphere;
         if(mp->transparent) {
            glEnable(GL_POLYGON_STIPPLE);
            glPolygonStipple(halftone1);
         }
      }
      else if(mp->ball_and_stick){
         factor = element[at->ord].rvdw * vdW_factor;
         if(quality == 0 ) sphere_obj = fine_sphere;
         else if(quality == 1) sphere_obj = medium_sphere;
         else sphere_obj = rough_sphere;
      }
      else if(mp->sticks){
//         if(!picking_mode && !bit.high_qual) return;
         if(!picking_mode && (quality == 2)) return;
         factor = 0.5*bond_factor;
         if(quality == 0) sphere_obj = medium_sphere;
         else sphere_obj = rough_sphere;
      }

      if(mp->wire_model) { /* picking in wire - model ! */
         factor = element[at->ord].rvdw * 0.1;
         sphere_obj = ugly_sphere;
      }

/* to be fixed
      lmbind(MATERIAL, at->ord + 1 + shadow_offset);

*/
      if(bit.texture && element[at->ord].texture) 
         enable_texture(element[at->ord].texture, element[at->ord].texenv,
                     element[at->ord].textype);

      glCallList(element_mat_base + at->ord);

      for(i=0; i<at->natoms; i++) {
         if(mp->sticks && !picking_mode &&
            (quality == 2) && at->atomList[i]->nbonds) continue;

         if(mp->residues) {
            if(!at->atomList[i]->main) continue;
         }

         if(mp==dynamics.molecule && dynamics.freeat && !bit.showfixed) {
             if(at->atomList[i]->fixed) continue;
         }

         glPushMatrix();

         glTranslatef(at->atomList[i]->coord[0],
                   at->atomList[i]->coord[1],
                   at->atomList[i]->coord[2]);

/*** draw a hemisphere, setting the rotation to 0 ***/
         glGetFloatv(GL_MODELVIEW_MATRIX, &m[0][0]);
         m[0][0] = m[1][1] = m[2][2] = scale_factor*factor;
         m[0][1] = m[0][2] = m[1][2] = m[1][0] = m[2][0] = m[2][1] = 0.0;
         glLoadMatrixf( &m[0][0]);

/* to be fixed
         if(at->atomList[i]->coordination != 0.0) {
            long cValue;
            long red, blue;

            blue = (at->atomList[i]->coordination-2.5)*200.0/1.5;
            red  = 200 - blue;
            if(red < 0) red = 0;
            cValue = 0xff002000 | red | (blue << 16);
            lmcolor(LMC_DIFFUSE);
            cpack(cValue);
            lmcolor(LMC_COLOR);
         }
*/

         glCallList(sphere_obj);
         glPopMatrix();
      }

      if(bit.texture && element[at->ord].texture) disable_texture();
   }
   if(mp->transparent) glDisable(GL_POLYGON_STIPPLE);

   if(picking_mode) {
      for(ap = mp->firstatom; ap; ap = ap->next) {
          if(ap->picked) draw_atom(mp, ap);
      }
   }
}

/* draw a sphere at the coordinates of the atom */
void draw_atom(Mol *mp, register AtoM *atom)
{
   float factor;
   Matrix m;
   GLuint sphere_obj;

   if(mp->suppress_H){
      if(atom->ord == 1) return;
   }
   if(mp==dynamics.molecule && dynamics.freeat && !bit.showfixed) {
      if(atom->fixed) return;
   }
   if(mp->spacefill) {
      factor = element[atom->ord].rvdw * 1.01;
      if(quality == 0) sphere_obj = fine_sphere;
      else if(quality == 1) sphere_obj = medium_sphere;
      else sphere_obj = medium_sphere;
      if(mp->transparent) {
         glEnable(GL_POLYGON_STIPPLE);
         glPolygonStipple(halftone1);
      }
   }
   else if(mp->ball_and_stick){
      factor = element[atom->ord].rvdw * vdW_factor * 1.01;
      if(quality == 0) sphere_obj = fine_sphere;
      else if(quality == 1) sphere_obj = medium_sphere;
      else sphere_obj = rough_sphere;
   }
   else if(mp->sticks){
      if(!picking_mode && atom->nbonds && (quality == 2)) return;
      factor = 0.5*bond_factor * 1.4;
      if(quality == 0) sphere_obj = medium_sphere;
      else sphere_obj = rough_sphere;
   }

   if(mp->wire_model) { /* picking in wire - model ! */
      factor = element[atom->ord].rvdw * 0.1;
      sphere_obj = ugly_sphere;
   }

   glPushMatrix();

   glTranslatef(atom->coord[0], atom->coord[1], atom->coord[2]);

/*** draw a hemisphere, setting the rotation to 0 ***/
   glGetFloatv(GL_MODELVIEW_MATRIX, &m[0][0]);
   m[0][0] = m[1][1] = m[2][2] = scale_factor*factor;
   m[0][1] = m[0][2] = m[1][2] = m[1][0] = m[2][0] = m[2][1] = 0.0;
   glLoadMatrixf( &m[0][0]);

/* to be fixed
   if(atom->picked) lmbind(MATERIAL, 11 + shadow_offset);
   else lmbind(MATERIAL, atom->ord + 1 + shadow_offset);
*/
   if(atom->picked) {
     glCallList(pick_mat);
   }
   else  glCallList(element_mat_base + atom->ord);

   if(bit.texture && element[atom->ord].texture) 
         enable_texture(element[atom->ord].texture, element[atom->ord].texenv,
          element[atom->ord].textype);

   glCallList(sphere_obj);

   if(bit.texture && element[atom->ord].texture) disable_texture();

   glPopMatrix();
   if(mp->transparent) glDisable(GL_POLYGON_STIPPLE);
}

void draw_stobj(Mol *mp, register Bond *bond)
{
/* to be fixed
   if(bond->picked) {
      lmbind(MATERIAL, 11 + shadow_offset);
      callobj(stick_object);
      return;
   }
*/

/* to be fixed
   if(!mp->bond_col) lmbind(MATERIAL, bond->a->ord + 1 + shadow_offset);
*/
   if(!mp->bond_col) glCallList(element_mat_base  + bond->a->ord);
   if(quality != 0) glCallList(rough_stick_object1);
   else glCallList(stick_object1);
/* to be fixed
   if(!mp->bond_col) lmbind(MATERIAL, bond->b->ord + 1 + shadow_offset);
*/
   if(!mp->bond_col) glCallList(element_mat_base + bond->b->ord);
   if(quality != 0) glCallList(rough_stick_object2);
   else glCallList(stick_object2);
   
}


void draw_bond(Mol *mp, register Bond *bond)
{
    Matrix bondMatrix;

   if(mp->suppress_H){
      if(bond->a->ord == 1 || bond->b->ord == 1) return;
   }
   if(mp->residues) {
      if(!bond->a->main || !bond->b->main) return;
   }
   if(mp==dynamics.molecule && dynamics.freeat && !bit.showfixed) {
      if(bond->a->fixed || bond->b->fixed) return;
   }

   glPushMatrix();

   getBondMatrix(bond, bondMatrix);
   glMultMatrixf(&bondMatrix[0][0]);

   if(mp->ball_and_stick && mp->multiple_bonds) {
      switch(bond->type){
         case SINGLEBOND : draw_stobj(mp, bond);         break;
         case DOUBLEBOND : glPushMatrix();
                           glTranslatef(0.13, 0, 0);
                           draw_stobj(mp, bond);
                           glTranslatef(-0.26, 0, 0);
                           draw_stobj(mp, bond);
                           glPopMatrix();
                           break;
         case TRIPLEBOND : glPushMatrix();
                           draw_stobj(mp, bond);
                           glTranslatef(0.18, 0, 0);
                           draw_stobj(mp, bond);
                           glTranslatef(-0.36, 0, 0);
                           draw_stobj(mp, bond);
                           glPopMatrix();
                           break;
         case H_BOND     : if(mp->h_bond) {
                              glCallList(bond_mat);
                              glCallList(dot_stick_object);
                           }
                           break;
      }
   }
   else if (bond->type == H_BOND) {
      if (mp->h_bond) {
         glCallList(bond_mat);
         glCallList(dot_stick_object);
      }
   }
   else draw_stobj(mp, bond);


   glPopMatrix();

}

void draw_color_stick(Mol *mp, register Bond *bond)
{
   Matrix bondMatrix;

   if(mp->suppress_H){
      if(bond->a->ord == 1 || bond->b->ord == 1) return;
   }
   if(mp->residues) {
      if(!bond->a->main || !bond->b->main) return;
   }
   if(mp==dynamics.molecule && dynamics.freeat && !bit.showfixed) {
      if(bond->a->fixed || bond->b->fixed) return;
   }

   glPushMatrix();

   getBondMatrix(bond, bondMatrix);
   glMultMatrixf(&bondMatrix[0][0]);

   if(bond->picked) {
      glCallList(pick_mat);
/* to be fixed
      lmbind(MATERIAL, 11 + shadow_offset);
*/
      if(quality != 0) glCallList(rough_stick_object);
      else glCallList(stick_object);
      glPopMatrix();
      return;
   }

/* to be fixed
   lmbind(MATERIAL, bond->a->ord + 1 + shadow_offset);
*/
   if(bond->type == H_BOND) {
      if(mp->h_bond) {
         glCallList(bond_mat);
         glCallList(dot_stick_object);
      }
   }
   else {
      glCallList(element_mat_base + bond->a->ord);
      if(quality != 0) glCallList(rough_stick_object1);
      else glCallList(stick_object1);
/* to be fixed
   lmbind(MATERIAL, bond->b->ord + 1 + shadow_offset);
*/
      glCallList(element_mat_base + bond->b->ord);
      if(quality != 0) glCallList(rough_stick_object2);
      else glCallList(stick_object2);
   }

   glPopMatrix();

}

void getBondMatrix(register Bond *bond, Matrix matrix)
{
   float l, v[3];
   float alpha, beta;
   static float konst = -180./M_PI;
   extern void vsub(float *src1, float *src2, float *dst);
   extern float dist(float *a, float *b);

   l = dist(bond->a->coord, bond->b->coord);
   vsub(bond->b->coord, bond->a->coord, v);

   alpha = konst*acos(v[2]/l);
   if(v[1] < 0) alpha = -alpha;

   if(v[1]) beta = konst*atan(v[0]/v[1]);
   else if(v[0] < 0) beta = 90.0;
   else              beta = -90.0;

   glPushMatrix();
   glLoadIdentity();
   glTranslatef(bond->a->coord[0], bond->a->coord[1], bond->a->coord[2]);
   glRotatef(beta, 0.0, 0.0, 1.0);
   glRotatef(alpha, 1.0, 0.0, 0.0);
   glScalef(1., 1., l);

   glGetFloatv(GL_MODELVIEW_MATRIX, &matrix[0][0]);
   glPopMatrix();
}

void draw_wirebox(void)
{
   float vert[8][3];

   vert[0][0] = box.x1; vert[0][1] = box.y1; vert[0][2] = box.z1;
   vert[1][0] = box.x2; vert[1][1] = box.y1; vert[1][2] = box.z1;
   vert[2][0] = box.x2; vert[2][1] = box.y2; vert[2][2] = box.z1;
   vert[3][0] = box.x1; vert[3][1] = box.y2; vert[3][2] = box.z1;
   vert[4][0] = box.x1; vert[4][1] = box.y1; vert[4][2] = box.z2;
   vert[5][0] = box.x2; vert[5][1] = box.y1; vert[5][2] = box.z2;
   vert[6][0] = box.x2; vert[6][1] = box.y2; vert[6][2] = box.z2;
   vert[7][0] = box.x1; vert[7][1] = box.y2; vert[7][2] = box.z2;

   set_color(0xff0000ff);
   glBegin(GL_LINE_LOOP);
   glVertex3fv(vert[0]); glVertex3fv(vert[3]); glVertex3fv(vert[2]); glVertex3fv(vert[1]);
   glEnd();
   glBegin(GL_LINE_LOOP);
   glVertex3fv(vert[4]); glVertex3fv(vert[5]); glVertex3fv(vert[6]); glVertex3fv(vert[7]);
   glEnd();
   glBegin(GL_LINES);
   glVertex3fv(vert[0]); glVertex3fv(vert[4]);
   glEnd();
   glBegin(GL_LINES);
   glVertex3fv(vert[1]); glVertex3fv(vert[5]);
   glEnd();
   glBegin(GL_LINES);
   glVertex3fv(vert[2]); glVertex3fv(vert[6]);
   glEnd();
   glBegin(GL_LINES);
   glVertex3fv(vert[3]); glVertex3fv(vert[7]);
   glEnd();
}

void draw_solidbox(void)
{
   static float norm[][3] = {{0, 0, -1}, {0, -1, 0}, {1, 0, 0},
                             {0, 1, 0}, {-1, 0, 0}, {0, 0, 1}};
   static short side[][4] = {{0, 3, 2, 1}, {0, 1, 5, 4}, {1, 2, 6, 5},
                             {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7} };
   short i;
   float vert[8][3];

   vert[0][0] = box.x1; vert[0][1] = box.y1; vert[0][2] = box.z1;
   vert[1][0] = box.x2; vert[1][1] = box.y1; vert[1][2] = box.z1;
   vert[2][0] = box.x2; vert[2][1] = box.y2; vert[2][2] = box.z1;
   vert[3][0] = box.x1; vert[3][1] = box.y2; vert[3][2] = box.z1;
   vert[4][0] = box.x1; vert[4][1] = box.y1; vert[4][2] = box.z2;
   vert[5][0] = box.x2; vert[5][1] = box.y1; vert[5][2] = box.z2;
   vert[6][0] = box.x2; vert[6][1] = box.y2; vert[6][2] = box.z2;
   vert[7][0] = box.x1; vert[7][1] = box.y2; vert[7][2] = box.z2;

   glEnable(GL_POLYGON_STIPPLE);
   glPolygonStipple(halftone1);
   for(i=0; i<6; i++){
      glBegin(GL_POLYGON);
      glNormal3fv(norm[i]);
      glVertex3fv(vert[side[i][0]]);
      glVertex3fv(vert[side[i][1]]);
      glVertex3fv(vert[side[i][2]]);
      glVertex3fv(vert[side[i][3]]);
      glEnd();
   }
   glDisable(GL_POLYGON_STIPPLE);
}
