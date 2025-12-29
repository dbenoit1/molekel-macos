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


#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "general.h"
#include "glutwin.h"
#include "trackball.h"
#include "pick.h"
#include "geometry.h"
#include "utils.h"
#include "manip.h"
#include "drawing.h"
#include "drawmol.h"
#include "maininterf.h"
#include "material.h"
#include "cubeInterpol.h"
#include "menu.h"
#include "calcdens.h"
#include "drawprop.h"
#include "chooseinterf.h"
#include "objects.h"
#include "action.h"
#include "colorinterf.h"
#include "texture.h"
#include "textureinterf.h"

int pick_key = -1;

void pick_action(int key, int x, int y)
{
   static AtoM *ap = NULL, *apd1 = NULL, *apd2 = NULL, *apa1 = NULL, *apa2 = NULL,\
                 *apa3 =  NULL, *apt1 = NULL, *apt2 = NULL, *apt3 = NULL, *apt4 = NULL;
   static Bond *bp = NULL;
   static Surface *sp = NULL;
   static Mol *molp = NULL;
   static Mol *mpsame = NULL, *m1 = NULL, *m2 = NULL;
   static AtoM **a1, **a2;
   static int divisor = 0, divisor_center = 0, na1 = 0, na2 = 0;
   static float xcoor = 0.0, ycoor = 0.0, zcoor = 0.0;
   static float mvec[3] = {0.0, 0.0, 0.0};
   Mol *mp;
   Matrix m;

   switch (key) {
      case PICK_MOL:
         if(!(molp = (Mol *)my_pick(PICK_MOL, x, y))) return;
         actualsurf = actualmol->firstsurf;
         molp = NULL;
         pick_action(DONE, 1, 1);
      break;
      case PICK_SURF:
         if(!(sp = (Surface *)my_pick(PICK_SURF, x, y))) return;
         actualsurf = sp;
         sp = NULL;
         pick_action(DONE, 1, 1);
      break;
      case MEP:
         if(!(sp = (Surface *)my_pick(PICK_SURF, x, y))) return;
         mainglui_top->enable();
         glutSetWindow(mainwin);
         mep_dot_surface(sp);
         sp = NULL;
         pick_action(DONE, 1, 1);
      break;
      case GRID_VALUES:
         if(!(sp = (Surface *)my_pick(PICK_SURF, x, y))) return;
         gridValues(sp);
         sp = NULL;
         pick_action(DONE, 1, 1);
      break;
      case ATOMDATA:
         if(!(ap = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
         measure_atom(ap);
         ap = NULL;
      break;
      case DISTANCE:
         if(!apd1) {
            if(!(apd1 = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
            mpsame = actualmol;
            return;
         }
         if(!apd2) {
            if(!(apd2 = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
            if(apd1 == apd2) {
               apd2 = NULL;
               return;
            }
            if(actualmol == mpsame) measure_distance(apd1, apd2);
            else measure_intermolecular_distance(apd1, mpsame, apd2, actualmol);
            apd1 = apd2 = NULL;
         }
       break;
       case VALENCE:
         if(!apa1) {
            if(!(apa1 = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
            mpsame = actualmol;
            return;
         }
         if(!apa2) {
            if(!(apa2 = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
            if(apa1 == apa2) {
               apa2 = NULL;
               return;
            }
            if(actualmol != mpsame) {
               apa2->picked = 0;
               apa2 = NULL;
               logprint("atom not in the same molecule");
               update_logs();
               glutPostRedisplay();
               return;
            }
            return;
         }
         if(!apa3) {
         if(!(apa3 = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
            if(apa2 == apa3) {
               apa3 = NULL;
               return;
            }
            if(actualmol != mpsame) {
               apa3->picked = 0;
               apa3 = NULL;
               logprint("atom not in the same molecule");
               update_logs();
               glutPostRedisplay();
               return;
            }   
            measure_valence(apa1, apa2, apa3);
            apa1 = apa2 = apa3 = NULL;
         }
       break;
       case DIHEDRAL:
         if(!apt1) {
            if(!(apt1 = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
            mpsame = actualmol;
            return;
         }
         if(!apt2) {
            if(!(apt2 = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
            if(apt1 == apt2) {
               apt2 = NULL;
               return;
            }
            if(actualmol != mpsame) {
               apt2->picked = 0;
               apt2 = NULL;
               logprint("atom not in the same molecule");
               update_logs();
               glutPostRedisplay();
               return;
            }
            return;
         }
         if(!apt3) {
            if(!(apt3 = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
            if(apt2 == apt3) {
               apt3 = NULL;
               return;
            }
            if(actualmol != mpsame) {
               apt3->picked = 0;
               apt3 = NULL;
               logprint("atom not in the same molecule");
               update_logs();
               glutPostRedisplay();
               return;
            }
            return;
         }
         if(!apt4) {
         if(!(apt4 = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
            if(apt3 == apt4) {
               apt4 = NULL;
               return;
            }
            if(actualmol != mpsame) {
               apt4->picked = 0;
               apt4 = NULL;
               logprint("atom not in the same molecule");
               update_logs();
               glutPostRedisplay();
               return;
            }   
            measure_dihedral(apt1, apt2, apt3, apt4);
            apt1 = apt2 = apt3 = apt4 = NULL;
         }
       break;
       case ADD_DUMMY:
          if(!(ap = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
          if(!mpsame) mpsame = actualmol;
          else if (mpsame != actualmol) {
             ap->picked = 0;
             ap = NULL;
             actualmol = mpsame;
             logprint("atom not in the same molecule");
             update_logs();
             glutPostRedisplay();
             return;
          } 
          xcoor += ap->coord[0];
          ycoor += ap->coord[1];
          zcoor += ap->coord[2];
          divisor++;
       break;
       case MATCH:
          if(!a1){
             if((a1 = (AtoM **)malloc(sizeof(AtoM *))) == NULL){
                showinfobox("match:\ncan't allocate atom-pointer 1");
                return;
             }
          }
          if(!a2){
             if((a2 = (AtoM **)malloc(sizeof(AtoM *))) == NULL){
                showinfobox("match:\ncan't allocate atom-pointer 2");
                return;
             }
          }
          if(!(ap = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
          if(!m1) m1 = actualmol;
          else if(actualmol != m1 && !m2) m2 = actualmol;
          if(actualmol == m1) { /* add atom to list a1 */
             na1++;
             if((a1 = (AtoM **)realloc(a1, na1 * sizeof(AtoM *))) == NULL){
                showinfobox("match:\ncan't reallocate atom-pointer 1");
                free(a1); free(a2);
                a1 = a2 = NULL;
                m1 = m2 = NULL;
                na1 = na2 = 0;
                return;
             }
             a1[na1-1] = ap;
          }
          else if(actualmol == m2) { /* add atom to list a2 */
             na2++;
             if((a2 = (AtoM **)realloc(a2, na2 * sizeof(AtoM *))) == NULL){
                showinfobox("match:\ncan't reallocate atom-pointer 1");
                free(a1); free(a2);
                a1 = a2 = NULL;
                m1 = m2 = NULL;
                na1 = na2 = 0;
                return;
             }
             a2[na2-1] = ap;
          }
          else {
             showinfobox("match only two molecules at once");
             ap->picked = 0;
             glutPostRedisplay();
             return;
          }
       break;
       case CENTER:
          if(!(ap = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
          if(!mpsame) mpsame = actualmol;
          else if (mpsame != actualmol) {
             ap->picked = 0;
             ap = NULL;
             actualmol = mpsame;
             logprint("atom not in the same molecule");
             update_logs();
             glutPostRedisplay();
             return;
          } 
          mvec[0] -= ap->coord[0];
          mvec[1] -= ap->coord[1];
          mvec[2] -= ap->coord[2];
          divisor_center++;
       break;
       case ADD_BOND:
         if(!apd1) {
            if(!(apd1 = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
            mpsame = actualmol;
            return;
         }
         if(!apd2) {
            if(!(apd2 = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
            if(apd1 == apd2 || actualmol != mpsame) {
               if(actualmol != mpsame) {
                  apd2->picked = 0;
                  logprint("atom not in the same molecule");
                  update_logs();
                  glutPostRedisplay();
               }
               apd2 = NULL;
               return;
            }
            add_bond(apd1, apd2);
            apd1->picked = apd2->picked = 0;
            apd1 = apd2 = NULL;
            glutPostRedisplay();
         }
       break;
       case DEL_BOND:
         if(!(bp = (Bond *)my_pick(PICK_BOND, x, y))) return;
         rem_bond_from_list(bp);
         bp = NULL;
         glutPostRedisplay();
       break;
       case SINGLEBOND:
         if(!(bp = (Bond *)my_pick(PICK_BOND, x, y))) return;
         if(actualmol && !actualmol->ball_and_stick) Set_ball_and_stick();
         modify_bond_type(SINGLEBOND, bp);
         bp = NULL;
         glutPostRedisplay();
       break;
       case DOUBLEBOND:
         if(!(bp = (Bond *)my_pick(PICK_BOND, x, y))) return;
         if(actualmol && !actualmol->ball_and_stick) Set_ball_and_stick();
         modify_bond_type(DOUBLEBOND, bp);
         bp = NULL;
         glutPostRedisplay();
       break;
       case TRIPLEBOND:
         if(!(bp = (Bond *)my_pick(PICK_BOND, x, y))) return;
         if(actualmol && !actualmol->ball_and_stick) Set_ball_and_stick();
         modify_bond_type(TRIPLEBOND, bp);
         bp = NULL;
         glutPostRedisplay();
       break;
       case H_BOND:
         if(!(bp = (Bond *)my_pick(PICK_BOND, x, y))) return;
         modify_bond_type(H_BOND, bp);
         bp = NULL;
         glutPostRedisplay();
       break;
       case ATOMCOLOR:
         if(!(ap = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
         assign_color_to_atom(ap);
         ap = NULL;
         pick_action(DONE, 1, 1);
       break;
       case SURFCOLOR:
         if(!(sp = (Surface *)my_pick(PICK_SURF, x, y))) return;
         actualsurf = sp;
         assign_color_to_actualsurf();
         sp = NULL;
         pick_action(DONE, 1, 1);
       break;
       case GET_ATOMCOLOR:
         if(!(ap = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
         get_color_from_atom(ap);
         ap = NULL;
         pick_action(DONE, 1, 1);
       break;
       case GET_SURFCOLOR:
         if(!(sp = (Surface *)my_pick(PICK_SURF, x, y))) return;
         actualsurf = sp;
         get_color_from_actualsurf();
         sp = NULL;
         pick_action(DONE, 1, 1);
       break;
       case SURF_TEXTURE:
          if(!(sp = (Surface *)my_pick(PICK_SURF, x, y))) return;
          actualsurf = sp;
          if(assign_no_tex) {
             actualsurf->texture = NULL;
          }
          else {
             glutSetWindow(mainwin);
             set_texture_range();
             mod_texture((GLenum *)texprop_live);
             actualsurf->texture = actual_texture;
             actualsurf->texenv = (GLenum)texenv_live;
             actualsurf->textype = TEX_MAP;
             glutSetWindow(textureinterfwin);
             glutPostRedisplay();
          }
          glutSetWindow(mainwin);
          sp = NULL;
          pick_action(DONE, 1, 1);
       break;
       case SURF_REFLECT:
          if(!(sp = (Surface *)my_pick(PICK_SURF, x, y))) return;
          actualsurf = sp;
          if(assign_no_tex) {
             actualsurf->texture = NULL;
          }
          else {
             glutSetWindow(mainwin);
             set_texture_range();
             mod_texture((GLenum *)texprop_live);
             actualsurf->texture = actual_texture;
             actualsurf->texenv = (GLenum)texenv_live;
             actualsurf->textype = TEX_REFLECT;
             glutSetWindow(textureinterfwin);
             glutPostRedisplay();
          }
          glutSetWindow(mainwin);
          sp = NULL;
          pick_action(DONE, 1, 1);
       break;
       case SURF_PHONG:
          if(!(sp = (Surface *)my_pick(PICK_SURF, x, y))) return;
          actualsurf = sp;
          if(assign_no_tex) {
             actualsurf->texture = NULL;
          }
          else {
             glutSetWindow(mainwin);
             set_texture_range();
             mod_texture((GLenum *)texprop_live);
             actualsurf->texture = actual_texture;
             actualsurf->texenv = (GLenum)texenv_live;
             actualsurf->textype = TEX_PHONG;
             glutSetWindow(textureinterfwin);
             glutPostRedisplay();
          }
          glutSetWindow(mainwin);
          sp = NULL;
          pick_action(DONE, 1, 1);
       break;
       case ATOM_TEXTURE:
          if(!(ap = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
          if(assign_no_tex) element[ap->ord].texture = NULL;
          else {
             glutSetWindow(mainwin);
             mod_texture((GLenum *)texprop_live);
             element[ap->ord].texture = actual_texture;
             element[ap->ord].texenv = (GLenum)texenv_live;
             element[ap->ord].textype = TEX_MAP;
             glutSetWindow(textureinterfwin);
             glutPostRedisplay();
          }
          glutSetWindow(mainwin);
          ap = NULL;
          pick_action(DONE, 1, 1);
       break;
       case ATOM_REFLECT:
          if(!(ap = (AtoM *)my_pick(PICK_ATOM, x, y))) return;
          if(assign_no_tex) element[ap->ord].texture = NULL;
          else {
             glutSetWindow(mainwin);
             mod_texture((GLenum *)texprop_live);
             element[ap->ord].texture = actual_texture;
             element[ap->ord].texenv = (GLenum)texenv_live;
             element[ap->ord].textype = TEX_REFLECT;
             glutSetWindow(textureinterfwin);
             glutPostRedisplay();
          }
          glutSetWindow(mainwin);
          ap = NULL;
          pick_action(DONE, 1, 1);
       break;
       case DONE:
/* finish ADD DUMMY */
          if(divisor > 1){
              xcoor /= divisor;
              ycoor /= divisor;
              zcoor /= divisor;
              ap = add_atom(0, xcoor, ycoor, zcoor);
              updateAtomList(ap);
          }

/* finish CENTER */
          if(divisor_center){
             mpsame->centervec[0] = mvec[0]/divisor_center;
             mpsame->centervec[1] = mvec[1]/divisor_center;
             mpsame->centervec[2] = mvec[2]/divisor_center;
          }

/* finish MATCH */
          if(a1 && a2){
             if(fit_atoms(MIN(na1, na2), a1, a2, m, m1, m2)) {
                move_fittedmol(m1, m2, m);
             }
             else showinfobox("matching error");
          }

/* exit picking mode and reset variables */
          picking_mode = 0;
          for(mp=firstmol; mp; mp=mp->next){
             for(ap = mp->firstatom; ap; ap = ap->next) {
                ap->picked = 0;
             }
          }
          ap = apd1 = apd2 = apa1 = apa2 = apa3 = apt1 = apt2 = apt3 = apt4 = NULL;
          mpsame = NULL;
          xcoor = ycoor = zcoor = divisor = divisor_center = 0;
          free(a1); free(a2);
          a1 = a2 = NULL;
          m1 = m2 = NULL;
          na1 = na2 = 0;
          mvec[0] = mvec[1] = mvec[2] = 0.0;

          action_key = 0;
          pick_key = -1;
          mainglui_top->enable();
          glutMouseFunc(mouse);
          glutMotionFunc(motion);
          logprint("");
          logprint("Picking mode done");
          update_logs();
 
          glutPostRedisplay();
       break;
   }
   update_interface_flags();
}

static float getRadius(Mol *mp, AtoM *ap)
{
   if(mp->spacefill)      return element[ap->ord].rvdw;
   if(mp->ball_and_stick) return element[ap->ord].rvdw * vdW_factor;
   if(mp->sticks)         return 0.5*bond_factor;
   if(mp->wire_model)     return 0.1;
   return -1;
}


static void getWorldSpaceRay(long ix, long iy, float rayStart[3], float rayEnd[3])
{
   Matrix projmat, inversemat;
   float a[4], b[4];

   glMatrixMode(GL_PROJECTION);
   glGetFloatv(GL_PROJECTION_MATRIX, &projmat[0][0]);
   glMatrixMode(GL_MODELVIEW);

   invert4x4(projmat, inversemat);

   a[0] = (ix*2.0 + 1.0)/(float)xsize - 1.0;
   a[1] = (iy*2.0 + 1.0)/(float)ysize - 1.0;
   a[2] = -1.0;
   a[3] = 1.0;

   vecxmat(a, inversemat, b);

   rayStart[0] = b[0]/b[3];
   rayStart[1] = b[1]/b[3];
   rayStart[2] = b[2]/b[3];

   a[2] =  1.0;

   vecxmat(a, inversemat, b);

   rayEnd[0] = b[0]/b[3];
   rayEnd[1] = b[1]/b[3];
   rayEnd[2] = b[2]/b[3];
}


static void getObjectSpaceRay(float wsStart[3], float wsEnd[3],
                              float osStart[3], float osEnd[3])
{
   Matrix modelMat, inverseMat;
   float a[4], b[4];

   glGetFloatv(GL_MODELVIEW_MATRIX, &modelMat[0][0]);
   invert4x4(modelMat, inverseMat);

   a[0] = wsStart[0]; a[1] = wsStart[1]; a[2] = wsStart[2]; a[3] = 1.0;
   vecxmat(a, inverseMat, b);
   osStart[0] = b[0]; osStart[1] = b[1]; osStart[2] = b[2];

   a[0] = wsEnd[0]; a[1] = wsEnd[1]; a[2] = wsEnd[2]; a[3] = 1.0;
   vecxmat(a, inverseMat, b);
   osEnd[0] = b[0]; osEnd[1] = b[1]; osEnd[2] = b[2];
}


#define DOT(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])


static unsigned long pickAtom(Mol *mp, float start[3], float end[3], float *ratio)
{
   float d[3], dd;
   float newRatio, r;
   AtoM *ap, *candidate;

   d[0] = end[0] - start[0];
   d[1] = end[1] - start[1];
   d[2] = end[2] - start[2];

   dd = DOT(d, d);

   candidate = 0;
   newRatio = 1.0;

   for(ap = mp->firstatom; ap; ap=ap->next) {
      float q[3], qq, dq, determinant;
      float ratio1;

      if(mp->suppress_H && (ap->ord == 1)) continue;
      if(mp->residues && !ap->main) continue;

      q[0] = start[0] - ap->coord[0];
      q[1] = start[1] - ap->coord[1];
      q[2] = start[2] - ap->coord[2];

      dq = DOT(d, q);
      qq = DOT(q, q);

      r = getRadius(mp, ap);

      determinant = dq*dq - dd*(qq - r*r);
      if(determinant < 0) continue;

      ratio1 = (-dq - sqrt(determinant))/dd;
      if((ratio1 > 0.0) && (ratio1 < newRatio)) {
         candidate = ap;
         newRatio = ratio1;
      }
   }

   if(candidate) {
      *ratio = newRatio;
   }

   return (unsigned long)candidate;
}

static unsigned long pickBond(Mol *mp, float start[3], float end[3], float *ratio)
{
   Bond *candidate, *bp;
   float r, newRatio;
   float rayStart[4], rayEnd[4];
   extern float bond_factor;

   candidate = 0;
   newRatio = 1.0;

   if(mp->spacefill) return 0;
   if(mp->ball_and_stick || mp->sticks) r = 0.5*bond_factor;
   else r = 0.1;

   rayStart[0] = start[0]; rayStart[1] = start[1]; 
   rayStart[2] = start[2]; rayStart[3] =  1.0;
   rayEnd[0] = end[0]; rayEnd[1] = end[1]; rayEnd[2] = end[2]; rayEnd[3] = 1.0;

   for(bp = mp->firstbond; bp; bp = bp->next) {
      Matrix inverse, bondMatrix;
      float startP[4], endP[4];
      float p[3], cylOffset[3];
      int nBonds, i;

      if(mp->suppress_H){
         if(bp->a->ord == 1 || bp->b->ord == 1) continue;
      }
      if(mp->residues) {
         if(!bp->a->main || !bp->b->main) continue;
      }

      getBondMatrix(bp, bondMatrix);
      invert4x4(bondMatrix, inverse);
      vecxmat(rayStart, inverse, startP);
      vecxmat(rayEnd, inverse, endP);

      /* the cylinder is now around the positive z-axis and goes from 0 to 1 */
      p[0] = startP[0]; p[1] = startP[1]; p[2] = 0;

      if(mp->ball_and_stick && mp->multiple_bonds){
         switch(bp->type) {
            case SINGLEBOND: nBonds = 1; cylOffset[0] = 0; break;
            case DOUBLEBOND: nBonds = 2; cylOffset[0] = -.13; cylOffset[1] = .13; break;
            case TRIPLEBOND: nBonds = 3;
                cylOffset[0] = -.18; cylOffset[1] = 0; cylOffset[2] = .18; break;
            case H_BOND: nBonds = 1; cylOffset[0] = 0; break;
         }
      }
      else {
         nBonds = 1; cylOffset[0] = 0;
      }

      for(i=0; i<nBonds; i++) {
         float q[3], d[3], startPrime[3], endPrime[3], dirPrime[3];
         float dd, dq, qq, determinant;
         float alpha1, test1;

         /* q = cylinderCenter - p */
         q[0] = cylOffset[i] - p[0]; q[1] = -p[1]; q[2] = -p[2];
         d[0] = endP[0] - startP[0]; d[1] = endP[1] - startP[1]; d[2] = 0;
         dd = DOT(d, d);
         if(dd == 0.0) continue;   /* the ray is parallel to the bond */

         dq = DOT(d, q);
         qq = DOT(q, q);
         determinant = dq*dq - dd*(qq - r*r);
         if(determinant < 0) continue;
            /* the ray passes too far from the bondline */

         alpha1 = (dq - sqrt(determinant))/dd;

         startPrime[0] = startP[0]; startPrime[1] = startP[1];
         startPrime[2] = startP[2];
         endPrime[0] = endP[0]; endPrime[1] = endP[1]; endPrime[2] = endP[2];
         dirPrime[0] = endPrime[0] - startPrime[0];
         dirPrime[1] = endPrime[1] - startPrime[1];
         dirPrime[2] = endPrime[2] - startPrime[2];

         test1 = startPrime[2] + alpha1*dirPrime[2];
         if((test1 < 0.0) || (test1 > 1.0)) continue; 
            /* the ray passes outside the ends of the bond */

         if((alpha1 > 0.0) && (alpha1 < newRatio)) {
            candidate = bp;
            newRatio = alpha1;
         }
      }
   }

   if(candidate) {
      *ratio = newRatio;
   }

   return (unsigned long)candidate;
}

void getClipplaneEquation(float eq[4], Surface *sp)
{
   Matrix m;
   float a[3], b[3];

   build_rotmatrix(m, sp->clipplane_rvec);

/* the transformed starting point of the plane-normal (originally 0, 0, 0)*/
   a[0] = sp->clipplane_tvec[0];
   a[1] = sp->clipplane_tvec[1];
   a[2] = sp->clipplane_tvec[2];

/* the transformed ending point of the plane-normal (originally 0, 0, -1)*/
   b[0] = a[0] - m[2][0];
   b[1] = a[1] - m[2][1];
   b[2] = a[2] - m[2][2];

/* the transformed plane-normal */
   eq[0] = b[0] - a[0];
   eq[1] = b[1] - a[1];
   eq[2] = b[2] - a[2];

/* the offset by backsubstitution of point a */
   eq[3] = -eq[0]*a[0] - eq[1]*a[1] - eq[2]*a[2];
}

static unsigned long pickSurface(Mol *mp, float p[3],
       float end[3], float *ratio, short key, Triangle **tri, float *value)
{
   Surface *candidate, *sp;
   Triangle *hit, *tp;
   float newRatio, alpha, beta, gamma;
   float *a, ab[3], ac[3], d[3], q[3], result[3];
   float *fp;
   float plane[4];
   float mat[3][3], inverse[3][3];
   int i;

   candidate = 0;
   hit = 0;
   newRatio = 1.0;

   for(sp=mp->firstsurf; sp; sp=sp->next) {
      if((key == PICK_TRI) && !sp->val) continue;
      if(sp->surf_clip) getClipplaneEquation(plane, sp);

      for(i=0, tp=sp->tri; i<sp->ntri; i++, tp++) {
         a  = sp->dot[tp->p1].v;
         ab[0] = sp->dot[tp->p2].v[0] - a[0];
         ab[1] = sp->dot[tp->p2].v[1] - a[1];
         ab[2] = sp->dot[tp->p2].v[2] - a[2];
         ac[0] = sp->dot[tp->p3].v[0] - a[0];
         ac[1] = sp->dot[tp->p3].v[1] - a[1];
         ac[2] = sp->dot[tp->p3].v[2] - a[2];
         /* triangle : a + beta*ab + gamma*ac
            beta, gamma in [0,1], beta + gamma <= 1 */

         d[0] = end[0] - p[0];
         d[1] = end[1] - p[1];
         d[2] = end[2] - p[2];
         q[0] = a[0] - p[0];
         q[1] = a[1] - p[1];
         q[2] = a[2] - p[2];

         fp=mat[0];
         *fp++ =  d[0];  *fp++ =  d[1];  *fp++ = d[2];
         *fp++ = ab[0];  *fp++ = ab[1];  *fp++ = ab[2];
         *fp++ = ac[0];  *fp++ = ac[1];  *fp++ = ac[2];

         invert3x3(mat, inverse);
         vec3mat(q, inverse, result);

         alpha =  result[0];
         beta  = -result[1];
         gamma = -result[2];

         if((beta < 0.0) || (beta > 1.0)) continue;
         if((gamma < 0.0) || (gamma > 1.0)) continue;
         if((beta+gamma) > 1.0) continue;
         if((alpha > 0.0) && (alpha < newRatio)) {

/*
// if cuttingplane: now check whether the point (p + alpha*d) 
// lies on the visible side of the cutting plane
*/
            if(sp->surf_clip) {
               float test[3];
               test[0] = p[0] + alpha*d[0];
               test[1] = p[1] + alpha*d[1];
               test[2] = p[2] + alpha*d[2];

               if((test[0]*plane[0] + test[1]*plane[1] +
                   test[2]*plane[2] + plane[3]) < 0) continue;
            }

/* and perhaps: get texture coord and see whether the point is visible...*/

            candidate = sp;
            hit = tp;
            newRatio  = alpha;
            if(key == PICK_TRI) {
               *value = (1.0 - beta - gamma)*sp->val[tp->p1] +
                        beta*sp->val[tp->p2] + gamma*sp->val[tp->p3];
            }
         }
      }
   }

   if(candidate) {
      *ratio = newRatio;
      *tri = hit;
   }

   return (unsigned long)candidate;
}

void highLight(unsigned long objName, short key, Mol *mp,
                      Triangle *trip)
{
   AtoM *ap;
   Bond *bp;
   Surface *sp;

   glPushMatrix();
   global_move();
   individual_move(mp);
   glDrawBuffer(GL_FRONT);

   switch(key){
      case PICK_ATOM :
         ap = (AtoM*)objName;
         ap->picked = 1;
         draw_atom(mp, ap); 
         break;

      case PICK_BOND : 
         bp = (Bond*)objName;
         bp->picked = 1; 
         draw_color_stick(mp, bp); 
         bp->picked = 0; 
         break;

/* to be fixed
      case PICK_BOX  : 
         break;

      case PICK_TRI  : 
*/
      case PICK_MOL  :
         for(ap=mp->firstatom; ap; ap=ap->next) {
            ap->picked = 1;
            draw_atom(mp, ap);
            ap->picked = 0;
         }
         if(mp->firstsurf) {
            glCallList(pick_mat);
            for(sp=mp->firstsurf; sp; sp=sp->next) {
                draw_highlight_surf(sp);
            }
         }
         break;
      case PICK_SURF : 
         sp = (Surface*)objName;
         glCallList(pick_mat);
         draw_highlight_surf(sp);
         break;

/* to be fixed
      case PICK_SEGMENT :
         break;
*/
   }

   glFlush();
   glDrawBuffer(GL_BACK);
   glPopMatrix();
}

unsigned long my_pick(int key, int x, int y)
{
   Mol *mp;
   Triangle *tri = NULL;
   float rIn[3], rOut[3], nearest, value;
   unsigned long objname;
   short textureflag;
   int xpixel, ypixel;

   picking_mode = 1;

   textureflag = bit.texture;
   bit.texture = 0;

   xpixel =  x ;
   ypixel =  ysize - 1 - y ;

   init_persp();

   getWorldSpaceRay(xpixel, ypixel, rIn, rOut);

   objname = 0;
   nearest = 1.0;

   glPushMatrix();
   global_move();

   for(mp=firstmol; mp; mp=mp->next){
      float start[3], end[3], ratio;
      unsigned long newObject;

      glPushMatrix();
      individual_move(mp);
      getObjectSpaceRay(rIn, rOut, start, end);
      newObject = 0; ratio = 1.0;

      switch(key){
         case PICK_ATOM : newObject = pickAtom(mp, start, end, &ratio);
                          break;
         case PICK_BOND : newObject = pickBond(mp, start, end, &ratio);
                          break;
         case PICK_MOL  : if(!(newObject = pickAtom(mp, start, end, &ratio))) {
                           if(!(newObject = pickBond(mp, start, end, &ratio))) {
                             newObject = pickSurface(mp, start, end,
                                          &ratio, key, &tri, &value);
                           }
                          }
                          break;
/* to be fixed
         case PICK_BOX  : newObject = pickBox(mp, start, end, &ratio);
                          break;
         case PICK_TRI  : 
*/
         case PICK_SURF : newObject = pickSurface(mp, start, end,
                                          &ratio, key, &tri, &value);
                          break;
/* to be fixed
         case PICK_SEGMENT :break;
*/
      }
      glPopMatrix();

      if(newObject && (ratio < nearest)) {
         objname = newObject;
         nearest = ratio;
         actualmol = mp;
         if(key == PICK_MOL) objname = (unsigned long)actualmol;
      }
   }

   glPopMatrix();

/*
   if(key == PICK_TRI) {
      char str[80];
      sprintf(str, " picked triangle : %f", value);
      printf("%s\n", str);
      screenprint(str);
   }
*/

   if(objname) {
      highLight(objname, key, actualmol, tri);
   }

//   picking_mode = 0;
   bit.texture = textureflag;

   return objname;
}


void pick_surface(int key)
{
   Surface *surf;
   if(!(surf = actualmol->firstsurf)){
      action_key = 0;
      mainglui_top->enable();
      showinfobox("Load or generate a molecular surface first!");
      return;
   }
   pick_key = key;
   if(surf->next == NULL) { 
     switch (key) {
         case GRID_VALUES:
            gridValues(surf);
         break;
         case MEP:
            mainglui_top->enable();
            glutSetWindow(mainwin);
            mep_dot_surface(surf);
         break;
     } 
     action_key = 0;
   }
   else {
      glutSetWindow(mainwin);
      showinfobox("Pick a surface!");
      enter_pick(key, "Pick surface");
   }
}

void center(void)
{
   enter_pick(CENTER, "Pick atoms");
}


void enter_pick(int value, char *string)
{
   picking_mode = 1;
   pick_key = value;
   glutMouseFunc(picking_mouse);
   glutMotionFunc(picking_motion);
   logprint("");
   logprint("Entering picking mode");
   logprint("   return with middle mouse");
   logprint("   or \"Done picking\"");
   logprint(string);
   update_logs();
}
