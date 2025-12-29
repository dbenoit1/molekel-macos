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
#include "geometry.h"
#include "maininterf.h"
#include "utils.h"
#include "manip.h"

void printAtomData(AtoM *ap);

float valence_angle(float *a, float *b, float *c)
{
   float v[3], w[3], dot_product, lv, lw;
   register int i;

   for(i=0; i<3; i++){
      v[i] = a[i] - b[i];
      w[i] = c[i] - b[i];
   }
   dot_product = v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
   lv = fsqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
   lw = fsqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

   return facos(dot_product/lv/lw)*180./M_PI;
}

float dihedral_angle(float *a, float *b, float *c, float *d)
{
   float u[3], v[3], w[3], n1[3], n2[3], angle;
   static float o[3] = {0, 0, 0};
   register short i;

   for(i=0; i<3; i++){
      u[i] = b[i] - a[i];
      v[i] = c[i] - b[i];
      w[i] = d[i] - c[i];
   }
   cross_product(u, v, n1);
   cross_product(v, w, n2);

   angle = valence_angle(n1, o, n2);

   if(triple_product(n1, n2, v) < 0) angle *= -1.0;

   return angle;
}


void measure_atom(AtoM *ap)
{
   printAtomData(ap);
   ap->picked = 0;
}

void measure_distance(AtoM *ap1, AtoM *ap2)
{
   float value;
   char text[40];

   add_distance(ap1, ap2);
   value = dist(ap1->coord, ap2->coord);
   sprintf(text, "distance %2s%2d - %2s%2d = %8.3f",
      element[ap1->ord].symbol, ap1->name,
      element[ap2->ord].symbol, ap2->name, 
       value);
   printf("%s\n", text);
   sprintf(text, "distance %2s%2d - %2s%2d =",
      element[ap1->ord].symbol, ap1->name,
      element[ap2->ord].symbol, ap2->name);
   logprint(text);
   sprintf(text, "           %8.3f", value);
   logprint(text);
   update_logs();

   ap1->picked = ap2->picked = 0;
   glutPostRedisplay();
}


void measure_intermolecular_distance(AtoM *ap1, Mol *mp1, AtoM *ap2, Mol *mp2)
{
   add_intermolecular_distance(ap1, mp1, ap2, mp2);
   logprint("intermolecular distance");
   update_logs();

   ap1->picked = ap2->picked = 0;
   glutPostRedisplay();
}


void measure_alldist(void)
{
   Bond *bp;
   float value;
   char text[40];

   for(bp = actualmol->firstbond; bp; bp = bp->next) {
      add_distance(bp->a, bp->b);
      value = dist(bp->a->coord,  bp->b->coord);
      sprintf(text, "distance %2s%2d - %2s%2d = %8.3f",
            element[bp->a->ord].symbol, bp->a->name,
            element[bp->b->ord].symbol, bp->b->name,
              value);
      printf("%s\n", text);
      sprintf(text, "distance %2s%2d - %2s%2d =",
            element[bp->a->ord].symbol, bp->a->name,
            element[bp->b->ord].symbol, bp->b->name);
      logprint(text);
      sprintf(text, "           %8.3f", value);
      logprint(text);
      update_logs();
   }
   glutPostRedisplay();
}


void measure_valence(AtoM *ap1, AtoM *ap2, AtoM *ap3)
{
   float value;
   char text[40];

   add_angle(ap1, ap2, ap3);
   value = valence_angle(ap1->coord, ap2->coord, ap3->coord);
   sprintf(text, "angle %2s%2d - %2s%2d - %2s%2d = %8.3f",
      element[ap1->ord].symbol, ap1->name,
      element[ap2->ord].symbol, ap2->name, 
      element[ap3->ord].symbol, ap3->name, value);
   printf("%s\n", text);
   sprintf(text, "angle %2s%2d - %2s%2d - %2s%2d =",
      element[ap1->ord].symbol, ap1->name,
      element[ap2->ord].symbol, ap2->name, 
      element[ap3->ord].symbol, ap3->name);
   logprint(text);
   sprintf(text, "         %8.3f", value);
   logprint(text);
   update_logs();

   ap1->picked = ap2->picked = ap3->picked = 0;
   glutPostRedisplay();
}


void measure_dihedral(AtoM *ap1, AtoM *ap2, AtoM *ap3, AtoM *ap4)
{
   float value;
   char text[50];

   add_torsion(ap1, ap2, ap3, ap4);
   value = dihedral_angle(ap1->coord, ap2->coord,
      ap3->coord, ap4->coord);
   sprintf(text,
       "dihedral %2s%2d - %2s%2d - %2s%2d - %2s%2d = %8.3f",
       element[ap1->ord].symbol, ap1->name,
       element[ap2->ord].symbol, ap2->name, 
       element[ap3->ord].symbol, ap3->name,
       element[ap4->ord].symbol, ap4->name, value);
   printf("%s\n", text);
   logprint("dihedral");
   sprintf(text,
       "  %2s%2d - %2s%2d - %2s%2d - %2s%2d =",
       element[ap1->ord].symbol, ap1->name,
       element[ap2->ord].symbol, ap2->name, 
       element[ap3->ord].symbol, ap3->name,
       element[ap4->ord].symbol, ap4->name);
   logprint(text);
   sprintf(text, "             %8.3f", value);
   logprint(text);
   update_logs();
   ap1->picked = ap2->picked = ap3->picked = ap4->picked = 0;
   glutPostRedisplay();
}

void printGaussianShell(Shell *sp)
{
   Gauss *gauss;
   static char *shell[] = {"",
                           "S",  /* 1 */
                           "",   /* 2 */
                           "P",  /* 3 */
                           "SP", /* 4 */
                           "D",  /* 5 */
                           "D",  /* 6 */
                           "F",  /* 7 */
                           "",   /* 8 */
                           "",   /* 9 */
                           "F"}; /* 10 */

   printf("          %s Shell, scale %f\n",
                        shell[sp->n_base], sp->scale_factor);

   for(gauss = sp->firstgauss; gauss; gauss = gauss->next){
      printf("                         %12.5f %8.4f %8.4f\n",
              gauss->exponent, gauss->coeff, gauss->coeff2);
   }
}



void printValenceShell(Valence_shell *vp)
{
   printf("           ns, np, nd = %d %d %d\n", vp->ns, vp->np, vp->nd);
   printf("           exp(s, p, d) = %f %f %f\n", vp->exps, vp->expp, vp->expd);
   printf("           coul(s, p, d) = %f %f %f\n", vp->couls, vp->coulp, vp->could);
}



void printSlaterShell(Slater *sp)
{
   printf("           %d%s, exp = %f\n", sp->n, sp->type, sp->exponent);
}


void printAtomData(AtoM *ap)
{
   Shell *sp;
   Slater *sl;
   char text[40];

   sprintf(text, "picked atom : %s %d\n", element[ap->ord].symbol, ap->name);
   logprint(text);
   update_logs();

   printf(" Picked Atom : %s %d\n", element[ap->ord].symbol, ap->name);
   printf("      coordinates          %f %f %f\n",
                          ap->coord[0], ap->coord[1], ap->coord[2]);
   printf("      van der Waals radius %f\n", element[ap->ord].rvdw);
   printf("      covalent radius      %f\n", element[ap->ord].rvdw);

   if(ap->firstshell) {
      printf("         gaussian type basis set:\n");
      for(sp = ap->firstshell; sp; sp = sp->next) {
         printGaussianShell(sp);
      }
   }

   if(ap->valence) {
      printf("         semiempirical basis set:\n");
      printValenceShell(ap->valence);
   }

   if(ap->firstslater) {
      printf("         slater type basis set:\n");
      for(sl = ap->firstslater; sl; sl = sl->next) {
         printSlaterShell(sl);
      }
   }

   printf("\n");
}


void removeAtomType(int ordinal)
{
   AtomType *at, *previous;

   for(at=actualmol->firstAtomType; at; at=at->next) {
      if(at->ord == ordinal){
         if(at == actualmol->firstAtomType)
            actualmol->firstAtomType = actualmol->firstAtomType->next;
         else
            previous->next = at->next;

         if(at->atomList) free(at->atomList);
         free(at);
         return;
      }
      previous = at;
   }
}

void remove_dummies(void)
{
   AtoM *ap, *prev;
   Bond *bp, *dp;

   if(!actualmol) return;

   for(bp=actualmol->firstbond; bp;){
      if(bp->a->ord == 0 || bp->b->ord == 0) {
         dp=bp;
         bp=bp->next;
         rem_bond_from_list(dp);
      }
      else bp=bp->next;
   }

   removeAtomType(0);

   for(ap=actualmol->firstatom; ap; ap=ap->next){
      if(ap->ord == 0){
         if(ap == actualmol->firstatom) actualmol->firstatom = ap->next;
         else {
            prev->next = ap->next;
         }
         if(!ap->next) actualmol->lastatom = prev;
         free(ap);
         actualmol->natoms--;
         ap = prev;
      }
      else {
         prev = ap;
      }

   }

   glutPostRedisplay();
}


#define TOLERANCE   0.001
#define MAXIT     400
int fit_atoms(int nat, AtoM **a, AtoM **b, Matrix rt, Mol *molA, Mol *molB)
/* from MOLCAD - pff aug '93 */
{
   int i, j, k, ict, ix, iy, iz, iflag;
   float xt, yt, zt, xx;
   float aa[3][3], da[3], db[3], cgm[3], cgf[3];
   float gam, gam2, sig, sig2, bb, cc, sg;

   if(!nat) return 0;

   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) aa[i][j] = rt[i][j] = 0;
      rt[i][i] = 1.;
   }

   get_center(nat, cgm, a);
   get_center(nat, cgf, b);

   molA->centervec[0] = -cgm[0];
   molA->centervec[1] = -cgm[1];
   molA->centervec[2] = -cgm[2];
   molB->centervec[0] = -cgf[0];
   molB->centervec[1] = -cgf[1];
   molB->centervec[2] = -cgf[2];

   for(i=0; i<nat; i++){
      da[0] = a[i]->coord[0] - cgm[0];
      da[1] = a[i]->coord[1] - cgm[1];
      da[2] = a[i]->coord[2] - cgm[2];
      db[0] = b[i]->coord[0] - cgf[0];
      db[1] = b[i]->coord[1] - cgf[1];
      db[2] = b[i]->coord[2] - cgf[2];
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

/** print the distances between the matched atoms *
   {
      int i;
      float u[4], v[4];
      extern void vecxmat(float a[4], Matrix m, float b[4]);

      printf("\n cgm = %7.3f, %7.3f, %7.3f\n cgf = %7.3f, %7.3f, %7.3f\n",
             cgm[0], cgm[1], cgm[2], cgf[0], cgf[1], cgf[2]);
      for(i=0; i<nat; i++){
         u[0] = a[i]->coord[0];
         u[1] = a[i]->coord[1];
         u[2] = a[i]->coord[2];
         u[3] = 1.0;
         vecxmat(u, rt, v);
         printf("dist a%-3d- b%-3d : %6.3f   (%7.3f,%7.3f,%7.3f)\n",
            a[i]->name, b[i]->name, dist(v, b[i]->coord),
            v[0]-b[i]->coord[0], v[1]-b[i]->coord[1], v[2]-b[i]->coord[2]);
      }
   }
**********************/
   return 1;
}


void move_fittedmol(Mol *mo, Mol *mp, Matrix  m)
{
   float e[4];

   /* the rotataion */
   get_rotation(m, e);

   mo->rvec[0] = e[0];
   mo->rvec[1] = e[1];
   mo->rvec[2] = e[2];
   mo->rvec[3] = e[3];
   mp->rvec[0] = mp->rvec[1] = mp->rvec[2] = mp->rvec[3] = 0;

   /* the translation */

   mo->tvec[0] = mo->tvec[1] = mo->tvec[2] = 0;
   mp->tvec[0] = mp->tvec[1] = mp->tvec[2] = 0;

   bit.manip_all = 1;
   mainglui_top->sync_live(); 
}
