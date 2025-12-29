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


/* utility functions */

#include "main.h"
#include "molekel.h"
#include "glutwin.h"
#include "constant.h"
#include "utils.h"
#include "general.h"
#include "drawing.h"
#include "trackball.h"
#include "manip.h"
#include "maininterf.h"
#include "surfaceinterf.h"
#include "box.h"
#include "macu.h"
#include "geometry.h"

#define DET3x3(a,b,c) (a[0]*b[1]*c[2]+b[0]*c[1]*a[2]+c[0]*a[1]*b[2]\
                      -a[2]*b[1]*c[0]-b[2]*c[1]*a[0]-c[2]*a[1]*b[0])
#define DIST(a,b) sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+\
                       (a[2]-b[2])*(a[2]-b[2]))

#define FULLVALENCY(a,b) (\
   (((a)==1)&&(b))\
   ||(((a)==6)&&((b)>=4))\
)

static float inv, cubeSize;
static Atomcube ***acube;
static int nx, ny, nz;
/*----------------------------------------*/

void vecxmat(float a[4], Matrix m, float b[4])
{
   register int i;

   for(i=0; i<4; i++){
      b[i] = a[0]*m[0][i] + a[1]*m[1][i] +
             a[2]*m[2][i] + a[3]*m[3][i];
   }
}


float determinant(float *a, float *b, float *c)
{
   return (  a[0]*b[1]*c[2] + b[0]*c[1]*a[2] + c[0]*a[1]*b[2]
           - a[2]*b[1]*c[0] - b[2]*c[1]*a[0] - c[2]*a[1]*b[0] );
}


void invert4x4(Matrix m, Matrix n)
{
   double det, det_1;
   short i, j;

   det = DET3x3(m[0], m[1], m[2]) * m[3][3] -
         DET3x3(m[0], m[1], m[3]) * m[2][3] +
         DET3x3(m[0], m[2], m[3]) * m[1][3] -
         DET3x3(m[1], m[2], m[3]) * m[0][3];
   if(det == 0){
      printf("singular matrix!\n");
      return;
   }
   det_1 = 1.0/det;

   for(i=0; i<4; i++){
      for(j=0; j<4; j++){
         n[i][j] = det_1 * adjunct4x4(m, i, j);
      }
   }
}

double adjunct4x4(Matrix m, short k, short l)
{
   float a[3][3];
   register short i, j, u, v;
   double x;

   for(i=u=0; i<4; i++){
      if(i == l) continue;
      for(j=v=0; j<4; j++){
         if(j == k) continue;
         a[u][v] = m[i][j];
         v++;
      }
      u++;
   }

   x = DET3x3(a[0], a[1], a[2]);

   return (k+l)%2 ? -x: x;
}

void set_scale_factor(void)
{
   Matrix m;

   glGetFloatv(GL_MODELVIEW_MATRIX, &m[0][0]);
   scale_factor = determinant(m[0], m[1], m[2]);
   scale_factor = pow(scale_factor, 0.3333333333333);
}

float dot_product(float *a, float *b)
{
   return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

void cross_product(float *a, float *b, float *c)
/* taken from connect.c */
{
   c[0] = a[1]*b[2] - a[2]*b[1];
   c[1] = a[2]*b[0] - a[0]*b[2];
   c[2] = a[0]*b[1] - a[1]*b[0];
}

float triple_product(float *a, float *b, float *c)
{
   float ab[3];

   cross_product(a, b, ab);
   return dot_product(ab, c);
}

void transpose(Matrix m)
{
   register int i, j;
   Matrix temp;

   memcpy(temp, m, 16*sizeof(float));
   for(i=0; i<4; i++){
      for(j=0; j<4; j++){
         m[i][j] = temp[j][i];
      }
   }
}

float length(float *a)
{
   return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}


void get_center(int n, float *c, AtoM **ap)
{
   int i;
   float x, y, z;

   if(!n) return;
   x = y = z = 0;

   for(i=0; i<n; i++){
      x += ap[i]->coord[0];
      y += ap[i]->coord[1];
      z += ap[i]->coord[2];
   }

   c[0] = x/n;   c[1] = y/n;   c[2] = z/n;
}


void clipplane_on(Surface *sp)
{
   Matrix m, n;

   glPushMatrix();
   if(clipplane_moving){
      glGetFloatv(GL_MODELVIEW_MATRIX, &n[0][0]);
      glLoadIdentity();
   }
   glTranslatef(sp->clipplane_tvec[0], sp->clipplane_tvec[1],
             sp->clipplane_tvec[2]);
   build_rotmatrix(m, sp->clipplane_rvec);
   glMultMatrixf(&m[0][0]);
   if(clipplane_moving) glMultMatrixf(&n[0][0]);
   glClipPlane(GL_CLIP_PLANE0, clipping_plane);
   glEnable(GL_CLIP_PLANE0);
   glPopMatrix();
}


void clipplane_off(void)
{
   glDisable(GL_CLIP_PLANE0);
}


void set_clipplane_matrix(Surface *sp, Mol *mp, float r[4], float t[4])
{
Matrix rotMat, totalMat, invMat;
float rotVec[4], transVec[4];

         glPushMatrix();
         glLoadIdentity();
         individual_move(mp);
         glTranslatef(sp->clipplane_tvec[0],\
            sp->clipplane_tvec[1], sp->clipplane_tvec[2]);
         build_rotmatrix(rotMat, sp->clipplane_rvec);
         glMultMatrixf(&rotMat[0][0]);
         glGetFloatv(GL_MODELVIEW_MATRIX, &totalMat[0][0]);

         glPopMatrix();

         get_translation(totalMat, transVec);
         get_rotation(totalMat, rotVec);

         add_eulers(r, rotVec, rotVec);
         vadd(t, transVec, transVec);
         transVec[3] = 0;
 
         glPushMatrix();
         glLoadIdentity();
         individual_move(mp);
         glGetFloatv(GL_MODELVIEW_MATRIX, &totalMat[0][0]);
         invert4x4(totalMat, invMat);
         glLoadMatrixf(&invMat[0][0]);
         glTranslatef(transVec[0], transVec[1], transVec[2]);
         build_rotmatrix(rotMat, rotVec);
         glMultMatrixf(&rotMat[0][0]);
         glGetFloatv(GL_MODELVIEW_MATRIX, &totalMat[0][0]);
         glPopMatrix();

         get_translation(totalMat, sp->clipplane_tvec);
         get_rotation(totalMat, sp->clipplane_rvec);

}

double adjunct3x3(float m[3][3], int k, int l)
{
   float a[2][2];
   register int i, j, u, v;
   double x;

   for(i=u=0; i<3; i++){
      if(i == l) continue;
      for(j=v=0; j<3; j++){
         if(j == k) continue;
         a[u][v] = m[i][j];
         v++;
      }
      u++;
   }

   x = a[0][0]*a[1][1] - a[1][0]*a[0][1];

   return (k+l)%2 ? -x: x;
}

void invert3x3(float a[3][3], float b[3][3])
{
   double det, det_1;
   int i, j;

   det = a[0][0]*a[1][1]*a[2][2] + a[1][0]*a[2][1]*a[0][2] + 
         a[2][0]*a[0][1]*a[1][2] - a[0][2]*a[1][1]*a[2][0] -
         a[1][2]*a[2][1]*a[0][0] - a[2][2]*a[0][1]*a[1][0];

   if(det == 0) {
      printf("singular matrix!\n");
      return;
   }
   det_1 = 1.0/det;

   for(i=0; i<3; i++){
      for(j=0; j<3; j++){
         b[i][j] = det_1 * adjunct3x3(a, i, j);
      }
   }
}


void vec3mat(float a[3], float m[3][3], float b[3])
{
   register int i;

   for(i=0; i<3; i++){
      b[i] = a[0]*m[0][i] + a[1]*m[1][i] + a[2]*m[2][i];
   }
}


/* 
 * test with rotation panel -> does not work properly and to slow
 *
void set_clipplane_matrix(Mol *mp)
{
Matrix clip_rot_mat, rotMat, totalMat, invMat;
float clip_rot_vec[4], clip_rot_vec_cp[4], rotVec[4], transVec[4];
static float clip_rot_vec_last[] = { 0.0, 0.0, 0.0, 1.0};

         rot->get_float_array_val(&clip_rot_mat[0][0]);
         get_rotation(clip_rot_mat, clip_rot_vec);
         get_rotation(clip_rot_mat, clip_rot_vec_cp);
printf("%f %f %f %f\n", clip_rot_mat[0][0], clip_rot_mat[1][0], clip_rot_mat[2][0], clip_rot_mat[3][0]);
printf("%f %f %f %f\n", clip_rot_mat[0][1], clip_rot_mat[1][1], clip_rot_mat[2][1], clip_rot_mat[3][1]);
printf("%f %f %f %f\n", clip_rot_mat[0][2], clip_rot_mat[1][2], clip_rot_mat[2][2], clip_rot_mat[3][2]);
printf("%f %f %f %f\n\n", clip_rot_mat[0][3], clip_rot_mat[1][3], clip_rot_mat[2][3], clip_rot_mat[3][3]);
printf("%f %f %f %f\n\n", clip_rot_vec[0], clip_rot_vec[1], clip_rot_vec[2], clip_rot_vec[3]);
         add_eulers(clip_rot_vec_last, clip_rot_vec, clip_rot_vec);

         glPushMatrix();
         glLoadIdentity();
         individual_move(mp);
         glTranslatef(mp->clipplane_tvec[0],\
            mp->clipplane_tvec[1], mp->clipplane_tvec[2]);
         build_rotmatrix(rotMat, mp->clipplane_rvec);
         glMultMatrixf(&rotMat[0][0]);
         glGetFloatv(GL_MODELVIEW_MATRIX, &totalMat[0][0]);

         glPopMatrix();

//         get_translation(totalMat, transVec);
         get_rotation(totalMat, rotVec);

         add_eulers(clip_rot_vec, rotVec, rotVec);
//         vadd(t, transVec, transVec);
//         transVec[3] = 0;
 
         glPushMatrix();
         glLoadIdentity();
         individual_move(mp);
         glGetFloatv(GL_MODELVIEW_MATRIX, &totalMat[0][0]);
         invert4x4(totalMat, invMat);
         glLoadMatrixf(&invMat[0][0]);
//         glTranslatef(transVec[0], transVec[1], transVec[2]);
         build_rotmatrix(rotMat, rotVec);
         glMultMatrixf(&rotMat[0][0]);
         glGetFloatv(GL_MODELVIEW_MATRIX, &totalMat[0][0]);
         glPopMatrix();

//         get_translation(totalMat, mp->clipplane_tvec);
         get_rotation(totalMat, mp->clipplane_rvec);
         invert_euler(clip_rot_vec_cp, clip_rot_vec_last);

         glutSetWindow(mainwin);
         glutPostRedisplay();
}

 */


int get_ordinal(char *sym)
{
   int num, len;
   register short i;
   char line[100];

   if(isdigit(sym[0]) || sym[0] == ' ' || sym[0] == '.') sym++;
   len = strlen(sym);
   if(len > 1 && (isdigit(sym[1]) || sym[1] == '.')) len = 1;
   if(len > 2) len = 2;


   num = -1;
   for(i=0; i<NLIST; i++){
      if(len != strlen(element[i].symbol)) continue;
      if(!strncasecmp(sym, element[i].symbol, len)){
         num = i;
         break;
      }
   }
   if(num < 0){
      sprintf(line, "Unknown atom type %s,", sym);
      logprint(line);
      sprintf(line, "   set as hydrogen!", sym);
      logprint(line);
      num = 1;
   }

   return num;
}





void get_conect(char *line)
{
   register AtoM *i;
   AtoM *a, *b;
   int na, nb;

   (void)sscanf(line+6, "%d %d", &na, &nb);

   a = b = NULL;
   for(i = actualmol->firstatom; i; i = i->next){
      if(!a){
         if(i->name == na) a = i;
      }
      else {
         if(i->name == nb) { b = i; break;}
      }
   }

   if(a && b) add_bond(a, b);
}


/************** bonds ******************/

float getMaxCoval(void)
{
   AtomType *at;
   float maxCoval, tmp;

   maxCoval = 0;
   for(at = actualmol->firstAtomType; at; at = at->next) {
      tmp = element[at->ord].coval;
      if(tmp > maxCoval) maxCoval = tmp;
   }
   return maxCoval;
}



void create_bonds(void)
{
   AtoM *ap;

   for(ap = actualmol->firstatom; ap; ap = ap->next) updateAtomList(ap);
   
   createBonds(0);
   create_Hbonds();
}

void create_Hbonds(void)
{
   Bond *bp, *tbp;

   if(actualmol->firstbond) {
      for(bp = actualmol->firstbond; bp; ){
         tbp = bp->next;
         if(bp->type == H_BOND) rem_bond_from_list(bp);
         bp = tbp;
      }
   }
   
   createBonds(H_BOND);
}

void createBonds(int key)
{

   if(key == H_BOND) cubeSize = max_h_bond;
   else cubeSize = 2.0*getMaxCoval();
   if(!cubeSize) return;

   inv = 1./cubeSize;

   create_box();

   box.x1 += 1.4; box.y1 += 1.4; box.z1 += 1.4;
   box.x2 -= 1.4; box.y2 -= 1.4; box.z2 -= 1.4;

   nx = (int)fceil((box.x2 - box.x1)*inv);
   ny = (int)fceil((box.y2 - box.y1)*inv);
   nz = (int)fceil((box.z2 - box.z1)*inv);

   if((acube = (Atomcube***) alloc_3D(nx, ny, nz, sizeof(Atomcube))) == NULL){
      logprint("too many cubes, retrying");
      inv = 0.5/cubeSize;
      nx = (int)fceil((box.x2 - box.x1)*inv);
      ny = (int)fceil((box.y2 - box.y1)*inv);
      nz = (int)fceil((box.z2 - box.z1)*inv);
      if((acube = (Atomcube***) alloc_3D(nx, ny, nz, sizeof(Atomcube))) == NULL) {
	    logprint("still too many cubes,");
	    logprint("- third attempt");
	    inv = 0.5/cubeSize;
	    nx = (int)fceil((box.x2 - box.x1)*inv);
	    ny = (int)fceil((box.y2 - box.y1)*inv);
	    nz = (int)fceil((box.z2 - box.z1)*inv);
	    if((acube = (Atomcube***) alloc_3D(nx, ny, nz, sizeof(Atomcube))) == NULL){
		showinfobox("can't allocate enough memory for the Atomcubes");
		return;
	    }
      }
   }

   do_atomcubes(key);

   do_connections(key);

   free_atomcubes();

   box.x1 -= 1.4; box.y1 -= 1.4; box.z1 -= 1.4;
   box.x2 += 1.4; box.y2 += 1.4; box.z2 += 1.4;
}


void do_atomcubes(int key)
{
   register AtoM *ap;
   register int i;
   int ix, iy, iz, j;
   register Atomcube *ac;

   ac = acube[0][0];
   for(i=0; i<nx*ny*nz; ac++, i++) ac->na = 0;

   for(ap=actualmol->firstatom, i=0; ap; ap=ap->next, i++){
      if(key == H_BOND) {
         switch (ap->ord) {
            case 1:
            case 7:
            case 8:
            break;
            default: continue;
         }
      }
      ix = (int)ffloor((ap->coord[0] - box.x1)*inv);
      iy = (int)ffloor((ap->coord[1] - box.y1)*inv);
      iz = (int)ffloor((ap->coord[2] - box.z1)*inv);
      ac = acube[iz][iy] + ix;

      if(!ac->na){
         if((ac->a = (AtoM**) malloc(sizeof(AtoM *))) == NULL){
            showinfobox("no place for ac->a");
            return;
         }
      }
      else {
         if((ac->a = (AtoM**) realloc(ac->a, (ac->na+1)*sizeof(AtoM *))) == NULL){
            showinfobox("no realloc place for ac->a");
            return;
         }
      }

      j = ac->na++;
      ac->a[j] = ap;
   }
}


/* is it woth adding N, O, F, Cl, Br...? */

void do_connections(int key)
{
   register AtoM *a, *b, *c, *d, *e;
   Atomcube *ac;
   Bond *bp;
   int j, jx, jy, jz;
   int ix, iy, iz;
   int jx1, jx2, jy1, jy2, jz1, jz2;
   int flag;
   float angle;

   for(a = actualmol->firstatom; a; a = a->next){
      if(key == H_BOND) {
         switch (a->ord) {
            case 1:
            case 7:
            case 8:
            break;
            default: continue;
         }
      }
      else {
         if(FULLVALENCY(a->ord, a->nbonds)) continue;
      }

      ix = (int)ffloor((a->coord[0] - box.x1)*inv);
      iy = (int)ffloor((a->coord[1] - box.y1)*inv);
      iz = (int)ffloor((a->coord[2] - box.z1)*inv);

      jx1 = MAX(ix - 1, 0);
      jx2 = MIN(ix + 2, nx);
      jy1 = MAX(iy - 1, 0);
      jy2 = MIN(iy + 2, ny);
      jz1 = MAX(iz - 1, 0);
      jz2 = MIN(iz + 2, nz);

      for(jz = jz1; jz < jz2; jz++){
         for(jy = jy1; jy < jy2; jy++){
            ac = acube[jz][jy]+jx1;
            for(jx = jx1; jx < jx2; jx++, ac++){
               for(j=0; j<ac->na; j++){
                  b = ac->a[j];
                  if(a->name >= b->name) continue;
                  if(key == H_BOND) {
                     if(a->ord > 1 && b->ord > 1) continue;
                     if(a->ord == b->ord) continue;
                     if(DIST(a->coord, b->coord) <= max_h_bond
                         && DIST(a->coord, b->coord) > element[a->ord].coval + element[b->ord].coval){
                        flag = 1;
                        if(a->ord == 1) {
                          c = a; d = b;
                        }
                        else {
                          c = b; d = a;
                        }
                        if(d->nbonds && actualmol->firstbond) {
                           for(bp = actualmol->firstbond; bp; bp = bp->next){
                              if(bp->type != H_BOND && (d == bp->a || d == bp->b)) {
                                 e = (d == bp->a) ? bp->b : bp->a;
                                 angle = valence_angle(c->coord, d->coord, e->coord);
                                 if(angle <= min_h_angle) {
                                    flag = 0;
                                    break;
                                 }
                              }
                           }
                        }
                        if(flag) {
                           actualmol->nbonds++;
                           actualmol->n_h_bonds++;
                           add_bond(a, b);
                           actualmol->lastbond->type = H_BOND;
                           a->nbonds--;
                           b->nbonds--;
                        }
                     }
                  }
                  else {
                     if(FULLVALENCY(b->ord, b->nbonds)) continue;
                     if(DIST(a->coord, b->coord) <= 
                        element[a->ord].coval + element[b->ord].coval){
                        actualmol->nbonds++;
                        add_bond(a, b);
                     }
                  }
               }
            }
         }
      }
   }
}




void free_atomcubes(void)
{
   register Atomcube *ac;
   register int i;

   for(i=0, ac=acube[0][0]; i<nx*ny*nz; i++, ac++){
      if(ac->na) free(ac->a);
   }
   free(acube[0][0]);
   free(acube[0]);
   free(acube);
}



/* try to find multiple-bonds : to each bond with two unsaturated atoms
                        a bond is added (in the order of the bond-list) */
void find_multiplebonds(void)
{
   register Bond *bond;
   for(bond = actualmol->firstbond; bond; bond = bond->next){
      if(bond->a->nbonds < element[bond->a->ord].nbond - 1 &&
         bond->b->nbonds < element[bond->b->ord].nbond - 1){
         bond->type = TRIPLEBOND;
         bond->a->nbonds += 2;
         bond->b->nbonds += 2;
      }
      else if(bond->a->nbonds < element[bond->a->ord].nbond &&
              bond->b->nbonds < element[bond->b->ord].nbond){
         bond->type = DOUBLEBOND;
         bond->a->nbonds++;
         bond->b->nbonds++;
      }
   }
}


void computeOccupations(Mol *mp)  /* only depending on nr. of electrons! */
{
   int i;

   if(!mp->alphaOrbital) return;

   for(i=1; i<=mp->nAlpha; i++) {
      if(i < mp->firstOrbital) continue;
      mp->alphaOrbital[i - mp->firstOrbital].occ = 1;
   }

   if(mp->alphaBeta) {
      if(!mp->betaOrbital) return;
      for(i=1; i<=mp->nBeta; i++) {
         if(i < mp->firstOrbital) continue;
         mp->betaOrbital[i - mp->firstOrbital].occ = 1;
      }
   }

   else {
      for(i=1; i<=mp->nBeta; i++) {
         if(i < mp->firstOrbital) continue;
         mp->alphaOrbital[i - mp->firstOrbital].occ += 1;
      }
   }
}

void norm1(Shell *sp, double fac, double ex)
{
   register Gauss *i, *j;
   double sum, norm;

   sum = 0;
   for(i=sp->firstgauss; i; i=i->next){
      for(j=sp->firstgauss; j; j=j->next){
         sum += i->coeff * j->coeff * pow(i->exponent + j->exponent, ex);
      }
   }
   sum *= pow(M_PI, 1.5) * fac;
   norm = 1.0/sqrt(sum);

   for(i=sp->firstgauss; i; i=i->next) i->coeff *= norm;
}



void norm2(Shell *sp, double fac, double ex)
{
   register Gauss *i, *j;
   double sum, norm;

   sum = 0;
   for(i=sp->firstgauss; i; i=i->next){
      for(j=sp->firstgauss; j; j=j->next){
         sum += i->coeff2 * j->coeff2 * pow(i->exponent + j->exponent, ex);
      }
   }
   sum *= pow(M_PI, 1.5) * fac;
   norm = 1.0/sqrt(sum);

   for(i=sp->firstgauss; i; i=i->next) i->coeff2 *= norm;
}


void normalize_amoss(void)
/* normalize contracted gaussians */
{
   Amoss_basis *ap;
   Shell *sp;

   for(ap=actualmol->firstamoss; ap; ap=ap->next){
      for(sp=ap->firstshell; sp; sp=sp->next){
         switch(sp->n_base){
            case  1 : norm1(sp, 1.0, -1.5); break;
            case  3 : norm1(sp, 0.5, -2.5); break;
            case  4 : norm1(sp, 1.0, -1.5);
                      norm2(sp, 0.5, -2.5); break;
            case  6 : norm1(sp, 0.25, -3.5); break;
            case 10 : norm1(sp, 0.125, -4.5); break;
         }
      }
   }
}



void normalize_gaussians(void)
/* normalize contracted gaussians */
{
   AtoM  *ap;
   Shell *sp;

   for(ap=actualmol->firstatom; ap; ap=ap->next){
      for(sp=ap->firstshell; sp; sp=sp->next){
         switch(sp->n_base){
            case  1 : norm1(sp, 1.0, -1.5); break;
            case  3 : norm1(sp, 0.5, -2.5); break;
            case  4 : norm1(sp, 1.0, -1.5);
                      norm2(sp, 0.5, -2.5); break;
            case  5 :
            case  6 : norm1(sp, 0.25, -3.5); break;
            case  7 :
            case 10 : norm1(sp, 0.125, -4.5); break;
         }
      }
   }
}

/***
          CONVERSION OF A DOUBLE ZETA S.T.O.
          TO SINGLE ZETA S.T.O.
		(translated from FORTRAN February 93 - pff)
***/


static double cc1, cc2, zz1, zz2;
static int nn;

static double fov(double z)
{
   double t1, t2;

   t1 = 2.0*sqrt(z*zz1)/(z+zz1);
   t2 = 2.0*sqrt(z*zz2)/(z+zz2);

   return (-cc1*pow(t1, nn+nn+1) - cc2*pow(t2, nn+nn+1));
}

static double golden(double ax, double bx, double cx, double tol)
{
   static double r = 0.61803399, c = 1.0;
   double x0, x1, x2, x3, f0, f1, f2, f3;

   x0 = ax;
   x3 = cx;

   if(fabs(cx-bx) > fabs(bx-ax)){
      x1 = bx;
      x2 = bx + c*(cx-bx);
   }
   else {
      x2 = bx;
      x1 = bx + c*(bx-ax);
   }

   f1 = fov(x1);
   f2 = fov(x2);

   while(fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))){
      if(f2 < f1){
         x0 = x1;
         x1 = x2;
         x2 = r*x1 + c*x3;
         f0 = f1;
         f1 = f2;
         f2 = fov(x2);
      }
      else {
         x3 = x2;
         x2 = x1;
         x1 = r*x2 + c*x0;
         f3 = f2;
         f2 = f1; 
         f1 = fov(x1);
      }
   }

   if(f1 < f2) return x1;
   else        return x2;
}

double dszeta(int n, double c1, double z1, double c2, double z2)
{
   double s1, s2;

   if(z1 == z2) return z1;
   if(c1 == 0)  return z2;
   if(c2 == 0)  return z1;

   s1 = MIN(z1, z2);
   s2 = MAX(z1, z2);
   cc1 = c1;   cc2 = c2;
   zz1 = z1;   zz2 = z2;
   nn  = n;

   return golden(s1, 0.5*(s1+s2), s2, 0.0001);
}

void free_frequencies(void)
{
   Vibration *vib;
   int n, i;

   if(!actualmol->vibration) return;
   n = actualmol->n_frequencies%3 ? actualmol->n_frequencies/3 + 1 : actualmol->n_frequencies/3;
   for(vib=actualmol->vibration, i=0; i<n; i++, vib++){
      if(vib->coord) free(vib->coord);
   }
   free(actualmol->vibration);
   actualmol->vibration = NULL;
   actualmol->n_frequencies = 0;
}
