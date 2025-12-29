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


/* lecture of own format .mkl */
#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "readmkl.h"
#include "general.h"
#include "maininterf.h"
#include "utils.h"
#include "box.h"

#define RHF          1
#define ROHF         2
#define UHF          4



static int read_atomic_coordinates(char *name);
static char *find_string(char *s);
static int addTrajectoryStep(void);
static int read_atomic_charges(void);
static int read_charge(void);
static int read_basis_set(void);
static int read_coefficients(void);
static int read_occupations(void);
static int read_frequencies(void);
static int read_dipole(void);

static FILE *fp;
static char line[256];
static long previous_line = 0, preprevious = 0;


void read_mkl(char *name)
{
   if((fp = fopen(name, "r")) == NULL){
      sprintf(line, "can't open file\n%s !", name);
      showinfobox(line);
      return;
   }

   if(!read_atomic_coordinates(name)){
      showinfobox("can't read the atomic coordinates");
      fclose(fp);
      return;
   }

   if(!actualmol->natoms){
      sprintf(line, "No atoms in %s!", name);
      showinfobox(line);
      fclose(fp);
      return;
   }

   create_bonds();
   find_multiplebonds();
//   create_box();
   new_mole(name);

   rewind(fp);
   if(!read_atomic_charges()){
      showinfobox("unable to read the atomic charges");
      actualmol->charges = 0;
   }

   rewind(fp);
   free_frequencies();
   if(read_frequencies()){
      logprint("frequencies present");
   }

   rewind(fp);
   if(read_dipole()) {
      logprint("dipole moment present");
   }

   rewind(fp);
   if(!read_charge()){
      showinfobox("can't read the charge and multiplicity");
      fclose(fp);
      return;
   }

   rewind(fp);
   if(!read_basis_set()){
      showinfobox("can't read the basis-set");
      fclose(fp);
      return;
   }

   rewind(fp);
   if(!read_coefficients()){
      showinfobox("can't read the MO-coefficients");
      fclose(fp);
      return;
   }

   if(!read_occupations()){
      logprint("No occupations read");
      logprint("Calculated occupations used");
      computeOccupations(actualmol);
   }


   fclose(fp);
   update_logs();
}


static char *find_string(char *s)
{
   previous_line = ftell(fp);
   do {
      if(!fgets(line, 255, fp)) return NULL;
      if(strstr(line, s)) return line;
      preprevious = previous_line;
      previous_line = ftell(fp);
   } while (1);
}


static int read_atomic_coordinates(char *name)
{
   long fpos;
   float x, y, z;
   int ord, natoms;

   free_dyna();

   if(!find_string("$COORD")) {
      return 0;
   }

   do {
      fpos = ftell(fp);
      natoms = addTrajectoryStep();
      } while(strncmp(line, "$END", 4)!= 0);

      add_mol(name);
      dynamics.molecule = actualmol;
      dynamics.current = dynamics.ntotalsteps - 1;

   fseek(fp, fpos, SEEK_SET);

   fgets(line, 255, fp);
   do {
      if(sscanf(line, "%d %f %f %f", &ord, &x, &y, &z) != 4) return 0;
      if(ord >= 0) add_atom(ord, x, y, z);
      fgets(line, 255, fp);
   } while(strncmp(line, "$END", 4) != 0);

   return 1;
}

static int addTrajectoryStep(void)
{
   long fpos, index, i;
   float x, y, z;
   static int natoms;
   
   index = dynamics.ntotalsteps++;

   if(!dynamics.trajectory){
      if((dynamics.trajectory = (Vector **)malloc(sizeof(Vector *))) == NULL){
         fprintf(stderr, "can't allocate dyna pointer\n");
         return 0;
      }

      fpos = ftell(fp);
      fgets(line, 255, fp);  
      natoms = 0;

/* count nr of atoms */
      do {
         natoms++;
         fgets(line, 255, fp);
      } while(strncmp(line, "$$", 2) != 0 && strncmp(line, "$END", 4) != 0);
      fseek(fp, fpos, SEEK_SET);

   }
   else {
      if((dynamics.trajectory = (Vector **)realloc(dynamics.trajectory,
          dynamics.ntotalsteps*sizeof(Vector *))) == NULL){
         fprintf(stderr, "can't reallocate dyna pointer\n");
         dynamics.ntotalsteps--;
         free_dyna();
         return 0;
      }
   }

   if((dynamics.trajectory[index] = (Vector *)malloc(natoms*sizeof(Vector))) == NULL){
      fprintf(stderr, "can't allocate timestep dyna[%d]\n", index);
      dynamics.ntotalsteps--;
      free_dyna();
      return 0;
   }


   fgets(line, 255, fp);
   i = 0;
   do {
      if(sscanf(line, "%*d %f %f %f", &x, &y, &z) != 3) return 0;
      dynamics.trajectory[index][i].x = x;
      dynamics.trajectory[index][i].y = y;
      dynamics.trajectory[index][i].z = z;
      i++;
      fgets(line, 255, fp);
   } while(strncmp(line, "$$", 2) != 0 && strncmp(line, "$END", 4) != 0);

   return natoms;

}


static int read_atomic_charges(void)
{
   AtoM *ap;

   if(!find_string("$CHARGES")) return 0;
   
   for(ap=actualmol->firstatom; ap; ap = ap->next) {
      if(!fgets(line, 255, fp)) return 0;
      if(sscanf(line, "%f", &ap->charge) != 1) return 0;
   }

   actualmol->charges = 1;
   return 1;
}


static int read_charge(void)
{
   AtoM *ap;

   if(!find_string("$CHAR_MULT")) return 0;

   if(!fgets(line, 255, fp)) return 0;
   if(sscanf(line, "%f %d", &actualmol->charge, &actualmol->multiplicity) != 2) return 0;
   
   for(ap=actualmol->firstatom; ap; ap = ap->next) {
      actualmol->nElectrons += ap->ord;
   }
   actualmol->nElectrons -= actualmol->charge;
   actualmol->nBeta = (actualmol->nElectrons - actualmol->multiplicity +1 )/2;
   actualmol->nAlpha = actualmol->nElectrons - actualmol->nBeta;
   sprintf(line, "Charge: %2.0f", actualmol->charge);
   logprint(line);
   sprintf(line, "Multiplicity: %d", actualmol->multiplicity);
   logprint(line);
   sprintf(line, "Alpha Electrons: %d", actualmol->nAlpha);
   logprint(line);
   sprintf(line, "Beta Electrons: %d", actualmol->nBeta);
   logprint(line);
   return 1;
}


static int read_basis_set(void)
{
   AtoM *ap;
   Gauss *gp;
   Shell *sp;
   double s_coeff, p_coeff, d_coeff, f_coeff, alpha, norm;
   char type[5], cp[256];

   if(!find_string("$BASIS")) return 0;

   for(ap = actualmol->firstatom; ap; ap = ap->next){
      fgets(line, 255, fp);
      do {
         if(!(sp = add_shell(ap))) return 0;
         if(sscanf(line, "%hd %s %f", &sp->n_base, type, &sp->scale_factor) != 3 ) return 0;
         actualmol->nBasisFunctions += sp->n_base;
         if(!fgets(line, 255, fp)) return 0;
         sscanf(line, "%*d %s %*f", cp);
         while(!isalpha(*cp) && *line != 36) {
            gp = add_gauss(sp);
            if(strcmp(type, "S") == 0 || strcmp(type, "s") == 0) {
               sscanf(line, "%lf %lf", &alpha, &s_coeff);
               gp->exponent = alpha * sp->scale_factor * sp->scale_factor;
               alpha = gp->exponent;
               norm  = pow(2.0 * alpha / M_PI, 0.75);
               gp->coeff  = s_coeff * norm;
            }
            else if(strcmp(type, "SP") == 0 || strcmp(type, "sp") == 0) {
               sscanf(line, "%lf %lf %lf", &alpha, &s_coeff, &p_coeff);
               gp->exponent = alpha * sp->scale_factor * sp->scale_factor;
               alpha = gp->exponent;
               norm  = pow(2.0 * alpha / M_PI, 0.75);
               gp->coeff  = s_coeff * norm;
               norm  = pow(128.0 * pow(alpha, 5) / pow(M_PI, 3), 0.25);
               gp->coeff2 = p_coeff * norm;
            }
            else if(strcmp(type, "P") == 0 || strcmp(type, "p") == 0){
               sscanf(line, "%lf %lf", &alpha, &p_coeff);
               gp->exponent = alpha * sp->scale_factor * sp->scale_factor;
               alpha = gp->exponent;
               norm  = pow(128.0 * pow(alpha, 5) / pow(M_PI, 3), 0.25);
               gp->coeff  = p_coeff * norm;
            }
            else if(strcmp(type, "D") == 0 || strcmp(type, "d") == 0){
               sscanf(line, "%lf %lf", &alpha, &d_coeff);
               gp->exponent = alpha * sp->scale_factor * sp->scale_factor;
               alpha = gp->exponent;
               norm = pow(2048. * pow(alpha, 7) / pow(M_PI, 3), .25);
               gp->coeff  = d_coeff * norm;
            }
            else if(strcmp(type, "F") == 0 || strcmp(type, "f") == 0){
               sscanf(line, "%lf %lf", &alpha, &f_coeff);
               gp->exponent = alpha * sp->scale_factor * sp->scale_factor;
               alpha = gp->exponent;
               norm = pow(32768. * pow(alpha, 9) / pow(M_PI, 3), .25);
               gp->coeff = f_coeff * norm;
            }
            else {
               showinfobox("No coefficients in gaussian primitive");
               return 0;
            }

            if(!fgets(line, 255, fp)) return 0;
            sscanf(line, "%*d %s %*f", cp);
         }
      } while(*line != 36);
   }

   normalize_gaussians();

   return 1;

}

static int read_coefficients(void)
{
   long fpos;
   int nblocks, opl;
   MolecularOrbital *alphaOrb, *betaOrb;
   register int i, j;
   char str[256];
   short n[6];

   if(!find_string("$COEFF_ALPHA")) return 0;
   fpos = ftell(fp);

   rewind(fp);

   if(find_string("$COEFF_BETA")) {
      actualmol->alphaBeta = 1;
   }

   fseek(fp, fpos, SEEK_SET);

   actualmol->nMolecularOrbitals = actualmol->nBasisFunctions;

   actualmol->lastOrbital = actualmol->firstOrbital + actualmol->nMolecularOrbitals;

   if((alphaOrb = allocOrbital(actualmol->nBasisFunctions, actualmol->nBasisFunctions, GAUSS_ORB)) == NULL){
      showinfobox("can't allocate the MO-structures");
      return 0;
   }

   if(actualmol->alphaBeta) {
      if((betaOrb = allocOrbital(actualmol->nBasisFunctions, actualmol->nBasisFunctions, GAUSS_ORB)) == NULL){
         showinfobox("can't allocate the MO-structures");
         return 0;
      }
   }

   nblocks = (int)fceil(actualmol->nBasisFunctions / 5.0);
   j = 0;
   while(j < nblocks){
      fgets(line, 255, fp);
      if(strncmp(line, "$END", 4) == 0) break;
      opl = sscanf(line, "%s %s %s %s %s", str, str, str, str, str);
      n[1] = n[2] = n[3] = n[4] = n[5] = 0;
      while(opl--) n[opl+1] = opl;
      sscanf(line, "%s %s %s %s %s",\
          alphaOrb[n[1]+j*5].type, alphaOrb[n[2]+j*5].type, alphaOrb[n[3]+j*5].type,
          alphaOrb[n[4]+j*5].type, alphaOrb[n[5]+j*5].type);
      fgets(line, 255, fp);
      sscanf(line, "%lf %lf %lf %lf %lf", &alphaOrb[n[1]+j*5].eigenvalue,
          &alphaOrb[n[2]+j*5].eigenvalue, &alphaOrb[n[3]+j*5].eigenvalue,
          &alphaOrb[n[4]+j*5].eigenvalue, &alphaOrb[n[5]+j*5].eigenvalue);
      for(i=0; i<actualmol->nBasisFunctions; i++){
         fgets(line, 255, fp);
         sscanf(line, "%lf %lf %lf %lf %lf", &alphaOrb[n[1]+j*5].coefficient[i],
         &alphaOrb[n[2]+j*5].coefficient[i], &alphaOrb[n[3]+j*5].coefficient[i],
         &alphaOrb[n[4]+j*5].coefficient[i], &alphaOrb[n[5]+j*5].coefficient[i]);
      }
      j++;
   }

   if(!actualmol->alphaBeta) {
      actualmol->alphaOrbital = alphaOrb;
      return 1;
   }

   rewind(fp);

   if(find_string("$COEFF_BETA")) {
      actualmol->alphaBeta = 1;
   }

   j = 0;
   while(j < nblocks){
      fgets(line, 255, fp);
      if(strncmp(line, "$END", 4) == 0) break;
      opl = sscanf(line, "%s %s %s %s %s", str, str, str, str, str);
      n[1] = n[2] = n[3] = n[4] = n[5] = 0;
      while(opl--) n[opl+1] = opl;
      sscanf(line, "%s %s %s %s %s",\
          betaOrb[n[1]+j*5].type, betaOrb[n[2]+j*5].type, betaOrb[n[3]+j*5].type,
          betaOrb[n[4]+j*5].type, betaOrb[n[5]+j*5].type);
      fgets(line, 255, fp);
      sscanf(line, "%lf %lf %lf %lf %lf", &betaOrb[n[1]+j*5].eigenvalue,
          &betaOrb[n[2]+j*5].eigenvalue, &betaOrb[n[3]+j*5].eigenvalue,
          &betaOrb[n[4]+j*5].eigenvalue, &betaOrb[n[5]+j*5].eigenvalue);
      for(i=0; i<actualmol->nBasisFunctions; i++){
         fgets(line, 255, fp);
         sscanf(line, "%lf %lf %lf %lf %lf", &betaOrb[n[1]+j*5].coefficient[i],
         &betaOrb[n[2]+j*5].coefficient[i], &betaOrb[n[3]+j*5].coefficient[i],
         &betaOrb[n[4]+j*5].coefficient[i], &betaOrb[n[5]+j*5].coefficient[i]);
      }
      j++;
   }

   actualmol->alphaOrbital = alphaOrb;
   actualmol->betaOrbital = betaOrb;
   return 1;

}


static int read_occupations(void)
{
   int nblocks, opl;
   register int j;
   float flo;
   short n[6];
 
   if(!actualmol->alphaOrbital) return 0;
   rewind(fp);

   if(!find_string("$OCC_ALPHA")) return 0;

   nblocks = actualmol->nBasisFunctions / 5;
   j = 0;
   while(j <= nblocks){
         fgets(line, 255, fp);
         if(strncmp(line, "$END", 4) == 0) break;
         opl = sscanf(line, "%f %f %f %f %f", &flo, &flo, &flo, &flo, &flo);
         n[1] = n[2] = n[3] = n[4] = n[5] = 0;
         while(opl--) n[opl+1] = opl;
         sscanf(line, "%f %f %f %f %f", &actualmol->alphaOrbital[n[1]+j*5].occ,
         &actualmol->alphaOrbital[n[2]+j*5].occ, &actualmol->alphaOrbital[n[3]+j*5].occ,
         &actualmol->alphaOrbital[n[4]+j*5].occ, &actualmol->alphaOrbital[n[5]+j*5].occ);
         j++;
   }

   if(actualmol->alphaBeta) {
      rewind(fp);
      if(!find_string("$OCC_BETA")) return 0;
         j = 0;
         while(j <= nblocks){
               fgets(line, 255, fp);
               if(strncmp(line, "$END", 4) == 0) break;
               opl = sscanf(line, "%f %f %f %f %f", &flo, &flo, &flo, &flo, &flo);
               n[1] = n[2] = n[3] = n[4] = n[5] = 0;
               while(opl--) n[opl+1] = opl;
               sscanf(line, "%f %f %f %f %f", &actualmol->betaOrbital[n[1]+j*5].occ,
               &actualmol->betaOrbital[n[2]+j*5].occ, &actualmol->betaOrbital[n[3]+j*5].occ,
               &actualmol->betaOrbital[n[4]+j*5].occ, &actualmol->betaOrbital[n[5]+j*5].occ);
               j++;
         }

   }

   return 1;
}


static int read_frequencies(void)
{
   long fpos;
   int count = 0;
   int n_freq;
   register short i, j;
   AtoM *ap;
   Vibration *vib;

   if(!find_string("$FREQ")) return 0;
   fpos = ftell(fp);
   fgets(line, 255, fp);
   do {
      fgets(line, 255, fp);
      count++;
   } while(strncmp(line, "$END", 4) != 0);
   n_freq = 3 * (count / (actualmol->natoms + 2));

   if((actualmol->vibration = (Vibration *)malloc(n_freq*sizeof(Vibration))) == NULL){
      showinfobox("can't realloc vibrational frequency");
      return 0;
   }
   for(i=0, vib=actualmol->vibration; i<n_freq; i++, vib++){
      if((vib->coord = (Vector *)malloc(actualmol->natoms*sizeof(Vector))) == NULL){
         showinfobox("can't allocate vibration");
         return 0;
      }
   }

   fseek(fp, fpos, SEEK_SET);
   for(i=0, vib=actualmol->vibration; i<n_freq/3; i++, vib += 3) {
      if(!fgets(line, 255, fp)) return 0;
      actualmol->n_frequencies += sscanf(line, "%s %s %s", vib->type, (vib+1)->type, (vib+2)->type);
      if(!fgets(line, 255, fp)) return 0;
      sscanf(line, "%f %f %f", &vib->frequency, &(vib+1)->frequency, &(vib+2)->frequency);
      for(ap=actualmol->firstatom, j=0; ap; ap=ap->next, j++){
         if(!fgets(line, 255, fp)) return 0;
         sscanf(line, "%f %f %f %f %f %f %f %f %f", &vib->coord[j].x,
            &vib->coord[j].y, &vib->coord[j].z, &(vib+1)->coord[j].x,
            &(vib+1)->coord[j].y, &(vib+1)->coord[j].z, &(vib+2)->coord[j].x,
            &(vib+2)->coord[j].y, &(vib+2)->coord[j].z);
      }
   }

   return 1;
}


static int read_dipole(void) {

   float x, y, z;

   if(!find_string("$DIPOLE")) return 0;
   fgets(line, 255, fp);
   sscanf(line, "%f %f %f", &x, &y, &z);
   actualmol->dipole = add_dipole(x, y, z);

   return 1;
}
