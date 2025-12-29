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


/* lecture of GAMESS output */
#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "readgamess.h"
#include "general.h"
#include "utils.h"
#include "maininterf.h"
#include "box.h"



static int read_atomic_coordinates(void);
static char *find_string(char *s);
static int read_charge(void);
static int read_basis_set(void);
static int read_eigenvectors(void);
static int read_frequencies(void);
static int read_dipole(void);
static int addGMTrajectoryStep(void);

static FILE *fp;
static char line[256];
static int nblocks;
static long previous_line = 0, preprevious = 0;


/**** lecture of GAMESS output ****/

void print_amoss_basis_set(void)
{
   AtoM *basis;
   Shell *shell;
   Gauss *gauss;

   for(basis = actualmol->firstatom; basis; basis = basis->next){
      printf("Atom %s\n", element[basis->ord].symbol);
      for(shell = basis->firstshell; shell; shell = shell->next){
         printf("          Shell (%d), scale %f\n",
            shell->n_base, shell->scale_factor);
         for(gauss = shell->firstgauss; gauss; gauss = gauss->next){
            printf("                         %12.5f %8.4f %8.4f\n",
                     gauss->exponent, gauss->coeff, gauss->coeff2);
         }
      }
   }
}

void read_gamess(char *file)
{

   if((fp = fopen(file, "r")) == NULL){
      sprintf(line, "Can't open file\n%s !", file);
      showinfobox(line);
      return;
   }

   if(!find_string("GAMESS")) {
      showinfobox("read_gamess : GAMESS (US) output only!");
      fclose(fp);
      return;
   }

   add_mol(file);

   if(!read_basis_set()){
      showinfobox("Can't read the basis-set!");
      rewind(fp);
   }

   if(!read_atomic_coordinates()){
      showinfobox("Can't read the atomic coordinates!");
      fclose(fp);
      free_mol(actualmol);
      return;
   }

   if(!actualmol->natoms){
      sprintf(line, "No atoms in %s!", file);
      showinfobox(line);
      free_mol(actualmol);
      fclose(fp);
      return;
   }

   create_bonds();
   find_multiplebonds();
//   create_box();
   new_mole(file);

/*
   print_amoss_basis_set();
*/

   if(!read_eigenvectors()){
      logprint("can't read the eigenvalues!");
   }
   computeOccupations(actualmol);

/*
   print_coefficients();
*/

   if(!read_charge()){
      logprint("unable to read");
      logprint("   the atomic charges!");
      actualmol->charges = 0;
   }

   if(read_frequencies()){
      logprint("frequencies present");
   }

   if(read_dipole()) {
      logprint("dipole moment present");
   }

   fclose(fp);

   if(maininterfwin) update_logs();
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




static Shell *get_basis(char *basis)
{
   register Amoss_basis *bp;

   for(bp = actualmol->firstamoss; bp; bp = bp->next){
      if(!strcmp(basis, bp->basis)) return bp->firstshell;
   }
   return NULL;
}


static int read_atomic_coordinates(void)
{
   AtoM *ap;
   long fpos;
   float x, y, z;
   char basis_symbol[20], *format;
   float ord;
   unsigned angst = 0;
   int natoms;

   free_dyna();

   rewind(fp);

   if(find_string("COORDINATES (BOHR)")) {
      do {
         fpos = ftell(fp);
      } while(find_string("COORDINATES (BOHR)"));

      fseek(fp, fpos, SEEK_SET);
      if(!find_string("CHARGE         X                   Y ")){
         fseek(fp, fpos, SEEK_SET);
         if(!find_string("ATOM     ZNUC       X             Y")) return 0;
         fgets(line, 255, fp);
         fgets(line, 255, fp);
         format = "%*d%s %f %f %f %f";
      }
      else format = "%s %f %f %f %f";
      fpos = ftell(fp);
   }
   else return 0;

   rewind(fp);
   if(find_string("COORDINATES OF ALL ATOMS ARE (ANGS)")) {
      do {
         angst = 1;
         fpos = ftell(fp);
         natoms = addGMTrajectoryStep();
      } while(find_string("COORDINATES OF ALL ATOMS ARE (ANGS)"));
      fseek(fp, fpos, SEEK_SET);
      fgets(line, 255, fp);
      fgets(line, 255, fp);
      fpos = ftell(fp);
   }

   dynamics.molecule = actualmol;
   dynamics.current = dynamics.ntotalsteps - 1;

   fseek(fp, fpos, SEEK_SET);
   while(1) {
      if(!fgets(line, 255, fp)) return 0;
      if(sscanf(line, format, basis_symbol, &ord, &x, &y, &z)
         != 5) break;
      if(angst) ap = add_atom((int)ord, x, y, z);
      else ap = add_atom((int)ord, x*BOHR, y*BOHR, z*BOHR);
      if(actualmol->firstamoss) {
         if(!(ap->firstshell = get_basis(basis_symbol))) {
            return 0;
         }
      }
   }

   return 1;
}




static int addGMTrajectoryStep(void)
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
/* count nr of atoms */
      fpos = ftell(fp);
      fgets(line, 255, fp);
      fgets(line, 255, fp);
      fgets(line, 255, fp);
      natoms = 0;
      do {
         natoms++;
         if(!fgets(line, 255, fp)) break;
      } while(strlen(line) > 1);
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
      fpos = ftell(fp);
   }

   if((dynamics.trajectory[index] = (Vector *)malloc(natoms*sizeof(Vector))) == NULL){
      fprintf(stderr, "can't allocate timestep dyna[%d]\n", index);
      dynamics.ntotalsteps--;
      free_dyna();
      return 0;
   }

   fgets(line, 255, fp);
   fgets(line, 255, fp);
   fgets(line, 255, fp);

   i = 0;
   do {
      if(sscanf(line, "%*s %*f %f %f %f", &x, &y, &z) != 3) return 0;
      dynamics.trajectory[index][i].x = x;
      dynamics.trajectory[index][i].y = y;
      dynamics.trajectory[index][i].z = z;
      if(!fgets(line, 255, fp)) break;
      i++;
   } while(strlen(line) > 1);



   fseek(fp, fpos, SEEK_SET);
   return natoms;
}

static int read_basis_set(void)
{
   Amoss_basis *ap;
   Gauss *gp;
   Shell *sp;
   int i1, i2;
   float f1, f2, f3;
   double s_coeff, p_coeff, d_coeff, f_coeff, alpha, norm;
   int new_gamess = 1;

   if(!find_string("ATOMIC BASIS SET")) return 0;
   if(!find_string("CONTRACTED PRIMITIVE FUNCTIONS")) return 0;
   if(!find_string("CONTRACTION COEFFICIENTS")) return 0;
   fgets(line, 255, fp); /*  blank line */

   while (!strstr(line, "TOTAL NUMBER OF")) {

      while(fgets(line, 255, fp)) {

         if(line[0] == '\n' || line[1] == '\n') continue;

         if(strlen(line) < 16) { /* new basis set symbol */
            if(!(ap = add_amoss(actualmol))) return 0;
            sscanf(line, "%s", ap->basis);
            break;
         }

         if(strchr(line, '(')){
            new_gamess = 0;
         }

         if(new_gamess) {
            if(sscanf(line, "%d %*s %d %f %f",
                           &i1, &i2, &f1, &f2) != 4) break;
         }
         else {
            if(sscanf(line, "%d %*s %d %f %f %*s %f",
                           &i1, &i2, &f1, &f2, &f3) != 5) break;
         }

         if(!(sp = add_shell_to_amoss(ap))) return 0;
         sp->scale_factor = 1.0;

         while(line[0] != '\n' && line[1] != '\n') {
            gp = add_gauss(sp);
            gp->coeff2 = 0;

            if(strchr(line, 'S')){
               if(new_gamess) {
                  sscanf(line, "%*d %*s %*d %lf %lf",
                     &gp->exponent, &s_coeff);
                  norm  = pow(2.0 * gp->exponent / M_PI, 0.75);
                  gp->coeff  = s_coeff * norm;
               }
               else {
                  sscanf(line, "%*d %*s %*d %lf %lf %*s %*f",
                     &gp->exponent, &gp->coeff);
               }
               sp->n_base = 1;
            }
            else if(strchr(line, 'P')){
               if(new_gamess) {
                  sscanf(line, "%*d %*s %*d %lf %lf",
                     &gp->exponent, &p_coeff);
                  norm  = pow(128.0 * pow(gp->exponent, 5) / pow(M_PI, 3), 0.25);
                  gp->coeff  = p_coeff * norm;
               }
               else {
                  sscanf(line, "%*d %*s %*d %lf %lf %*s %*f",
                     &gp->exponent, &gp->coeff);
               }
               sp->n_base = 3;
            }
            else if(strchr(line, 'L')){
/*
               sscanf(line, "%*d %*s %*d %lf %lf %*s %*f",
                  &gp->exponent, &gp->coeff);
               sscanf(line+54, "%lf", &gp->coeff2);
*/
               if(new_gamess) {
                  sscanf(line, "%*d %*s %*d %lf %lf %lf",
                     &gp->exponent, &s_coeff, &p_coeff);
                  norm  = pow(2.0 * gp->exponent / M_PI, 0.75);
                  gp->coeff  = s_coeff * norm;
                  norm  = pow(128.0 * pow(gp->exponent, 5) / pow(M_PI, 3), 0.25);
                  gp->coeff2 = p_coeff * norm;
               }
               else {
                  if(sscanf(line, "%*d %*s %*d %lf %lf %*s %*f%*s %lf %*s %*f",
                     &gp->exponent, &gp->coeff, &gp->coeff2) == 2) {
                       fgets(line, 255, fp);
                       sscanf(line, "%lf", &gp->coeff2);
                  }
               }
               sp->n_base = 4;
            }
            else if(strchr(line, 'D')){
               if(new_gamess) {
                  sscanf(line, "%*d %*s %*d %lf %lf",
                     &gp->exponent, &d_coeff);
                  norm = pow(2048. * pow(gp->exponent, 7) / pow(M_PI, 3), .25);
                  gp->coeff  = d_coeff * norm;
// 4.4.01 STP: the line below gives the same number, that was present in
// the gamess output before. But what I need is above.
//                  gp->coeff  = d_coeff * norm / sqrt(3);
               }
               else {
                  sscanf(line, "%*d %*s %*d %lf %lf %*s %*f",
                     &gp->exponent, &d_coeff);
                     gp->coeff  = d_coeff * 1.7320508;
// 4.4.01 STP: d_coeff seems not to be what was thought it to be:
// d_coeff is:  CDINP * (2048 EX^7/PI^3)^1/4 *(1/3)^1/2
// need to correct with factor 3^1/2
//                     gp->coeff  = d_coeff;
               }
               sp->n_base = 6;
            }
            else if(strchr(line, 'F')){
               if(new_gamess) {
                  sscanf(line, "%*d %*s %*d %lf %lf",
                     &gp->exponent, &f_coeff);
                  norm = pow(32768. * pow(gp->exponent, 9) / pow(M_PI, 3), .25);
                  gp->coeff  = f_coeff * norm;
// 4.4.01 STP: the line below gives the same number, that was present in
// the gamess output before. But what I need is above.
//                  gp->coeff  = f_coeff * norm / sqrt(15);
               }
               else {
                  sscanf(line, "%*d %*s %*d %lf %lf %*s %*f",
                     &gp->exponent, &f_coeff);
                     gp->coeff = f_coeff * 3.8729833;
// 4.4.01 STP: f_coeff seems not to be what was thought it to be:
// f_coeff is: CFINP * (32768 EX^9/PI^3)^1/4 *(1/15)^1/2
// need to correct with factor 15^1/2
//                     gp->coeff = f_coeff;
               }
               sp->n_base = 10;
            }
            else {
               logprint("can't read");
               logprint("   gaussian primitives");
               return 0;
            }
            if(!fgets(line, 255, fp)) return 0;
         }                                /* end of gaussian primitives */
      }                                   /* end of shell */
   }                                      /* end of basis-set */

   normalize_gaussians();

   return 1; 
}


static int read_eigenvectors(void)
{
   long fpos, itest;
   MolecularOrbital *mo, *mbeta;
   int i1, i2, i3, i4, i5, i6, i7, i8, i9, i0, nnn, n_mo;
   int incr = 38;
   register int i, j;
   char *eigenstr = "EIGENVECTORS",
               *alphastr = "----- ALPHA SET -----",
               *betastr  = "----- BETA SET -----";

   char *pkey;
   unsigned equgeo = 0;

   rewind(fp);

   if(!find_string("TOTAL NUMBER OF BASIS FUNCTIONS")) {
      incr = 47;
      rewind(fp);
      if(!find_string("NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS")) return 0;
   }
   sscanf(line + incr, "%d", &actualmol->nBasisFunctions);
   if(!find_string("NUMBER OF OCCUPIED ORBITALS (ALPHA)")) return 0;
   sscanf(line + incr, "%d", &actualmol->nAlpha);
   if(!find_string("NUMBER OF OCCUPIED ORBITALS (BETA )")) return 0;
   sscanf(line + incr, "%d", &actualmol->nBeta);

   if(find_string("EQUILIBRIUM GEOMETRY LOCATED")) {
      alphastr = "**** ALPHA SET ****";
      betastr = "**** BETA SET ****";
      eigenstr = "MOLECULAR ORBITALS\n";
      equgeo = 1;
   }

   rewind(fp);

   if(find_string(alphastr)) {
      actualmol->alphaBeta = 1;
      pkey = alphastr;
   }
   else {
      actualmol->alphaBeta = 0;
      pkey = eigenstr;
   }
   rewind(fp);

   if(!find_string(eigenstr)) {
      rewind(fp);
      pkey = " ORBITALS\n";
      if(!find_string(pkey)) return 0;
   }

   do {
      fpos = ftell(fp);
   } while(find_string(pkey));
   fseek(fp, fpos, SEEK_SET);

   if(actualmol->alphaBeta){
      if(!equgeo) if(!find_string(eigenstr))  return 0;
   }
   nblocks = (actualmol->nBasisFunctions - 1)/5 + 1;

   nnn = actualmol->nBasisFunctions;


   if((mo = allocOrbital(nnn, nnn, GAMESS_ORB)) == NULL){
      logprint("can't allocate the MO-structures");
      return 0;
   }
 
   fseek(fp, fpos, SEEK_SET);

   j = nblocks;
   n_mo = 0;

   i1 = i2 = i3 = i4 = i5 = i6 = i7 = i8 = i9 = i0 = 0;

   if(!(equgeo && actualmol->alphaBeta)) fgets(line, 255, fp); /* ------------ */

   while(j--){
      fgets(line, 255, fp); /* blank */
      fgets(line, 255, fp); /* MO-numbers */
      itest = n_mo + 1;
      n_mo += sscanf(line, " %d %d %d %d %d %d %d %d %d %d",
         &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &i0);
      if(i1 != itest) break;

      fgets(line, 255, fp); /* Eigenvalues */
      sscanf(line, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &mo[--i1].eigenvalue, 
         &mo[--i2].eigenvalue, &mo[--i3].eigenvalue, &mo[--i4].eigenvalue, 
         &mo[--i5].eigenvalue, &mo[--i6].eigenvalue, &mo[--i7].eigenvalue,
         &mo[--i8].eigenvalue, &mo[--i9].eigenvalue, &mo[--i0].eigenvalue);

      fgets(line, 255, fp); /* Symmetries */

      for(i=0; i<actualmol->nBasisFunctions; i++){
         fgets(line, 255, fp);
         sscanf(line+15, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",  &mo[i1].coefficient[i], 
            &mo[i2].coefficient[i], &mo[i3].coefficient[i], &mo[i4].coefficient[i],
            &mo[i5].coefficient[i], &mo[i6].coefficient[i], &mo[i7].coefficient[i], 
            &mo[i8].coefficient[i], &mo[i9].coefficient[i], &mo[i0].coefficient[i]);
      }
   }
   actualmol->nMolecularOrbitals = n_mo;
   actualmol->alphaOrbital = mo;

   if(!actualmol->alphaBeta) return 1;

   if((mbeta = allocOrbital(nnn, nnn, GAMESS_ORB)) == NULL){
      logprint("can't allocate the MO-structures");
      return 0;
   }

   if(!find_string(betastr)) return 0;
   if(!equgeo) if(!find_string(eigenstr)) return 0;

   j = nblocks;
   n_mo = 0;
   i1 = i2 = i3 = i4 = i5 = i6 = i7 = i8 = i9 = i0 = 0;

   if(!equgeo) fgets(line, 255, fp); /* ------------ */

   while(j--){
      fgets(line, 255, fp); /* blank */
      fgets(line, 255, fp); /* MO-numbers */
      itest = n_mo + 1;
      n_mo += sscanf(line, " %d %d %d %d %d %d %d %d %d %d",
         &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &i0);
      if(i1 != itest) break;

      fgets(line, 255, fp); /* Eigenvalues */
      sscanf(line, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &mbeta[--i1].eigenvalue, 
         &mbeta[--i2].eigenvalue, &mbeta[--i3].eigenvalue, &mbeta[--i4].eigenvalue, 
         &mbeta[--i5].eigenvalue, &mbeta[--i6].eigenvalue, &mbeta[--i7].eigenvalue,
         &mbeta[--i8].eigenvalue, &mbeta[--i9].eigenvalue, &mbeta[--i0].eigenvalue);

      fgets(line, 255, fp); /* Symmetries */

      for(i=0; i<nnn; i++){
         fgets(line, 255, fp);
         sscanf(line+15, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",  &mbeta[i1].coefficient[i], 
            &mbeta[i2].coefficient[i], &mbeta[i3].coefficient[i], &mbeta[i4].coefficient[i],
            &mbeta[i5].coefficient[i], &mbeta[i6].coefficient[i], &mbeta[i7].coefficient[i], 
            &mbeta[i8].coefficient[i], &mbeta[i9].coefficient[i], &mbeta[i0].coefficient[i]);
      }
   }

   actualmol->nMolecularOrbitals = n_mo;
   actualmol->alphaOrbital = mo;
   actualmol->betaOrbital = mbeta;

   return 1;
}






static int read_charge(void)
{
   long fpos;
   float ch;
   AtoM *ap;

   rewind(fp);
   if(!find_string("TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS")) return 0;
   do {
      fpos = ftell(fp);
   } while(find_string("TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS"));
   fseek(fp, fpos, SEEK_SET);

   fgets(line, 255, fp);
   for(ap = actualmol->firstatom; ap; ap = ap->next){
      if(!fgets(line, 255, fp)) return 0;
      if(!sscanf(line, "%*d%*s%*f%f", &ch)) return 0;
      ap->charge = ch;
   }
   actualmol->charges = 1;
   return 1;
}




static int read_frequencies(void)
{
   long fpos;
   int n_freq;
   register short i, j;
   AtoM *ap;
   Vibration *vib;
   char *iptr = NULL;

   rewind(fp);

   if(!find_string("FREQUENCIES IN CM**-1")) return 0;
   fpos = ftell(fp);

   n_freq = 0;
   while(find_string(" FREQUENCY: ")) n_freq += 5;
   if(!n_freq) return 0;

   if((actualmol->vibration = (Vibration *)malloc(n_freq*sizeof(Vibration))) == NULL){
      logprint("can't realloc");
      logprint("   vibrational frequency");
      return 0;
   }
   for(i=0, vib=actualmol->vibration; i<n_freq; i++, vib++){
      if((vib->coord = (Vector *)malloc(actualmol->natoms*sizeof(Vector))) == NULL){
         logprint("can't allocate vibration");
         return 0;
      }
   }

   fseek(fp, fpos, SEEK_SET);
   for(i=0, vib=actualmol->vibration; i<n_freq/5; i++, vib += 5) {
      find_string(" FREQUENCY: ");
/* replace all I for imaginary freq with a blank */
      iptr = strchr(line, 'I');
      while (iptr) {
         *iptr = ' ';
         iptr = strchr(line, 'I');
      }
      actualmol->n_frequencies += sscanf(line+20, "%f%f%f%f%f",
         &vib->frequency, &(vib+1)->frequency, &(vib+2)->frequency,
         &(vib+3)->frequency, &(vib+4)->frequency);

      vib->type[0] = (vib+1)->type[0] = (vib+2)->type[0] = 0;
      (vib+3)->type[0] = (vib+4)->type[0] = 0;

      fgets(line, 255, fp);
      /* newer gamess versions print also the reduced mass -> additional fgets */
      if (strstr(line, "REDUCED MASS")) {
         fgets(line, 255, fp); /* either empty or with Intensities */
      }
      if(strlen(line) > 1) fgets(line, 255, fp);

      for(ap=actualmol->firstatom, j=0; ap; ap=ap->next, j++){
         if(!fgets(line, 255, fp)) return 0;
         if(line[19] != 'X')    return 0;
         sscanf(line+20, "%f%f%f%f%f",
            &vib->coord[j].x, &(vib+1)->coord[j].x, &(vib+2)->coord[j].x,
            &(vib+3)->coord[j].x, &(vib+4)->coord[j].x);

         if(!fgets(line, 255, fp)) return 0;
         if(line[19] != 'Y')    return 0;
         sscanf(line+20, "%f%f%f%f%f",
            &vib->coord[j].y, &(vib+1)->coord[j].y, &(vib+2)->coord[j].y,
            &(vib+3)->coord[j].y, &(vib+4)->coord[j].y);

         if(!fgets(line, 255, fp)) return 0;
         if(line[19] != 'Z')    return 0;
         sscanf(line+20, "%f%f%f%f%f",
            &vib->coord[j].z, &(vib+1)->coord[j].z, &(vib+2)->coord[j].z,
            &(vib+3)->coord[j].z, &(vib+4)->coord[j].z);

      }
   }


   return 1;
}

static int read_dipole(void)
{
   long fpos;
   float x, y, z;

   rewind(fp);
   if(!find_string("(DEBYE)")) return 0;
   fpos = ftell(fp);
   while(find_string("(DEBYE)")) fpos = ftell(fp);
   fseek(fp, fpos, SEEK_SET);
   fgets(line, 255, fp);

   sscanf(line, "%f %f %f", &x, &y, &z);
   actualmol->dipole = add_dipole(x, y, z);

   return 1;
}
