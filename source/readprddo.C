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


/* lecture of PRDDO-output */
#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "readprddo.h"
#include "general.h"
#include "maininterf.h"
#include "utils.h"
#include "box.h"

static int read_atomic_parameters(void);
static int read_atomic_coordinates(void);
static int read_energy_levels(void);
static int read_energy_levels2(void);
static int read_energy_levels3(void);
static int read_energy_levels4(void);
static int read_wavefunctions(void);
static int read_wavefunctions2(void);
static int read_occupations(void);
static int read_occupations2(void);
static int read_atomic_charges(void);
static void prddo_norm(void);
static char *find_string(char *s);

static FILE *fp;
static char line[256];
static long previous_line = 0, preprevious = 0;


/**** lecture of PRDDO output ****/

void read_prddo(char *file)
{
   char str[20];

   if((fp = fopen(file, "r")) == NULL){
      sprintf(line, "can't open file\n%s !", file);
      showinfobox(line);
      return;
   }

   if(!find_string("<<<<<< P R D D O / M >>>>>>")) {
      rewind(fp);
      if(!find_string("******P R D D O******")) {
         showinfobox("read PRDDO : PRDDO output only !");
         fclose(fp);
         return;
      }
   }

   add_mol(file);

   if(!read_atomic_parameters()){
      showinfobox("can't read the atomic parameters");
      fclose(fp);
      return;
   }

   if(!actualmol->natoms){
      sprintf(line, "No atoms in %s!", file);
      showinfobox(line);
      return;
   }

   create_bonds();
   find_multiplebonds();
//   create_box();
   new_mole(file);

   if(!read_energy_levels()){
      showinfobox("can't read the energy levels");
   }

   if(!read_wavefunctions()){
      showinfobox("can't read the wave_functions");
      rewind(fp);
   }
   else prddo_norm();

   if(!read_atomic_charges()){
      showinfobox("can't read the atomic charges");
   }

   fclose(fp);

   sprintf(str, "%d PRDDO wavefunctions", actualmol->nMolecularOrbitals);
   logprint(str);
   update_logs();

}




static char *find_string(char *s)
{
   previous_line = ftell(fp);
   for(;;) {
      if(!fgets(line, 255, fp)) return NULL;
      if(strstr(line, s)) return line;
      preprevious = previous_line;
      previous_line = ftell(fp);
   }
}



static char *getNextLine(void)
{
   char *retVal;

   retVal = fgets(line, 255, fp);
   if(retVal) {
      if(line[0] == '1') { /* page break for printer */
         (void)fgets(line, 255, fp); /* title */
         (void)fgets(line, 255, fp); /* blank */
         (void)fgets(line, 255, fp); /* blank */
         retVal = fgets(line, 255, fp); /* new data line */
      }
   }
   return retVal;
}



static int read_atomic_parameters(void)
{
   long fpos;
   AtoM *ap;
   Slater *sp;
   float x, y, z;
   float dum, c1, exp2, c2;
   int ord, offset, i;
   char str[30];

   if(!find_string("ATOM  NAME/CHG       X         Y         Z")) return 0;
   getNextLine();


   while(getNextLine()) {

      if(line[0] == '\n') break;

      for(i=strlen(line); i<256; i++) line[i] = 0;

      if(sscanf(line+18, "%f%f%f", &x, &y, &z) == 3) {
         sscanf(line+14, "%d", &ord);
         ap = add_atom(ord, x, y, z);
      }

      if(ord == 0) continue;
      offset = 49;
      if(!sscanf(line+offset, "%s", str)) return 0;
      offset = strstr(line, str) - line;    /* for D orbitals on new lines */
      sp = add_slater(ap);
      if(!isdigit(str[0])) return 0;
      sp->n = str[0] - '0';
      if(str[1] != 'S' && str[1] != 'P' && str[1] != 'D') return 0;
      sp->type[0] = str[1]; sp->type[1] = 0;
      if(!sscanf(line+offset+7, "%f%f", &sp->exponent, &c1)) return 0;
      if(sp->type[0] == 'D') { /* new version */
         fpos = ftell(fp);
         getNextLine();
         if(sscanf(line, "%f%f%f", &exp2, &c2, &dum) != 2)
            fseek(fp, fpos, SEEK_SET);
         else {
            sp->exponent = dszeta(sp->n, c1, sp->exponent, c2, exp2);
            continue;
         }
      }

      offset = 69;
      if(!line[offset]) continue;
      if(!sscanf(line+offset, "%s", str)) continue;
      sp = add_slater(ap);
      if(!isdigit(str[0])) return 0;
      sp->n = str[0] - '0';
      if(str[1] != 'S' && str[1] != 'P' && str[1] != 'D') return 0;
      sp->type[0] = str[1]; sp->type[1] = 0;
      if(!sscanf(line+offset+9, "%f", &sp->exponent)) return 0;

      offset = 92;
      if(!line[offset] || line[offset] == '\n') continue;
      if(!sscanf(line+offset, "%s", str)) continue;
      sp = add_slater(ap);
      if(!isdigit(str[0])) return 0;
      sp->n = str[0] - '0';
      if(str[1] != 'S' && str[1] != 'P' && str[1] != 'D') return 0;
      sp->type[0] = str[1]; sp->type[1] = 0;
      if(!sscanf(line+offset+9, "%f", &sp->exponent)) return 0;
      sscanf(line+109, "%f", &c1);
      if(sp->type[0] == 'D') { /* new version */
         fpos = ftell(fp);
         getNextLine();
         if(sscanf(line, "%f%f%f", &exp2, &c2, &dum) != 2)
            fseek(fp, fpos, SEEK_SET);
         else {
            sp->exponent = dszeta(sp->n, c1, sp->exponent, c2, exp2);
         }
      }
   }

   if(!find_string("BASIS FUNCTIONS")) return 0;
   sscanf(line + 11, "%d", &actualmol->nBasisFunctions);
   actualmol->nMolecularOrbitals = actualmol->nBasisFunctions;

   if(!find_string("COORDINATES ABOVE ARE IN")) return 0;
   if(strstr(line, "ATOMIC UNITS")) {
      for(ap = actualmol->firstatom; ap; ap = ap->next){
         ap->coord[0] *= BOHR;
         ap->coord[1] *= BOHR;
         ap->coord[2] *= BOHR;
      }
   }

   if(!read_atomic_coordinates()) rewind(fp);

   return 1;
}





static int read_atomic_coordinates(void)
{
   long fpos;
   int nbas, ord;
   float ford, x, y, z;
   AtoM *ap;

   if(find_string("ATOMIC COORDINATES FOR THIS CALCULATION")) {
      do {
         fpos = ftell(fp);
      } while (find_string("ATOMIC COORDINATES FOR THIS CALCULATION"));
   }
   else return 0;
   fseek(fp, fpos, SEEK_SET);

   if(!find_string("X         Y         Z")) return 0;

   ap = actualmol->firstatom;
   while(getNextLine()) {
      if(line[0] == '\n') continue;
      if(sscanf(line, "%d%f%f%f%f", &nbas, &ford, &x, &y, &z) != 5) break;
      if(!ap) return 0;               /* too many lines in file */
      ord = ford;
      if(ap->ord != ord) return 0;    /* incoherent file */
      ap->coord[0] = x*BOHR; ap->coord[1] = y*BOHR; ap->coord[2] = z*BOHR;
      ap = ap->next;
   }
   if(ap) return 0;                   /* not enough lines in file */

   return 1;
}

static int read_energy_levels(void)
{
   MolecularOrbital *alphaOrb;
   long fpos;
   int level, i, n_mo, n_electrons;
   float nel;
   char *cp;

   if(!find_string("EIGENVALUES (AU)......")) {
      return read_energy_levels2();
   }
   getNextLine();
   fpos = ftell(fp);
   n_mo = actualmol->nMolecularOrbitals;

   if((alphaOrb = allocOrbital(n_mo, 0, PRDDO_ORB)) == NULL){
      showinfobox("can't allocate the MO-structures");
      return 0;
   }

   level = 0;
   n_electrons = 0;
   while(getNextLine()) {
      if(line[0] == '\n') break;
      cp = line;
      for(i=0; cp=strchr(cp, '='); i++, level++) { /* x columns */
         if(level == n_mo) break;
         sscanf(++cp, " %lf%*c%f", &alphaOrb[level].eigenvalue, &nel);
         alphaOrb[level].occ = nel;
         n_electrons += nel;
      }
   }

   actualmol->nElectrons = n_electrons;
   actualmol->alphaOrbital = alphaOrb;

   return 1;
}


static int read_energy_levels2(void)
{
   MolecularOrbital *prddo_orb;
   int i, offset, n_electrons;
   float nel;
   char *cp;

   rewind(fp);
   if(!find_string(" EIGENVALUES (AU)")) {
      rewind(fp);
      return read_energy_levels4();
   }

   do getNextLine(); while(line[0] == '\n');
   if(strstr(line, "MO       ALPHA SPIN    BETA SPIN")) 
      return read_energy_levels3();

   if((prddo_orb = allocOrbital(actualmol->nMolecularOrbitals, 0, PRDDO_ORB)) == NULL){
      showinfobox("can't allocate the MO-structures");
      return 0;
   }

   n_electrons = 0;
   for(i=0; i<actualmol->nMolecularOrbitals; i++) {
      if((cp = strchr(line, '=')) == NULL) return 0;
      offset = cp - line;
      sscanf(++cp, " %lf", &prddo_orb[i].eigenvalue);
      if((cp = strchr(cp, '(')) == NULL) return 0;
      sscanf(++cp, "%f", &nel);
      prddo_orb[i].occ = nel;
      n_electrons += nel;
      getNextLine();
   }

   actualmol->nElectrons = n_electrons;
   actualmol->alphaOrbital = prddo_orb;

   return 1;
}



static int read_energy_levels3(void)
{
   int i, n_orbs, dum;
   MolecularOrbital *a_orb, *b_orb;

   n_orbs = actualmol->nMolecularOrbitals;
   if((a_orb = allocOrbital(n_orbs, 0, PRDDO_ORB)) == NULL){
      showinfobox("can't allocate the alpha MO-structures");
      return 0;
   }
   if((b_orb = allocOrbital(n_orbs, 0, PRDDO_ORB)) == NULL){
      showinfobox("can't allocate the beta MO-structures");
      return 0;
   }

   getNextLine();

   for(i=0; i<n_orbs; i++) {
      getNextLine();
      sscanf(line, "%d%lf%lf", &dum,
             &a_orb[i].eigenvalue, &b_orb[i].eigenvalue);

      if(dum != i+1) return 0;
   }

   actualmol->alphaOrbital = a_orb;
   actualmol->betaOrbital  = b_orb;
   return 1;
}


static int read_energy_levels4(void)
{
   int i, n_orbs, level, nel, n_alpha = 0, n_beta = 0;
   MolecularOrbital *a_orb, *b_orb;
   char *cp;

   rewind(fp);
   if(!find_string(" ALPHA OCCUPATIONS")) {
      return 0;
   }

   n_orbs = actualmol->nMolecularOrbitals;
   if((a_orb = allocOrbital(n_orbs, 0, PRDDO_ORB)) == NULL){
      showinfobox("can't allocate the alpha MO-structures");
      return 0;
   }
   if((b_orb = allocOrbital(n_orbs, 0, PRDDO_ORB)) == NULL){
      showinfobox("can't allocate the beta MO-structures");
      return 0;
   }

   fgets(line, 255, fp);
   for(i=0; i<n_orbs; i++) {
      sscanf(line+i, "%1d", &nel);
      a_orb[i].occ = nel;
      n_alpha += nel;
   }

   if(!find_string(" BETA OCCUPATIONS")) {
      return 0;
   }

   fgets(line, 255, fp);
   for(i=0; i<n_orbs; i++) {
      sscanf(line+i, "%1d", &nel);
      b_orb[i].occ = nel;
      n_beta += nel;
   }

   if(!find_string(" EIGENVALUES IN AU FOR ALPHA SPIN")) {
      return 0;
   }

   getNextLine();

   level = 0;
   while(getNextLine()) {
      if(line[0] == '\n') break;
      cp = line;
      for(i=0; cp=strchr(cp, '='); i++, level++) { /* x columns */
         if(level == n_orbs) break;
         sscanf(++cp, " %lf", &a_orb[level].eigenvalue);
      }
   }

   if(!find_string(" EIGENVALUES IN AU FOR BETA SPIN")) {
      return 0;
   }

   getNextLine();

   level = 0;
   while(getNextLine()) {
      if(line[0] == '\n') break;
      cp = line;
      for(i=0; cp=strchr(cp, '='); i++, level++) { /* x columns */
         if(level == n_orbs) break;
         sscanf(++cp, " %lf", &b_orb[level].eigenvalue);
      }
   }
   actualmol->alphaOrbital = a_orb;
   actualmol->betaOrbital  = b_orb;
   actualmol->nAlpha = n_alpha;
   actualmol->nBeta = n_beta;
   actualmol->nElectrons = n_alpha + n_beta;
   return 1;

}


static int read_wavefunctions(void)
{
   MolecularOrbital *prddo_orb;
   int n, nblocks, dum, n_mo;
   register int i, j;

   if(!find_string("CANONICAL MOS (BY COLUMNS, OVER SLATER AOS)")) {
      return read_wavefunctions2();
   }

   n_mo = actualmol->nMolecularOrbitals;

   if(!actualmol->alphaOrbital) {
      if((prddo_orb = allocOrbital(n_mo, 0, PRDDO_ORB)) == NULL){
         showinfobox("can't allocate the MO-structures");
         return 0;
      }
   }
   else prddo_orb = actualmol->alphaOrbital;

   nblocks = n_mo/12;
   if(n_mo%12) nblocks++;

   for(i=0; i<n_mo; i++){
      if((prddo_orb[i].coefficient = (double *)malloc(n_mo*sizeof(double))) == NULL){
         showinfobox("can't allocate the MO-structures");
         return 0;
      }
   }

   for(j=n=0; j<nblocks; j++, n+=12){
      do getNextLine(); while (line[0] == '\n');
      getNextLine(); /* empty line (one space!) */

      for(i=0; i<n_mo; i++){
         getNextLine();
         sscanf(line, " %d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
            &dum,
            &prddo_orb[n].coefficient[i],    &prddo_orb[n+1].coefficient[i], 
            &prddo_orb[n+2].coefficient[i],  &prddo_orb[n+3].coefficient[i], 
            &prddo_orb[n+4].coefficient[i],  &prddo_orb[n+5].coefficient[i], 
            &prddo_orb[n+6].coefficient[i],  &prddo_orb[n+7].coefficient[i], 
            &prddo_orb[n+8].coefficient[i],  &prddo_orb[n+9].coefficient[i], 
            &prddo_orb[n+10].coefficient[i], &prddo_orb[n+11].coefficient[i]);
         if(dum != i+1) { fprintf(stderr, "out of sync 1\n"); return 0; }
      }
/*
      if(!getNextLine() && j == nblocks-1) return 1;
      if(line[0] != '\n'){ fprintf(stderr, "out of sync 2\n"); return 0; }
*/
   }

   rewind(fp);
   actualmol->alphaOrbital = prddo_orb;
   if(!actualmol->nElectrons) read_occupations();

   return 1;
}





static int read_wavefunctions2(void)
{
   MolecularOrbital *prddo_orb;
   int n, nblocks, dum, norbs;
   register int i, j;

   rewind(fp);

   if(!find_string("EIGENVECTORS (BY COLUMNS)")) return 0;
   if(strstr(line, "FOR ALPHA-SPIN MOS")) actualmol->alphaBeta = 1;
   else actualmol->alphaBeta = 0;

   norbs = actualmol->nMolecularOrbitals;

   prddo_orb = actualmol->alphaOrbital;
   if(!prddo_orb) {
      if((prddo_orb = allocOrbital(norbs, 0, PRDDO_ORB)) == NULL){
         showinfobox("can't allocate the MO-structures");
         return 0;
      }
   }
   for(i=0; i<norbs; i++) {
      if((prddo_orb[i].coefficient = (double *)malloc(norbs*sizeof(double))) == NULL) {
         showinfobox("can't allocate the MO-structures");
         return 0;
      }
   }

   nblocks = actualmol->nMolecularOrbitals/12;
   if(norbs%12) nblocks++;

   for(j=n=0; j<nblocks; j++, n+=12){
      do getNextLine(); while (line[0] == '\n');
      getNextLine(); /* empty line (one space!) */

      for(i=0; i<norbs; i++){
         getNextLine();
         sscanf(line, " %d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
            &dum,
            &prddo_orb[n].coefficient[i],    &prddo_orb[n+1].coefficient[i], 
            &prddo_orb[n+2].coefficient[i],  &prddo_orb[n+3].coefficient[i], 
            &prddo_orb[n+4].coefficient[i],  &prddo_orb[n+5].coefficient[i], 
            &prddo_orb[n+6].coefficient[i],  &prddo_orb[n+7].coefficient[i], 
            &prddo_orb[n+8].coefficient[i],  &prddo_orb[n+9].coefficient[i], 
            &prddo_orb[n+10].coefficient[i], &prddo_orb[n+11].coefficient[i]);
         if(dum != i+1) { fprintf(stderr, "out of sync 1\n"); return 0; }
      }

      if(!getNextLine() && j == nblocks-1) return 1;
/*
      if(line[0] != '\n'){ fprintf(stderr, "out of sync 2\n"); return 0; }
*/
   }

   rewind(fp);
   actualmol->alphaOrbital = prddo_orb;

   if(!actualmol->alphaBeta) {
      if(!actualmol->nElectrons) read_occupations();
      return 1;
   }

   if(!find_string("EIGENVECTORS (BY COLUMNS) FOR BETA -SPIN MOS")) return 0;

   prddo_orb = actualmol->betaOrbital;
   if(!prddo_orb) {
      if((prddo_orb = allocOrbital(norbs, norbs, PRDDO_ORB)) == NULL){
         showinfobox("can't allocate the MO-structures");
         return 0;
      }
   }
   for(i=0; i<norbs; i++) {
      if((prddo_orb[i].coefficient = (double *)malloc(norbs*sizeof(double))) == NULL) {
         showinfobox("can't allocate the MO-structures");
         return 0;
      }
   }

   for(j=0, n=0; j<nblocks; j++, n+=12){
      do getNextLine(); while (line[0] == '\n');
      getNextLine(); /* empty line (one space!) */

      for(i=0; i<norbs; i++){
         getNextLine();
         sscanf(line, " %d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
            &dum,
            &prddo_orb[n].coefficient[i],    &prddo_orb[n+1].coefficient[i], 
            &prddo_orb[n+2].coefficient[i],  &prddo_orb[n+3].coefficient[i], 
            &prddo_orb[n+4].coefficient[i],  &prddo_orb[n+5].coefficient[i], 
            &prddo_orb[n+6].coefficient[i],  &prddo_orb[n+7].coefficient[i], 
            &prddo_orb[n+8].coefficient[i],  &prddo_orb[n+9].coefficient[i], 
            &prddo_orb[n+10].coefficient[i], &prddo_orb[n+11].coefficient[i]);
         if(dum != i+1) { fprintf(stderr, "out of sync 1\n"); return 0; }
      }

      if(!getNextLine() && j == nblocks-1) return 1;
/*
      if(line[0] != '\n'){ fprintf(stderr, "out of sync 2\n"); return 0; }
*/
   }
   actualmol->betaOrbital = prddo_orb;

   if(!actualmol->nElectrons) read_occupations();

   return 1;
}



static int read_occupations(void)
{
   char *cp;
   int i;

   rewind(fp);
   if(!find_string(" OCCUPATIONS ")) {
      rewind(fp);
      return read_occupations2();
   }

   for(cp=line, i=0; *cp; cp++) {
      switch(*cp) {
         case '0' :
            i++;
            break;
         case '1' :
            actualmol->alphaOrbital[i].occ = 1;
            actualmol->nElectrons++;
            i++;
            break;
         case '2' :
            actualmol->alphaOrbital[i].occ = 2;
            actualmol->nElectrons += 2;
            i++;
            break;
      }
   }
   actualmol->nAlpha = actualmol->nElectrons;

   if(!actualmol->alphaBeta) return 1;

   if(!find_string("BETA OCCUPATIONS")) {
      rewind(fp);
      return 0;
   }

   for(cp=line, i=0; *cp; cp++) {
      switch(*cp) {
         case '0' :
            i++;
            break;
         case '1' :
            actualmol->betaOrbital[i].occ = 1;
            actualmol->nElectrons++;
            i++;
            break;
         case '2' :
            actualmol->betaOrbital[i].occ = 2;
            actualmol->nElectrons += 2;
            i++;
            break;
      }
   }
   actualmol->nBeta = actualmol->nElectrons - actualmol->nAlpha;

   return 1;
}



static int read_occupations2(void)
{
   char *cp;
   int norbs, ncharge, i;

   if(!find_string("NUMBER OF DOUBLY OCCUPIED ORBITALS")) {
      rewind(fp);
      return 0;
   }

   cp = strstr(line, "ORBITALS");
   sscanf(cp+8, "%d", &norbs);
   cp = strstr(line, "IONIC CHARGE");
   if(!cp) return 0;
   sscanf(cp+12, "%d", &ncharge);
   actualmol->charge = ncharge;

   actualmol->nElectrons = 2*norbs;

   for(i=0; i<actualmol->nMolecularOrbitals; i++) {
      actualmol->alphaOrbital[i].occ = 2*(i<norbs);
   }

   return 1;
}





static int read_atomic_charges(void)
{
   int nitems, num;
   float valency; /* not used */
   AtoM *ap;

   actualmol->charges = 0;
   if(find_string("ATOM      ATOM CHG")) {
      getNextLine();
      getNextLine();

      for(ap=actualmol->firstatom; ap; ap=ap->next){
         if(ap->ord == 0) continue;
         getNextLine();
         nitems = sscanf(line, "%d%f%f", &num, &ap->charge, &valency);
         if(nitems != 3) {
            showinfobox("read_atomic_charges: not enough entries!");
           return 0;
         }
      }
      actualmol->charges = 1;
   }
   else {
      rewind(fp);
   }

   if(find_string("ATOMIC CHARGES")) {
      int i, dum;
      char *cp;

      getNextLine();  /* empty line */
      getNextLine();
      if(strchr(line, '=')) { /* format with several entries per line */
         cp = line;
         for(ap=actualmol->firstatom, i=0; ap; ap=ap->next, i++){
            cp = strchr(cp, '=');
            if(!cp) {
               getNextLine();
               cp = strchr(line, '=');
               if(!cp) return 0;
            }
            nitems = sscanf(++cp, "%f", &ap->charge);

            if(nitems != 1) {
               showinfobox("read_atomic_charges: not enough entries!");
               return 0;
            }
         }
      }
      else {    /* format with one entry per line */
         for(ap=actualmol->firstatom, i=0; ap; ap=ap->next, i++){
            if(ap->ord != 0) {
               nitems = sscanf(line, "%d%f", &dum, &ap->charge);
               if(dum != ap->name) return 0;
               if(nitems != 2) return 0;
            }
            getNextLine();
         }
      }

      actualmol->charges = 1;
   }

   return actualmol->charges;
}

#define POW3(x)  ((x)*(x)*(x))
#define POW5(x)  ((x)*(x)*(x)*(x)*(x))
#define POW7(x)  ((x)*(x)*(x)*(x)*(x)*(x)*(x))
#define POW9(x)  ((x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x))
#define POW11(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x))


static void prddo_norm(void)
{
   AtoM *ap;
   Slater *vp;

   for(ap = actualmol->firstatom; ap; ap = ap->next){
      for(vp = ap->firstslater; vp; vp = vp->next) {
 
      switch(vp->type[0]) {
         case 'S' :
            switch(vp->n) {
               case 1: vp->norm[0] = sqrt(POW3(vp->exponent) * M_1_PI); break;
               case 2: vp->norm[0] = sqrt(POW5(vp->exponent) * M_1_PI / 3.); break;
               case 3: vp->norm[0] = sqrt(POW7(vp->exponent) * M_2_PI / 45.); break;
               case 4: vp->norm[0] = sqrt(POW9(vp->exponent) * M_1_PI / 315.); break;
               case 5: vp->norm[0] = sqrt(POW11(vp->exponent) * M_2_PI / 14175.); break;
            }
            break;
         case 'P' :
            switch(vp->n) {
               case 2: vp->norm[1] = sqrt(POW5(vp->exponent) * M_1_PI); break;
               case 3: vp->norm[1] = sqrt(POW7(vp->exponent) * M_2_PI / 15.); break;
               case 4: vp->norm[1] = sqrt(POW9(vp->exponent) * M_1_PI / 105.); break;
               case 5: vp->norm[1] = sqrt(POW11(vp->exponent) * M_2_PI / 4725.); break;
            }
            break;
         case 'D' :
            switch(vp->n) {
               case 3: vp->norm[2] = sqrt(POW7(vp->exponent) * M_1_PI / 6.);
                       vp->norm[3] = sqrt(POW7(vp->exponent) * M_1_PI / 18.);
                       vp->norm[4] = sqrt(POW7(vp->exponent) * M_2_PI / 3.);
                       break;
               case 4: vp->norm[2] = sqrt(POW9(vp->exponent) * M_1_PI / 84.);
                       vp->norm[3] = sqrt(POW9(vp->exponent) * M_1_PI / 252.);
                       vp->norm[4] = sqrt(POW9(vp->exponent) * M_1_PI / 21.);
                       break;
               case 5: vp->norm[2] = sqrt(POW11(vp->exponent) * M_1_PI / 1890.);
                       vp->norm[3] = sqrt(POW11(vp->exponent) * M_1_PI / 5670.);
                       vp->norm[4] = sqrt(POW11(vp->exponent) * M_2_PI / 945.);
                       break;
            }
            break;
      }
      }
   }

   return;
}

