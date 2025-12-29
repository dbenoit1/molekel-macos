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


/* reading GAUSSIAN output G92 - G98 should work */
#include "main.h"
#include "constant.h"
#include "molekel.h"
#include "general.h"
#include "maininterf.h"
#include "readgauss.h"
#include "box.h"
#include "utils.h"

#define RHF          1
#define ROHF         2
#define UHF          4

void read_gauss(char *name);
int read_atomic_coordinates(char *file);
char *find_string(char *s);
int read_charge(void);
int read_basis_set(void);
int  read_coefficients(void);
int read_density_matrix(void);
float **read_trimat(void);
void read_coeffs(char *s, double *v1, double *v2, double *v3,
                                double *v4, double *v5);
double get_value(char *s, int n);
void all_uppercase(char *s);
int read_frequencies(void);
int read_dipole(void);
void print_frequencies(void);
int read_atomic_charges(void);
int addTrajectoryStep(void);

void print_basis_set(void);
void print_coefficients(void);
void print_density_matrix(void);

int n_primitive_gaussians;

static FILE *fp;
static char line[256];
static int nblocks, d_type, f_type;
static long previous_line = 0, preprevious = 0, preprepre = 0;
unsigned short flagG98 = 1;


/**** lecture of gaussian output ****/

void read_gauss(char *name)
{
   unsigned long position;
   int basisread = 1;

   if((fp = fopen(name, "r")) == NULL){
      sprintf(line, "read_gauss : can't open %s\n", name);
      showinfobox(line);
      return;
   }
   if(!find_string("Gaussian")) {
      showinfobox("read_gauss : Gaussian output only !");
      fclose(fp);
      return;
   }

   if(!find_string("Gaussian 98")) flagG98 = 0;
   else {
      flagG98 = 1;
      logprint("G98");
   }
   rewind (fp);

   if(!read_atomic_coordinates(name)){
      showinfobox("can't read the atomic coordinates");
      fclose(fp);
      return;
   }

   if(!actualmol->natoms){
      sprintf(line, "No atoms in\n%s !", name);
      showinfobox(line);
      return;
   }

   create_bonds();
   find_multiplebonds();
//   create_box();
   new_mole(name);


   if(!read_atomic_charges()){
      logprint("unable to read");
      logprint("the atomic charges");
      update_logs();
      actualmol->charges = 0;
   }

   rewind(fp);

   if(!find_string("----")) return;
   fgets(line, 255, fp);
   if(line[1] != '#'){
      while(find_string("----")) {
         fgets(line, 255, fp);
         if(line[1] == '#') break;
      }
   }
   all_uppercase(line);
   logprint(line);
   if(strlen(line) > 29) logprint(line+29);
   if(strlen(line) > 58) logprint(line+58);

   if(!strstr(line, "GFPRINT") || !strstr(line, "POP")) {
      logprint("");
      logprint("for visualizing orbitals");
      logprint("and densities, you must");
      logprint("run the GAUSSIAN-job");
      logprint("with the keywords");
      logprint("GFPRINT,");
      logprint("  POP = FULL or POP = REGULAR");
      logprint("");
   }

   find_string("-\n");
   if(!find_string("-\n")) {
      logprint("can't read the title");
   }   
   else {
      fgets(line, 255, fp);
      logprint("title:");
      logprint(line);     /* the title */
   }
   update_logs();

   position = ftell(fp);
   if(!read_charge()){
      logprint("can't read the charge");
      logprint("  and multiplicity");
      fseek(fp, position, SEEK_SET);
   }

   position = ftell(fp);
   if(!read_basis_set()){
      logprint("can't read the basis-set");
      fseek(fp, position, SEEK_SET);
      basisread = 0;
   }

/*
   print_basis_set();
*/

   position = ftell(fp);
   if(basisread) {
      if(!read_coefficients()){
         logprint("can't read");
         logprint("  the MO-coefficients");
         fseek(fp, position, SEEK_SET);
      }
      computeOccupations(actualmol);


/*
   print_coefficients();
*/
      position = ftell(fp);
      if(!read_density_matrix()){
         logprint("can't read");
         logprint("  the density matrix");
         fseek(fp, position, SEEK_SET);
      }
/*
   print_density_matrix();
*/
   }

   free_frequencies();
   if(read_frequencies()){
      logprint("frequencies present");
/*
      print_frequencies();
*/
   }

   if(read_dipole()) {
      logprint("dipole moment present");
   }

   fclose(fp);
   if(maininterfwin) update_logs();

}


char *find_string(char *s)
{
   previous_line = ftell(fp);
   do {
      if(!fgets(line, 255, fp)) return NULL;
      if(strstr(line, s)) return line;
      preprepre = preprevious;
      preprevious = previous_line;
      previous_line = ftell(fp);
   } while (1);
}


int read_atomic_coordinates(char *file)
{
   long fpos;
   float x, y, z;
   int ord, natoms;

   free_dyna();

   if(find_string("Standard orientation")) {
      do {
         fpos = ftell(fp);
         natoms = addTrajectoryStep();
         } while(find_string("Standard orientation"));
   }
   else {
      rewind(fp);
      if(find_string("Z-Matrix orientation")) {
         do {
            fpos = ftell(fp);
            natoms = addTrajectoryStep();
         } while(find_string("Z-Matrix orientation"));
      }
      else {
           rewind(fp);
           if(find_string("Input orientation")) {
                do {
                 fpos = ftell(fp);
                 natoms = addTrajectoryStep();
                } while(find_string("Input orientation"));
           }
      }
   }

      add_mol(file);
      dynamics.molecule = actualmol;
      dynamics.current = dynamics.ntotalsteps - 1;
 
   fseek(fp, fpos, SEEK_SET);

   if(!find_string("Coordinates (Angstroms)")) return 0;
   if(!find_string("-----")) return 0;

   fgets(line, 255, fp);
   do {
      if(flagG98) {
      if(sscanf(line, "%*d %d %*d %f %f %f", &ord, &x, &y, &z) != 4) return 0;
      }
      else {
      if(sscanf(line, "%*d %d %f %f %f", &ord, &x, &y, &z) != 4) return 0;
      }
      if(ord >= 0) add_atom(ord, x, y, z);
      fgets(line, 255, fp);
   } while(!strstr(line, "------"));
   

/*
   {

   AtoM *ap;
   FILE *pdb;

   if((pdb = fopen("temp.pdb", "w")) == NULL){
      screenprint("can't open temp.pdb for writing");
      return 0;
   }

   for(ap = actualmol->firstatom; ap; ap = ap->next){
      fprintf(pdb, "ATOM      1  %2s  pff     1    %8.3f%8.3f%8.3f\n",
              element[ap->ord].symbol, ap->coord[0], ap->coord[1], ap->coord[2]);
   }

   fclose(pdb);

   }
*/
   return 1;
}




int addTrajectoryStep(void)
{
   long fpos, index, i;
   float x, y, z;
   static int natoms;
   char str[30];

   index = dynamics.ntotalsteps++;

   if(!dynamics.trajectory){
      if((dynamics.trajectory = (Vector**) malloc(sizeof(Vector *))) == NULL){
         showinfobox("can't allocate dyna pointer\n");
         return 0;
      }

/* count nr of atoms */
      fpos = ftell(fp);
      if(!find_string("Coordinates (Angstroms)")) return 0;
      if(!find_string("-----")) return 0;
      fgets(line, 255, fp);
      natoms = 0;
      do {
         natoms++;
         if(!fgets(line, 255, fp)) break;
      } while(!strstr(line, "------"));
      fseek(fp, fpos, SEEK_SET);

   }
   else {
      if((dynamics.trajectory = (Vector**) realloc(dynamics.trajectory,
          dynamics.ntotalsteps*sizeof(Vector *))) == NULL){
         showinfobox("can't reallocate dyna pointer\n");
         dynamics.ntotalsteps--;
         free_dyna();
         return 0;
      }
      fpos = ftell(fp);
   }

   if((dynamics.trajectory[index] = (Vector*) malloc(natoms*sizeof(Vector))) == NULL){
      sprintf(str, "can't allocate timestep dyna[%d]\n", index);
      showinfobox(str);
      dynamics.ntotalsteps--;
      free_dyna();
      return 0;
   }

   if(!find_string("Coordinates (Angstroms)")) return 0;
   if(!find_string("-----")) return 0;

   fgets(line, 255, fp);
   i = 0;
   do {
      if(flagG98) {
      if(sscanf(line, "%*d %*d %*d %f %f %f", &x, &y, &z) != 3) return 0;
      } 
      else {
      if(sscanf(line, "%*d %*d %f %f %f", &x, &y, &z) != 3) return 0;
      }
      dynamics.trajectory[index][i].x = x;
      dynamics.trajectory[index][i].y = y;
      dynamics.trajectory[index][i].z = z;
      if(!fgets(line, 255, fp)) break;
      i++;
   } while(!strstr(line, "------"));

   fseek(fp, fpos, SEEK_SET);
   return natoms;
}




int read_charge(void)
{
   char str[30];
   if(!find_string("Multiplicity =")) return 0;

   sscanf(line, "%*s %*c%f %*s %*c%d", &actualmol->charge, &actualmol->multiplicity);
   sprintf(str, "charge: %3.1f, multiplicity: %2d", actualmol->charge, actualmol->multiplicity);
   logprint(str);

   return 1;
}



int read_basis_set(void)
{
   AtoM *ap;
   Gauss *gp;
   Shell *sp;
   double s_coeff, p_coeff, d_coeff, f_coeff, alpha, norm;
   char *cp;
   char tmpline[256];
   long fpos;

   if(!find_string(" Basis read")) {
      rewind(fp);
      if(!find_string(" basis")) return 0;

   }
   fpos = ftell(fp);
   fgets(tmpline, 255, fp);
   if(!strstr(tmpline, " ******")) {
      rewind(fp);
       if(!find_string(" basis")) return 0;
   }
   else fseek(fp, fpos, SEEK_SET);
   
   if(!strstr(line, "General basis") && 
      !strstr(line, "Standard basis") &&
      !strstr(line, "from chk"))
      return 0;
   logprint(line);
   if(strlen(line) > 29) logprint(line+29);

   if(!(cp = strrchr(line, '('))) return 0;
   if(!(cp = strchr(cp, 'D'))) return 0;
   cp--;
   if(*cp == '5') d_type = 5;
   else if(*cp == '6') d_type = 6;
   else {
      showinfobox("weird d-orbital");
      return 0;
   }
   if(!(cp = strchr(cp, 'F'))) return 0;
   cp--;
   if(*cp == '7') f_type = 7;
   else if(*cp == '0') f_type = 10;
   else {
      showinfobox("weird f-orbital");
      return 0;
   }

   if(!find_string("GAUSSIAN FUNCTIONS")) return 0;
   if(!find_string("*    ATOM   X-COORD  Y-COORD  Z-COORD   *  NUMBER     TYPE\
     FACTOR *  EXPONENT    S-COEF      P-COEF      D-COEF      F-COEF   *"))
      return 0;
   if(!find_string("****")) return 0;

   for(ap = actualmol->firstatom; ap; ap = ap->next){
      fgets(line, 255, fp); /* atom - symbol, x, y, z */
      fgets(line, 255, fp);
      do {
         if(!(sp = add_shell(ap))) return 0;
/*         sscanf(line+52, "%s", sp->type); */
         sscanf(line+63, "%f", &sp->scale_factor);
         fgets(line, 255, fp);
         while(strstr(line, "                                             \
                       ")) {
            gp = add_gauss(sp);
            alpha = 0.0; s_coeff = 0.0; p_coeff = 0.0;
            d_coeff = 0.0; f_coeff = 0.0;
            read_coeffs(line+2, &alpha, &s_coeff, &p_coeff,
                                               &d_coeff, &f_coeff);
            gp->exponent = alpha * sp->scale_factor * sp->scale_factor;
            alpha = gp->exponent;
            if(s_coeff){
               norm  = pow(2.0 * alpha / M_PI, 0.75);
               gp->coeff  = s_coeff * norm;
               gp->coeff2 = 0.0;
               sp->n_base = 1;
               if(p_coeff){
                  norm  = pow(128.0 * pow(alpha, 5) / pow(M_PI, 3), 0.25);
                  gp->coeff2 = p_coeff * norm;
                  sp->n_base = 4;
               }
            }
            else if(p_coeff){
               norm  = pow(128.0 * pow(alpha, 5) / pow(M_PI, 3), 0.25);
               gp->coeff  = p_coeff * norm;
               gp->coeff2  = 0.0;
               sp->n_base = 3;
            }
            else if(d_coeff){
               norm = pow(2048. * pow(alpha, 7) / pow(M_PI, 3), .25);
               gp->coeff  = d_coeff * norm;
               gp->coeff2  = 0.0;
               sp->n_base = d_type;
            }
            else if(f_coeff){
               norm = pow(32768. * pow(alpha, 9) / pow(M_PI, 3), .25);
               gp->coeff = f_coeff * norm;
               gp->coeff2  = 0.0;
               sp->n_base = f_type;
            }
            else {
               showinfobox("No coefficients in gaussian primitive");
               return 0;
            }

            if(!fgets(line, 255, fp)) return 0;
         }                                /* end of gaussian primitives */
      } while(!strstr(line, "-----") && 
              !strstr(line, "*****"));    /* end of atom */
   }                                      /* end of basis-set */

   if(!strstr(line, "****")) return 0;

   normalize_gaussians();

   return 1; 
}



void read_coeffs(char *s, double *v1, double *v2, double *v3,
                                double *v4, double *v5)
/* read FORTRAN-type double-precision numbers */
{
   while(*s == ' ') s++;

   *v1 = get_value(s, 12);    s += 12;
   *v2 = get_value(s, 12);    s += 12;
   *v3 = get_value(s, 12);    s += 12;
   *v4 = get_value(s, 12);    s += 12;
   *v5 = get_value(s, 12);    s += 12;
}



double get_value(char *s, int n)
{
   char digits[20];
   register char *letter;
   double v;

   strncpy(digits, s, n);
   digits[n] = 0;
   for(letter=digits; *letter; letter++){
      if(*letter == 'D') *letter = 'e';
   }

   if(sscanf(digits, "%le", &v) != 1){
      showinfobox("problems reading the coefficients...");
      return 0;
   }
   return v;
}



void print_basis_set(void)
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



void all_uppercase(char *s)
{
   while(*s) {
      *s = toupper(*s);
      s++;
   }
}



int read_coefficients(void)
{
   long fpos, denspos, coefficientsPos;
   int i1, i2, i3, i4, i5, n_mo, norbs;
   int firstOrb, firstPass;
   register int i, j;
   static char *keystr = "Alpha Molecular Orbital Coefficients", *pkey;
   MolecularOrbital *alphaOrb, *betaOrb;

   if(!find_string("primitive gaussians")) return 0;
   sscanf(line, "%d %*s %*s %d", &actualmol->nBasisFunctions, &n_primitive_gaussians);
   sprintf(line, "%3d basis functions", 
          actualmol->nBasisFunctions);
   logprint(line);
   sprintf(line, "  and %3d primitive gaussians", 
          n_primitive_gaussians);
   logprint(line);
   fgets(line, 255, fp);
   sscanf(line, "%d %*s %*s %d", &actualmol->nAlpha, &actualmol->nBeta);
   logprint(line);

   pkey = strchr(keystr, 'O');
   if(!find_string(pkey)) return 0;
   if(strstr(line, keystr)) {
      actualmol->alphaBeta = 1;
      pkey = keystr;
   }

   do {
      fpos = ftell(fp);
   } while(find_string(pkey));
   fseek(fp, fpos, SEEK_SET);
   coefficientsPos = fpos;

   find_string("DENSITY MATRIX.");
   denspos = ftell(fp);
   fseek(fp, fpos, SEEK_SET);

   nblocks = 0;
   i1 = i2 = i3 = i4 = i5 = 0;
   n_mo = firstPass = 0;
   while(find_string("EIGENVALUES")){
      if(ftell(fp) >= denspos) break;
      nblocks++;
      
 	if(actualmol->firstOrbital == 0) {
	    fseek(fp, coefficientsPos, SEEK_SET);
	    fgets(line, 255, fp); /* contains the orbital numbers */
	}
	else {
	    fseek(fp, preprevious, SEEK_SET);
	    fgets(line, 255, fp);   /* contains the orbital numbers or the symmetry */
	    if(strchr(line, '(') || strchr(line, 'O') || strchr(line, 'V')) {
	    /* contains symmetry symbols */
		fseek(fp, preprepre, SEEK_SET);
	        fgets(line, 255, fp); /* contains the orbital numbers */
	    }
	}

      n_mo += sscanf(line, " %d %d %d %d %d", &i1, &i2, &i3, &i4, &i5);
      if(!firstPass) actualmol->firstOrbital = firstPass = i1;
      find_string("EIGENVALUES"); /* "EIGENVALUES" string */
   }
   actualmol->lastOrbital = actualmol->firstOrbital + n_mo;

   norbs = (actualmol->alphaBeta ? n_mo/2 : n_mo);

   if((alphaOrb = allocOrbital(norbs, actualmol->nBasisFunctions, GAUSS_ORB)) == NULL){
      showinfobox("can't allocate the MO-structures");
      return 0;
   }

   if(actualmol->alphaBeta) {
      if((betaOrb = allocOrbital(norbs, actualmol->nBasisFunctions, GAUSS_ORB)) == NULL){
         showinfobox("can't allocate the MO-structures");
         return 0;
      }
   }

   fseek(fp, fpos, SEEK_SET);

   firstOrb = actualmol->firstOrbital;
   j = nblocks;
   if(actualmol->alphaBeta) j = nblocks/2;
   while(j--){
      fgets(line, 255, fp);
      sscanf(line, " %d %d %d %d %d", &i1, &i2, &i3, &i4, &i5);
      find_string("EIGENVALUES");
      sscanf(line+20, "%lf %lf %lf %lf %lf", &alphaOrb[i1-firstOrb].eigenvalue, 
          &alphaOrb[i2-firstOrb].eigenvalue, &alphaOrb[i3-firstOrb].eigenvalue,
          &alphaOrb[i4-firstOrb].eigenvalue, &alphaOrb[i5-firstOrb].eigenvalue);
      for(i=0; i<actualmol->nBasisFunctions; i++){
         fgets(line, 255, fp);
         sscanf(line+20, "%lf %lf %lf %lf %lf", &alphaOrb[i1-firstOrb].coefficient[i], 
         &alphaOrb[i2-firstOrb].coefficient[i], &alphaOrb[i3-firstOrb].coefficient[i],
         &alphaOrb[i4-firstOrb].coefficient[i], &alphaOrb[i5-firstOrb].coefficient[i]);
      }
   }

   if(!actualmol->alphaBeta) {
      actualmol->alphaOrbital = alphaOrb;
      actualmol->nMolecularOrbitals = n_mo;
      return 1;
   }

   if(!find_string("Beta Molecular Orbital Coefficients")) return 0;
   j = nblocks/2;
   while(j--){
      fgets(line, 255, fp);
      sscanf(line, " %d %d %d %d %d", &i1, &i2, &i3, &i4, &i5);
      find_string("EIGENVALUES");
      sscanf(line+20, "%lf %lf %lf %lf %lf", &betaOrb[i1-firstOrb].eigenvalue, 
           &betaOrb[i2-firstOrb].eigenvalue, &betaOrb[i3-firstOrb].eigenvalue,
           &betaOrb[i4-firstOrb].eigenvalue, &betaOrb[i5-firstOrb].eigenvalue);
      for(i=0; i<actualmol->nBasisFunctions; i++){
         fgets(line, 255, fp);
         sscanf(line+20, "%lf %lf %lf %lf %lf", &betaOrb[i1-firstOrb].coefficient[i], 
          &betaOrb[i2-firstOrb].coefficient[i], &betaOrb[i3-firstOrb].coefficient[i],
          &betaOrb[i4-firstOrb].coefficient[i], &betaOrb[i5-firstOrb].coefficient[i]);
      }
   }

   actualmol->alphaOrbital = alphaOrb;
   actualmol->betaOrbital = betaOrb;
   actualmol->nMolecularOrbitals = norbs;
   return 1;
}


/* reads the alpha- (and beta-) matrices */
int read_density_matrix(void)
{
   if(!find_string("DENSITY MATRIX.")) return 0;

   if(actualmol->alphaDensity) {
      free(actualmol->alphaDensity[0]);
      free(actualmol->alphaDensity);
   }
   if(actualmol->betaDensity) {
      free(actualmol->betaDensity[0]);
      free(actualmol->betaDensity);
   }
   actualmol->alphaDensity = actualmol->betaDensity = NULL;

   if(!strstr(line, "ALPHA")){                 /* only one matrix */
      if((actualmol->alphaDensity = read_trimat()) == NULL){
         showinfobox("can't allocate density matrix");
         return 0;
      }
      return 1;
   }

   /* two matrices */
   if((actualmol->alphaDensity = read_trimat()) == NULL){
      showinfobox("can't allocate alpha density matrix");
      return 0;
   }
   if(!find_string("BETA DENSITY MATRIX.")) return 0;
   if((actualmol->betaDensity = read_trimat()) == NULL){
      showinfobox("can't allocate beta density matrix");
      return 0;
   }

   return 1;
}



float **read_trimat(void)
{
   register short i, j;
   float **trimat;

   if((trimat = (float**) alloc_trimat(actualmol->nBasisFunctions, sizeof(float))) == NULL){
      showinfobox("can't allocate density matrix");
      return NULL;
   }

   for(i=0; i<=((actualmol->nBasisFunctions-1)/5); i++){
      fgets(line, 255, fp);
      for(j=5*i; j<actualmol->nBasisFunctions; j++){
         fgets(line, 255, fp);
         sscanf(line + 20, "%f %f %f %f %f", trimat[j]+5*i, trimat[j]+5*i+1,
                          trimat[j]+5*i+2, trimat[j]+5*i+3, trimat[j]+5*i+4);
      }
   }

   return trimat;
}




void print_density_matrix(void)
{
   register short i, j;

   for(i=0; i<actualmol->nBasisFunctions; i++){
      for(j=0; j<=i; j++) printf("%8.5f", actualmol->alphaDensity[i][j]);
      printf("\n");
   }
   if(actualmol->betaDensity){
      for(i=0; i<actualmol->nBasisFunctions; i++){
         for(j=0; j<=i; j++) printf("%8.5f", actualmol->betaDensity[i][j]);
         printf("\n");
      }
   }
}



void print_coefficients(void)
{
   register int i, j;

   printf("\n   MO  eigenvalue  occupation\n\n");
   for(i=0; i<actualmol->nMolecularOrbitals; i++){
      printf(" %3d  %10.5f      %1f", i+1,
            actualmol->alphaOrbital[i].eigenvalue,
            actualmol->alphaOrbital[i].occ);

      for(j=0; j<actualmol->nBasisFunctions; j++){
         if(!(j%10)) printf("\n");
         printf(" %8.5f", actualmol->alphaOrbital[i].coefficient[j]);
      }
      printf("\n");

   }
   printf("\n");
}


int read_frequencies(void)
{
   long fpos;
   int n_freq;
   register short i, j;
   AtoM *ap;
   Vibration *vib;

   rewind(fp);

   if(!find_string("Harmonic frequencies (cm**-1)")) return 0;

   fpos = ftell(fp);
   if(!find_string(" Frequencies ---")) fseek(fp, fpos, SEEK_SET);
   else fpos = ftell(fp);
   n_freq = 0;
   while(find_string(" Frequencies -- ")) n_freq += 3;

   if((actualmol->vibration = (Vibration*) malloc(n_freq*sizeof(Vibration))) == NULL){
      showinfobox("can't realloc vibrational frequency");
      return 0;
   }
   for(i=0, vib=actualmol->vibration; i<n_freq; i++, vib++){
      if((vib->coord = (Vector*) malloc(actualmol->natoms*sizeof(Vector))) == NULL){
         showinfobox("can't allocate vibration");
         return 0;
      }
   }

   fseek(fp, fpos, SEEK_SET);
   for(i=0, vib=actualmol->vibration; i<n_freq/3; i++, vib += 3) {
      find_string(" Frequencies -- ");
      actualmol->n_frequencies += sscanf(line+15, "%f %f %f", &vib->frequency,
         &(vib+1)->frequency, &(vib+2)->frequency);
      fseek(fp, preprevious, SEEK_SET);
      fgets(line, 255, fp);
      sscanf(line, "%s %s %s", vib->type, (vib+1)->type, (vib+2)->type);
      if(!find_string("Atom AN      X      Y      Z        X      Y      Z"))
         return 0;
      for(ap=actualmol->firstatom, j=0; ap; ap=ap->next, j++){
         if(!fgets(line, 255, fp)) return 0;
         sscanf(line, "%*d %*d %f %f %f %f %f %f %f %f %f", &vib->coord[j].x,
            &vib->coord[j].y, &vib->coord[j].z, &(vib+1)->coord[j].x,
            &(vib+1)->coord[j].y, &(vib+1)->coord[j].z, &(vib+2)->coord[j].x,
            &(vib+2)->coord[j].y, &(vib+2)->coord[j].z);
      }
   }


   return 1;
}

int read_dipole(void)
{
   long fpos;
   float x, y, z;

   rewind(fp);
   if(!find_string("Dipole moment")) return 0;
   fpos = ftell(fp);
   while(find_string("Dipole moment")) fpos = ftell(fp);
   fseek(fp, fpos, SEEK_SET);
   fgets(line, 255, fp);

   sscanf(line, "%*s %f %*s %f %*s %f", &x, &y, &z);
   actualmol->dipole = add_dipole(x, y, z);

   return 1;
}

void print_frequencies(void)
{
   register short i, j;
   AtoM *ap;
   Vibration *vp;

   if(!actualmol->vibration) return;

   printf("%d vibrational modes\n", actualmol->n_frequencies);

   for(i=0, vp=actualmol->vibration; i<actualmol->n_frequencies; i++, vp++){
      printf("%3d %10.3f \n", i+1, vp->frequency);
      for(ap=actualmol->firstatom, j=0; ap; ap=ap->next, j++){
         printf("                %5.2f%5.2f%5.2f\n", vp->coord[j].x,
            vp->coord[j].y, vp->coord[j].z);
      }
   }  
}



int read_atomic_charges(void)
{
   long total, fitted, choice;
   AtoM *ap;

   total = fitted = 0;
   rewind(fp);
   if(find_string("Total atomic charges")) {
      do {
         total = ftell(fp);
      } while(find_string("Total atomic charges"));
   }
   rewind(fp);
   if(find_string("Charges from ESP fit")) {
      do {
         fitted = ftell(fp);
      } while(find_string("Charges from ESP fit"));
   }

   if(total && fitted){   /* select which set to read */
/* to be fixed
      choice = selection(&xselection, &yselection,
          "Choose one of the two sets of atomic charges", 2, list);
*/
      choice = 1;
      if(choice < 0) return 0;
      if(choice == 0) {
         fitted = 0;
         logprint("Reading the total atomic charges");
      }
      else if(choice == 1) {
         total = 0;
         logprint("Reading the atomic charges");
         logprint("   from ESP fit");
      }
   }

   if(total){
      fseek(fp, total, SEEK_SET);
      fgets(line, 255, fp);

      for(ap=actualmol->firstatom; ap; ap = ap->next) {
         if(!fgets(line, 255, fp)) return 0;
         if(sscanf(line, "%*d %*s %f", &ap->charge) != 1) return 0;
/*
      printf("Atom %s, charge = %f\n", element[ap->ord].symbol, ap->charge);
*/
      }
   }
   else if(fitted){
      fseek(fp, fitted, SEEK_SET);
      fgets(line, 255, fp);
      fgets(line, 255, fp);

      for(ap=actualmol->firstatom; ap; ap = ap->next) {
         if(!fgets(line, 255, fp)) return 0;
         if(sscanf(line, "%*d %*s %f", &ap->charge) != 1) return 0;
      }
   }
   else return 0;

   actualmol->charges = 1;
   return 1;
}


