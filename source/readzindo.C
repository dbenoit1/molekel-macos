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


/* lecture of ZINDO output */
#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "readzindo.h"
#include "general.h"
#include "maininterf.h"
#include "utils.h"
#include "box.h"

static int read_atomic_coordinates(char *file);
static char *find_string(char *s);
static int read_charge(void);
static int read_basis_set(void);
static int read_coefficients(void);
static int read_eigenvalues(void);
static int fill_mos(int nnn, MolecularOrbital *m);
static int read_occupations(Mol *mp);

static FILE *fp;
static char line[256];
static long previous_line = 0, preprevious = 0;


/**** lecture of ZINDO output ****/

void read_zindo(char *file)
{
   int nel, mult;
   unsigned long position;

   if((fp = fopen(file, "r")) == NULL){
      sprintf(line, "can't open file\n%s !", file);
      showinfobox(line);
      return;
   }

   if(!find_string("PROGRAM DEFAULTS")) {
      showinfobox("read_zindo : ZINDO output only !");
      fclose(fp);
      return;
   }

   rewind(fp);
   if(!find_string("USER DEFINED VARIABLES")) {
      showinfobox("Can't find the USER DEFINED VARIABLES");
      fclose(fp);
      return;
   }
   if(!find_string(" NEL ")) {
      showinfobox("Can't find the number of electrons");
   }
   else {
       char *cp;
       
       cp = strstr(line, " NEL ");
       sscanf(cp+5, "%d", &nel);
       cp = strstr(line, " MULT ");
       if(!cp) {
	    showinfobox("Can't find the multiplicity");
       }
       else {
	   sscanf(cp+6, "%d", &mult);
       }
   }
    

   if(!read_atomic_coordinates(file)){
      showinfobox("can't read the atomic coordinates");
      fclose(fp);
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

    actualmol->nElectrons = nel;
    actualmol->nBeta = nel/2;
    actualmol->nAlpha = nel - nel/2;
    actualmol->multiplicity = mult;

   if(!read_charge()){
      showinfobox("unable to read the atomic charges");
      actualmol->charges = 0;
   }

   if(!read_basis_set()){
      showinfobox("can't read the basis-set");
   }

    position = ftell(fp);
   if(!read_eigenvalues()){
      showinfobox("can't read the eigenvalues");
      
   }

    fseek(fp, position, SEEK_SET);
   if(!read_coefficients()){
      showinfobox("can't read the MO-coefficients");
      fseek(fp, position, SEEK_SET);
   }

    rewind(fp);
   if(!read_occupations(actualmol)){
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




static int read_atomic_coordinates(char *file)
{
   float x, y, z, charge;
   char atom_symbol[5];
   AtoM *ap;

    if(find_string("SUMMARY OF POPULATION STUDY IN THE ZDO BASIS")) {
	if(!find_string("ATOM   TYPE    FORMAL     TOTAL")) return 0;
	if(!find_string("CHARGE     CHARGE        X          Y          Z")) return 0;
	fgets(line, 255, fp);

	add_mol(file);

	while(1) {
	    if(!fgets(line, 255, fp)) return 0;
	    if(sscanf(line, "%*d %s %f %*f %f %f %f",
/*			    atom_symbol, &charge, &x, &y, &z) */
			    atom_symbol, &charge, &z, &x, &y)
		!= 5) break;
	    ap = add_atom(get_ordinal(atom_symbol), x, y, z);
	    ap->charge = charge;
	}
	actualmol->charges = 1;
    }
    else {
	return 0;
    }
    
   return 1;
}




static int read_charge(void)
/* Mulliken charge (if present) overwrites atomic charge */
{
   float ch;
   AtoM *ap;

   if(!find_string("ATOM   TYPE    FORMAL     TOTAL")) return 0;
   if(!find_string("CHARGE     CHARGE")) return 0;
   if(strstr(line, "X          Y          Z")) return 0;

    fgets(line, 255, fp);
   for(ap = actualmol->firstatom; ap; ap = ap->next){
      if(!fgets(line, 255, fp)) return 0;
      if(!sscanf(line, "%*d%*s%f", &ch)) return 0;
      ap->charge = ch;
   }
   actualmol->charges = 1;
   return 1;
}




static int read_basis_set(void)
{
   AtoM *ap;
   Slater *sp;
   float exponent;
   char sym[10], type;
   int index, nFunctions, prin;

   rewind(fp);
   if(!find_string("0  A.O.     ATOM  SYM. TYPE     PRIN.QUANT.NO.     VALENCE TYPE ")) return 0;
   fgets(line, 255, fp); /*  blank */
   fgets(line, 255, fp); /*  blank */

    ap = actualmol->firstatom;
    type = 0;
    nFunctions = 0;
    while(fgets(line, 255, fp)) {
	if(sscanf(line, "%*d %d %s %*d %d %*d %f",
			&index, sym, &prin, &exponent) != 4) break;
	nFunctions++;
	if(index > ap->name) {
	    ap = ap->next;
	    if(!ap) return 0;
	    type = 0;
	}
	if(sym[0] == type) continue;
	type = sym[0];
	
	sp = add_slater(ap);
        if(sym[0] != 'S' && sym[0] != 'P' && sym[0] != 'D') return 0;
	sp->type[0] = sym[0]; sp->type[1] = 0;
	sp->n = prin;
	sp->exponent = (double)exponent;
    }
    if(ap->next) return 0;
    
    actualmol->nBasisFunctions = nFunctions;

    return 1; 
}





static int read_eigenvalues(void)
{
   MolecularOrbital *mo;
   int i, k, nitems;
   int index[12];
   long fpos;
   double eig[12];
   static char *pkey = "0EIGENVALUES";

   if(!find_string(pkey)) return 0;
    actualmol->alphaBeta = 0;

   do {
      fpos = ftell(fp);
   } while(find_string(pkey));
   fseek(fp, fpos, SEEK_SET);

    nitems = 0;
 /* first scan through just for counting the eigenvalues */
    while(fgets(line, 255, fp)) {
	if(line[0] != ' ' && line[0] != '\n') break;
	if(strstr(line, " NO.")) {
	    nitems += sscanf(line, "%*s%d%d%d%d%d%d%d%d%d%d%d%d", 
		index, index+1, index+2, index+3, index+4, index+5,
		index+6, index+7, index+8, index+9, index+10, index+11);
	    if(!actualmol->firstOrbital) actualmol->firstOrbital = index[0];
	}
    }
    actualmol->lastOrbital = nitems;
    actualmol->nMolecularOrbitals = 
	actualmol->lastOrbital - actualmol->firstOrbital + 1;
   
   fseek(fp, fpos, SEEK_SET);

   if((mo = allocOrbital(actualmol->nMolecularOrbitals, 0, ZINDO_ORB)) == NULL){
      showinfobox("can't allocate the ZINDO MO-structures");
      return 0;
   }

    actualmol->alphaOrbital = mo;
 /* read the eigenvalues */
    while(fgets(line, 255, fp)) {
	if(line[0] == '\n') continue;
	if(line[0] != ' ') break;
	if(strstr(line, " NO.")) {
	    sscanf(line, "%*s%d", &i);
	    i -= actualmol->firstOrbital;
	}
	else if(strstr(line, " SYMMETRY")) continue;
	else if(strstr(line, " VALUE(A.U.)")) {
	    nitems = sscanf(line, "%*s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
		eig, eig+1, eig+2, eig+3, eig+4, eig+5, 
		eig+6, eig+7, eig+8, eig+9, eig+10, eig+11);
	    for(k=0; k<nitems; k++) {
		mo[i+k].eigenvalue = eig[k];
	    }
	}
    }

   if(!actualmol->alphaBeta) return 1;

   return 0;
}






static int read_coefficients(void)
{
   MolecularOrbital *mo, *mbeta;
   long fpos;
   int nnn;
   register int i;
   static char *pkey = "0EIGENVECTOR MATRIX";

   nnn   = actualmol->nMolecularOrbitals;
   mo    = actualmol->alphaOrbital;
   mbeta = actualmol->betaOrbital;

   if(!find_string(pkey)) return 0;

   do {
      fpos = ftell(fp);
   } while(find_string(pkey));
   fseek(fp, fpos, SEEK_SET);

   for(i=0; i<nnn; i++){
      if((mo[i].coefficient = (double *)malloc(nnn*sizeof(double))) == NULL){
         showinfobox("can't allocate the MO-coeffs");
         return 0;
      }
   }

    if(!fill_mos(nnn, mo)) return 0;
    
   if(!actualmol->alphaBeta) return 1;
   else return 0;
}



static int fill_mos(int nnn, MolecularOrbital *m)
{
    int i, j, k, nitems;
    double coeff[12];
    
    while(1) {
	if(!fgets(line, 255, fp)) return 0;
	if(line[0] == '\n') continue;
	if(line[0] != ' ') break;
	if(strstr(line, "M.O.")) {
	    sscanf(line+20, "%d", &i);
	    i -= actualmol->firstOrbital;
	}
	else {
	    sscanf(line, "%d",  &j);
	    j--;
	    nitems = sscanf(line+20, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
		coeff, coeff+1, coeff+2, coeff+3, coeff+4, coeff+5, 
		coeff+6, coeff+7, coeff+8, coeff+9, coeff+10, coeff+11);
	    
	    for(k=0; k<nitems; k++) {
		m[i+k].coefficient[j] = coeff[k];
	    }
	}
    }
    actualmol->nBasisFunctions = j+1;
    
    return 1;
}

static int read_occupations(Mol *mp)
{
   int i, k, nitems;
   float Occ[8];


   if(mp->alphaBeta) {
      showinfobox("Only excited states\nof closed shells supported!\nCalculated occupations used.\n");
      return 0;
   }

   if(!mp->alphaOrbital) return 0; 

   if(!find_string("STATe")) {
      showinfobox("No occupations specified to read in output");
      return 0;
   }

   logprint("Reading occupations of");
   logprint(strpbrk(line, "S"));
   fgets(line, 255, fp);
   for(i=1; i<=mp->nMolecularOrbitals; i++) {
      if(i < mp->firstOrbital) continue;
       while(1) {
           if(!fgets(line, 255, fp)) return 0;
           if(line[0] == '\n') break;
           nitems = sscanf(line, "%*d %f%*s %*d %f%*s %*d %f%*s\
               %*d %f%*s %*d %f%*s %*d %f%*s %*d %f%*s %*d %f%*s",
               Occ, Occ+1, Occ+2, Occ+3, Occ+4, Occ+5, Occ+6, Occ+7);
           for(k=0; k<nitems; k++) {
                mp->alphaOrbital[i+k - mp->firstOrbital].occ = Occ[k];
           }
           i += k;
       }
   }

   if(i - 1 -  mp->firstOrbital != mp->nMolecularOrbitals) {
      showinfobox("No. of MO's does not match");
      return 0;
   }

   return 1;
}
