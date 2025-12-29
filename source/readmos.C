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


/* lecture of MOS output (Baumann,  ETHZ) */

#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "readmos.h"
#include "general.h"
#include "maininterf.h"
#include "utils.h"
#include "box.h"

static int read_atomic_coordinates(char *file);
static char *find_string(char *s);
static int readVars();
static int read_basis_set(void);
static int read_coefficients(void);

static FILE *fp;
static char line[256];
static long previous_line = 0, preprevious = 0;
static int nAtoms, nOrbs, nElecs;


#include "exponents.h"

/**** lecture of Mos output ****/

void read_mos(char *file)
{
   if((fp = fopen(file, "r")) == NULL){
      sprintf(line, "can't open file\n%s !", file);
      showinfobox(line);
      return;
   }


    if(!readVars()) {
      showinfobox("read_mos : MOS output only !");
      fclose(fp);
      return;
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

   if(!read_basis_set()){
      showinfobox("can't read the basis-set");
   }

   if(!read_coefficients()){
      showinfobox("can't read the MO-coefficients");
      return;
   }

   computeOccupations(actualmol);

   fclose(fp);
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



static int readVars()
{
    char str[20];
    if(!find_string("NATOMS=")) return 0;
    if(!sscanf(line+7, "%d", &nAtoms)) return 0;
    if(nAtoms <= 0) return 0;
    
    if(!find_string("NORBS=")) return 0;
    if(!sscanf(line+6, "%d", &nOrbs)) return 0;
    if(nOrbs <= 0) return 0;

    if(!find_string("NELECS=")) return 0;
    if(!sscanf(line+7, "%d", &nElecs)) return 0;
    if(nElecs <= 0) return 0;

    sprintf(str, "%d atoms", nAtoms);
    logprint(str);
    sprintf(str, "%d orbitals", nOrbs);
    logprint(str);
    sprintf(str, "%d electons", nElecs);
    logprint(str);
    update_logs();

    return 1;
}


static int read_atomic_coordinates(char *file)
{
   long fpos;
   float x, y, z;
   int ord, i;

   if(find_string("END")) {
      do {
         fpos = ftell(fp);
      } while(find_string("END"));
   }
   else return 0;

   fseek(fp, fpos, SEEK_SET);

   add_mol(file);

   for(i=0; i<nAtoms; i++) {
      if(!fgets(line, 255, fp)) return 0;
      if(sscanf(line, "%d %f %f %f", &ord, &z, &x, &y) != 4) return 0;
/*      if(sscanf(line, "%d %f %f %f", &ord, &x, &y, &z) != 4) return 0; */
      add_atom(ord, x, y, z);
   }

   return 1;
}


#define POW3(x)  ((x)*(x)*(x))
#define POW5(x)  ((x)*(x)*(x)*(x)*(x))
#define POW7(x)  ((x)*(x)*(x)*(x)*(x)*(x)*(x))
#define POW9(x)  ((x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x))
#define POW11(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x))


static int read_basis_set(void)
{
    AtoM *ap;
    Slater *sp;
    int principalQuantumNumber;

    for(ap = actualmol->firstatom; ap; ap = ap->next){

	if(ap->ord <= 2) principalQuantumNumber = 1;
	else if(ap->ord <= 10) principalQuantumNumber = 2;
	else if(ap->ord <= 18) principalQuantumNumber = 3;
	else if(ap->ord <= 36) principalQuantumNumber = 4;
	else if(ap->ord <= 54) principalQuantumNumber = 5;
/*
	else if(ap->ord <= 86) principalQuantumNumber = 6;
	else principalQuantumNumber = 7;
*/
	else {
	    showinfobox("atoms heavier than Xenon are not supported\n");
	    return 0;
	}

    /* S - shell */
	if(!(sp = add_slater(ap))) return 0;
	sp->n = principalQuantumNumber;
	sp->type[0] = 'S'; sp->type[1] = 0;
	sp->exponent = Zs[ap->ord];
	switch(principalQuantumNumber) {
	    case 1: sp->norm[0] = sqrt(POW3(sp->exponent) * M_1_PI); break;
	    case 2: sp->norm[0] = sqrt(POW5(sp->exponent) * M_1_PI / 3.); break;
	    case 3: sp->norm[0] = sqrt(POW7(sp->exponent) * M_2_PI / 45.); break;
	    case 4: sp->norm[0] = sqrt(POW9(sp->exponent) * M_1_PI / 315.); break;
	    case 5: sp->norm[0] = sqrt(POW11(sp->exponent) * M_2_PI / 14175.); break;
	}

	if(ap->ord <= 2) continue;

    /* P - shell */
	if(!(sp = add_slater(ap))) return 0;
	sp->n = principalQuantumNumber;
	sp->type[0] = 'P'; sp->type[1] = 0;
	sp->exponent = Zp[ap->ord];
	switch(principalQuantumNumber) {
	    case 2: sp->norm[1] = sqrt(POW5(sp->exponent) * M_1_PI); break;
	    case 3: sp->norm[1] = sqrt(POW7(sp->exponent) * M_2_PI / 15.); break;
	    case 4: sp->norm[1] = sqrt(POW9(sp->exponent) * M_1_PI / 105.); break;
	    case 5: sp->norm[1] = sqrt(POW11(sp->exponent) * M_2_PI / 4725.); break;
	}
    }

    return 1; 
}


static int read_coefficients(void)
{
    MolecularOrbital *mo;
    int index;
    float eigenValue;
    register int i, j;

    actualmol->alphaBeta = 0;

    if((mo = allocOrbital(nOrbs, nOrbs, MOS_ORB)) == NULL){
	showinfobox("can't allocate the MO-structures");
	return 0;
    }

    actualmol->nMolecularOrbitals = nOrbs;
    actualmol->nBasisFunctions = nOrbs;
    actualmol->alphaOrbital = mo;
    actualmol->betaOrbital = NULL;
    actualmol->firstOrbital = 1;
    actualmol->nElectrons = nElecs;
    actualmol->nBeta = nElecs/2;
    actualmol->nAlpha = nElecs - actualmol->nBeta;
    actualmol->alphaDensity = actualmol->betaDensity = NULL;

    for(i=0; i<nOrbs; i++) {
	if(!fgets(line, 255, fp)) return 0;
	if(sscanf(line, "%d%f", &index, &eigenValue) != 2 || index != i+1) {
	    fprintf(stderr, "unexpected line %s", line);
	    return 0;
	}
	mo[i].eigenvalue = eigenValue;
	
	for(j=0; j<nOrbs; j++) {
	    if(!fscanf(fp, "%lf", &mo[i].coefficient[j])) return 0;
        }
	fgets(line, 255, fp);
    }

    return 1;
}

