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


/* lecture of HONDO output */
#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "readhondo.h"
#include "general.h"
#include "utils.h"
#include "maininterf.h"
#include "box.h"
#include "browser.h"


static int read_atomic_coordinates(void);
static char *find_string(char *s);
static int read_charge(void);
static int read_basis_set(void);
static int read_coefficients(void);
static int read_eigenvalues(void);
static int fill_mos(int nnn, int nbase, MolecularOrbital *m);

static FILE *fp;
static char line[256];
static long previous_line = 0, preprevious = 0;


/**** lecture of HONDO output ****/

void read_hondo(char *file)
{

   if((fp = fopen(file, "r")) == NULL){
      sprintf(line, "can't open file\n%s !", file);
      showinfobox(line);
      return;
   }

   if(!find_string("HONDO")) {
      showinfobox("read_hondo : HONDO output only!");
      fclose(fp);
      return;
   }

   add_mol(file);

   if(!read_atomic_coordinates()){
      showinfobox("can't read the atomic coordinates");
      free_mol(actualmol);
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

   if(!read_charge()){
      showinfobox("unable to read the atomic charges");
      actualmol->charges = 0;
   }


   if(!read_basis_set()){
      showinfobox("can't read the basis-set");
   }

/*
   print_basis_set();
*/

   fclose(fp);

   logprint("");
   logprint("select corresponding");
   logprint("   HONDO PUNCH file!");
   if(maininterfwin) update_logs();
   file_select(LOAD_HONDO_PUNCH, pun_ext);
}


void read_hondo_punch(char *file)
{

   if((fp = fopen(file, "r")) == NULL){
      sprintf(line, "can't open file\n%s !", file);
      showinfobox(line);
      return;
   }

   if(!read_eigenvalues()){
      showinfobox("can't read the eigenvalues");
      return;
   }

   rewind(fp);

   if(!read_coefficients()){
      showinfobox("can't read the MO-coefficients");
      return;
   }

   computeOccupations(actualmol);

/*
   print_coefficients();
*/

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




static int read_atomic_coordinates(void)
{
   long fpos;
   float x, y, z, ford;
   int ord;

    if(find_string("COORDINATES (BOHR)")) {
	do {
	    fpos = ftell(fp);
	} while(find_string("COORDINATES (BOHR)"));

	fseek(fp, fpos, SEEK_SET);
	if(!find_string("ATOM     ZNUC       X             Y             Z")) return 0;
	if(!find_string("-----")) return 0;

	fgets(line, 255, fp);
	while(1) {
	    if(!fgets(line, 255, fp)) return 0;
	    if(sscanf(line, "%*d %*s %d %f %f %f", &ord, &x, &y, &z)
		!= 4) break;
	    add_atom(ord, -x*BOHR, y*BOHR, z*BOHR);
	}
	
	return 1;
   }
   else {
	rewind(fp);
	if(find_string("ATOM   ATOMIC                COORDINATES")) {
	    do {
		fpos = ftell(fp);
	    } while(find_string("ATOM   ATOMIC                COORDINATES"));
	}
	else return 0;

	fseek(fp, fpos, SEEK_SET);
	if(!find_string("-----")) return 0;

	fgets(line, 255, fp);
	while(1) {
	    if(!fgets(line, 255, fp)) return 0;
	    if(strstr(line, "---------")) break;
	    if(sscanf(line, "%*s %f %f %f %f", &ford, &x, &y, &z)
		!= 4) continue;
	    add_atom((int)ford, -x*BOHR, y*BOHR, z*BOHR);
	}
   }

   return 1;
}




static int read_charge(void)
{
   long fpos;
   float ch;
   AtoM *ap;

   if(!find_string("MULLIKEN CHARGE")) return 0;
   do {
      fpos = ftell(fp);
   } while(find_string("MULLIKEN CHARGE"));
   fseek(fp, fpos, SEEK_SET);

   for(ap = actualmol->firstatom; ap; ap = ap->next){
      if(!fgets(line, 255, fp)) return 0;
      if(!sscanf(line, "%*d%*s%*f%*f%*f%f", &ch)) return 0;
      ap->charge = ap->ord - ch;
   }
   actualmol->charges = 1;
   return 1;
}




static int read_basis_set(void)
{
   AtoM *ap;
   Gauss *gp;
   Shell *sp;
   char sym[5], ptgrp[6], symtmp[5];
   int i1, i2;
   float f1, f2, f3;
   long fpos;

   rewind(fp);
   if(!find_string("THE POINT GROUP OF THE MOLECULE IS")) return 0;
   sscanf(line+39, "%s", ptgrp);
   if(!find_string("MOLECULAR BASIS SET")) return 0;
   if(!find_string("CONTRACTED PRIMITIVE FUNCTIONS")) return 0;
   if(!find_string("CONTRACTION COEFFICIENTS")) return 0;
   fgets(line, 255, fp); /*  "UNNORM.        NORM." */
   fgets(line, 255, fp); /*  blank */
   fgets(line, 255, fp); /*  first atom symbol */
   sscanf(line, "%s", sym);
   strcpy(symtmp, sym); /* We are now at the possition of first symbol
                           which is copied to symtmp */
   if(strcmp(ptgrp, "C1") != 0) {
      logprint("WARNING:");
      logprint("   Reading only first");
      logprint("   occurrance of basis set.");
      logprint("   Assigning this same");
      logprint("   basis set to each atom");
      logprint("   of the same type!");
      fpos = ftell(fp);
      for(ap = actualmol->firstatom; ap; ap = ap->next){
         fseek(fp, fpos, SEEK_SET); /* back to first symbol */
         strcpy(sym, symtmp);       /* set symbol to to initial value */
         if(get_ordinal(sym) != ap->ord) {
            while(fgets(line, 255, fp)) {
               if(line[0] == '\n') continue;
               if(strlen(line) < 16) { /* new atom symbol */
                  sscanf(line, "%s", sym);
                  if(get_ordinal(sym) != ap->ord) continue;
                  else break;
                  }
               continue;
            }
         }
         while(fgets(line, 255, fp)) {
            if(line[0] == '\n') continue;
            if(strlen(line) < 16) { /* new atom symbol */
               break;
            }
            if(sscanf(line, "%d %*s %d %f %f %*s %f",
                             &i1, &i2, &f1, &f2, &f3) != 5) break;
  
            if(!(sp = add_shell(ap))) return 0;
  
            while(line[0] != '\n') {
               gp = add_gauss(sp);

               if(strchr(line, 'S')){
                  sscanf(line, "%*d %*s %*d %lf %*f %*s %lf",
                     &gp->exponent, &gp->coeff);
                  sp->n_base = 1;
                  if(strchr(line, 'P')){
                     sscanf(line+82, "%*f %*s %lf", &gp->coeff2);
                     sp->n_base = 4;
                  }
               }
               else if(strchr(line, 'P')){
                  sscanf(line, "%*d %*s %*d %lf %*f %*s %lf",
                     &gp->exponent, &gp->coeff);
                  sp->n_base = 3;
               }
               else if(strchr(line, 'D')){
                  sscanf(line, "%*d %*s %*d %lf %*f %*s %lf",
                     &gp->exponent, &gp->coeff);
                  sp->n_base = 6; /* 5 or 6 ??? */
               }
               else {
                  showinfobox("can't read gaussian primitives");
                  return 0;
               }
               if(!fgets(line, 255, fp)) return 0;
            }                                /* end of gaussian primitives */
         }                                   /* end of shell */
         fseek(fp, fpos, SEEK_SET); 
      }                                      /* end of loop over atoms */
   } 
   else {
   for(ap = actualmol->firstatom; ap; ap = ap->next){
     if(get_ordinal(sym) != ap->ord) {
        showinfobox("read HONDO : incoherent output-file");
        return 0;
     } 

      while(fgets(line, 255, fp)) {
         if(line[0] == '\n') {
            continue;
         }
         if(strlen(line) < 16) { /* new atom symbol */
            sscanf(line, "%s", sym);
            break;
         }
         if(sscanf(line, "%d %*s %d %f %f %*s %f",
                           &i1, &i2, &f1, &f2, &f3) != 5) break;

         if(!(sp = add_shell(ap))) return 0;

         while(line[0] != '\n') {
            gp = add_gauss(sp);

            if(strchr(line, 'S')){
               sscanf(line, "%*d %*s %*d %lf %*f %*s %lf",
                  &gp->exponent, &gp->coeff);
               sp->n_base = 1;
               if(strchr(line, 'P')){
                  sscanf(line+82, "%*f %*s %lf", &gp->coeff2);
                  sp->n_base = 4;
               }
            }
            else if(strchr(line, 'P')){
               sscanf(line, "%*d %*s %*d %lf %*f %*s %lf",
                  &gp->exponent, &gp->coeff);
               sp->n_base = 3;
            }
            else if(strchr(line, 'D')){
               sscanf(line, "%*d %*s %*d %lf %*f %*s %lf",
                  &gp->exponent, &gp->coeff);
               sp->n_base = 6; /* 5 or 6 ??? */
            }
            else {
               showinfobox("can't read gaussian primitives");
               return 0;
            }
            if(!fgets(line, 255, fp)) return 0;
         }                                /* end of gaussian primitives */
      }                                   /* end of shell */
   }                                      /* end of basis-set */
   }

    if(!find_string("TOTAL NUMBER OF BASIS FUNCTIONS")) return 0;
    sscanf(line+38, "%d", &actualmol->nBasisFunctions);

    normalize_gaussians();

   return 1; 
}





static int read_eigenvalues(void)
{
   MolecularOrbital *mo, *mbeta;
   int nocc, nsing, nvirt, i, nnn;
   long fpos;
   static char *str1 = "$EIG", *str2 = "$EIGA", *pkey;

   if(!find_string(str1)) return 0;
   if(strstr(line, str2)) {
      actualmol->alphaBeta = 1;    /* UHF */
      pkey = str2;
   }
   else {
      actualmol->alphaBeta = 0;
      pkey = str1;
   }

   do {
      fpos = ftell(fp);
   } while(find_string(pkey));
   fseek(fp, fpos, SEEK_SET);

   fgets(line, 255, fp);
   if(sscanf(line, "%d%d%d", &nsing, &nocc, &nvirt) != 3) return 0;
   if(nsing > nocc) { /* inverse */
	int tmp;
	tmp = nsing;
	nsing = nocc;
	nocc = tmp;
   }
   actualmol->nAlpha = nocc + nsing;
   actualmol->nBeta = nocc;
   nnn = nocc + nsing + nvirt;
   actualmol->nMolecularOrbitals = nnn;

   if((mo = allocOrbital(nnn, 0, HONDO_ORB)) == NULL){
      showinfobox("can't allocate the MO-structures");
      return 0;
   }

   actualmol->alphaOrbital = mo;
   i=0;
   while(fgets(line, 255, fp)){
      if(strstr(line, "$END")) break;
      sscanf(line, "%*s %lf %*s %lf %*s %lf %*s %lf", &mo[i].eigenvalue,
         &mo[i+1].eigenvalue, &mo[i+2].eigenvalue, &mo[i+3].eigenvalue);
      i += 4;
   }

   if(!actualmol->alphaBeta) return 1;

   if(!find_string("$EIGB")) return 0;
   fgets(line, 255, fp);
   if(sscanf(line, "%d%d%d", &nocc, &nsing, &nvirt) != 3) return 0;

   if((mbeta = allocOrbital(nnn, 0, HONDO_ORB)) == NULL){
      showinfobox("can't allocate the MO-structures");
      return 0;
   }

   actualmol->betaOrbital = mbeta;
   i=0;
   while(fgets(line, 255, fp)){
      if(strstr(line, "$END")) break;
      sscanf(line, "%*s %lf %*s %lf %*s %lf %*s %lf", &mbeta[i].eigenvalue,
         &mbeta[i+1].eigenvalue, &mbeta[i+2].eigenvalue,
         &mbeta[i+3].eigenvalue);
      i += 4;
   }

   return 1;
}

static int read_coefficients(void)
{
   MolecularOrbital *mo, *mbeta;
   long fpos;
   int nnn, nbase;
   register int i;
   static char *str1 = "$VEC", *str2 = "$VECA", *pkey;

   if(actualmol->alphaBeta) { pkey = str2; }
   else                     { pkey = str1; }
   nbase = actualmol->nBasisFunctions;
   nnn   = actualmol->nMolecularOrbitals;
   mo    = actualmol->alphaOrbital;
   mbeta = actualmol->betaOrbital;

   if(!find_string(pkey)) return 0;

   do {
      fpos = ftell(fp);
   } while(find_string(pkey));
   fseek(fp, fpos, SEEK_SET);

   for(i=0; i<nnn; i++){
      if((mo[i].coefficient = (double *)malloc(nbase*sizeof(double))) == NULL){
         showinfobox("can't allocate the MO-coeffs");
         return 0;
      }
   }

    if(!fill_mos(nnn, nbase, mo)) {
	rewind(fp);
	find_string(pkey);
	if(!fill_mos(nnn, nbase, mo)) return 0;
    }

   if(!actualmol->alphaBeta) return 1;

   if(!find_string("$VECB")) return 0;

   for(i=0; i<nnn; i++){
      if((mbeta[i].coefficient = (double *)malloc(nbase*sizeof(double))) == NULL){
         showinfobox("can't allocate the beta MO-coeffs");
         return 0;
      }
   }

    if(!fill_mos(nnn, nbase, mbeta)) return 0;
    
    return 1;
}



static int fill_mos(int nnn, int nbase, MolecularOrbital *m)
{
    register int i, j;
        
    for(i=0; i<nnn; i++) {
	for(j=0; j<nbase; j+=5) {
	    if(!fgets(line, 255, fp)) return 0;
	    if(strstr(line, "$END")) return 0;
	    sscanf(line, "%*d%*d%lf%lf%lf%lf%lf", &m[i].coefficient[j],
                 &m[i].coefficient[j+1], &m[i].coefficient[j+2],
                 &m[i].coefficient[j+3], &m[i].coefficient[j+4]);
	}
    }
    
    return 1;
}

