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


/* lecture of ADF-output */
#include <dirent.h>
#include <sys/stat.h>
#ifndef WIN32
#include <sys/param.h>
#endif
#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "readadf.h"
#include "general.h"
#include "utils.h"
#include "maininterf.h"
#include "box.h"
#include "browser.h"
#include "calcdens.h"

#ifdef WIN32
#define MAXPATHLEN 512
#endif

static int read_atomic_parameters(char *name);
static int read_atomic_coordinates(void);
static int read_energy_levels(void);
static int read_atomic_charges(void);
static void read_freq(void);
static int read_dipole(void);
static char *find_string(char *s);
static int compare_orbitals(const void *a, const void *b);
static Basis *get_basis(char *basis);

char tape21file[256];
static FILE *fp;
static char line[256];
static long previous_line = 0, preprevious = 0;


/**** lecture of ADF output ****/

void read_adf(char *file)
{

   if((fp = fopen(file, "r")) == NULL){
      sprintf(line, "can't open file\n%s !", file);
      showinfobox(line);
      return;
   }

   if(!find_string("Amsterdam Density Functional  (ADF)")) {
      rewind(fp);
      if(!find_string(">>>> ADF")) {
         showinfobox("read ADF : ADF output only !");
         fclose(fp);
         return;
      }
   }

   if(!read_atomic_parameters(file)){
      showinfobox("can't read the atomic parameters");
      fclose(fp);
      return;
   }

   rewind(fp);
   if(!read_atomic_charges()){
      showinfobox("can't read the atomic charges");
      logprint("can't read the atomic charges");
      update_logs();
/*
      fclose(fp);
      return;
*/
   }

   rewind(fp);
   free_frequencies();
   read_freq();

   if(read_dipole()) {
      logprint("dipole moment present");
   }
   rewind(fp);
   if(!read_energy_levels()){
      showinfobox("can't read the energy levels");
      fclose(fp);
      return;
   }

   fclose(fp);

   logprint("");
   logprint("select corresponding");
   logprint("   ADF tape21 file!");
   if(maininterfwin) update_logs();
   file_select(LOAD_T21, t21_ext);
}

void select_t21_file(char *file)
{
   char *cpt;

   strcpy(tape21file, file);
   logprint("");
   if(!(cpt = strrchr(file, '/')))
      cpt = (char *)file;
   else cpt++;
   logprint(cpt);
   logprint("   as tape21 file selected");
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




static int read_atomic_parameters(char *name)
{
   add_mol(name);
   rewind(fp);

   if(!read_atomic_coordinates()){
      showinfobox("can't read the atomic coordinates");
      fclose(fp);
      free_mol(actualmol);
      return 0;
   }

   if(!actualmol->natoms){
      sprintf(line, " No atoms in %s!", name);
      showinfobox(line);
      free_mol(actualmol);
      fclose(fp);
      return 0;
   }

   create_bonds();
   find_multiplebonds();
//   create_box();
   new_mole(name);

   return 1;
}

int addOptStep(void)
{
   long index, fpos, i;
   static int natoms;
   float x, y, z;
   char str[30];

   index = dynamics.ntotalsteps++;

   if(!dynamics.trajectory){
      if((dynamics.trajectory = (Vector**) malloc(sizeof(Vector *))) == NULL){
         showinfobox("can't allocate dyna pointer\n");
         return 0;
      }

      fpos = ftell(fp);
      natoms = 0;
      while(1) {
         if(!fgets(line, 255, fp)) return 0;
         if(sscanf(line, "%*d%*s%*f%*f%*f%f%f%f", &x, &y, &z)
                  != 3) break;
         natoms++;
      }
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

   i=0;
   while(1) {
      if(!fgets(line, 255, fp)) return 0;
      if(sscanf(line, "%*d%*s%*f%*f%*f%f%f%f", &x, &y, &z)
                != 3) break;
      dynamics.trajectory[index][i].x = x;
      dynamics.trajectory[index][i].y = y;
      dynamics.trajectory[index][i].z = z;
      i++;
   }

   return natoms;
}

static int read_atomic_coordinates(void)
{
   long fpos;
   short ord;
   int natoms;
   float x, y, z;
   AtoM *ap;
   Basis *bp;
   char basis_symbol[100], *format_string;

   free_dyna();

   fpos = ftell(fp);
   if(find_string("Atom   X        Y        Z")) {
      format_string = "%s%f%f%f\0";
      if(!strstr(line, "(Angstr.)")) {
         showinfobox("only Angstroms are supported now\n");
         
         return 0;
      }
      do {
         fpos = ftell(fp);
      } while(find_string("Atom   X        Y        Z"));
   }
   else {
      fseek(fp, fpos, SEEK_SET);
      format_string = "%*d%s%f%f%f\0";

      if(find_string("=====                            X Y Z")) {
         fgets(line, 255, fp);
         if(!strstr(line, "(Angstrom)")) {
            showinfobox("only Angstroms are supported now\n");
            return 0;
         }
         fgets(line, 255, fp);
         fpos = ftell(fp);
         while(find_string("=====                            X Y Z")) {
            fgets(line, 255, fp);
            fgets(line, 255, fp);
            fpos = ftell(fp);
         }
      }
      fseek(fp, fpos, SEEK_SET);
/*      if(find_string("X           Y           Z")) { */
/* pattern enlarged: adf99 prints in the log part (end of the file)
 * X Y Z with the exact same spacing
 * as it did since ever in the output part */
      if(find_string("Coordinates (Cartesian)")) {
         format_string = "%*d%s%*f%*f%*f%f%f%f\0";
         fseek(fp, fpos, SEEK_SET);
/*         while(find_string("X           Y           Z")) { */
         while(find_string("Coordinates (Cartesian)")) {
            fgets(line, 255, fp);
            fgets(line, 255, fp);
            fgets(line, 255, fp);
            fgets(line, 255, fp);
            fgets(line, 255, fp);
            fpos = ftell(fp);
            natoms = addOptStep();
         }
      }
   }

   dynamics.molecule = actualmol;
   dynamics.current = dynamics.ntotalsteps - 1;

   fseek(fp, fpos, SEEK_SET);

   while(1) {
      if(!fgets(line, 255, fp)) return 0;
      if(sscanf(line, format_string, basis_symbol, &x, &y, &z)
               != 4) return 1;
      if(actualmol->basisset){
         if(basis_symbol[0] == 'X' || basis_symbol[0] == 'x') ord = 0;
         else {
            if(!(bp = get_basis(basis_symbol))) return 0;
            ord = bp->ord;
         }
      }
      else {
         ord = get_ordinal(basis_symbol);
      }
      ap = add_atom(ord, x, y, z);
/*      if(ord) ap->firstslater = bp->firstslater; */
   }
}

int isEmpty(char *str)
{
    for(; *str; str++) if(!isspace(*str)) return 0;
    return 1;
}


static int read_energy_levels(void)
{

   typedef struct REPRES { char *degen; /* all symbols per line in one array */
                           char base_symb[20];
                           unsigned short n_symb;
                           struct REPRES *next;
                         } Repres; /* Representation as printed by ADF */
                                   /* (Subspecies) line in character table */

   long fpos;
   int nlevels, n_alpha, n_beta = 0, n_alpha_el, n_beta_el = 0, i, j;
   MolecularOrbital *alphaOrb = NULL, *betaOrb = NULL;
   float nel;
   double eigen;
   int orbnumber;
   char *cp, str[20],  symmetry[20], printline[30];
   Repres *first_repres = NULL, *last_repres = NULL, *temp;

/* read character table as printed in ADF */
/* memory dynamic -> most flexible */
/* symbol required to extract orbital with densf */

   rewind(fp);
   if(!find_string("Subspecies")) return 0;
    do {
        fpos = ftell(fp);
    } while(find_string("Subspecies"));
    fseek(fp, fpos, SEEK_SET);
    fgets(line, 255, fp);
    fgets(line, 255, fp);
    while(!isEmpty(line)) {
       if((temp = (Repres *)malloc(sizeof(Repres))) == NULL) {
          printf("can't allocate enough memory\n");
          exit(-1);
       }
       temp->next = NULL;
       cp = line; /* count number of symbols per line */
       cp = strstr(cp, "  ");
       temp->n_symb = 1;
       while(cp != NULL) {
         cp = strstr(cp+2, "  ");
         temp->n_symb++;
       }
       if((temp->degen = (char *)malloc(temp->n_symb*20*sizeof(char))) == NULL) {
           printf(" can't allocate enough memory\n");
           exit(-1);
       } 
       cp = line; /* determine base symbol */
       strncpy(str, cp, strcspn(cp, ":\0"));
       sscanf(str, "%s", temp->base_symb);
       strncpy(str,"", 20); /* fill str with \0 */
       sscanf(cp, "%s", temp->degen); /* fill array with symbols */
       cp = strstr(cp, "  ");
       for(i=1; i < temp->n_symb; i++) {
         sscanf(cp, "%s", temp->degen+i*20);
         cp = strstr(cp+2, "  ");
       }
       if(!first_repres) first_repres = temp;
       else last_repres->next = temp;
       last_repres = temp;
       fgets(line, 255, fp);
    }

    
   rewind(fp);
   if(find_string("SPIN 1")) actualmol->alphaBeta = 1; /* restricted or unrestricted */

   rewind(fp);
   if(actualmol->alphaBeta) {
     while(find_string("  ***** SPIN 1")) {
        fgets(line, 255, fp);
        fgets(line, 255, fp);
//        if(!strstr(line,"Orbital Energies, per Irrep and Spin:") == NULL) fpos = ftell(fp);
//        if(!strstr(line,"Orbital Energies, per Irrep and Spin:") == NULL) logprint("spin 1");
        if(strstr(line,"Orbital Energies, per Irrep and Spin:")) fpos = ftell(fp);
        if(strstr(line,"Orbital Energies, per Irrep and Spin:")) logprint("spin 1");
        else continue;
     }
   }
   else {
      if(!find_string("Orbital Energies, per Irrep and Spin:")) return 0;
      do {
         fpos = ftell(fp);
      } while(find_string("Orbital Energies, per Irrep and Spin:"));
   }

   fseek(fp, fpos, SEEK_SET);
   fgets(line, 255, fp);
   fgets(line, 255, fp);
   fgets(line, 255, fp);
   fpos = ftell(fp);

   n_alpha = 0; /* count levels */
   while(1) {
      if(!fgets(line, 255, fp)) return 0;
      if(isEmpty(line)) break;
      if(strlen(line) < 16) {
        if(sscanf(line, "%s", symmetry) != 1) return 0;
        for(temp = first_repres; temp; temp = temp->next) { /* determine representation */
           if(!strcmp(temp->base_symb, symmetry)) break;
        }
      }
// in case of mixed output file (create and compute jobs) temp needs a test
// in order not to crash
      if(temp) if(strlen(line) >= 16) n_alpha += temp->n_symb;  /* count */
   }

   nlevels = n_alpha;
   if((alphaOrb = allocOrbital(n_alpha, 0, ADF_ORB_A)) == NULL){ /* changed from ADF_ORB */
      showinfobox("can't allocate the MO-structures");
      return 0;
   }

   fseek(fp, fpos, SEEK_SET); /* read levels */
   n_alpha_el = 0; line[0] = '\0';
   for(i = 0; i<n_alpha;){
      fgets(line, 255, fp);
      if(isEmpty(line)) break;
      if(strlen(line) < 16) {
	 if(sscanf(line, "%s", symmetry) != 1) return 0;
         for(temp = first_repres; temp; temp = temp->next) { /* determine representation */
            if(!strcmp(temp->base_symb, symmetry)) break;
         }
	 continue;
      }
      if(sscanf(line, "%d%f%*f%lf",
         &orbnumber, &nel, &eigen) != 3) return 0;
      for(j = 0; j < temp->n_symb; j++) { /* read and doublicate */
         alphaOrb[i].number = orbnumber;
         alphaOrb[i].occ = nel/temp->n_symb;
         alphaOrb[i].eigenvalue = eigen;
         strncpy(alphaOrb[i].type, temp->degen+j*20, 7);
         n_alpha_el += (int)alphaOrb[i].occ;
         i++;
      }
   }

   if(actualmol->alphaBeta) { /* the same for beta electrons */
     if(!find_string("  ***** SPIN 2")) return 0;
     if(!find_string("Orbital Energies, per Irrep and Spin:")) return 0;
     fgets(line, 255, fp);
     fgets(line, 255, fp);
     fgets(line, 255, fp);
     fpos = ftell(fp);
     n_beta = 0;
     while(1) {
        if(!fgets(line, 255, fp)) return 0;
        if(isEmpty(line)) break;
        if(strlen(line) < 16) {
          if(sscanf(line, "%s", symmetry) != 1) return 0;
          for(temp = first_repres; temp; temp = temp->next) {
             if(!strcmp(temp->base_symb, symmetry)) break;
          }
        }
        if(strlen(line) >= 16) n_beta += temp->n_symb;
     }

     nlevels = (n_alpha >= n_beta ? n_alpha : n_beta); /* free_mo requires equal number of */
                                              /* alpha and beta levels, take biger value */
                                     /* results in dummy levels at the end of the list */
     if((betaOrb = allocOrbital(nlevels, 0, ADF_ORB_B)) == NULL){
        showinfobox("can't allocate the MO-structures");
        return 0;
     }

     fseek(fp, fpos, SEEK_SET);
     n_beta_el = 0; line[0] = '\0';
     for(i = 0; i<n_beta;){
        fgets(line, 255, fp);
        if(isEmpty(line)) break;
        if(strlen(line) < 16) {
           if(sscanf(line, "%s", symmetry) != 1) return 0;
           for(temp = first_repres; temp; temp = temp->next) {
              if(!strcmp(temp->base_symb, symmetry)) break;
           }
           continue;
        }
        if(sscanf(line, "%d%f%*f%lf",
           &orbnumber, &nel, &eigen) != 3) return 0;
          for(j = 0; j < temp->n_symb; j++) {
           betaOrb[i].number = orbnumber;
           betaOrb[i].occ = nel/temp->n_symb;
           betaOrb[i].eigenvalue = eigen;
           strncpy(betaOrb[i].type, temp->degen+j*20, 7);
           n_beta_el += (int)betaOrb[i].occ;
           i++;
        }
     }
   }


   actualmol->alphaOrbital = alphaOrb;
   if(actualmol->alphaBeta) actualmol->betaOrbital = betaOrb;
   actualmol->nMolecularOrbitals = nlevels;
   actualmol->nElectrons = n_alpha_el + n_beta_el;

   logprint("");
   sprintf(printline, "%d energy levels", n_alpha + n_beta);
   logprint(printline);
   sprintf(printline, "%d electrons", actualmol->nElectrons);
   logprint(printline);

   qsort(actualmol->alphaOrbital,
         n_alpha,
         sizeof(MolecularOrbital),
         compare_orbitals);
   qsort(actualmol->betaOrbital,
         n_beta,
         sizeof(MolecularOrbital),
         compare_orbitals);

   return 1;
}

static int read_atomic_charges(void)
{
   long fpos;
   char sym[3];
   int tmp, nitems;
   AtoM *ap;

   if(!find_string("Atom              Charge")) return 0;
    do {
        fpos = ftell(fp);
    } while(find_string("Atom              Charge"));
    fseek(fp, fpos, SEEK_SET);

   fgets(line, 255, fp);

   for(ap=actualmol->firstatom; ap; ap=ap->next){
      fgets(line, 255, fp);
      if(sscanf(line, "%d", &tmp) != 1) fgets(line, 255, fp); /* unrestricted: 2 lines per atom */
      nitems = sscanf(line, "%*d%s%f", sym, &ap->charge);
      if(ap->ord != get_ordinal(sym)) {
         showinfobox("read_atomic_charges: inconsistent atom numbering!");
         return 0;
      }
   }

   actualmol->charges = 1;
   return 1;
}

static void read_freq(void)
{
   long fpos;
   int i, j, n, n_freq = 0;
   float f1, f2, f3;
   Vibration *vib;
   AtoM *ap;


   if(!find_string("Vibrations and Normal Modes")) return;
   for (i = 0; i < 6; i++) fgets(line, 255, fp);
   fpos = ftell(fp);
   fgets(line, 255, fp);
   n = sscanf(line, "%f %f %f", &f1, &f2, &f3);
   if(n > 0) n_freq += n;
   while(!isEmpty(line)) {
      while(!isEmpty(line)) {
         fgets(line, 255, fp);
      }
      fgets(line, 255, fp);
      fgets(line, 255, fp);
      n = sscanf(line, "%f %f %f", &f1, &f2, &f3);
      if(n > 0) n_freq += n;
   }

   if((actualmol->vibration = (Vibration*) malloc(n_freq*sizeof(Vibration))) == NULL){
      showinfobox("can't realloc vibrational frequency");
      return;
   }
   for(i=0, vib=actualmol->vibration; i<n_freq; i++, vib++){
      if((vib->coord = (Vector*) malloc(actualmol->natoms*sizeof(Vector))) == NULL){
         showinfobox("can't allocate vibration");
         return;
      }
   }

   fseek(fp, fpos, SEEK_SET);
   fgets(line, 255, fp);
   for(i=0, vib=actualmol->vibration; i<n_freq/3; i++, vib += 3) {
      actualmol->n_frequencies += sscanf(line, "%f %f %f", &vib->frequency,
         &(vib+1)->frequency, &(vib+2)->frequency);
      fgets(line, 255, fp);
      for(ap=actualmol->firstatom, j=0; ap; ap=ap->next, j++){
         fgets(line, 255, fp);
         sscanf(line, "%*s %f %f %f %f %f %f %f %f %f", &vib->coord[j].x,
            &vib->coord[j].y, &vib->coord[j].z, &(vib+1)->coord[j].x,
            &(vib+1)->coord[j].y, &(vib+1)->coord[j].z, &(vib+2)->coord[j].x,
            &(vib+2)->coord[j].y, &(vib+2)->coord[j].z);
      }
      fgets(line, 255, fp);
      fgets(line, 255, fp);
      fgets(line, 255, fp);
   }

   logprint("frequencies present");
}

int read_dipole(void)
{
   long fpos;
   float x, y, z;

   rewind(fp);
   if(!find_string("Dipole Moment")) return 0;
   fpos = ftell(fp);
   while(find_string("Dipole Moment")) fpos = ftell(fp);
   fseek(fp, fpos, SEEK_SET);
   if(!find_string("Vector   :")) return 0;

   sscanf(line, "%*s %*s %f %f %f", &x, &y, &z);
   actualmol->dipole = add_dipole(x, y, z);

   return 1;
}

Basis *get_basis(char *basis)
{
   register Basis *bp;

   for(bp = actualmol->basisset; bp; bp = bp->next){
      if(!strcmp(basis, bp->basisname)) return bp;
   }
   return NULL;
}

int compare_orbitals(const void *aa, const void *bb)
{
    MolecularOrbital *a, *b;
    
    a = (MolecularOrbital *)aa;
    b = (MolecularOrbital *)bb;

   if(a->eigenvalue < b->eigenvalue) return -1;
   if(a->eigenvalue > b->eigenvalue) return  1;

   return strcmp(a->type, b->type);
}

/**************************************
 * 
 *  reading of TAPE41
 * 
 */
 


int getIntValue(char *str, int *ival)
{
    if(!find_string(str)) {
	fprintf(stderr, "can't find %s\n", str);
	return 0;
    }
    fgets(line, 255, fp);
    fgets(line, 255, fp);
    sscanf(line, "%d", ival);
    return 1;
}

 
void tape41ToMacu(char *tape41, char *macu, char *density)
{
    long fpos, fpos2;
    FILE *fout;
    char *adfbin, command[512], asciiFile[256], *spin = NULL;
    float dim[6], x, y, z, value, value2;
    int n[3], len, i, j;

    if((fout = fopen(macu, "w")) == NULL){
	fprintf(stderr, "can't open .macu file\n");
	sprintf(timestring, "can't open .macu file\n");
	return;
    }

    if((adfbin = getenv("ADFBIN")) == NULL) {
        fclose(fout);
        unlink(macu);
	fprintf(stderr, " can't find ADFBIN variable\n");
	sprintf(timestring, " can't find ADFBIN variable\n");
	return;
    }
    sprintf(asciiFile, "%s.ascii", tape41);
    sprintf(command, "%s/dmpkf %s > %s\n", adfbin, tape41, asciiFile);
    system(command);
    
    if((fp = fopen(asciiFile, "r")) == NULL) {
        fclose(fout);
        unlink(macu);
	fprintf(stderr, "can't open %s\n", command);
	sprintf(timestring, "can't open %s\n", command);
	return;
    }

    if(!find_string("Start_point")) {
        fclose(fout);
        unlink(macu);
	fprintf(stderr, "can't find start-point\n");
	sprintf(timestring, "can't find start-point\n");
	return;
    }
    fgets(line, 255, fp);
    fgets(line, 255, fp);
    sscanf(line, "%f%f%f", dim, dim+2, dim+4);

    if(!getIntValue("nr of points x", n)) return;
    if(!getIntValue("nr of points y", n+1)) return;
    if(!getIntValue("nr of points z", n+2)) return;

    if(!find_string("x-vector")) {
        fclose(fout);
        unlink(macu);
	fprintf(stderr, "can't find x-vector\n");
	sprintf(timestring, "can't find x-vector\n");
	return;
    }
    fgets(line, 255, fp);
    fgets(line, 255, fp);
    sscanf(line, "%f", &x);

    if(!find_string("y-vector")) {
        fclose(fout);
        unlink(macu);
	fprintf(stderr, "can't find y-vector\n");
	sprintf(timestring, "can't find y-vector\n");
	return;
    }
    fgets(line, 255, fp);
    fgets(line, 255, fp);
    sscanf(line, "%*f%f", &y);

    if(!find_string("z-vector")) {
        fclose(fout);
        unlink(macu);
	fprintf(stderr, "can't find z-vector\n");
	sprintf(timestring, "can't find z-vector\n");
	return;
    }
    fgets(line, 255, fp);
    fgets(line, 255, fp);
    sscanf(line, "%*f%*f%f", &z);

    dim[1] = dim[0] + (n[0]-1)*x;
    dim[3] = dim[2] + (n[1]-1)*y;
    dim[5] = dim[4] + (n[2]-1)*z;
    for(i=0; i<6; i++) dim[i] *= BOHR;

    len = 36;
    fwrite(&len, sizeof(int), 1, fout);
    fwrite(dim, sizeof(float), 6, fout);
    fwrite(n, sizeof(int), 3, fout);
    fwrite(&len, sizeof(int), 1, fout);

    if(strstr(density, "Spin")) spin = strcpy(density, "Density");
    if(!find_string(density)) {
        fclose(fout);
        unlink(macu);
	fprintf(stderr, "can't find Density\n");
	sprintf(timestring, "can't find Density\n");
	return;
    }

    if(strstr(density, "Density")) {
        if(strstr(line, "Density_A")) { /* in unrestricted calc. density is */
                                        /* printed for spin A and spin b seperately */
                                        /* add values for density, substract for */
                                        /* spin density */
                                        /* read A and B at the same time, jump with */
                                        /* two pointers. Maybe not very efficient but */
                                        /* straight forward to implement */

          fgets(line, 255, fp);
          fpos = ftell(fp);
          if(!find_string("Density_B")) {
              fclose(fout);
              unlink(macu);
              fprintf(stderr, "can't find Density\n");
              sprintf(timestring, "can't find Density\n");
              return;
          }
          fgets(line, 255, fp);
          fpos2 = ftell(fp);
          len = 4*n[0]*n[1]; 

          for(i=0; i<n[2]; i++) {
              fwrite(&len, sizeof(int), 1, fout);
              for(j=0; j<n[0]*n[1]; j++) {
                  fseek(fp, fpos, SEEK_SET);
                  if(!fscanf(fp, "%f", &value)) {
                      fclose(fout);
                      unlink(macu);
                      fprintf(stderr, "can't read value %d in plane %d\n", j, i);
                      sprintf(timestring, "can't read value %d in plane %d\n", j, i);
                      return;
                  }
                  fpos = ftell(fp);
                  fseek(fp, fpos2, SEEK_SET);
                  if(!fscanf(fp, "%f", &value2)) {
                      fclose(fout);
                      unlink(macu);
                      fprintf(stderr, "can't read value %d in plane %d\n", j, i);
                      sprintf(timestring, "can't read value %d in plane %d\n", j, i);
                      return;
                  }
                  fpos2 = ftell(fp);
                  if(spin) { /* subtract -> spin */
                  value -= value2;
                  fwrite(&value, sizeof(float), 1, fout);
                  }
                  else { /* add -> density */
                  value += value2;
                  fwrite(&value, sizeof(float), 1, fout);
                  }
              }
              fwrite(&len, sizeof(int), 1, fout);
          }

          fclose(fout); /* clean exit from tape41tomacu if successful */
          fclose(fp);
          unlink(asciiFile);
          return;
        }
        fgets(line, 255, fp);
    }
    else {
        do {
            fpos = ftell(fp);
        } while(find_string(strcat(density," ")));
//        } while(find_string(density));
        fseek(fp, fpos, SEEK_SET);
        fgets(line, 255, fp);
        fgets(line, 255, fp);
    }

    len = 4*n[0]*n[1];

    for(i=0; i<n[2]; i++) {
	fwrite(&len, sizeof(int), 1, fout);
	for(j=0; j<n[0]*n[1]; j++) {
	    if(!fscanf(fp, "%f", &value)) {
                fclose(fout);
                unlink(macu);
		fprintf(stderr, "can't read value %d in plane %d\n", j, i);
		sprintf(timestring, "can't read value %d in plane %d\n", j, i);
		return;
	    }
	    fwrite(&value, sizeof(float), 1, fout);
	}
	fwrite(&len, sizeof(int), 1, fout);
    }

    fclose(fout);
    fclose(fp);
 
    unlink(asciiFile);
}

/* routines to start "densf" to compute grids externally */
/*
 * not neaded anymore
static int getLabel(MolecularOrbital *mp)
{
    int label;
    MolecularOrbital *orb;

    label = 1;

    if(mp->flag == ADF_ORB_A) { 
       for(orb=actualmol->alphaOrbital; orb<mp; orb++) {
   	if(!strcmp(orb->type, mp->type)) label++;
       }
    }
    else {
           for(orb=actualmol->betaOrbital; orb<mp; orb++) {
        if(!strcmp(orb->type, mp->type)) label++;
       }
    }
    return label;
}
*/

#ifdef WIN32
void executeAdfUtilities(char *macuName, float *dim, int *nc, int key)
{
    char *adfbin, command[1024], keyStr[512], spin[20], density[20];
//    char currentDirectory[MAXPATHLEN];
    char *currentDirectory;
    char printline[100], adf_alpha_string[6];
    char *s;
    extern char tape21file[256];
    struct stat buf;
    FILE *fp;
    
    strcpy(adf_alpha_string, "ALPHA");
    if((s = getenv("MOLEKEL_ADF_ALPHA_KEY")) != NULL) {
      strcpy(adf_alpha_string, s);
    }

    sprintf(timestring, " Done\n");

    unlink(macuName);

/* on WIN32 create just the inputfile in C:/TEMP */


    fp = fopen("C:\\TEMP\\input", "w");
    if(!fp) {
       fprintf(stderr,"can't open input file for densf C:\\TEMP\\input\n");
       sprintf(timestring,"can't open input file for densf C:\\TEMP\\input\n");
       return;
    }

    if(key == EL_DENS) {
	sprintf(keyStr, "DENSITY scf");
        strcpy(density, "Density"); //density = "Density";
    }
    else if(key == CALC_ORB) {
	if(molOrb->flag == ADF_ORB_A) strcpy(spin, adf_alpha_string);
        else strcpy(spin ,"BETA");
	sprintf(keyStr, "ORBITALS SCF\n%s\n%s %d\nEND", spin,
//      molOrb->type, getLabel(molOrb));
	molOrb->type, molOrb->number);
        strcpy(density, "SCF_"); //density = "SCF_";
        strcat(density, molOrb->type);
        if(actualmol->alphaBeta) {
          if(molOrb->flag == ADF_ORB_A) strcat(density, "_A");
          else strcat(density, "_B");
        }
    }
    else if(key == SPIN_DENS) { /* spin denisty for unrestricted */
        sprintf(keyStr, "DENSITY scf");
        strcpy(density, "Spin"); //density = "Spin";
    }

/* The next line makes no sence, but it helps to get rid of a "A"
   that confuses the next print statement. So lets print nonsense first. */

//printf("%s\n", density);

    if((getenv("MOLEKEL_ADF2x")) == NULL) {
       fprintf(fp, "\
UNITS\n\
LENGTH Bohr\n\
END\n");
    }

    fprintf(fp, "\
GRID\n\
%f, %f, %f\n\
%d, %d, %d\n\
1.0, 0.0, 0.0, %f\n\
0.0, 1.0, 0.0, %f\n\
0.0, 0.0, 1.0, %f\n\
END\n\
%s\n\
END INPUT\n",
    dim[0]*_1_BOHR, dim[2]*_1_BOHR, dim[4]*_1_BOHR,
    nc[0], nc[1], nc[2], 
    (dim[1]-dim[0])*_1_BOHR, (dim[3]-dim[2])*_1_BOHR,
    (dim[5]-dim[4])*_1_BOHR, keyStr);
    fclose(fp);

    sprintf(timestring, "\
                the input has been saved in C:\\TEMP\\input");
}
#else

void executeAdfUtilities(char *macuName, float *dim, int *nc, int key)
{
    char *adfbin, command[1024], keyStr[512], spin[20], density[20];
//    char currentDirectory[MAXPATHLEN];
    char *currentDirectory;
    char printline[100], adf_alpha_string[6];
    char *s;
    extern char tape21file[256];
    struct stat buf;
    FILE *fp;
    
//  the keyword for alpha changed in ADF2002.02 from ALFA to ALPHA!
//  keyword can be set to the old value for use with an old ADF version
//  by setting the environmental variable MOLEKEL_ADF_ALPHA_KEY to "ALFA"
    strcpy(adf_alpha_string, "ALPHA");
    if((s = getenv("MOLEKEL_ADF_ALPHA_KEY")) != NULL) {
      strcpy(adf_alpha_string, s);
    }

    sprintf(timestring, " Done\n");

    unlink(macuName);

//    if(!getwd(currentDirectory)) {
    if(!(currentDirectory = getcwd(NULL, MAXPATHLEN + 1))) {
        sprintf(printline, "getcwd error: %s", currentDirectory);
	fprintf(stderr, printline);
	strcpy(timestring, printline);
	return;
    }

/* densf looks for a file named TAPE21 
   in the current working directory! 
   So let's do it in /usr/tmp and make a symb. link to the t21 file
*/

    if(stat("/usr/tmp/TAPE21", &buf) == 0) {
	fprintf(stderr,"Can't start the grid computation:\n\
	/usr/tmp/TAPE21 already exists -\n\
	make sure no other densf-job is running,\n\
	then remove the file /usr/tmp/TAPE21\n\
	before restarting the grid computation.\n");
	sprintf(timestring,"\
        can't start the grid computation:\n\
        /usr/tmp/TAPE21 already exists -\n\
        make sure no other densf-job is running,\n\
        then remove the file /usr/tmp/TAPE21\n\
        before restarting the grid computation.\n");
        return;
    }

    if(chdir("/usr/tmp") == -1) {
	fprintf(stderr,"can't change dir to /usr/tmp\n");
	sprintf(timestring,"can't change dir to /usr/tmp\n");
	return;
    }
    sprintf(command, "ln -s %s TAPE21\n", tape21file);
    system(command);

    fp = fopen("/usr/tmp/input", "w");
    if(!fp) {
       fprintf(stderr,"can't open input file for densf /usr/tmp/input\n");
       sprintf(timestring,"can't open input file for densf /usr/tmp/input\n");
       return;
    }

    if(key == EL_DENS) {
	sprintf(keyStr, "DENSITY scf");
        strcpy(density, "Density"); //density = "Density";
    }
    else if(key == CALC_ORB) {
	if(molOrb->flag == ADF_ORB_A) strcpy(spin, adf_alpha_string); // alpha or beta orbital
        else strcpy(spin ,"BETA");
	sprintf(keyStr, "ORBITALS SCF\n%s\n%s %d\nEND", spin,
//        molOrb->type, getLabel(molOrb));
	molOrb->type, molOrb->number);
        strcpy(density, "SCF_"); //density = "SCF_";
        strcat(density, molOrb->type);
        if(actualmol->alphaBeta) {
          if(molOrb->flag == ADF_ORB_A) strcat(density, "_A");
          else strcat(density, "_B");
        }
    }
    else if(key == SPIN_DENS) { /* spin denisty for unrestricted */
        sprintf(keyStr, "DENSITY scf");
        strcpy(density, "Spin"); //density = "Spin";
    }

/* The next line makes no sence, but it helps to get rid of a "A"
   that confuses the next print statement. So lets print nonsense first. */

//printf("%s\n", density);

    if((getenv("MOLEKEL_ADF2x")) == NULL) {
       fprintf(fp, "\
UNITS\n\
LENGTH Bohr\n\
END\n");
    }

    fprintf(fp, "\
GRID\n\
%f, %f, %f\n\
%d, %d, %d\n\
1.0, 0.0, 0.0, %f\n\
0.0, 1.0, 0.0, %f\n\
0.0, 0.0, 1.0, %f\n\
END\n\
%s\n\
END INPUT\n",
    dim[0]*_1_BOHR, dim[2]*_1_BOHR, dim[4]*_1_BOHR,
    nc[0], nc[1], nc[2], 
    (dim[1]-dim[0])*_1_BOHR, (dim[3]-dim[2])*_1_BOHR,
    (dim[5]-dim[4])*_1_BOHR, keyStr);
    fclose(fp);

    sprintf(command, "cat < /usr/tmp/input\n");
    printf(command);
    system(command);
    printf("\n");

    if((adfbin = getenv("ADFBIN")) == NULL) {
	fprintf(stderr,"can't find ADFBIN variable\n");
        unlink("/usr/tmp/TAPE21");
        sprintf(command, "mv /usr/tmp/input %s/input.adf", currentDirectory);
        sprintf(printline, "%s/input.adf", currentDirectory);
        if(system(command) != 0) {
        sprintf(timestring, "\
                can't find ADFBIN variable\n\
                set ADFBIN to point to the directory\n\
                containing \"densf\"\n\
                the input has been saved in /usr/tmp/input");
        }
        else {
        sprintf(timestring, "\
                can't find ADFBIN variable\n\
                set ADFBIN to point to the directory\n\
                containing \"densf\"\n\
                the input has been saved in %s", printline);
        }
//        unlink("/usr/tmp/input");
	return;
    }

    sprintf(command, "%s/densf < /usr/tmp/input\n", adfbin);
    printf(command);
    system(command);

/** remove the symbolic link **/
    unlink("/usr/tmp/TAPE21");

/** remove the input file **/
    unlink("/usr/tmp/input");
    unlink("/usr/tmp/logfile");
    
/* the result is in /usr/tmp/TAPE41, read the data and write into 
    the .macu file
*/
    tape41ToMacu("/usr/tmp/TAPE41", macuName, density);
    unlink("/usr/tmp/TAPE41");
    unlink("/usr/tmp/TAPE41.ascii");

    if(chdir(currentDirectory) == -1) {
        sprintf(printline, "can't change back to %s", currentDirectory);
	fprintf(stderr, printline);
	strcpy(timestring, printline);
	return;
    }
}

#endif

