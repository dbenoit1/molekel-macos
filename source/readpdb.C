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


/*
   read a Brookhaven Protein Data Bank (PDB) file
   (extracted from general.h, improved by using a cubing
    algorithm for bond-generation)
*/
#ifndef WIN32
#include <sys/param.h>
#endif

#include "constant.h"
#include "main.h"
#include "molekel.h"
#include "general.h"
#include "readpdb.h"
#include "macu.h"
#include "drawing.h"
#include "maininterf.h"
#include "box.h"
#include "utils.h"

/*----------------------------------------*/
/* global variables */


/*----------------------------------------*/
/* local variables */

/*----------------------------------------*/

/* read the ".pdb"-file, pass the lines containing atomic coordinates */
void read_pdb(char *file)
{
   char line[100];
   FILE *fp;

   if((fp = fopen(file, "r")) == NULL){
      sprintf(line, "Can't open file\n%s !", file);
      showinfobox(line);
      return;
   }

   (void)add_mol(file);

   while(fgets(line, 99, fp)){
      if(!strncmp(line, "ATOM", 4))   decipher(line);
      if(!strncmp(line, "HETATM", 6)) decipher(line);
      if(!strncmp(line, "CONECT", 6)) get_conect(line);
   }
   (void)fclose(fp);

   if(!actualmol->natoms){
      sprintf(line, "No atoms in\n%s !", file);
      showinfobox(line);
      free_mol(actualmol);
      return;
   }

   create_bonds();
   find_multiplebonds();
 
   new_mole(file);

   if(maininterfwin) update_logs();
}




/* read the items in a line from a pdb_file */
void decipher(char *line)
{
   int num;
   char sym[3], hetflag, residue[20];
   float x, y, z, coord;
   AtoM *ap;

   (void)sscanf(line, "%*s %d", &num);
   (void)sscanf(line+12, "%c%c%c %s", sym, sym+1, &hetflag, residue);
   sym[2] = 0;
   coord = 0;
   (void)sscanf(line+30, "%f%f%f%f", &x, &y, &z, &coord);
/*   (void)printf("%3d %-2s %8.3f %8.3f %8.3f\n", num, sym, x, y, z); */

   num = get_ordinal(sym);

   ap = add_atom(num, x, y, z);

/*
   ap->coordination = coord;
*/
   
   if(ap->ord == 6 || ap->ord == 7) {
      if(hetflag == 'A' || hetflag == ' ') ap->main = 1;
   }
   if(ap->ord == 16 && !strcmp(residue, "CYS")) ap->main = 1;
}





void read_spartan (char *file)
{
   char line[100];
   FILE *fp;
   int ord;
   float x, y, z;

   if((fp = fopen(file, "r")) == NULL){
      sprintf(line, "Can't open file\n%s !", file);
      showinfobox(line);
      return;
   }

#ifdef REDUCED_VERSION
   if(firstmol) reduced_version();
   free_mol(firstmol);
#endif

   (void)add_mol(file);

   while(fgets(line, 99, fp)){
      if(sscanf(line, "%d%f%f%f", &ord, &x, &y, &z) == 4)
         break;
   }
   do {
      if(!strncmp(line, "ENDCART", 7)) break;
      (void)sscanf(line, "%d%f%f%f", &ord, &x, &y, &z);
      (void)add_atom(ord, x, y, z);
   } while(fgets(line, 99, fp));

   (void)fclose(fp);

   if(!actualmol->natoms){
      (void)sprintf(line, "No atoms in:\n %s !", file);
      showinfobox(line);
      return;
   }

//   create_box();
   create_bonds();
   find_multiplebonds();
 
   new_mole(file);
}



/* read multipel PDB */


void readmultipdb(char *file)
{
   static char line[100];
   static long fpos;
   int natoms = 0;

   FILE *fp;

   if((fp = fopen(file, "r")) == NULL){
      sprintf(line, "Can't open file\n%s !", file);
      showinfobox(line);
      return;
   }

   free_dyna();

   while(fgets(line, 99, fp)){
      if(!strncmp(line, "MODEL", 5)) {
         fpos = ftell(fp);
         natoms = addPDBstep(fp, fpos);
         if(!natoms) break;
      }
   }

   if(!natoms) {
      rewind(fp);
      while(fgets(line, 99, fp)){
          if((!strncmp(line, "ATOM", 4)) || (!strncmp(line, "HETATM", 6))) {
#ifdef WIN32
/* line contains also additional \r not counted in strlen */
              fpos = ftell(fp) - strlen(line) - 1;
#else
              fpos = ftell(fp) - strlen(line);
#endif
              natoms = addPDBstep(fp, fpos);
              if(!natoms) break;
          }
      }
   }
      
   if(!natoms) {
      sprintf(line, "Error in:\n%s !", file);
      showinfobox(line);
      fclose(fp);
      return;
   }

   add_mol(file);
   dynamics.molecule = actualmol;
   dynamics.current = dynamics.ntotalsteps - 1;

   fseek(fp, fpos, SEEK_SET);

   while(fgets(line, 99, fp)){
      if(!strncmp(line, "ATOM", 4))   decipher(line);
      if(!strncmp(line, "HETATM", 6)) decipher(line);
      if(!strncmp(line, "CONECT", 6)) get_conect(line);
   }

   fclose(fp);

   if(!actualmol->natoms){
      sprintf(line, "No atoms in\n%s !", file);
      showinfobox(line);
      return;
   }

   create_bonds();
   find_multiplebonds();

   new_mole(file);

   if(maininterfwin) update_logs();
}


int addPDBstep(FILE *fp, long fpos)
{
   long index, i; 
   float x, y, z;
   char line[100];
   static int natoms;

   index = dynamics.ntotalsteps++;

   fseek(fp, fpos, SEEK_SET);

   if(!dynamics.trajectory){
      if((dynamics.trajectory = (Vector**) malloc(sizeof(Vector *))) == NULL){
         showinfobox("can't allocate dyna pointer\n");
         return 0;
      }

/* count nr of atoms */
      natoms = 0;
      while(fgets(line, 99, fp)){
          if((!strncmp(line, "ATOM", 4)) || (!strncmp(line, "HETATM", 6))) natoms++;
          else if((!strncmp(line, "TER", 3)) || (!strncmp(line, "SIGATM", 6)) || \
                  (!strncmp(line, "ANISOU", 6)) || (!strncmp(line, "SIGUIJ", 6))) continue;
          else break;
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
   }

   if((dynamics.trajectory[index] = (Vector*) malloc(natoms*sizeof(Vector))) == NULL){
      sprintf(line, "can't allocate timestep dyna[%d]\n", index);
      showinfobox(line);
      dynamics.ntotalsteps--;
      free_dyna();
      return 0;
   }

   i = 0;
   while(fgets(line, 99, fp)) {
      if((!strncmp(line, "TER", 3)) || (!strncmp(line, "SIGATM", 6)) || \
          (!strncmp(line, "ANISOU", 6)) || (!strncmp(line, "SIGUIJ", 6))) continue;
      if(!((!strncmp(line, "ATOM", 4)) || (!strncmp(line, "HETATM", 6)))) break;
          if(sscanf(line+30, "%f %f %f", &x, &y, &z) != 3) return 0;
          dynamics.trajectory[index][i].x = x;
          dynamics.trajectory[index][i].y = y;
          dynamics.trajectory[index][i].z = z;
          i++;
   }

   if(i != natoms) natoms = 0;
   return natoms;

}

