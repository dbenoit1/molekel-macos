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


/* read any xyz coordinate format */
#include "main.h"
#include "constant.h"
#include "molekel.h"
#include "readxyz.h"
#include "general.h"
#include "utils.h"
#include "chooseinterf.h"
#include "maininterf.h"
#include "browser.h"



void read_xyz(char *file)
{
   char line[100];
   char symb[100];
   int nitems, ord, ordinal = 0;
   float x, y, z;
   FILE *fp;

   if((fp = fopen(file, "r")) == NULL){
      sprintf(line, "Can't open file\n%s !", file);
      showinfobox(line);
      return;
   }

   add_mol(file);

   while(fgets(line, 99, fp)) {
      switch(xyzformat) {
         case 0:
            nitems = sscanf(line, "%s %f %f %f", symb, &x, &y, &z);
            break;
         case 1:
            nitems = sscanf(line, "%d %f %f %f", &ord, &x, &y, &z);
            ordinal = 1;
            break;
         case 2:
            nitems = sscanf(line, "%f %f %f %s", &x, &y, &z, symb);
            break;
         case 3:
            nitems = sscanf(line, "%f %f %f %d", &x, &y, &z, &ord);
            ordinal = 1;
            break;
      }
      if(nitems == 4) {
         if(!ordinal) ord = get_ordinal(symb);
         add_atom(ord, x, y, z);
      }
   }
   fclose(fp);

   if(!actualmol->natoms){
      sprintf(line, "No atoms in %s !", file);
      showinfobox(line);
      return;
   }

   create_bonds();
   find_multiplebonds();
 
   new_mole(file);

   if(maininterfwin) update_logs();
}

/* readmultixyz: read xyz file with multiple structures (bable -oxyz)*/

void addXYZstep(FILE *fp, int natoms)
{
   long index, i; 
   float x, y, z;
   char line[100], symb[100];

   index = dynamics.ntotalsteps++;

   if(!dynamics.trajectory){
      if((dynamics.trajectory = (Vector**) malloc(sizeof(Vector *))) == NULL){
         showinfobox("can't allocate dyna pointer\n");
         return;
      }
   }
   else {
      if((dynamics.trajectory = (Vector**) realloc(dynamics.trajectory,
          dynamics.ntotalsteps*sizeof(Vector *))) == NULL){
          showinfobox("can't reallocate dyna pointer\n");
         dynamics.ntotalsteps--;
         free_dyna();
         return;
      }
   }

   if((dynamics.trajectory[index] = (Vector*) malloc(natoms*sizeof(Vector))) == NULL){
      sprintf(line, "can't allocate timestep dyna[%d]\n", index);
      showinfobox(line);
      dynamics.ntotalsteps--;
      free_dyna();
      return;
   }
   for(i=0; i<natoms; i++) {
      fgets(line, 99, fp);
      sscanf(line, "%s %f %f %f", symb, &x, &y, &z);
      dynamics.trajectory[index][i].x = x;
      dynamics.trajectory[index][i].y = y;
      dynamics.trajectory[index][i].z = z;
   }
}

void readmultixyz(char *file)
{
   static char name[200], line[100];
   static long fpos;
   int natoms = 0, ord, i;
   char symb[100];
   float x, y, z;

   FILE *fp;

   if((fp = fopen(file, "r")) == NULL){
      sprintf(line, "Can't open file\n%s !", file);
      showinfobox(line);
      return;
   }

   fgets(line, 99, fp);
   if((sscanf(line, "%d", &natoms)) != 1){
      sprintf(line, "number of atoms expexted in first line!");
      showinfobox(line);
   }
   rewind(fp);

   free_dyna();

   while(fgets(line, 99, fp)){
      if(line[0] == '\n') continue;
      sscanf(line, "%d", &i);
      if(i != natoms) {
         sprintf(line, "Error in:\n%s !", file);
         showinfobox(line);
         fclose(fp);
         return;
      }
      fgets(line, 99, fp);
      fpos = ftell(fp);
      addXYZstep(fp, natoms);
   }

   add_mol(name);
   dynamics.molecule = actualmol;
   dynamics.current = dynamics.ntotalsteps - 1;

   fseek(fp, fpos, SEEK_SET);

   for(i=0; i<natoms; i++) {
      fgets(line, 99, fp);
      sscanf(line, "%s %f %f %f", symb, &x, &y, &z);
      ord = get_ordinal(symb);
      add_atom(ord, x, y, z);
   }      

   fclose(fp);

   if(!actualmol->natoms){
      sprintf(line, "No atoms in\n%s !", file);
      showinfobox(line);
      return;
   }

   create_bonds();
   find_multiplebonds();

   new_mole(name);

   if(maininterfwin) update_logs();
}
