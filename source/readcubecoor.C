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


#include "main.h"
#include "molekel.h"
#include "general.h"
#include "readpdb.h"
#include "maininterf.h"
#include "utils.h"


void read_gcubecoor(char *file)
{
   int natoms;
   register int i;
   char line[132];
   FILE *fp;

   int ord;
   float x, y, z;

   if((fp = fopen(file, "r")) == NULL){
      sprintf(line, "Can't open file\n%s !", file);
      showinfobox(line);
      return;
   }


   fgets(line, 132, fp);
   fgets(line, 132, fp);
   fgets(line, 132, fp);
   sscanf(line, "%d%*f%*f%*f", &natoms);
   if( natoms < 0 ) natoms=-natoms;
   fgets(line, 132, fp);
   fgets(line, 132, fp);
   fgets(line, 132, fp);

   if(natoms < 1 ){
      sprintf(line, "No atoms in\n%s !", file);
      showinfobox(line);
      return;
   }

   add_mol(file);

   i=0;
   while (i<natoms) {
     i=i+1;
     fgets(line, 132, fp);
     sscanf(line, "%d%*f%f%f%f", &ord, &x, &y, &z);
     add_atom(ord, x * BOHR, y * BOHR, z * BOHR);
   }

   fclose(fp);

   create_bonds();
   find_multiplebonds();
   new_mole(file);

   logprint("only coordinates read");
   logprint("read grid using \"Surface\"");
   if(maininterfwin) update_logs();
}

