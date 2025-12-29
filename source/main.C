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


#include "version.h"
#include "main.h"
#include "maininterf.h"
#include "molekel.h"
#include "general.h"
#include "glutwin.h"
#include "readpdb.h"
#include "macu.h"
#include "constant.h"
#include "browser.h"
/*----------------------------------------*/

/* MAIN PROG */
int main(int argc, char **argv)
{
   char *argfname1 = NULL;
   char *argfname2 = NULL;
   int argcp = argc;
   int key = 0;
   char str[30], infostring[400];

   read_table();

/* handle command line arguments */
   while(--argcp > 0){
     if(argv[argcp][0] != '-') {
        if(!argfname1) argfname1 = argv[argcp];
/*
        else if(!argfname2) {
           argfname2 = argfname1;
           argfname1 = argv[argcp];
        }
*/
        else {
           printf("molekel: can only read one input file from the command line!\n");
           printf("Try `molekel -help' for more information.\n");
           return 1;
        }
     } else {
        if(strcmp(argv[argcp], "-display") == 0 || strcmp(argv[argcp], "-geometry") == 0) ;
        else if(strcmp(argv[argcp], "-iconic") == 0 || strcmp(argv[argcp], "-indirect") == 0) ;
        else if(strcmp(argv[argcp], "-direct") == 0 || strcmp(argv[argcp], "-gldebug") == 0) ;
        else if(strcmp(argv[argcp], "-sync") == 0) ;
        else if(strcmp(argv[argcp], "--help") == 0 || strcmp(argv[argcp], "-help") == 0 || strcmp(argv[argcp], "-h") == 0) {
           printf("Usage: molekel [glutInit options] [OPTION] [file type FILE]\n");
           printf("Advanced Interactive 3D-Graphics for Molecular Sciences\n\n");
           printf("glutInit options:   see the glut manual (-geometry is not recognized)\n\n");
           printf("OPTION:   \n");
           printf("  -h, -help, --help   display this help and exit\n");
           printf("  -version            output version information and exit\n\n");
           printf("file type: spezifies the format of the following file\n");
           printf("(only one file can be read from the command line at start up)\n");
           printf("recognized file types:\n");
           printf("  -xyz, -sxyz   xyz coordinate file (symbol x y z)\n");
           printf("  -dxyz         xyz coordinate file (ordinal x y z)\n");
           printf("  -xyzs         xyz coordinate file (x y z symbol)\n");
           printf("  -xyzd         xyz coordinate file (x y z ordinal)\n");
           printf("  -mxyz         xyz coordinate file containing multiple structures\n");
           printf("  -pdb          PDF file\n");
           printf("  -mpdb         multiple structure PDB file\n");
           printf("  -gau          gaussian log file\n");
           printf("  -gcb          gaussian cube file (reads coordinates only at start up)\n");
           printf("  -gam          gamess log file\n");
           printf("  -hon          hondo log file\n");
           printf("  -adf          adf log file\n");
           printf("  -mos          mos log file\n");
           printf("  -prd          prddo log file\n");
           printf("  -zin          zindo log file\n");
           printf("  -mld          molden file\n");
           printf("  -mldf         molden frequency file\n");
           printf("  -mldg         molden geometry file\n");
           printf("  -mkl          molekel format file\n");
           printf("\nfor more detailed info see the manual available at http://www.cscs.ch/molekel/\n");
           return 0;
        } else if(strcmp(argv[argcp], "-version") == 0) {
           printf("MOLEKEL, Version %s, Date: %s\n", MOLEKEL_VERSION, MOLEKEL_VERSION_DATE);
           return 0;
        } else if(argfname1) {
           if(strcmp(argv[argcp], "-xyz") == 0 || strcmp(argv[argcp], "-sxyz") == 0) {
              key = LOAD_XYZ;
              xyzformat = 0;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-dxyz") == 0) {
              key = LOAD_XYZ;
              xyzformat = 1;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-xyzs") == 0) {
              key = LOAD_XYZ;
              xyzformat = 2;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-xyzd") == 0) {
              key = LOAD_XYZ;
              xyzformat = 3;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-mxyz") == 0) {
              key = LOAD_MULTI_XYZ;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-pdb") == 0) {
              key = LOAD_PDB;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-mpdb") == 0) {
              key = LOAD_MULTI_PDB;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-gau") == 0) {
              key = LOAD_GAUSS;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-gcb") == 0) {
              key = LOAD_GCUBE;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-gam") == 0) {
              key = LOAD_GAMESS;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-hon") == 0) {
              key = LOAD_HONDO;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-adf") == 0) {
              key = LOAD_ADF;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-mos") == 0) {
              key = LOAD_MOS;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-prd") == 0) {
              key = LOAD_PRDDO;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-zin") == 0) {
              key = LOAD_ZINDO;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-mld") == 0) {
              key = LOAD_MOLDEN;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-mldf") == 0) {
              key = LOAD_MOLDEN_FREQ;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-mldg") == 0) {
              key = LOAD_MOLDEN_GEOM;
              strcpy(filename,argfname1);
           } else if(strcmp(argv[argcp], "-mkl") == 0) {
              key = LOAD_MKL;
              strcpy(filename,argfname1);
           } else {
              printf("molekel: invalid option: %s\nTry `molekel -help' for more information.\n", argv[argcp]);
              return 1;
           }
        } else {
           printf("molekel: no file name spezified!\nTry `molekel -help' for more information.\n");
           return 1;
        }
     }
   }
   if(argfname1 && key == 0) {
      printf("molekel: file type not specified!\nTry `molekel -help' for more information.\n");
      return 1;
   }

    glutwin(&argc, argv);
//    makefilewin();
    makebrowsergui();
//    makeorbwin();
    makeorbgui();
//    makefreqwin();
    makefreqgui();
    maininterf();
    logprint("Welcome to MOLEKEL!");
    sprintf(str, "Version %s", MOLEKEL_VERSION);
    logprint(str);
    logprint("");
    sprintf(str, "Your right mouse button");
    logprint(str);
    sprintf(str, "is your menu!");
    logprint(str);

    if(argfname1) timer_cb(key);

/* enter main program loop */
    sprintf(infostring, "MOLEKEL, Version %s, Date: %s,\n",\
       MOLEKEL_VERSION, MOLEKEL_VERSION_DATE);
    strcat(infostring, "by Stefan Portmann, Copyright (C) 2002 CSCS/ETHZ\n");
    strcat(infostring, "(original IRIX GL implementation, concept and data structure\n");
    strcat(infostring, "by Peter F. Fluekiger, CSCS/UNI Geneva)\n");
    strcat(infostring, "MOLEKEL comes with ABSOLUTELY NO WARRANTY; for details see menu \"Warranty\".\n");
    strcat(infostring, "The binary code is available free of charge,\nbut is not in the public domain.\n");
    strcat(infostring, "See menu \"License\" for details on conditions and restrictions.\n");
    strcat(infostring, "Info: http://www.cscs.ch/molekel/\n");
/*
    strcat(infostring, "This is free software, and you are welcome to redistribute it\n");
    strcat(infostring, "under certain conditions; menu \"License\" for details.\n");
    strcat(infostring, "Info: http://igc.ethz.ch/molekel/\n");
*/
    printf("%s\n", infostring);
    printf("OpenGL Version: %s\n", glGetString(GL_VERSION));
    glGetError();
    glutMainLoop();
    return 0;
}

