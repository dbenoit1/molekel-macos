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


extern MolecularOrbital *molOrb;
extern int did_sigset, pipein, pipeout;
extern int child_action_key;
extern char timestring[];

void calcdens(int key);
void define_box(float *dim, int *ncub);
void sighandler(int sig);
void process_calc(char *s, float *dim, int *ncubes, int key);
double calc_point(float x, float y, float z);
double calculateSomo(float x, float y, float z);
double calculate_density(float x, float y, float z);
void calc_chi(float x, float y, float z);
int generate_density_matrix(int key);
void mep_dot_surface(Surface *surf);
void mep_dot_surface_no_pick(void);
double calc_mep(float x, float y, float z);
void spin_on_atoms(void);
double calc_prddo_point(float x, float y, float z);
double calc_prddo_density(float x, float y, float z);
double calc_prddo_spindensity(float x, float y, float z);
double calc_sltr_point(float x, float y, float z);
double calc_sltr_density(float x, float y, float z);
double calc_sltr_spindensity(float x, float y, float z);

