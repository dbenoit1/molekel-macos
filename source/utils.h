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


/* utility functions */
void vecxmat(float a[4], Matrix m, float b[4]);
float determinant(float *a, float *b, float *c);
void invert4x4(Matrix m, Matrix n);
double adjunct4x4(Matrix m, short k, short l);
float dot_product(float *a, float *b);
void cross_product(float *a, float *b, float *c);
float triple_product(float *a, float *b, float *c);
void transpose(Matrix m);
float length(float *a);
void get_center(int n, float *c, AtoM **ap);
void set_scale_factor(void);
void clipplane_on(Surface *sp);
void clipplane_off(void);
void set_clipplane_matrix(Surface *sp, Mol *mp, float r[4], float t[4]);
double adjunct3x3(float m[3][3], int k, int l);
void invert3x3(float a[3][3], float b[3][3]);
void vec3mat(float a[3], float m[3][3], float b[3]);
int get_ordinal(char *sym);
void get_conect(char *line);
float getMaxCoval(void);
void create_bonds(void);
void create_Hbonds(void);
void createBonds(int key);
void do_atomcubes(int key);
void do_connections(int key);
void free_atomcubes(void);
void find_multiplebonds(void);
void computeOccupations(Mol *mp);
void normalize_amoss(void);
void normalize_gaussians(void);
double dszeta(int n, double c1, double z1, double c2, double z2);
void free_frequencies(void);
