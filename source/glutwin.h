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


#include <GLUT/glut.h>
extern float x_offset, y_offset;
extern int xorig, yorig, xsize, ysize;
extern int mainwin, filewin, orbwin, freqwin;
extern GLdouble clipping_plane[];
extern GLfloat clearcolor[4];
extern GLfloat fogstart, fogend;

/*----------------------------------------*/

void glutwin(int *, char **);
void init_GL(void);
void init_persp(void);
void clip_persp(void);
void reshape(int w, int h);
void dummy_gluiMouseFunc(int button, int state, int x, int y);
void dummy_gluiKeyboardFunc(unsigned char key, int x, int y);
void dummy_gluiSpecialFunc(int key, int x, int y);

