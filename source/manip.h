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


/* glx window manipulations */

extern int picking_mode, spinning;

void mouse(int button, int state, int x, int y);
void picking_mouse(int button, int state, int x, int y);
void motion(int x, int y);
void picking_motion(int x, int y);
void spinmotion(int x, int y);
void keyboard(unsigned char key, int x, int y);
void spec_keyboard(int key, int x, int y);
void mouse_rotation(long *omx, long *omy, long *nmx, long *nmy);
void mouse_clip_rotation(long *omx, long *omy, long *nmx, long *nmy);
void idle_spin(void);
void idle_freq(void);
void idle_coord(void);
void spinmouse_rotation(long *omx, long *omy, long *nmx, long *nmy);
void mouse_translation(void);
void mouse_clip_translation(void);
void flip(char axis, float angle);
void u_to_world(long sx, long sy, float *wx, float *wy);
void new_trackballcenter(Mol *mp);
void switch_mol(void);
void switch_surf(void);
void get_translation(Matrix m, float *t);
void get_rotation(Matrix m, float *e);
void reset_globals(void);
void reset(void);
void center_of_mass();
void update_interface_flags(void);
