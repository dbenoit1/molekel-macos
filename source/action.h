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


/* event actions!
 * actions caused by the glutwin, menu or keyboard, start with a capital letter
*/

void Quit(void);
void Tgl_manip_all(void);
void Tgl_depth(void);
void Set_wire(void);
void Set_ball_and_stick(void);
void Set_stick(void);
void Set_spacefill(void);
void Load_pdb(void);
void Tgl_fullscreen(void);
void Tgl_spinning(void);
void Set_surf_transp(void);
void Set_chickenwire(void);
void Set_dot_surface(void);
void Set_flatshade(void);
void Tgl_idle(void);
void Set_glui_idle(void);
void Set_glut_idle(void);

