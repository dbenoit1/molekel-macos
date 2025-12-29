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


extern Texture *firsttexture;
extern Texture *actual_texture;
extern GLubyte *contour_image1;

extern Texture *bond_texture;
extern int bond_textype;
extern GLenum bond_texenv;

extern Texture *back_texture;
extern int back_textype;
extern GLenum back_texenv;

extern float start_texture, end_texture;

extern int assign_no_tex;
Texture *add_texture(GLubyte *image, int idepth, int iwidth, int iheight, GLenum *texprop);
void mod_texture(GLenum *texprop);
void init_texture(void);
void add_contour_textures(void);
void enable_texture(Texture *texture, GLenum texenv, int type);
void disable_texture(void);
void texture_up(void);
void texture_down(void);
void read_img(char *file);
void set_texture_range(void);
void get_pmat(void);
void myt3f(float a[3]);
