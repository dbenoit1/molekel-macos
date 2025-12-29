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


extern GLUI *coef_or_matglui, *overwrite_filegui, *job_prioritygui, *grid_or_dotglui;
extern GLUI *select_surfgui, *vmin_vmaxglui, *optionglui, *bondglui, *freqglui;
extern GLUI *playglui, *boxglui, *connollyglui, *labelglui, *renderglui, *dipoleglui;
extern GLUI_StaticText *freqtext, *play_no, *dipoletext;
extern GLUI_Button *coord_cancel;
extern int datasource;
extern int job_priority_live_var, surface_live_var, speed_live_var;
extern int vibrating, playing;

void select_coeff_or_matrix(void);
void select_overwrite_file(char *fname);
void select_job_priority(void);
void select_quit(void);
void select_grid_or_dot(void);
void select_surface(int key);
void set_vmin_vmax(Surface *surf);
void select_option(void);
void bond_interf(void);
void freq_interf(void);
void play_interf(void);
void box_interf(void);
void connolly_interf(void);
void triangglui_interf(void);
void label_interf(void);
void xyzformat_interf(void);
void t41content_interf(void);
void render_interf(void);
void dipole_interf(void);

