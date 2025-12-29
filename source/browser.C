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


/*
 * Copyright (c) 1993-1997, Silicon Graphics, Inc.
 * ALL RIGHTS RESERVED 
 * Permission to use, copy, modify, and distribute this software for 
 * any purpose and without fee is hereby granted, provided that the above
 * copyright notice appear in all copies and that both the copyright notice
 * and this permission notice appear in supporting documentation, and that 
 * the name of Silicon Graphics, Inc. not be used in advertising
 * or publicity pertaining to distribution of the software without specific,
 * written prior permission. 
 *
 * THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
 * AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
 * FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON
 * GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
 * SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
 * KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
 * LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
 * THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
 * POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
 * 
 * US Government Users Restricted Rights 
 * Use, duplication, or disclosure by the Government is subject to
 * restrictions set forth in FAR 52.227.19(c)(2) or subparagraph
 * (c)(1)(ii) of the Rights in Technical Data and Computer Software
 * clause at DFARS 252.227-7013 and/or in similar or successor
 * clauses in the FAR or the DOD or NASA FAR Supplement.
 * Unpublished-- rights reserved under the copyright laws of the
 * United States.  Contractor/manufacturer is Silicon Graphics,
 * Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.
 *
 * OpenGL(R) is a registered trademark of Silicon Graphics, Inc.
 */
#define GL_SILENCE_DEPRECATION
#include <GLUT/glut.h>
/*#include <mui/mui.h>*/
#include "mui/mui.h"
#include <dirent.h>
#include <string.h>
#include <sys/stat.h>
#include "main.h"
#include "molekel.h"
#include "browser.h"
#include "glutwin.h"
#include "general.h"
#include "constant.h"
#include "menu.h"
#include "readpdb.h"
#include "readgauss.h"
#include "maininterf.h"
#include "chooseinterf.h"
#include "surfaceinterf.h"
#include "calcdens.h"
#include "macu.h"
#include "readcubecoor.h"
#include "drawing.h"
#include "connect.h"
#include "output.h"
#include "snap.h"
#include "readxyz.h"
#include "readgamess.h"
#include "readhondo.h"
#include "readadf.h"
#include "texture.h"
#include "readmos.h"
#include "readprddo.h"
#include "readzindo.h"
#include "readmkl.h"
#include "textureinterf.h"
#include "readnbo.h"
#include "render.h"
#include "readmolden.h"
#include <tiffio.h>



#ifdef LINUX
#define NOFSTAT
#endif
#ifdef WIN32
#define MAXNAMLEN 512
#define NOFSTAT
#endif

/* settings from the original */
#define MINFILELIST 22

#ifndef WIN32
char SEPARATOR[] = "/";
#else
char SEPARATOR[] = "\\";
#endif

char	**filelist=NULL, **orblist, **freqlist;
char	err[80];
char	*dot = ".";
char	*dotdot = "..";
struct stat	d, dd;
struct dirent	*dir;

DIR	*file;
int	off;

int selectedfile = -1;
int cd(char *s);
void pwd(void);
void ls(void);

muiObject *tt, *tl, *vs, *vs_orb, *vs_freq, *l4, *trb1, *tl1, *tlf1;
muiObject *pl[10];


/*------------------------------*/

int winposition[] = {20, 80};
int file_select_key = 0;
char **file_select_pattern = NULL;
char filename[200];
int o_nlines, f_nlines, o_choice, f_choice;
Mol *molptr;

/*------------------------------*/

/* extensions */
char *pdb_ext[] = {".pdb", ".ent", NULL};
char *xyz_ext[] = {".xyz", NULL};
char *macu_ext[] = {".macu", NULL};
char *gcube_ext[] = {".cube", ".dat", NULL};
char *t21_ext[] = {".t21", NULL};
char *t41_ext[] = {".t41", ".ft41", NULL};
char *ms_ext[]  = {".ms", NULL};
char *out_ext[] = {".out", ".log", NULL};
char *pun_ext[] = {".pun", NULL};
char *mkl_ext[] = {".mkl", NULL};
char *molden_ext[] = {".molf", NULL};
char *sld_ext[] = {".sld", NULL};
char *eht_ext[] = {".eht", NULL};
char *ehtout_ext[] = {".ehtout", NULL};
char *qual_ext[] = {".qual", ".eint", ".nddo", ".mull", ".ct", ".rep", ".eorb", ".esol", ".eint2", NULL};
char *masq_ext[] = {".masq", NULL};
char *rgb_ext[] = {".rgb", ".bw", NULL};
char *tif_ext[] = {".tif", NULL};
char *spartan_ext[] = {".input", NULL};
char *mos_ext[] = {".out", ".mos", NULL};
char *nboorb_ext[] = {".orb", NULL};
char *ppm_ext[] = {".ppm", NULL};
char *jpeg_ext[] = {".jpg", NULL};
char *dotval_ext[] = {".dtv", NULL};

/*------------------------------*/

char *adjstr(char *dir)
{
   int len = strlen(dir);
   static char rstr[128];
   if(len < 52) return dir;
   else {
      sprintf(rstr, ".....%s", dir+(len - 47));
      return rstr;
   }
   
}

void update_interf_wread(void)
{
   char str[30];

   if(freqglui) {
      freqtext->set_text("no frequency selected");
      bit.show_freq_arrow = actualmol->show_freq_arrow;
      freqglui->sync_live();
   }
   if(planeglui) {
      planeglui->close();
      planeglui = NULL;
/*
      bit.cubeplanes = actualmol->cubeplanes;
      planeglui->sync_live();
*/
   }
   if(dipoleglui) {
      if(actualmol && actualmol->dipole) 
         sprintf(str, "dipole moment: %8.4f Debye", actualmol->dipole->absolute);
      else sprintf(str, "no dipole moment info");
      dipoletext->set_text(str);
      dipoleglui->sync_live();
   }
}

void timer_cb(int key)
/* used to force a return to glutMainLoop when reading
 * a larg file. If read_*** is called direclty in file_action
 * nothing happens on the screen, and this is confusing
*/
{
   green_note("reading/writing");
   switch(key) {
      case LOAD_XYZ:
         read_xyz(filename);
         update_interf_wread();
      break;
      case LOAD_MULTI_XYZ:
         readmultixyz(filename);
         update_interf_wread();
      break;
      case LOAD_PDB:
         read_pdb(filename);
         update_interf_wread();
      break;
      case LOAD_MULTI_PDB:
         readmultipdb(filename);
         update_interf_wread();
      break;
      case LOAD_GAUSS:
         read_gauss(filename);
         update_interf_wread();
      break;
      case LOAD_GAMESS:
         read_gamess(filename);
         update_interf_wread();
      break;
      case LOAD_HONDO:
         read_hondo(filename);
         update_interf_wread();
      break;
      case LOAD_HONDO_PUNCH:
         read_hondo_punch(filename);
         update_interf_wread();
      break;
      case LOAD_ADF:
         read_adf(filename);
         update_interf_wread();
      break;
      case LOAD_T21:
         select_t21_file(filename);
         update_interf_wread();
      break;
      case LOAD_MOS:
         read_mos(filename);
         update_interf_wread();
      break;
      case LOAD_PRDDO:
         read_prddo(filename);
         update_interf_wread();
      break;
      case LOAD_ZINDO:
         read_zindo(filename);
         update_interf_wread();
      break;
      case LOAD_MKL:
         read_mkl(filename);
         update_interf_wread();
      break;
      case LOAD_MOLDEN:
         read_molden(filename);
         update_interf_wread();
      break;
      case LOAD_MOLDEN_FREQ:
         read_molden_freq(filename);
         update_interf_wread();
      break;
      case LOAD_MOLDEN_GEOM:
         read_molden_geom(filename);
         update_interf_wread();
      break;
      case LOAD_DOTS:
         read_dots(filename);
      break;
      case LOAD_DOTVAL:
         read_dots_with_val(filename);
      break;
      case LOAD_SLD:
         read_sld(filename);
      break;
      case LOAD_GCUBE:
         read_gcubecoor(filename);
         update_interf_wread();
      break;
      case LOAD_NBO_ORB:
         read_nbo(filename);
      break;
      case LOAD_MACU:
      case ADD_MACU:
      case SUBST_MACU:
         switch (format_rbgrp->get_int_val()) {
            case 0:
               logprint(" it is a macu file");
               update_logs();
               readmacu(filename, key);
            break;
            case 1:
               logprint(" it is a gaussian cube file");
               update_logs();
               read_gcube(filename, key);
            break;
            case 2:
               logprint(" it is a adf tape41 file");
               update_logs();
               read_t41(filename, key, t41cont);
            break;
         }
/*
         if(!format_rbgrp->get_int_val()) {
            logprint(" it is a macu file");
            update_logs();
            readmacu(filename, key);
         }
         else {
            logprint(" it is a gaussian cube file");
            update_logs();
            read_gcube(filename, key);
         }
*/
      break;
      case LOAD_TEXTURE:
         read_img(filename);
      break;
      case WRITE_XYZ_ORIG:
         write_xyz_origorient(filename);
      break;
      case WRITE_XYZ_CURRENT:
         write_xyz_currentorient(filename);
      break;
      case WRITE_PDB_ORIG:
         write_pdb_origorient(filename);
      break;
      case WRITE_PDB_CURRENT:
         write_pdb_currentorient(filename);
      break;
      case WRITE_MS:
         write_ms(filename);
      break;
      case WRITE_DOTVAL:
         write_dots_with_value(filename);
      break;
      case WRITE_SLD:
         write_sld(filename);
      break;
      case WRITE_RGB:
         writergb_rle(filename, "molekel", imgx, imgy);
      break;
      case WRITE_TIFF:
         /*writetiff(filename, "molekel", imgx, imgy , "none");*/
         writetiff(filename, "molekel", imgx, imgy , 0);
      break;
      case WRITE_TIFF_LZW:
         /*writetiff(filename, "molekel", imgx, imgy ,"LZW");*/
         writetiff(filename, "molekel", imgx, imgy ,5);
      break;
      case WRITE_TIFF_PACK:
         /*writetiff(filename, "molekel", imgx, imgy ,"PackBits");*/
         writetiff(filename, "molekel", imgx, imgy ,32773);
      break;
      case WRITE_JPEG:
         writejpeg(filename, imgx, imgy, jpeg_qual);
      break;
      case RGB_CURR_IMG:
         render(0);
      break;
      case TIF_CURR_IMG:
         render(COMPRESSION_NONE);
      break;
      case TIFLZW_CURR_IMG:
         render(COMPRESSION_LZW);
      break;
      case TIFPACK_CURR_IMG:
         render(COMPRESSION_PACKBITS);
      break;
      case PPM_CURR_IMG:
         render(PPM);
      break;
      case JPEG_CURR_IMG:
         render(jpeg_qual);
      break;
      case RGB_STRC_IMG:
         generate_struct_imgs(0);
      break;
      case TIF_STRC_IMG:
         generate_struct_imgs(COMPRESSION_NONE);
      break;
      case TIFLZW_STRC_IMG:
         generate_struct_imgs(COMPRESSION_LZW);
      break;
      case TIFPACK_STRC_IMG:
         generate_struct_imgs(COMPRESSION_PACKBITS);
      break;
      case PPM_STRC_IMG:
         generate_struct_imgs(PPM);
      break;
      case JPEG_STRC_IMG:
         generate_struct_imgs(jpeg_qual);
      break;
      case RGB_FREQ_IMG:
         generate_freq_imgs(0);
      break;
      case TIF_FREQ_IMG:
         generate_freq_imgs(COMPRESSION_NONE);
      break;
      case TIFLZW_FREQ_IMG:
         generate_freq_imgs(COMPRESSION_LZW);
      break;
      case TIFPACK_FREQ_IMG:
         generate_freq_imgs(COMPRESSION_PACKBITS);
      break;
      case PPM_FREQ_IMG:
         generate_freq_imgs(PPM);
      break;
      case JPEG_FREQ_IMG:
         generate_freq_imgs(jpeg_qual);
      break;
   }  
#ifdef SUN
/* need a glutPostRedisplay on SUN for unknown reasons */
   glutPostRedisplay();
#endif
#ifdef LINUX
/* need a glutPostRedisplay on LINUX for unknown reasons */
   glutPostRedisplay();
#endif
#ifdef WIN32
/* need a glutPostRedisplay on WIN32 for unknown reasons */
   glutPostRedisplay();
#endif
}

int file_exists(char *fname)
{
  FILE *fp;

  fp = fopen(fname,"r");
  if(!fp) {
     return 0;
  }
  else {
     fclose(fp);
     return 1;
  }
}

void file_action(char *file)
{
   char  fname[200];

#ifdef WIN32
   int length;
   char bslsh[] = "\\";

   length = strlen(actual_directory);
   if(strcmp(&actual_directory[length-1], bslsh) != 0) {
      sprintf(fname, "%s%s%s", actual_directory, SEPARATOR, file);
   }
   else {
      sprintf(fname, "%s%s", actual_directory, file);
   }
#else
   sprintf(fname, "%s%s%s", actual_directory, SEPARATOR, file);
#endif

   switch(file_select_key) {
      case LOAD_XYZ:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading XYZ");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_MULTI_XYZ:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading multiple XYZ");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_PDB:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading PDB");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_MULTI_PDB:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading multiple PDB");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_GAUSS:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading GAUSSIAN LOG");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_GAMESS:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading GAMESS US LOG");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_HONDO:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading HONDO LOG");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_HONDO_PUNCH:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading HONDO PUNCH");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_ADF:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading ADF LOG");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_T21:
         strcpy(filename, fname);
         logprint("");
         logprint("Selecting T21");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_MOS:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading MOS");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_PRDDO:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading PRDDO");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_ZINDO:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading ZINDO");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_MKL:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading MKL");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_MOLDEN:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading MOLDEN FILE");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_MOLDEN_FREQ:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading MOLDEN FREQ FILE");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_MOLDEN_GEOM:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading MOLDEN GEOM FILE");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_NBO_ORB:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading nbo orbitals");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_DOTS:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading MS");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_DOTVAL:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading dots with values");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_SLD:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading SLD");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_GCUBE:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading GAUSSIAN CUBE");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case LOAD_MACU:
      case ADD_MACU:
      case SUBST_MACU:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading density file");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case EL_DENS:
      case CALC_ORB:
      case SPIN_DENS:
      case MEP:
         strcpy(filename, fname);
         if (file_exists(fname)) select_overwrite_file(fname);
         else select_job_priority();
      break;
      case LOAD_TEXTURE:
         strcpy(filename, fname);
         logprint("");
         logprint("Loading texture");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case  WRITE_XYZ_ORIG:
         strcpy(filename, fname);
         if (file_exists(fname)) select_overwrite_file(fname);
         else {
            logprint("");
            logprint("Saving XYZ (orig. orient.)");
            update_logs();
            glutTimerFunc(1, timer_cb, file_select_key);
         }
      break;
      case  WRITE_XYZ_CURRENT:
         strcpy(filename, fname);
         if (file_exists(fname)) select_overwrite_file(fname);
         else {
            logprint("");
            logprint("Saving XYZ (current orient.)");
            update_logs();
            glutTimerFunc(1, timer_cb, file_select_key);
         }
      break;
      case  WRITE_PDB_ORIG:
         strcpy(filename, fname);
         if (file_exists(fname)) select_overwrite_file(fname);
         else {
            logprint("");
            logprint("Saving PDB (orig. orient.)");
            update_logs();
            glutTimerFunc(1, timer_cb, file_select_key);
         }
      break;
      case  WRITE_PDB_CURRENT:
         strcpy(filename, fname);
         if (file_exists(fname)) select_overwrite_file(fname);
         else {
            logprint("");
            logprint("Saving PDB (current orient.)");
            update_logs();
            glutTimerFunc(1, timer_cb, file_select_key);
         }
      break;
      case  WRITE_MS:
         strcpy(filename, fname);
         if (file_exists(fname)) select_overwrite_file(fname);
         else {
            logprint("");
            logprint("Saving MS");
            update_logs();
            glutTimerFunc(1, timer_cb, file_select_key);
         }
      break;
      case  WRITE_DOTVAL:
         strcpy(filename, fname);
         if (file_exists(fname)) select_overwrite_file(fname);
         else {
            logprint("");
            logprint("Saving dots with value");
            update_logs();
            glutTimerFunc(1, timer_cb, file_select_key);
         }
      break;
      case WRITE_SLD:
         strcpy(filename, fname);
         if (file_exists(fname)) select_overwrite_file(fname);
         else {
            logprint("");
            logprint("Saving SLD");
            update_logs();
            glutTimerFunc(1, timer_cb, file_select_key);
         }
      break;
      case WRITE_RGB:
         strcpy(filename, fname);
         if (file_exists(fname)) select_overwrite_file(fname);
         else {
            logprint("");
            logprint("Saving RGB");
            update_logs();
            glutTimerFunc(1, timer_cb, file_select_key);
         }
      break;
      case WRITE_TIFF:
         strcpy(filename, fname);
         if (file_exists(fname)) select_overwrite_file(fname);
         else {
            logprint("");
            logprint("Saving TIFF");
            update_logs();
            glutTimerFunc(1, timer_cb, file_select_key);
         }
      break;
      case WRITE_TIFF_LZW:
         strcpy(filename, fname);
         if (file_exists(fname)) select_overwrite_file(fname);
         else {
            logprint("");
            logprint("Saving TIFF LZW");
            update_logs();
            glutTimerFunc(1, timer_cb, file_select_key);
         }
      break;
      case WRITE_TIFF_PACK:
         strcpy(filename, fname);
         if (file_exists(fname)) select_overwrite_file(fname);
         else {
            logprint("");
            logprint("Saving TIFF packbits");
            update_logs();
            glutTimerFunc(1, timer_cb, file_select_key);
         }
      break;
      case WRITE_JPEG:
         strcpy(filename, fname);
         if (file_exists(fname)) select_overwrite_file(fname);
         else {
            logprint("");
            logprint("Saving JPEG");
            update_logs();
            glutTimerFunc(1, timer_cb, file_select_key);
         }
      break;
      case TIF_CURR_IMG:
      case TIFLZW_CURR_IMG:
      case TIFPACK_CURR_IMG:
      case PPM_CURR_IMG:
      case RGB_CURR_IMG:
      case JPEG_CURR_IMG:
         strcpy(filename, fname);
         if (file_exists(fname)) select_overwrite_file(fname);
         else {
            logprint("");
            logprint("Rendering image");
            update_logs();
            glutTimerFunc(1, timer_cb, file_select_key);
         }
      break;
      case TIF_STRC_IMG:
      case TIFLZW_STRC_IMG:
      case TIFPACK_STRC_IMG:
      case RGB_STRC_IMG:
      case PPM_STRC_IMG:
      case JPEG_STRC_IMG:
         strcpy(filename, fname);
         logprint("");
         logprint("Rendering structure images");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      case TIF_FREQ_IMG:
      case TIFLZW_FREQ_IMG:
      case TIFPACK_FREQ_IMG:
      case RGB_FREQ_IMG:
      case PPM_FREQ_IMG:
      case JPEG_FREQ_IMG:
         strcpy(filename, fname);
         logprint("");
         logprint("Rendering freq images");
         update_logs();
         glutTimerFunc(1, timer_cb, file_select_key);
      break;
      default:
         printf("unknown case\n");
   }
   file_select_key = 0;
}

void	controltltop(muiObject *obj, enum muiReturnValue value)
{
    float sliderval;

    if ((value != MUI_SLIDER_RETURN) && (value != MUI_SLIDER_THUMB)) return;
    sliderval = muiGetVSVal(obj);
    muiSetTLTop(tl, sliderval);
}

void	handlefileselection(muiObject *obj, enum muiReturnValue value)
{
    char *fname;
    int len;

    if (value == MUI_TEXTLIST_RETURN_CONFIRM) {
	selectedfile = muiGetTLSelectedItem(obj);
	fname = filelist[selectedfile];
        if(!fname) return;
	len = strlen(fname);
	if (fname[len-1] == *SEPARATOR) {
	    fname[len-1] = 0;
	    cd(fname);
            muiClearTBString(tt);
	    return;
	} else {
            if(fname) muiSetTBString(tt, fname);
//            glutHideWindow();
            glutSetWindow(mainwin);
            if(filewin) glutDestroyWindow(filewin);
            filewin = NULL;
	    file_action(fname);
	}
    }
    if (value == MUI_TEXTLIST_RETURN) {
       selectedfile = muiGetTLSelectedItem(obj);
       fname = filelist[selectedfile];
       if(fname) muiSetTBString(tt, fname);
    }
    else return;
    selectedfile = muiGetTLSelectedItem(obj);
    muiSetVSValue(vs, 1.0);
    obj = 0;	/* for lint's sake */
}

void handleaccept(muiObject *obj, enum muiReturnValue value)
{
    char *fname;
    int len;

    if (value != MUI_BUTTON_PRESS) return;
    if (selectedfile == -1) return;
//    fname = filelist[selectedfile];
    fname = muiGetTBString(tt);
    len = strlen(fname);
    if (fname[len-1] == *SEPARATOR) {
	fname[len-1] = 0;
	cd(fname);
        muiClearTBString(tt);
	return;
    } else {
//        glutHideWindow();
        glutSetWindow(mainwin);
        if(filewin) glutDestroyWindow(filewin);
        filewin = NULL;
	file_action(fname);
    }
    obj = 0;	/* for lint's sake */
}

void handleoriginal(muiObject *obj, enum muiReturnValue value)
{
    if (value != MUI_BUTTON_PRESS) return;
    cd(starting_directory);
    obj = 0;	/* for lint's sake */
}

void handleupdir(muiObject *obj, enum muiReturnValue value)
{
    if (value != MUI_BUTTON_PRESS) return;
    cd("..");
    obj = 0;	/* for lint's sake */
}

void handlefilter(muiObject *obj, enum muiReturnValue value)
{
    if (value != MUI_BUTTON_PRESS) return;
    cd(".");
    obj = 0;    /* for lint's sake */
}

void handlecancel(muiObject *obj, enum muiReturnValue value)
{
    if (value != MUI_BUTTON_PRESS) return;
    action_key = 0;
//    glutHideWindow();
    glutSetWindow(mainwin);
    if(filewin) glutDestroyWindow(filewin);
    filewin = NULL;
    mainglui_top->enable();
    if(pixels) {
       free(pixels);
       pixels = NULL;
    }
    obj = 0;	/* for lint's sake */
}

void handletextbox(muiObject *obj, enum muiReturnValue value)
{
    char *s, *slash;

    if (value != MUI_TEXTBOX_RETURN) return;
    s = muiGetTBString(obj);
    if (0 == chdir(s)) {
	pwd();
	ls();
	muiSetTLStrings(tl, filelist);
	selectedfile = 0;
	muiChangeLabel(l4, adjstr(actual_directory));
	muiClearTBString(obj);
	return;
    }
    /* hack up the path, if any */
    slash = strrchr(s, *SEPARATOR);
    if (slash == 0) {
	slash = s-1;	/* to make filename == slash+1 */
    } else {
	if (*s == *SEPARATOR) { /* absolute path */
	    strncpy(actual_directory, s, slash-s);
	    actual_directory[slash-s] = 0;
	} else {
	    strcat(actual_directory, SEPARATOR);
	    strncat(actual_directory, s, slash-s);
	}
    }
     /* now filename == slash+1 */
//    glutHideWindow();
    glutSetWindow(mainwin);
    if(filewin) glutDestroyWindow(filewin);
    filewin = NULL;
    file_action(slash+1);
}

#define THUMBHEIGHT 20
#define ARROWSPACE 5

void makebrowsergui(void)
{
    muiObject *l1, *l2, *l3, *b1, *b2, *b3, *b4;
    int xmin, ymin, xmax, ymax;
    int i;

    muiNewUIList(1);
    l1 = muiNewBoldLabel(10, 475, "Directory:");
    muiAddToUIList(1, l1);
#ifndef WIN32
    l4 = muiNewLabel(80, 475, "./");
#else
    l4 = muiNewLabel(80, 475, ".\\");
#endif
    muiAddToUIList(1, l4);
    l2 = muiNewBoldLabel(10, 430, "Set directory:");
    muiAddToUIList(1, l2);
    b1 = muiNewButton(10, 100, 390, 415);
    muiLoadButton(b1, "Up");
    muiAddToUIList(1, b1);
    muiSetCallback(b1, handleupdir);
    b2 = muiNewButton(10, 100, 355, 380);
    muiLoadButton(b2, "Original");
    muiAddToUIList(1, b2);
    muiSetCallback(b2, handleoriginal);

    trb1 = muiNewTinyRadioButton(10, 320);
    muiAddToUIList(1, trb1);
    muiLoadButton(trb1, "Filter:");
    muiSetActive(trb1, 1);
    muiSetCallback(trb1, handlefilter);
    for(i = 0; i < 10; i++) {
       pl[i] = muiNewLabel(30, 305 - 15*i, "");
       muiAddToUIList(1, pl[i]);
    }
    tl = muiNewTextList(120, 80, 370, MINFILELIST);
    muiAddToUIList(1, tl);
    muiGetObjectSize(tl, &xmin, &ymin, &xmax, &ymax);
    vs = muiNewVSlider(xmax, ymin+2, ymax, 0, THUMBHEIGHT);
    muiSetVSValue(vs, 1.0);
    muiSetVSArrowDelta(vs, ARROWSPACE);
    muiAddToUIList(1, vs);
    tt = muiNewTextbox(120, 390, 40);
    muiSetActive(tt, 1);
    muiAddToUIList(1, tt);
    muiSetCallback(tt, handletextbox);
    l3 = muiNewBoldLabel(10, 50, "Open File/Path:");
    muiAddToUIList(1, l3);
    b3 = muiNewButton(130, 230, 9, 34);
    muiLoadButton(b3, "Accept");
    muiSetCallback(b3, handleaccept);
    muiAddToUIList(1, b3);
    b4 = muiNewButton(250, 350, 9, 34);
    muiLoadButton(b4, "Cancel");
    muiSetCallback(b4, handlecancel);
    muiAddToUIList(1, b4);
    muiSetCallback(vs, controltltop);
    muiSetCallback(tl, handlefileselection);
    
    cd(actual_directory);
}

void file_select(int key, char **pattern)
{
    int i;

    makefilewin();
    glutSetWindow(filewin);
    file_select_key = key;
    file_select_pattern = pattern;
    muiSetTBString(tt, pattern[0]);
    muiSetActiveUIList(1);
    for(i = 0; i < 10; i++) {
       if(pattern){
          if(pattern[i]) muiChangeLabel(pl[i], pattern[i]);
          else {
             pattern = NULL;
             muiChangeLabel(pl[i], "");
          }
       }
       else muiChangeLabel(pl[i], "");
    }
    ls();
    muiSetTLStrings(tl, filelist);
    glutShowWindow();
    glutPopWindow();
}

void errormsg(char *s)
{
    showinfobox(s);
}

void prname(void)
{
	actual_directory[0] = *SEPARATOR;
	if (off == 0)
		off = 1;
	actual_directory[off] = 0;
}

int dirlevels(char *s)
{
    int levels;

    for (levels = 0; *s; s++)
	if (*s == *SEPARATOR)
	    levels++;
    return(levels);
}

int cat(void)
{
	register int i, j;
	char *name = actual_directory + 1;	/* I love C */

	i = -1;
	while (dir->d_name[++i] != 0) 
	if ((off+i+2) > MAXNAMLEN - 1) {
		prname();
		return 1;
	}
	for(j=off+1; j>=0; --j)
		name[j+i+1] = name[j];
	off=i+off+1;
	name[i] = *SEPARATOR;
	for(--i; i>=0; --i)
		name[i] = dir->d_name[i];
	return 0;
}

#ifdef WIN32
void pwd(void)
{
   getcwd(actual_directory, 128);
}

#else
/* Standard modern version for macOS/Linux */
#include <unistd.h>
#include <string.h>

void pwd(void)
{
    if (getcwd(actual_directory, 1024) == NULL) {
        // Fallback or error message
        strncpy(actual_directory, ".", 128); 
    }
}
#endif
/* get the current working directory (the following 3 routines are from pwd.c) 
void pwd(void)
{

	for(off = 0;;) {
		if(stat(dot, &d) < 0) {
			showinfobox("pwd: cannot stat .!\n");
			exit(2);
		}
		if ((file = opendir(dotdot)) == NULL) {
			showinfobox("pwd: cannot open ..\n");
			exit(2);
		}
#ifdef NOFSTAT
                if(stat(dotdot, &dd) < 0) {
#else
		if(fstat(file->dd_fd, &dd) < 0) {
#endif
			showinfobox("pwd: cannot stat ..!\n");
			exit(2);
		}
		if(chdir(dotdot) < 0) {
			showinfobox("pwd: cannot chdir to ..\n");
			exit(2);
		}
		if(d.st_dev == dd.st_dev) {
			if(d.st_ino == dd.st_ino) {
				prname();
				chdir(actual_directory);
				return;
			}
			do
				if ((dir = readdir(file)) == NULL) {
					showinfobox("pwd: read error in ..\n");
					exit(2);
				}
			while (dir->d_ino != d.st_ino);
		}
		else do {
				if((dir = readdir(file)) == NULL) {
					showinfobox("pwd: read error in ..\n");
					exit(2);
				}
				stat(dir->d_name, &dd);
		} while(dd.st_ino != d.st_ino || dd.st_dev != d.st_dev);
		(void)closedir(file);
		if (cat()) {
			chdir(actual_directory);
			return;
		}
	}
}
#endif */

void freels(void)
{
    char **p;

    if(!filelist) return;
    p = filelist;
    while (*p != 0) {
	free(*p);
	*p = 0;
	p++;
    }
    free(filelist);
    filelist = NULL;
}

int mystrcmp(char **s1, char **s2)
{
    return strcmp(*s1,*s2);
}

void ls(void)
{
    DIR			*dirp;
    int			i = 0, k = 0, j = 0;
    int			len;
    struct dirent	*dir;
    struct stat		statbuf;
    char *cp;
    

    if ((dirp = opendir(actual_directory)) == NULL) {
	strcpy(actual_directory, starting_directory);
        if ((dirp = opendir(actual_directory)) == NULL) {
           logprint("bad directory!");
           update_logs();
	   return;
        }
    }
    freels();
    if((filelist = (char **)malloc(MINFILELIST*sizeof(char *))) == NULL){
       showinfobox("Can't allocate string-pointer");
       return;
    }
    for(j=0; j<MINFILELIST; j++) {
       filelist[j] = NULL;
    }
    chdir(actual_directory);
    while ((dir = readdir(dirp)) != NULL) {
	if (dir->d_name[0] == '.')
	    continue;
	stat(dir->d_name,&statbuf);
	len = strlen(dir->d_name) + 1 + (statbuf.st_mode & S_IFDIR? 1 : 0);
        if(muiGetActive(trb1) && file_select_pattern) {
           if(statbuf.st_mode & S_IFDIR) {
              filelist[i] = (char *)malloc(len);
              strcpy(filelist[i], dir->d_name);
              if (statbuf.st_mode & S_IFDIR) {
              filelist[i][len-2] = *SEPARATOR; filelist[i][len-1] = 0;
              }
              i++;
              if((filelist = (char **)realloc(filelist, (i+MINFILELIST)*sizeof(char *))) == NULL){
                 showinfobox("Can't allocate string-pointer");
                 return;
              }
           }
           for(k = 0; file_select_pattern[k]; k++) {
              cp = dir->d_name + (strlen(dir->d_name) - strlen(file_select_pattern[k]));
              if(!strcmp(cp, file_select_pattern[k])) {
                 filelist[i] = (char *)malloc(len);
                 strcpy(filelist[i], dir->d_name);
	         i++;
                 if((filelist = (char **)realloc(filelist, (i+MINFILELIST)*sizeof(char *))) == NULL){
                    showinfobox("Can't allocate string-pointer");
                    return;
                 }
              }
           }
        }
        else {
	   filelist[i] = (char *)malloc(len);
	   strcpy(filelist[i], dir->d_name);
	   if (statbuf.st_mode & S_IFDIR) {
	       filelist[i][len-2] = *SEPARATOR; filelist[i][len-1] = 0;
	   }
	   i++;
           if((filelist = (char **)realloc(filelist, (i+MINFILELIST)*sizeof(char *))) == NULL){
              showinfobox("Can't allocate string-pointer");
              return;
           }
        }
    }
    filelist[i] = 0;
    qsort(&filelist[0], i, sizeof (char *), (int (*)(const void *, const void *))mystrcmp);
    closedir(dirp);
}

int cd(char *s)
{
    char str[100];

    if(chdir(s) < 0) {
	sprintf(str,"cannot open:\n%s !",s);
        showinfobox(str);
	return -1;
    }

    pwd();
    ls();
    muiSetTLStrings(tl, filelist);
    muiChangeLabel(l4, adjstr(actual_directory));
    selectedfile = 0;
    return 0;
}

void makefilewin(void)
{
//    glutInitWindowPosition(glutGet(GLUT_SCREEN_WIDTH)/2 - 200,\
//             glutGet(GLUT_SCREEN_HEIGHT)/2 - 250);
//    glutInitWindowPosition(glutGet(GLUT_WINDOW_X) + xsize/2 - 50, \
//        glutGet(GLUT_WINDOW_Y) + ysize/2 - 120);
    if(filewin) {
       glutSetWindow(filewin);
       glutShowWindow();
       glutPopWindow();
    }
    else {
       glutInitWindowPosition(winposition[0], winposition[1]);
       glutInitWindowSize(400, 500);
       glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
       filewin = glutCreateWindow("file browser");
//       makebrowsergui();
       muiInit();
       muiAttachUIList(1);
//       glutHideWindow();
    }
}

/* --------------------------------------------- */
/* orbital selection */

void handle_orb_accept(muiObject *obj, enum muiReturnValue value)
{
    register short i;
    char str[100];

    if (value != MUI_BUTTON_PRESS) return;
    else {
      if(molptr->alphaBeta){
         sprintf(str, "you selected orbital %d ", o_choice/2 + molptr->firstOrbital);
         if(o_choice%2) {
            strcat(str, "beta\n");
            molOrb = molptr->betaOrbital + o_choice/2;
         }
         else {
            strcat(str, "alpha\n");
            molOrb = molptr->alphaOrbital + o_choice/2;
         }
      }
      else {
         sprintf(str, "you selected orbital %d\n", o_choice + molptr->firstOrbital);
         molOrb = molptr->alphaOrbital + o_choice;
      }
//      glutHideWindow();
      glutSetWindow(mainwin);
      if(orbwin) glutDestroyWindow(orbwin);
      orbwin = NULL;
      file_select(CALC_ORB, macu_ext);

      for(i=0; i <= o_nlines; i++) free(orblist[i]);
      free(orblist);

      logprint("");
      logprint(str);
      update_logs();

    }
    obj = 0;	/* for lint's sake */
}

void handle_orb_cancel(muiObject *obj, enum muiReturnValue value)
{
    register short i;

    if (value != MUI_BUTTON_PRESS) return;
    action_key = 0;
//    glutHideWindow();
    glutSetWindow(mainwin);
    if(orbwin) glutDestroyWindow(orbwin);
    orbwin = NULL;
    mainglui_top->enable();

    for(i=0; i <= o_nlines; i++) free(orblist[i]);
    free(orblist);

    obj = 0;	/* for lint's sake */
}

void handle_orb_selection(muiObject *obj, enum muiReturnValue value)
{
   register short i;
   char str[100];

   if (value == MUI_TEXTLIST_RETURN_CONFIRM) {
      o_choice = muiGetTLSelectedItem(obj);
      if(molptr->alphaBeta){
         if(o_choice < 0 || o_choice > 2*actualmol->nMolecularOrbitals-1) o_choice = 0;
         sprintf(str, "you selected orbital %d ", o_choice/2 + molptr->firstOrbital);
         if(o_choice%2) {
            strcat(str, "beta\n");
            molOrb = molptr->betaOrbital + o_choice/2;
         }
         else {
            strcat(str, "alpha\n");
            molOrb = molptr->alphaOrbital + o_choice/2;
         }
      }
      else {
         if(o_choice < 0 || o_choice > actualmol->nMolecularOrbitals) o_choice = 0;
         sprintf(str, "you selected orbital %d\n", o_choice + molptr->firstOrbital);
         molOrb = molptr->alphaOrbital + o_choice;
      }
//      glutHideWindow();
      glutSetWindow(mainwin);
      if(orbwin) glutDestroyWindow(orbwin);
      orbwin = NULL;
      file_select(CALC_ORB, macu_ext);

      for(i=0; i <= o_nlines; i++) free(orblist[i]);
      free(orblist);

      logprint("");
      logprint(str);
      update_logs();

      return;
   }
   if (value == MUI_TEXTLIST_RETURN) {
      o_choice = muiGetTLSelectedItem(obj);
      if(molptr->alphaBeta){
         if(o_choice < 0 || o_choice > 2*actualmol->nMolecularOrbitals-1) o_choice = 0;
      }
      else {
         if(o_choice < 0 || o_choice > actualmol->nMolecularOrbitals) o_choice = 0;
      }
      return;
   }
   else return;
}

void	control_tl1_top(muiObject *obj, enum muiReturnValue value)
{
    float sliderval;

    if ((value != MUI_SLIDER_RETURN) && (value != MUI_SLIDER_THUMB)) return;
    sliderval = muiGetVSVal(obj);
    muiSetTLTop(tl1, sliderval);
}

void makeorbgui(void)
{
    muiObject *l1, *l2, *b1, *b2;
    int xmin, ymin, xmax, ymax;

    muiNewUIList(2);
    l1 = muiNewBoldLabel(10, 475, "Select Orbital:");
    muiAddToUIList(2, l1);
    l2 = muiNewLabel(140, 475, "MO       Eigenvalue      Occupation");
    muiAddToUIList(2, l2);
    tl1 = muiNewTextList(120, 80, 370, 22);
    muiAddToUIList(2, tl1);
    muiGetObjectSize(tl1, &xmin, &ymin, &xmax, &ymax);
    vs_orb = muiNewVSlider(xmax, ymin+2, ymax, 0, THUMBHEIGHT);
    muiSetVSValue(vs_orb, 1.0);
    muiSetVSArrowDelta(vs_orb, ARROWSPACE);
    muiAddToUIList(2, vs_orb);
    b1 = muiNewButton(250, 350, 9, 34);
    muiLoadButton(b1, "Cancel");
    muiSetCallback(b1, handle_orb_cancel);
    muiAddToUIList(2, b1);
    b2 = muiNewButton(130, 230, 9, 34);
    muiLoadButton(b2, "Accept");
    muiSetCallback(b2, handle_orb_accept);
    muiAddToUIList(2, b2);
    muiSetCallback(vs_orb, control_tl1_top);
    muiSetCallback(tl1, handle_orb_selection);
}

void makeorbwin(void)
{
    if(orbwin) {
       glutSetWindow(orbwin);
       glutShowWindow();
       glutPopWindow();
    }
    else {
       glutInitWindowPosition(winposition[0], winposition[1]);
       glutInitWindowSize(400, 500);
       glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
       orbwin = glutCreateWindow("select orbital");
//       makeorbgui();
       muiInit();
       muiAttachUIList(2);
//       glutHideWindow();
    }
}

void select_mo(Mol *mp)
{
   register short i, j;
   char str[1000];

   makeorbwin();

   molptr = mp;
   o_nlines = mp->alphaBeta ? 2*mp->nMolecularOrbitals: mp->nMolecularOrbitals;

   if((orblist = (char **)malloc((o_nlines+1)*sizeof(char *))) == NULL){
      showinfobox("Can't allocate string-pointer");
      return;
   }

   if((orblist[o_nlines] =strdup(str)) == NULL){
      showinfobox("Can't allocate selector-string");
      return;
   } 
//   free(*orblist);
   *orblist[o_nlines] = 0;
   orblist[o_nlines] = 0;

   if(mp->alphaBeta){
      for(i=0, j=0; i<mp->nMolecularOrbitals; i++){
         sprintf(str, "%3d alpha   %10.5f          %5.3f  %s", i+mp->firstOrbital,
                 mp->alphaOrbital[i].eigenvalue,
                 mp->alphaOrbital[i].occ,
                 mp->alphaOrbital[i].type);
         if((orblist[j++] = strdup(str)) == NULL){
            showinfobox("Can't allocate selector-string");
            return;
         }
         sprintf(str, "    beta    %10.5f          %5.3f  %s",
                 mp->betaOrbital[i].eigenvalue,
                 mp->betaOrbital[i].occ,
                 mp->betaOrbital[i].type);
         if((orblist[j++] = strdup(str)) == NULL){
            showinfobox("Can't allocate selector-string");
            return;
         }
      }
   }
   else {
      for(i=0; i<mp->nMolecularOrbitals; i++){
         sprintf(str, "    %3d   %10.5f          %5.3f  %s", i+mp->firstOrbital,
               mp->alphaOrbital[i].eigenvalue,
                 mp->alphaOrbital[i].occ,
                 mp->alphaOrbital[i].type);
         if((orblist[i] = strdup(str)) == NULL){
            showinfobox("Can't allocate selector-string");
            return;
         }
      }
   }
   glutSetWindow(orbwin);
   muiSetTLStrings(tl1, orblist);
   glutShowWindow();
   glutPopWindow();
   mainglui_top->disable();
}


/* --------------------------------------------- */
/* frequency selection */

void handle_freq_accept(muiObject *obj, enum muiReturnValue value)
{
    register short i;
    char str[100];

    if (value != MUI_BUTTON_PRESS) return;
    else {
      glutSetWindow(mainwin);
      if(freqwin) glutDestroyWindow(freqwin);
      freqwin = NULL;
      actualvib = actualmol->freq_arrow = actualmol->vibration + f_choice;
      if(actualvib) {
         sprintf(str, "%5.0f/cm %s", actualvib->frequency, actualvib->type); 
      }
      else {
         actualvib = actualmol->freq_arrow = NULL;
         sprintf(str, "no frequency selected");
      }
      freqtext->set_text(str);
      sprintf(str, "freq %3d %5.0f/cm selected", f_choice+1, actualvib->frequency); 

      for(i=0; i <= f_nlines; i++) free(freqlist[i]);
//      free(freqlist);

      logprint("");
      logprint(str);
      update_logs();

    }
    obj = 0;	/* for lint's sake */
}

void handle_freq_cancel(muiObject *obj, enum muiReturnValue value)
{
    register short i;
    char str[100];

    if (value != MUI_BUTTON_PRESS) return;

    action_key = 0;
    actualvib = actualmol->freq_arrow = NULL;

    sprintf(str, "no frequency selected");
    freqtext->set_text(str);
//    glutHideWindow();
    glutSetWindow(mainwin);
    if(freqwin) glutDestroyWindow(freqwin);
    freqwin = NULL;
    for(i=0; i <= f_nlines; i++) free(freqlist[i]);
//   free(freqlist);

    obj = 0;	/* for lint's sake */
}

void handle_freq_selection(muiObject *obj, enum muiReturnValue value)
{
   register short i;
   char str[100];

   if (value == MUI_TEXTLIST_RETURN_CONFIRM) {
      f_choice = muiGetTLSelectedItem(obj);
      if(f_choice < 0 || f_choice > actualmol->n_frequencies) f_choice = 0;
      glutSetWindow(mainwin);
      if(freqwin) glutDestroyWindow(freqwin);
      freqwin = NULL;
      actualvib = actualmol->freq_arrow = actualmol->vibration + f_choice;
      if(actualvib) {
         sprintf(str, "%5.0f/cm %s", actualvib->frequency, actualvib->type); 
      }
      else {
         sprintf(str, "no frequency selected");
         actualvib = actualmol->freq_arrow = NULL;
      }
      freqtext->set_text(str);
      sprintf(str, "freq %3d %5.0f/cm selected", f_choice+1, actualvib->frequency); 

      for(i=0; i <= f_nlines; i++) free(freqlist[i]);
      free(freqlist);

      logprint("");
      logprint(str);
      update_logs();

      return;
   }
   if (value == MUI_TEXTLIST_RETURN) {
      f_choice = muiGetTLSelectedItem(obj);
      if(f_choice < 0 || f_choice > actualmol->n_frequencies) f_choice = 0;
      return;
   }
   else return;

}

void control_tlf1_top(muiObject *obj, enum muiReturnValue value)
{
    float sliderval;

    if ((value != MUI_SLIDER_RETURN) && (value != MUI_SLIDER_THUMB)) return;
    sliderval = muiGetVSVal(obj);
    muiSetTLTop(tlf1, sliderval);
}


void makefreqgui(void)
{
    muiObject *l1, *l2, *b1, *b2;
    int xmin, ymin, xmax, ymax;

    muiNewUIList(3);
    l1 = muiNewBoldLabel(10, 475, "Select Frequency:");
    muiAddToUIList(3, l1);
    l2 = muiNewLabel(147, 475, "No        1/cm        Symmetry");
    muiAddToUIList(3, l2);
    tlf1 = muiNewTextList(120, 80, 370, 22);
    muiAddToUIList(3, tlf1);
    muiGetObjectSize(tlf1, &xmin, &ymin, &xmax, &ymax);
    vs_freq = muiNewVSlider(xmax, ymin+2, ymax, 0, THUMBHEIGHT);
    muiSetVSValue(vs_freq, 1.0);
    muiSetVSArrowDelta(vs_freq, ARROWSPACE);
    muiAddToUIList(3, vs_freq);
    b1 = muiNewButton(250, 350, 9, 34);
    muiLoadButton(b1, "Cancel");
    muiSetCallback(b1, handle_freq_cancel);
    muiAddToUIList(3, b1);
    b2 = muiNewButton(130, 230, 9, 34);
    muiLoadButton(b2, "Accept");
    muiSetCallback(b2, handle_freq_accept);
    muiAddToUIList(3, b2);
    muiSetCallback(vs_freq, control_tlf1_top);
    muiSetCallback(tlf1, handle_freq_selection);
}

void makefreqwin(void)
{
    if(freqwin) {
       glutSetWindow(freqwin);
       glutShowWindow();
       glutPopWindow();
    }
    else {    
       glutInitWindowPosition(winposition[0], winposition[1]);
       glutInitWindowSize(400, 500);
       glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
       freqwin = glutCreateWindow("select frequency");
//       makefreqgui();
       muiInit();
       muiAttachUIList(3);
//       glutHideWindow();
    }
}


void select_vib(void)
{
   register short i;
   char str[1000];

   makefreqwin();

   f_nlines = actualmol->n_frequencies;
   if((freqlist = (char **)malloc((f_nlines+1)*sizeof(char *))) == NULL){
      showinfobox("Can't allocate string-pointer");
      return;
   }

   if((freqlist[f_nlines] =strdup(str)) == NULL){
      showinfobox("Can't allocate selector-string");
      return;
   } 
//   free(*freqlist);
   *freqlist[f_nlines] = 0;
   freqlist[f_nlines] = 0;

   for(i=0; i<actualmol->n_frequencies; i++){
      sprintf(str, "    %3d   %12.2f   %s", i+1, actualmol->vibration[i].frequency,
                                                 actualmol->vibration[i].type);
      if((freqlist[i] = strdup(str)) == NULL){
         showinfobox("Can't allocate selector-string");
         return;
      }
   }

   glutSetWindow(freqwin);
   muiSetTLStrings(tlf1, freqlist);
   glutShowWindow();
   glutPopWindow();
}
