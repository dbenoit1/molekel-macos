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
#include "molekel.h"
#include "constant.h"
#include "glutwin.h"
#include "general.h"
#include "menu.h"
#include "action.h"
#include "browser.h"
#include "maininterf.h"
#include "chooseinterf.h"
#include "surfaceinterf.h"
#include "manip.h"
#include "pick.h"
#include "geometry.h"
#include "macu.h"
#include "connolly.h"
#include "connect.h"
#include "snap.h"
#include "colorinterf.h"
#include "textureinterf.h"
#include "license.h"
#include "render.h"

int action_key = 0;
int main_menu_id, truncC_menu_id, truncF_menu_id, truncB_menu_id;

void make_main_menu(void)
{
   int load_menuID, compute_menuID, geom_menuID, animate_menuID, write_menuID;
   int snap_menuID;

   load_menuID = glutCreateMenu(load_menu);
   glutAddMenuEntry("xyz", LOAD_XYZ);
   glutAddMenuEntry("multiple structure xyz", LOAD_MULTI_XYZ);
   glutAddMenuEntry("pdb", LOAD_PDB);
   glutAddMenuEntry("multiple structure pdb", LOAD_MULTI_PDB);
   glutAddMenuEntry("--------", DUMMY);
   glutAddMenuEntry("gaussian log", LOAD_GAUSS);
   glutAddMenuEntry("gaussian cube", LOAD_GCUBE);
   glutAddMenuEntry("nbo orb", LOAD_NBO_ORB);
   glutAddMenuEntry("--------", DUMMY);
   glutAddMenuEntry("gamess (us) log", LOAD_GAMESS);
   glutAddMenuEntry("hondo", LOAD_HONDO);
   glutAddMenuEntry("adf", LOAD_ADF);
   glutAddMenuEntry("--------", DUMMY);
   glutAddMenuEntry("mos", LOAD_MOS);
   glutAddMenuEntry("prddo", LOAD_PRDDO);
   glutAddMenuEntry("zindo", LOAD_ZINDO);
   glutAddMenuEntry("--------", DUMMY);
   glutAddMenuEntry("molden format", LOAD_MOLDEN);
   glutAddMenuEntry("molden freq format", LOAD_MOLDEN_FREQ);
   glutAddMenuEntry("molden geom format", LOAD_MOLDEN_GEOM);
   glutAddMenuEntry("--------", DUMMY);
   glutAddMenuEntry("mkl", LOAD_MKL);

   compute_menuID = glutCreateMenu(compute_menu);
   glutAddMenuEntry("Orbital",  CALC_ORB);
   glutAddMenuEntry("El. density",  EL_DENS);
   glutAddMenuEntry("Spin density",  SPIN_DENS);
   glutAddMenuEntry("Spin on atom", ATOM_SPIN);
   glutAddMenuEntry("MEP",  MEP);
   glutAddMenuEntry("--------", DUMMY);
   glutAddMenuEntry("Fast surface", FAST_SURF);
   glutAddMenuEntry("Connolly", CONNOLLY);
   glutAddMenuEntry("--------", DUMMY);
   glutAddMenuEntry("Adjust box", BOX);

   geom_menuID = glutCreateMenu(geom_menu);
   glutAddMenuEntry("Atom", ATOMDATA);
   glutAddMenuEntry("Distance", DISTANCE);
   glutAddMenuEntry("All distances", ALLDIST);
   glutAddMenuEntry("Angle", VALENCE);
   glutAddMenuEntry("Dihedral", DIHEDRAL);
   glutAddMenuEntry("Remove values", REM_VALUES);
   glutAddMenuEntry("--------", DUMMY);
   glutAddMenuEntry("Add dummy", ADD_DUMMY);
   glutAddMenuEntry("Remove dummies", REM_DUMMY);
   glutAddMenuEntry("--------", DUMMY);
   glutAddMenuEntry("Match", MATCH);

   animate_menuID = glutCreateMenu(animate_menu);
   glutAddMenuEntry("Frequency", VIBRATION);
   glutAddMenuEntry("Series of coords", PLAY);

   write_menuID = glutCreateMenu(write_menu);
   glutAddMenuEntry("xyz (orig. orient.)", WRITE_XYZ_ORIG);
   glutAddMenuEntry("xyz (current orient.)", WRITE_XYZ_CURRENT);
   glutAddMenuEntry("pdb (orig. orient.)", WRITE_PDB_ORIG);
   glutAddMenuEntry("pdb (current orient.)", WRITE_PDB_CURRENT);
   glutAddMenuEntry("--------", DUMMY);
   glutAddMenuEntry("Surface (sld)", WRITE_SLD);
   glutAddMenuEntry("Surface dot (ms)", WRITE_MS);
   glutAddMenuEntry("Surface dot w. value", WRITE_DOTVAL);

   snap_menuID = glutCreateMenu(snap_menu);
   glutAddMenuEntry("RGB", WRITE_RGB);
   glutAddMenuEntry("TIFF", WRITE_TIFF);
   glutAddMenuEntry("TIFF LZW compressed", WRITE_TIFF_LZW);
   glutAddMenuEntry("TIFF Packbits", WRITE_TIFF_PACK);
   glutAddMenuEntry("JPEG", WRITE_JPEG);

   main_menu_id = glutCreateMenu(main_menu);
   glutAddSubMenu("Load", load_menuID);
   glutAddMenuEntry("--------", DUMMY);
   glutAddSubMenu("Geometry", geom_menuID);
   glutAddMenuEntry("Dipole moment", DIPOLE);
   glutAddMenuEntry("Labels", LABELS);
   glutAddMenuEntry("--------", DUMMY);
   glutAddSubMenu("Compute", compute_menuID);
   glutAddMenuEntry("Surface", SURF_MENU);
   glutAddMenuEntry("Frequency",  FREQ_ARROW);
   glutAddSubMenu("Animate", animate_menuID);
   glutAddMenuEntry("--------", DUMMY);
   glutAddMenuEntry("Color", COLOR_MENU);
   glutAddMenuEntry("Texture", TEXTURE_MENU);
   glutAddMenuEntry("--------", DUMMY);
   glutAddSubMenu("Write", write_menuID);
   glutAddSubMenu("Snapshot", snap_menuID);
   glutAddMenuEntry("Render", RENDER);
   glutAddMenuEntry("Pop maininterface", POP_MAININTERF);
   glutAddMenuEntry("--------", DUMMY);
   glutAddMenuEntry("Done picking", DONE);
   glutAddMenuEntry("--------", DUMMY);
   glutAddMenuEntry("Warranty", WARRANTY);
   glutAddMenuEntry("License", LICENSE);
   glutAddMenuEntry("--------", DUMMY);
   glutAddMenuEntry("Quit", QUIT);
   glutAttachMenu(GLUT_RIGHT_BUTTON);

   truncF_menu_id = glutCreateMenu(truncF_menu);
   glutAddMenuEntry("Frequency",  VIBRATION);
   glutAddMenuEntry("Pop maininterface", POP_MAININTERF);

   truncC_menu_id = glutCreateMenu(truncC_menu);
   glutAddMenuEntry("Series of coords", PLAY);
   glutAddMenuEntry("Pop maininterface", POP_MAININTERF);

   truncB_menu_id = glutCreateMenu(truncC_menu);
   glutAddMenuEntry("Adjust box", BOX);
   glutAddSubMenu("Compute", compute_menuID);
   glutAddMenuEntry("Pop maininterface", POP_MAININTERF);
}

void truncF_menu(int value)
{
   switch(value) {
      case VIBRATION:
         freq_interf();
      break;
      case POP_MAININTERF:
         glutSetWindow(maininterfwin);
         glutPopWindow();
      break;
   }
}

void truncC_menu(int value)
{
   switch(value) {
      case PLAY:
         play_interf();
      break;
      case POP_MAININTERF:
         glutSetWindow(maininterfwin);
         glutPopWindow();
      break;
   }
}

void truncB_menu(int value)
{
   switch(value) {
      case BOX:
         box_interf();
      break;
      case POP_MAININTERF:
         glutSetWindow(maininterfwin);
         glutPopWindow();
      break;
   }
}


void main_menu(int value)
{
   switch(value) {
      case LABELS:
         label_interf();
      break;
      case SURF_MENU:
         if(!surfaceglui) makesurfaceglui();
         else {
            surfaceglui->show();
            glutSetWindow(surfaceglui->get_glut_window_id());
            glutPopWindow();
         }
      break;
      case FREQ_ARROW:
         freq_interf();
      break;
      case POP_MAININTERF:
         glutSetWindow(maininterfwin);
         glutPopWindow();
      break;
      case DONE:
         if(actualmol) {
            pick_action(value, 1, 1); /* dummy screen coordinates 1,1 not needed */
         }
      break;
      case TEXTURE_MENU:
         textureinterf();
      break;
      case COLOR_MENU:
         colorinterf();
      break;
      case RENDER:
         render_interf();
      break;
      case DIPOLE:
         dipole_interf();
      break;
      case WARRANTY:
         printf("%s\n", warranty_text);
      break;
      case LICENSE:
         printf("MOLEKEL, Version %s, Date: %s,\n",\
            MOLEKEL_VERSION, MOLEKEL_VERSION_DATE);
         printf("%s\n", "by Stefan Portmann, Copyright (C) 2002 CSCS/ETHZ");
         printf("%s\n", "(original IRIX GL implementation, concept and data structure");
         printf("%s\n", "by Peter F. Fluekiger, CSCS/UNI Geneva)");
         printf("%s\n", license_text);
      break;
      case QUIT:
         Quit();
      break;
   }
   glutPostRedisplay();
}

void load_menu(int value)
{
   switch(value) {
      case LOAD_XYZ:
         action_key = value;
         xyzformat_interf();
      break;
      case LOAD_MULTI_XYZ:
         file_select(value, xyz_ext);
      break;
      case LOAD_PDB:
      case LOAD_MULTI_PDB:
         file_select(value, pdb_ext);
         break;
      case LOAD_GAUSS:
      case LOAD_GAMESS:
      case LOAD_HONDO:
      case LOAD_ADF:
      case LOAD_PRDDO:
      case LOAD_ZINDO:
         file_select(value, out_ext);
         break;
      case LOAD_MOS:
         file_select(value, mos_ext);
         break;
      case LOAD_MKL:
         file_select(value, mkl_ext);
         break;
      case LOAD_MOLDEN:
      case LOAD_MOLDEN_FREQ:
      case LOAD_MOLDEN_GEOM:
         file_select(value, molden_ext);
         break;
      case LOAD_GCUBE:
         file_select(value, gcube_ext);
         break;
      case LOAD_NBO_ORB:
         file_select(value, nboorb_ext);
         break;
   }
   glutPostRedisplay();
}

void compute_menu(int value)
{
   if(!firstmol) {
      logprint("No molecule loaded");
      update_logs();
      return;
   }

   switch(value) {
      case EL_DENS:
         if(!actualmol || !actualmol->alphaOrbital) {
            logprint("Load an output-file first");
            update_logs();
            return;
         }
         action_key = value;
         select_coeff_or_matrix();
      break;
      case SPIN_DENS:
      case ATOM_SPIN:
         if(!actualmol || !actualmol->alphaOrbital) {
            logprint("Load an output-file first");
            update_logs();
            return;
         }
         if(!(actualmol->betaDensity ||
               actualmol->alphaBeta ||
               (actualmol->nAlpha != actualmol->nBeta))) {
            logprint("No spin density");
            update_logs();
            return;
         }
         action_key = value;
         select_coeff_or_matrix();
      break;
      case CALC_ORB:
         if(!actualmol || !actualmol->alphaOrbital) {
            logprint("Load an output-file first");
            update_logs();
            return;
         }
         action_key = value;
         select_mo(actualmol);
      break;
      case MEP:
         if(!actualmol || !actualmol->alphaOrbital) {
            logprint("Load an output-file first");
            update_logs();
            return;
         }
         action_key = value;
         select_grid_or_dot();
      break;
      case CONNOLLY:
         connolly_interf();
      break;
      case FAST_SURF:
         fast_surf();
      break;
      case BOX:
         box_interf();
      break;
   }
   glutPostRedisplay();
}

void geom_menu(int value)
{
   if(!firstmol) {
      logprint("No molecule loaded");
      update_logs();
      return;
   }

   switch(value) {
       case ATOMDATA:
          enter_pick(value, "Show atomic data");
       break;
       case DISTANCE:
          enter_pick(value, "Measure distances");
       break;
       case VALENCE:
          enter_pick(value, "Measure angles");
       break;
       case DIHEDRAL:
          enter_pick(value, "Measure dihedral angles");
       break;
       case ALLDIST:
          logprint("Measure all distances");
          update_logs();
          measure_alldist();
       break;
       case REM_VALUES:
          logprint("Remove all values");
          update_logs();
          remove_values();
       break;
       case ADD_DUMMY:
          enter_pick(value, "Add dummy atom");
       break;
       case REM_DUMMY:
          remove_dummies();
          logprint("Remove dummy atoms");
          update_logs();
       break;
       case MATCH:
          if(!firstmol || !firstmol->next){
             showinfobox("need two molecules to make a match!");
             return;
          }
          enter_pick(value, "Match two molecules");
       break;
   }
   update_interface_flags();
   glutPostRedisplay();
}

void label_menu(int value)
{
   if(!firstmol) {
      logprint("No molecule loaded");
      update_logs();
      return;
   }

   switch(value) {
       case LABELS:
          if(actualmol) {actualmol->labels = !actualmol->labels; bit.labels = actualmol->labels;} 
       break;
       case ATM_CHAR:
          if(actualmol) {actualmol->atm_char = !actualmol->atm_char; bit.atm_char = actualmol->atm_char;}
       break;
       case ATM_SPIN:
          if(actualmol) {actualmol->atm_spin = !actualmol->atm_spin; bit.atm_spin = actualmol->atm_spin;}
       break;
   }
   glutPostRedisplay();
}

void animate_menu(int value)
{
   if(!firstmol) {
      logprint("No molecule loaded");
      update_logs();
      return;
   }
   switch(value) {
      case VIBRATION:
         freq_interf();
      break;
      case PLAY:
         play_interf();
      break;
   }
}

void write_menu(int value)
{
   switch(value) {
      case  WRITE_XYZ_ORIG:
      case  WRITE_XYZ_CURRENT:
         action_key = value;
         file_select(value, xyz_ext);
      break;
      case  WRITE_PDB_ORIG:
      case  WRITE_PDB_CURRENT:
         action_key = value;
         file_select(value, pdb_ext);
      break;
      case  WRITE_MS:
         action_key = value;
         file_select(value, ms_ext);
      break;
      case WRITE_SLD:
         action_key = value;
         file_select(value, sld_ext);
      break;
      case  WRITE_DOTVAL:
         action_key = value;
         file_select(value, dotval_ext);
      break;
   }
}

void snap_menu(int value)
{
   switch(value) {
      case WRITE_RGB:
         getpixels(0, 0, xsize, ysize);
         imgx = xsize;
         imgy = ysize;
         action_key = value;
         file_select(value, rgb_ext);
      break;
      case WRITE_TIFF:
      case WRITE_TIFF_LZW:
      case WRITE_TIFF_PACK:
         getpixels(0, 0, xsize, ysize);
         imgx = xsize;
         imgy = ysize;
         action_key = value;
         file_select(value, tif_ext);
      break;
      case WRITE_JPEG:
         getpixels(0, 0, xsize, ysize);
         imgx = xsize;
         imgy = ysize;
         action_key = value;
         file_select(value, jpeg_ext);
      break;
   }
}
