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


#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "general.h"
#include "render.h"
#include "glutwin.h"
#include "drawing.h"
#include "snap.h"
#include "maininterf.h"
#include "browser.h"
#include "material.h"
#include <tiffio.h>

/*------------------------------------------------------------*/
TRcontext *tr;
unsigned int output_width, output_height;
int save_w, save_h;
int frame_nbr = 1;

/*------------------------------------------------------------*/
void tile_persp(TRcontext *tr)
{
   float ratio;

   ratio = (float)output_width/output_height;

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   if(bit.persp){
      trPerspective(tr, 30.0,  ratio,  0.23, 20);
/* gluLookAt is in trBeginTile */
   }
   else {
      trOrtho(tr, -ratio, ratio, -1.0, 1.0, -3.5, 17);
   }
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}

void DrawScene(TRcontext *tr)
{
   glClearColor(0.0, 0.0, 0.0, 0.0);
   glClear(GL_DEPTH_BUFFER_BIT);
   draw_back();

   if(bit.depthcue) glEnable(GL_FOG);
   else             glDisable(GL_FOG);

   drawit2();
}

void tile_render(int compression) {

   FILE *f;
   char line[100];
   GLint tile_width, tile_height, tile_border = 0;
   GLint tile_width_max = 512, tile_height_max = 512;
   int tileCount = 0, i, more;

   GLubyte *buffer;
   GLubyte *tile;

   output_width = (output_width / 4) * 4;
   output_height = (output_height / 4) * 4;
   tile_width = output_width;
   tile_height = output_height;
   scale_line(output_width, output_height);

/*
too large tiles seem not to work!
   if(glutGet(GLUT_SCREEN_WIDTH)) tile_width_max = glutGet(GLUT_SCREEN_WIDTH);
   if(glutGet(GLUT_SCREEN_HEIGHT)) tile_height_max = glutGet(GLUT_SCREEN_HEIGHT);
*/
   

   while(tile_width > tile_width_max) {
      tile_border = 16;
      tile_width = tile_width / 2;
   }
   while(tile_height > tile_height_max) {
      tile_border = 16;
      tile_height = tile_height / 2;
   }

   bit.render = 1;

   if(compression == PPM) {
// -----
// write direct to ppm file
// taken from trdemo2 of the tr package
//
      tile =  (GLubyte *)malloc(tile_width * tile_height * 3 * sizeof(GLubyte));
      if (!tile) {
         showinfobox("Malloc of tile buffer failed!");
         return;
      }

      /* allocate buffer to hold a row of tiles */
      buffer = (GLubyte *)malloc(output_width * tile_height * 3 * sizeof(GLubyte));
      if (!buffer) {
         free(tile);
         showinfobox("Malloc of tile row buffer failed!");
         return;
      }

      /* Setup.  Each tile is TILE_WIDTH x TILE_HEIGHT pixels. */
      tr = trNew();
      trTileSize(tr, tile_width, tile_height, tile_border);
      trTileBuffer(tr, GL_RGB, GL_UNSIGNED_BYTE, tile);
      trImageSize(tr, output_width, output_height);
      trRowOrder(tr, TR_TOP_TO_BOTTOM);
      tile_persp(tr);


      /* Prepare ppm output file */
      f = fopen(filename, "w");
      if (!f) {
         sprintf(line, "Can't open %s", filename);
         showinfobox(line);
         free(tile);
         free(buffer);
         return;
      }
      fprintf(f,"P6\n");
      fprintf(f,"# ppm-file created by %s\n", "MOLEKEL");
      fprintf(f,"%i %i\n", output_width, output_height);
      fprintf(f,"255\n");
      fclose(f);
      f = fopen(filename, "ab");  /* now append binary data */
      if (!f) {
         sprintf(line, "Couldn't append to image file: %s", filename);
         showinfobox(line);
         free(tile);
         free(buffer);
         return;
      }

      /* just to be safe... */
      glPixelStorei(GL_PACK_ALIGNMENT, 1);

      /* Draw tiles */
      more = 1;
         while (more) {
         int curColumn;
         trBeginTile(tr);
         curColumn = trGet(tr, TR_CURRENT_COLUMN);
         DrawScene(tr);      /* draw our stuff here */
         more = trEndTile(tr);

         /* save tile into tile row buffer*/
         {
            int curTileWidth = trGet(tr, TR_CURRENT_TILE_WIDTH);
            int bytesPerImageRow = output_width*3*sizeof(GLubyte);
            int bytesPerTileRow = (tile_width-2*tile_border) * 3*sizeof(GLubyte);
            int xOffset = curColumn * bytesPerTileRow;
            int bytesPerCurrentTileRow = (curTileWidth-2*tile_border)*3*sizeof(GLubyte);
            int i;
            for (i=0;i<tile_height;i++) {
                memcpy(buffer + i*bytesPerImageRow + xOffset, /* Dest */
                tile + i*bytesPerTileRow,              /* Src */
                bytesPerCurrentTileRow);               /* Byte count*/
            }
         }
      
         if (curColumn == trGet(tr, TR_COLUMNS)-1) {
            /* write this buffered row of tiles to the file */
            int curTileHeight;
            int bytesPerImageRow = output_width*3*sizeof(GLubyte);
            int i;
            GLubyte *rowPtr;
#ifdef WIN32
// color order is different on WIN32: rgb -> gbr
            GLubyte tmpbyte[1];
            for (i=0;i<(output_width * tile_height * 3);i=i+3) {
                memcpy(tmpbyte, buffer + i*sizeof(GLubyte), 1);
                memmove(buffer + i*sizeof(GLubyte), buffer + i*sizeof(GLubyte) + 1, 1);
                memmove(buffer + i*sizeof(GLubyte) + 1, buffer + i*sizeof(GLubyte) + 2, 1);
                memcpy(buffer + i*sizeof(GLubyte) + 2, tmpbyte, 1);
            }
#endif
            curTileHeight = trGet(tr, TR_CURRENT_TILE_HEIGHT) - 2*tile_border;
            for (i=0;i<curTileHeight;i++) {
            /* Remember, OpenGL images are bottom to top.  Have to reverse. */
            rowPtr = buffer + (curTileHeight-1-i) * bytesPerImageRow;
            fwrite(rowPtr, 1, output_width*3, f);
	 }
      }

   }
   trDelete(tr);

   fclose(f);
   free(tile);
   free(buffer);
// -----
   }
   else {
      pixels = (GLubyte *) malloc(output_width * output_height * 3 * sizeof(GLubyte));
      if (!pixels) {
         showinfobox("Can't allocate memory for image");
         return;
      }

      tr = trNew();
      trTileSize(tr, tile_width, tile_height, tile_border);
      trImageSize(tr, output_width, output_height);
      trImageBuffer(tr, GL_RGB, GL_UNSIGNED_BYTE, pixels);
      tile_persp(tr);

      do {
         trBeginTile(tr);
         DrawScene(tr);
         tileCount++;
      } while (trEndTile(tr));
      sprintf(line, "%d tiles [%dx%d] drawn\n", tileCount, tile_width, tile_height);
      logprint(line);

      trDelete(tr);
   }

   bit.render = 0;

}

void render(int compression)
{
   char str[30];

   glutSetWindow(mainwin);
   glutPopWindow();
   save_w = glutGet(GLUT_WINDOW_WIDTH);
   save_h = glutGet(GLUT_WINDOW_HEIGHT);
   tile_render(compression);

   switch(compression) {
      case 0:
         writergb_rle(filename, "molekel", output_width, output_height);
      break;
      case PPM:
         // nothing to do
      break;
      case COMPRESSION_NONE:
      case COMPRESSION_LZW:
      case COMPRESSION_PACKBITS:
         writetiff(filename, "molekel", output_width, output_height, compression);
      break;
      default:
         writejpeg(filename, output_width, output_height, compression);
   }

   sprintf(str, "image [%dx%d] written", output_width, output_height);
   logprint(str);
   update_logs();
   glutSetWindow(mainwin);
   reshape(save_w, save_h);  
}


int gen_img_filename(int nbr)
{
   char basis[200] = "", tmp[200] = "";
   char *ext, *fnp;
   size_t i = 0;

   if((ext = strrchr(filename, '.')) == NULL) return 0;
   fnp = filename;
   while(fnp != ext) {
      i++;
      fnp++;
   }
   if((strlen(filename) - i) > 6) return 0;
   strncpy(basis, filename, i);
   sprintf(tmp, "%s.%03d%s", basis, nbr, ext);
   strcpy(filename, tmp);
   return 1;
}

void generate_struct_imgs(int compression)
{
   int i, j, w, h;
   char tmpfile[200] = "";
   j = frame_nbr;
   w = output_width;
   h = output_height;
   strcpy(tmpfile, filename);
   if(!dynamics.trajectory){
      showinfobox("no structures to animate!"); return;
   }
   if(actualmol != dynamics.molecule){
      showinfobox("wrong molecule!"); return;
   }
   if(!dynamics.end) dynamics.end = dynamics.ntotalsteps - 1;
   
   for(i = dynamics.start; i <= dynamics.end; i++) {
      output_width = w;
      output_height = h;
      dynamicsGoto(i);
      if(!gen_img_filename(j)){
         showinfobox("no valid filename extension!"); return;
      }
      render(compression);
      strcpy(filename, tmpfile);
      j++;
   }
}


void generate_freq_imgs(int compression)
{
   float factor, step = M_PI/10;
   register AtoM *ap;
   Vector *vw;
   int i, j, k, w, h;
   char tmpfile[200] = "";

   j = frame_nbr;
   w = output_width;
   h = output_height;
   strcpy(tmpfile, filename);

   if(actualvib) {
      if(!save_coords(actualmol)) return;
      for(k = 0; k < 21; k++) {
         output_width = w;
         output_height = h;
         factor = actualmol->sc_freq_ar * sin(freq_pos) / 2.5;
         freq_pos += step;
         for(ap = actualmol->firstatom, vw=coord_backup, i=0; ap; ap=ap->next, vw++, i++){
            ap->coord[0] = vw->x + actualvib->coord[i].x * factor;
            ap->coord[1] = vw->y + actualvib->coord[i].y * factor;
            ap->coord[2] = vw->z + actualvib->coord[i].z * factor;
         }
         if(!gen_img_filename(j)){
            showinfobox("no valid filename extension!"); return;
         }
         render(compression);
         strcpy(filename, tmpfile);
         j++;
      }
      freq_pos = 0;
      restore_coords(actualmol);
   }
   else {
      showinfobox("No vibration selected!");
      logprint("No vibration selected");
      logprint("   choose a vibration first");
      update_logs();
      return;
   }
}
