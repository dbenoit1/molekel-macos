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


/* ideas taken from writetiff of Mark J. Kilgard */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <GLUT/glut.h>
/*#include <GLUI/glui.h>*/
#include <GL/glui.h>
#ifdef __cplusplus      //include the jpeg library as a C file.
extern "C"{
#endif
/*#include<UnixImageIO/jpeglib.h>*/
#include <jpeglib.h>
#include <tiffio.h>
#ifdef __cplusplus
}
#endif



#include "snap.h"
#include "maininterf.h"

GLubyte *pixels = NULL;
int imgx, imgy;
int jpeg_qual = 75;

#ifdef LINUX
#define NORGB
#endif
#ifdef WIN32
#define NORGB
#endif

void getpixels( int xx, int yy, int width, int height)
{
  if((pixels = (GLubyte *) malloc(width * height * 3 * sizeof(GLubyte))) == NULL){
     showinfobox("Can't allocate memory for pixels");
     return;
  }
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadBuffer(GL_FRONT);
  glReadPixels(xx, yy, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);
  glReadBuffer(GL_BACK);
}

void writetiff(char *filename, char *description,
  int width, int height, int compression)
{
  TIFF *file;
  GLubyte *p;
  int i;
  char line[100];

  file = TIFFOpen(filename, "wb");
  if (file == NULL) {
    sprintf(line, "Can't open %s\n", filename);
    showinfobox(line);
    free(pixels);
    pixels = NULL;
    return;
  }

  TIFFSetField(file, TIFFTAG_IMAGEWIDTH, (uint32) width);
  TIFFSetField(file, TIFFTAG_IMAGELENGTH, (uint32) height);
  TIFFSetField(file, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(file, TIFFTAG_COMPRESSION, compression);
  TIFFSetField(file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(file, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFSetField(file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(file, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(file, TIFFTAG_IMAGEDESCRIPTION, description);
  p = pixels;
  for (i = height - 1; i >= 0; i--) {
    if (TIFFWriteScanline(file, p, i, 0) < 0) {
      free(pixels);
      pixels = NULL;
      TIFFClose(file);
      return;
    }
    p += width * sizeof(GLubyte) * 3;
  }
  free(pixels);
  pixels = NULL;
  TIFFClose(file);
  update_logs();
}


#ifndef NORGB
void writergb_rle(char *filename, char *description, int width, int height)
{
  GLushort *rrow, *grow, *brow;
  int  x, y;
  char line[100];
  IMAGE *image;
  FILE *dummy;

/* this dummy needed because iopen exits the program if it fails */
  if((dummy = fopen(filename, "w")) == NULL){
    sprintf(line, "Can't open %s\n", filename);
    showinfobox(line);
    free(pixels);
    pixels = NULL;
    return;
  }
  fclose(dummy);

  image = iopen(filename, "wb", RLE(1), 3, width, height, 3);
  if (image == NULL) {
    sprintf(line, "Can't open %s\n", filename);
    showinfobox(line);
    free(pixels);
    pixels = NULL;
    return;
  }
  if((rrow = (GLushort *) malloc(width * sizeof(GLushort))) == NULL) {
    showinfobox("Can't allocate memory for rrow");
    return;
  }
  if((grow = (GLushort *) malloc(width * sizeof(GLushort))) == NULL) {
    showinfobox("Can't allocate memory for grow");
    return;
  }
  if((brow = (GLushort *) malloc(width * sizeof(GLushort))) == NULL) {
    showinfobox("Can't allocate memory for brow");
    return;
  }

  for(y = 0; y < height; y++) {
     for(x = 0; x < width; x++) {
        rrow[x] = pixels[(y * width + x) * 3 + 0];
        grow[x] = pixels[(y * width + x) * 3 + 1];
        brow[x] = pixels[(y * width + x) * 3 + 2];
     }
     putrow(image, rrow, y, 0);
     putrow(image, grow, y, 1);
     putrow(image, brow, y, 2);
  }
  iclose(image);
  free(rrow);
  free(grow);
  free(brow);
  free(pixels);
  pixels = NULL;
  update_logs();
}
#else
/* writes uncompressed RGB */
void putbyte(FILE *outf, unsigned char val)
{
   unsigned char buf[1];
    
   buf[0] = val;
   fwrite(buf,1,1,outf);
}
     
void putshort(FILE *outf, unsigned short val)
{
   unsigned char buf[2];
    
   buf[0] = (val>>8);
   buf[1] = (val>>0);
   fwrite(buf,2,1,outf);
}
    
static int putlong(FILE *outf, unsigned long val)
{
   unsigned char buf[4];
    
   buf[0] = (val>>24);
   buf[1] = (val>>16);
   buf[2] = (val>>8);
   buf[3] = (val>>0);
   return fwrite(buf,4,1,outf);
}

void writergb(char *filename, char *description,
              int width, int height)
{
  FILE *file;
  GLubyte *row;
  int i, x, y;

  unsigned short xsize, ysize, zsize = 3;
  char name[80];
  int slot;

  file = fopen(filename, "wb");
  if (file == NULL) {
    return;
  }
  row = (GLubyte *) malloc(width * sizeof(GLubyte));
  xsize = (unsigned short) width;
  ysize = (unsigned short) height;
  strcpy(name, description);

  putshort(file,474);       /* MAGIC               */
  putbyte(file,0);          /* STORAGE is VERBATIM */
  putbyte(file,1);          /* BPC is 1            */
  putshort(file,3);         /* DIMENSION is 3      */
  putshort(file, xsize);    /* XSIZE               */
  putshort(file, ysize);    /* YSIZE               */
  putshort(file, zsize);    /* ZSIZE               */
  putlong(file,0);          /* PIXMIN is 0         */
  putlong(file,1);          /* PIXMAX is 1         */
  for(i=0; i<4; i++)        /* DUMMY 4 bytes       */
     putbyte(file,0);
  fwrite(name,80,1,file);   /* IMAGENAME           */
  putlong(file,0);          /* COLORMAP is 0       */
  for(i=0; i<404; i++)      /* DUMMY 404 bytes     */
     putbyte(file,0);

  for(y = 0; y < height; y++) {
     for(x = 0; x < width; x++) {
        row[x] = pixels[(y * width + x) * 3 + 0];
     }
     fwrite (row, 1, xsize, file);
  }
  for(y = 0; y < height; y++) {
     for(x = 0; x < width; x++) {
        row[x] = pixels[(y * width + x) * 3 + 1];
     }
     fwrite (row, 1, xsize, file);
  }
  for(y = 0; y < height; y++) {
     for(x = 0; x < width; x++) {
        row[x] = pixels[(y * width + x) * 3 + 2];
     }
     fwrite (row, 1, xsize, file);
  }

  free (row);
  
  fclose(file);
}

void writergb_rle(char *filename, char *description, int width, int height)
{
   writergb(filename, description, width, height);
   update_logs();
}
#endif


void writejpeg(char *filename, int width, int height, int compression)
{
  FILE *file;
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  GLubyte *p;
  int i;
  char line[100];

  file = fopen(filename, "wb");
  if (file == NULL) {
    sprintf(line, "Can't open %s\n", filename);
    showinfobox(line);
    free(pixels);
    pixels = NULL;
    return;
  }

  cinfo.err = jpeg_std_error (&jerr);
  jpeg_create_compress (&cinfo);
  jpeg_stdio_dest (&cinfo, file);

  cinfo.image_width = width;
  cinfo.image_height = height;
  cinfo.input_components = 3;
  cinfo.in_color_space = JCS_RGB;
  jpeg_set_defaults (&cinfo);
  jpeg_set_quality (&cinfo, compression, TRUE);
  jpeg_start_compress(&cinfo, TRUE);

   p = pixels;
   p += (height - 1) * width * sizeof(GLubyte) * 3;
   for (i = height - 1; i >= 0; i--) {
     jpeg_write_scanlines(&cinfo, &p, 1);
     p -= width * sizeof(GLubyte) * 3;
   } 

   jpeg_finish_compress (&cinfo);
   jpeg_destroy_compress (&cinfo);

   free(pixels);
   pixels = NULL;

   fclose(file);

   update_logs();
}
