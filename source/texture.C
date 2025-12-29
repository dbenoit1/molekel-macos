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
#include "texture.h"
#include "maininterf.h"
#include "textureinterf.h"
#include "sgi.h"

Texture *firsttexture = NULL;
Texture *actual_texture = NULL;
//GLenum texprop1[4] = {GL_REPEAT, GL_REPEAT, GL_LINEAR, GL_NEAREST_MIPMAP_NEAREST};
//GLenum texprop2[4] = {GL_CLAMP, GL_CLAMP, GL_POINT, GL_LINEAR_MIPMAP_NEAREST}; 
// MIPMAP does not work on default OpenGL on SGI's, Mesa 3.0 is o.k.
GLenum texprop1[4] = {GL_REPEAT, GL_REPEAT, GL_LINEAR, GL_LINEAR};
//GLenum texprop2[4] = {GL_CLAMP, GL_CLAMP, GL_POINT, GL_LINEAR}; 
//GL_REPEAT porduces better results with the default texture and property mapping with Mesa 3.0
GLenum texprop2[4] = {GL_CLAMP, GL_REPEAT, GL_LINEAR, GL_LINEAR}; 
GLubyte *contour_image1;
GLubyte *contour_image2;
GLubyte *contour_image3;
GLubyte *contour_image4;

Texture *bond_texture = NULL;
int bond_textype = 0;
GLenum bond_texenv = 0;

Texture *back_texture = NULL;
int back_textype = 0;
GLenum back_texenv = 0;

float start_texture = 0.0, end_texture = 1.0;

int assign_no_tex = 0;

void makeContourImage1(void)
{
   int i, linelen;
   GLubyte *cp;

   linelen = 512;

   if((contour_image1 = (GLubyte *)malloc(sizeof(GLubyte)*linelen*4)) == NULL){
      fprintf(stderr, "Can't allocate default texture\n");
      return;
   }

   cp = contour_image1;

/* first and last black line replaced by first and last color
 * better with flexible map range values vmin vmax
*/
   *cp++ = 0xff; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   for(i=1; i<51; i++) {
      *cp++ = 0xff; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff;
   } 
   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<103; i++) {
      *cp++ = 0xff; *cp++ = 0x5f;  *cp++ = 0x00;  *cp++ = 0xff;
   } 
   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<154; i++) {
      *cp++ = 0xff; *cp++ = 0xbf;  *cp++ = 0x00;  *cp++ = 0xff;
   } 
   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<205; i++) {
      *cp++ = 0xff; *cp++ = 0xff;  *cp++ = 0x00;  *cp++ = 0xff;
   } 
   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<256; i++) {
      *cp++ = 0x9f; *cp++ = 0xff;  *cp++ = 0x00;  *cp++ = 0xff;
   } 
   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<307; i++) {
      *cp++ = 0x00; *cp++ = 0xff;  *cp++ = 0x00;  *cp++ = 0xff;
   } 
   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<358; i++) {
      *cp++ = 0x00; *cp++ = 0xff;  *cp++ = 0x9f;  *cp++ = 0xff;
   } 
   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<409; i++) {
      *cp++ = 0x00; *cp++ = 0xff;  *cp++ = 0xff;  *cp++ = 0xff;
   } 
   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<460; i++) {
      *cp++ = 0x00; *cp++ = 0x7f;  *cp++ = 0xff;  *cp++ = 0xff;
   } 
   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<511; i++) {
      *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0xff;  *cp++ = 0xff;
   } 
   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0xff;  *cp++ = 0xff; 
}

void makeContourImages(void)
{
   int i, linelen;
   GLubyte *cp;

   linelen = 512;

   if((contour_image2 = (GLubyte *)malloc(sizeof(GLubyte)*linelen*4)) == NULL){
      fprintf(stderr, "Can't allocate default texture\n");
      return;
   }

   cp = contour_image2;

/* first and last black line replaced by first and last color
 * better with flexible map range values vmin vmax
*/
   *cp++ = 0xff; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   for(i=1; i<26; i++) {
      *cp++ = 0xff; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<51; i++) {
      *cp++ = 0xff; *cp++ = 0x2f;  *cp++ = 0x00;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<77; i++) {
      *cp++ = 0xff; *cp++ = 0x5f;  *cp++ = 0x00;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<102; i++) {
      *cp++ = 0xff; *cp++ = 0x8e;  *cp++ = 0x00;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<128; i++) {
      *cp++ = 0xff; *cp++ = 0xbf;  *cp++ = 0x00;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<153; i++) {
      *cp++ = 0xff; *cp++ = 0xee;  *cp++ = 0x00;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<179; i++) {
      *cp++ = 0xff; *cp++ = 0xff;  *cp++ = 0x00;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<204; i++) {
      *cp++ = 0xcf; *cp++ = 0xff;  *cp++ = 0x00;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<230; i++) {
      *cp++ = 0x9f; *cp++ = 0xff;  *cp++ = 0x00;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<255; i++) {
      *cp++ = 0x6f; *cp++ = 0xff;  *cp++ = 0x00;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<281; i++) {
      *cp++ = 0x00; *cp++ = 0xff;  *cp++ = 0x00;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<306; i++) {
      *cp++ = 0x00; *cp++ = 0xff;  *cp++ = 0x40;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<332; i++) {
      *cp++ = 0x00; *cp++ = 0xff;  *cp++ = 0x6f;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<357; i++) {
      *cp++ = 0x00; *cp++ = 0xff;  *cp++ = 0x9f;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<383; i++) {
      *cp++ = 0x00; *cp++ = 0xff;  *cp++ = 0xcf;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<408; i++) {
      *cp++ = 0x00; *cp++ = 0xff;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<434; i++) {
      *cp++ = 0x00; *cp++ = 0xbf;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<459; i++) {
      *cp++ = 0x00; *cp++ = 0x7f;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<485; i++) {
      *cp++ = 0x00; *cp++ = 0x3f;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<511; i++) {
      *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0xff;  *cp++ = 0xff;
   } 
   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0xff;  *cp++ = 0xff; 


   if((contour_image3 = (GLubyte *)malloc(sizeof(GLubyte)*linelen*4)) == NULL){
      fprintf(stderr, "Can't allocate default texture\n");
      return;
   }

   cp = contour_image3;

/* first and last black line replaced by first and last color
 * better with flexible map range values vmin vmax
*/
   *cp++ = 0xff; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   for(i=1; i<47; i++) {
      *cp++ = 0xff; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<94; i++) {
      *cp++ = 0xff; *cp++ = 0x33;  *cp++ = 0x33;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<140; i++) {
      *cp++ = 0xff; *cp++ = 0x66;  *cp++ = 0x66;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<186; i++) {
      *cp++ = 0xff; *cp++ = 0x99;  *cp++ = 0x99;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<232; i++) {
      *cp++ = 0xff; *cp++ = 0xcc;  *cp++ = 0xcc;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<279; i++) {
      *cp++ = 0xff; *cp++ = 0xff;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<325; i++) {
      *cp++ = 0xcc; *cp++ = 0xcc;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<371; i++) {
      *cp++ = 0x99; *cp++ = 0x99;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<417; i++) {
      *cp++ = 0x66; *cp++ = 0x66;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<464; i++) {
      *cp++ = 0x33; *cp++ = 0x33;  *cp++ = 0xff;  *cp++ = 0xff;
   } 
 
   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<511; i++) {
      *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0xff;  *cp++ = 0xff;
   } 
   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0xff;  *cp++ = 0xff;


   if((contour_image4 = (GLubyte *)malloc(sizeof(GLubyte)*linelen*4)) == NULL){
      fprintf(stderr, "Can't allocate default texture\n");
      return;
   }

   cp = contour_image4;

/* first and last black line replaced by first and last color
 * better with flexible map range values vmin vmax
*/
   *cp++ = 0xff; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   for(i=1; i<25; i++) {
      *cp++ = 0xff; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<50; i++) {
      *cp++ = 0xff; *cp++ = 0x19;  *cp++ = 0x19;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<75; i++) {
      *cp++ = 0xff; *cp++ = 0x32;  *cp++ = 0x32;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<99; i++) {
      *cp++ = 0xff; *cp++ = 0x4b;  *cp++ = 0x4b;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<123; i++) {
      *cp++ = 0xff; *cp++ = 0x64;  *cp++ = 0x64;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<147; i++) {
      *cp++ = 0xff; *cp++ = 0x7d;  *cp++ = 0x7d;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<171; i++) {
      *cp++ = 0xff; *cp++ = 0x96;  *cp++ = 0x96;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<195; i++) {
      *cp++ = 0xff; *cp++ = 0xaf;  *cp++ = 0xaf;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<219; i++) {
      *cp++ = 0xff; *cp++ = 0xc8;  *cp++ = 0xc8;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<243; i++) {
      *cp++ = 0xff; *cp++ = 0xe1;  *cp++ = 0xe1;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<268; i++) {
      *cp++ = 0xff; *cp++ = 0xff;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<292; i++) {
      *cp++ = 0xe1; *cp++ = 0xe1;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<316; i++) {
      *cp++ = 0xc8; *cp++ = 0xc8;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<340; i++) {
      *cp++ = 0xaf; *cp++ = 0xaf;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<364; i++) {
      *cp++ = 0x96; *cp++ = 0x96;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<388; i++) {
      *cp++ = 0x7d; *cp++ = 0x7d;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<412; i++) {
      *cp++ = 0x64; *cp++ = 0x64;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<436; i++) {
      *cp++ = 0x4b; *cp++ = 0x4b;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<460; i++) {
      *cp++ = 0x32; *cp++ = 0x32;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<485; i++) {
      *cp++ = 0x19; *cp++ = 0x19;  *cp++ = 0xff;  *cp++ = 0xff;
   } 

   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0x00;  *cp++ = 0xff; 
   i++;
   for(; i<511; i++) {
      *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0xff;  *cp++ = 0xff;
   } 
   *cp++ = 0x00; *cp++ = 0x00;  *cp++ = 0xff;  *cp++ = 0xff;

}


Texture *add_texture(GLubyte *image, int idepth, int iwidth, int iheight, GLenum *texprop)
/* texprop: GL_TEXTURE_WRAP_S, GL_TEXTURE_WRAP_T,
 *          GL_TEXTURE_MAG_FILTER, GL_TEXTURE_MIN_FILTER
*/
{
   Texture *texture, *tp;
   GLuint texName[2];

   if((texture = (Texture *)malloc(sizeof(Texture))) == NULL){
      showinfobox("can't allocate memory for new texture-struct\n");
      return NULL;
   }

#ifdef NOTEXBIND
   texName[0] = glGenLists(1);
   glNewList(texName[0], GL_COMPILE);
#else
   glGenTextures(1, texName);
   glBindTexture(GL_TEXTURE_2D, texName[0]);
#endif
/* read_sgi_add_alpha() returns just idepth 2 and 4 */
   gluBuild2DMipmaps(GL_TEXTURE_2D,  idepth, iwidth, iheight,
		      idepth == 4 ? GL_RGBA : GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, image);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, texprop[0]);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, texprop[1]);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, texprop[2]);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, texprop[3]);
#ifdef NOTEXBIND
   glEndList();
#else
   glBindTexture(GL_TEXTURE_2D, 0);
#endif

   texture->texname = texName[0];
   texture->image = image;
   texture->texprop[0] = texprop[0];
   texture->texprop[1] = texprop[1];
   texture->texprop[2] = texprop[2];
   texture->texprop[3] = texprop[3];
   texture->idepth = idepth;
   texture->iwidth = iwidth;
   texture->iheight = iheight;
   texture->next = NULL;

   if(!firsttexture) firsttexture = texture;
   else {
      for(tp=firsttexture; tp->next; tp=tp->next);
      tp->next = texture;
   }

   return texture;
}

void mod_texture(GLenum *texprop)
/* texprop: GL_TEXTURE_WRAP_S, GL_TEXTURE_WRAP_T,
 *          GL_TEXTURE_MAG_FILTER, GL_TEXTURE_MIN_FILTER
*/
{
#ifdef NOTEXBIND
   glNewList(actual_texture->texname, GL_COMPILE);
#else
   glBindTexture(GL_TEXTURE_2D, actual_texture->texname);
#endif
/* read_sgi_add_alpha() returns just idepth 2 and 4 */
   gluBuild2DMipmaps(GL_TEXTURE_2D,  actual_texture->idepth, actual_texture->iwidth, actual_texture->iheight,
		      actual_texture->idepth == 4 ? GL_RGBA : GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, actual_texture->image);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, texprop_live[0]);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, texprop_live[1]);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, texprop_live[2]);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, texprop_live[3]);
#ifdef NOTEXBIND
   glEndList();
#else
   glBindTexture(GL_TEXTURE_2D, 0);
#endif

   actual_texture->texprop[0] = texprop_live[0];
   actual_texture->texprop[1] = texprop_live[1];
   actual_texture->texprop[2] = texprop_live[2];
   actual_texture->texprop[3] = texprop_live[3];
}

void init_texture(void)
{
   makeContourImage1();
   actual_texture = add_texture(contour_image1, 4, 512, 1, texprop2); 
}

void add_contour_textures(void)
{
   makeContourImages();
   actual_texture = add_texture(contour_image2, 4, 512, 1, texprop2);
   actual_texture = add_texture(contour_image3, 4, 512, 1, texprop2);
   actual_texture = add_texture(contour_image4, 4, 512, 1, texprop2);
   actual_texture = firsttexture;
}


void enable_texture(Texture *texture, GLenum texenv, int type)
{
   static float plane1[] = {1, 0, 0, 0.5}, plane2[] = {0, 1, 0, 0.5};

   if(bit.texture){
/* to be fixed
      if((previous_texture == texnum) && (previous_textype == type)) return;
      previous_texture = texnum;
      previous_textype = type;
*/
      glEnable(GL_TEXTURE_2D);
#ifdef NOTEXBIND
      glCallList(texture->texname);
#else
      glBindTexture(GL_TEXTURE_2D, texture->texname);
#endif
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, texenv);
      switch(type) {
         case TEX_MAP :
            glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
            glTexGenfv(GL_S, GL_OBJECT_PLANE, plane1);
            glEnable(GL_TEXTURE_GEN_S);
            glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
            glTexGenfv(GL_T, GL_OBJECT_PLANE, plane2);
            glEnable(GL_TEXTURE_GEN_T);
            break;
         case TEX_REFLECT :
            glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
            glEnable(GL_TEXTURE_GEN_S);
            glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
            glEnable(GL_TEXTURE_GEN_T);
            break;
         case TEX_PHONG :
            glDisable(GL_TEXTURE_GEN_S);
            glDisable(GL_TEXTURE_GEN_T);
            break;
      }

      glAlphaFunc(GL_NOTEQUAL, 0);

   }
}

void disable_texture(void)
{
   glDisable(GL_TEXTURE_2D);
#ifndef NOTEXBIND
   glBindTexture(GL_TEXTURE_2D, 0);
#endif
   glDisable(GL_TEXTURE_GEN_S);
   glDisable(GL_TEXTURE_GEN_T);
   glAlphaFunc(GL_ALWAYS, 0);
}


void texture_down(void)
{
   if(actual_texture->next) actual_texture = actual_texture->next;
}


void texture_up(void)
{
   Texture *tp, *init;

   init = actual_texture;
   if(firsttexture->next) {
      for(tp=firsttexture; tp->next; tp=tp->next) {
         if(tp->next == init) {
            actual_texture = tp;
            return;
         }
      }
   }
   

}

void read_img(char *file)
{
   int iheight, iwidth, idepth;
   GLubyte *image;

   image = (GLubyte *)read_sgi_add_alpha(file, &iwidth, &iheight, &idepth);
   actual_texture = add_texture(image, idepth, iwidth, iheight, texprop1);
   glutSetWindow(textureinterfwin);
   glutPostRedisplay();
}

void set_texture_range(void)
{
   register int x, y, i;
   int firstpixel, lastpixel;
   GLubyte *cp;

//   if(start_texture == 0.0 && end_texture == 1.0) retun;

      firstpixel = (int)(actual_texture->iwidth * start_texture);
      lastpixel  = (int)(actual_texture->iwidth * end_texture);

      for(y=0; y<actual_texture->iheight; y++) {
         cp = actual_texture->image + y*actual_texture->iwidth*actual_texture->idepth + actual_texture->idepth - 1;
i = y*actual_texture->iwidth*actual_texture->idepth + actual_texture->idepth - 1;
         for(x=0; x<actual_texture->iwidth; x++, cp+=actual_texture->idepth, i+=actual_texture->idepth) {
            if(firstpixel <= lastpixel){
               if(x >=firstpixel && x <= lastpixel) *cp = 0xff;
               else                                 *cp = 0x00;
            }
            else {
               if(x >=firstpixel || x <= lastpixel) *cp = 0xff;
               else                                 *cp = 0x00;
            }
         }
      }
}

/**** functions for simulated phong-shading using textures ***/

static Matrix pmat;

void get_pmat(void)
{
   glGetFloatv(GL_MODELVIEW_MATRIX, &pmat[0][0]);
}


static void vecxmat(float a[3], float b[3])
/* rotation only */
{
   register int i;

   for(i=0; i<3; i++)
      b[i] = a[0]*pmat[0][i] + a[1]*pmat[1][i] + a[2]*pmat[2][i];
}




void myt3f(float a[3])
{
   float b[3], t[2], l;

   vecxmat(a, b);

   l = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
   if(l) {b[0]/=l; b[1]/=l;}

   t[0] = b[0]*0.5 + 0.5;
   t[1] = b[1]*0.5 + 0.5;

   glTexCoord2fv(t);
}
