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


double Zs[120]={0,
    /* 1 H */   0.967807,0,0,0.877439,0,
    /* 6 C */   1.565085,2.028094,3.796544,4.708555,0,
    /*11 Na*/   0,0.698552,1.702888,1.635075,2.017563,
    /*16 S */   1.891185,2.24621,0,0,0,
    /*21 Sc*/   0,0,0,0,0,
    /*26 Fe*/   0,0,0,0,1.819989,
    /*31 Ga*/   1.84704,2.237353,2.636177,2.828051,5.348457,
    /*36 Kr*/   0,0,0,0,0,
    /*41 Nb*/   0,0,0,0,0,
    /*46 Pd*/   0,0,1.679351,2.016116,2.373328,
    /*51 Sb*/   2.343039,4.165492,7.00113,0,0,
    /*56 Ba*/   0,0,0,0,0,
    /*61 Pm*/   0,0,0,0,0,
    /*66 Tb*/   0,0,0,0,0,
    /*71 Lu*/   0,0,0,0,0,
    /*76 Os*/   0,0,0,0,1.476885,
    /*81 Tl*/   6.867921,3.141289,4.916451,0,0,
    /*86 Rn*/   0,0,0,0,0,
    /*91 Pa*/   0,0,0,0,0,
    /*96 Cm*/   0,0,0,0,0,
   /*101 Md*/   0,0,0,0,0,0,0}; /* au */
/*alle PM3-Parameter neu */
double Zp[120]={0,
    /* 1 H */   0,0,0,1.508755,0,
    /* 6 C */   1.842345,2.313728,2.389402,2.491178,0,
    /*11 Na*/   0,1.483453,1.073629,1.313088,1.504732,
    /*16 S */   1.658972,2.15101,0,0,0,
    /*21 Sc*/   0,0,0,0,0,
    /*26 Fe*/   0,0,0,0,1.506922,
    /*31 Ga*/   0.839411,1.592432,1.703889,1.732536,2.12759,
    /*36 Kr*/   0,0,0,0,0,
    /*41 Nb*/   0,0,0,0,0,
    /*46 Pd*/   0,0,2.066412,1.44535,1.638233,
    /*51 Sb*/   1.899992,1.647555,2.454354,0,0,
    /*56 Ba*/   0,0,0,0,0,
    /*61 Pm*/   0,0,0,0,0,
    /*66 Tb*/   0,0,0,0,0,
    /*71 Lu*/   0,0,0,0,0,
    /*76 Os*/   0,0,0,0,2.479951,
    /*81 Tl*/   1.969445,1.892418,1.934935,0,0,
    /*86 Rn*/   0,0,0,0,0,
    /*91 Pa*/   0,0,0,0,0,
    /*96 Cm*/   0,0,0,0,0,
   /*101 Md*/   0,0,0,0,0,0,0}; /* au */
