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



typedef struct GAUSS { double exponent, coeff, coeff2;
                       struct GAUSS *next;
                     } Gauss;
typedef struct SHELL { short  n_base;
                       float  scale_factor;
                       Gauss  *firstgauss, *actualgauss;
                       struct SHELL *next;
                     } Shell;
typedef struct SLATER { short  n;
                        char   type[2];
                        unsigned a : 4;
                        unsigned b : 4;
                        unsigned c : 4;
                        unsigned d : 4;
                        float  exponent, norm[5];
                        struct SLATER *next;
                      } Slater;
typedef struct       { short ns, np, nd;
                       float exps, expp, expd;
                       float couls, coulp, could;
                       double *norm;
                     } Valence_shell;
typedef struct AMOSS { char center[3], basis[21];
                       Shell  *firstshell;
                       struct AMOSS *next;
                     } Amoss_basis;
typedef struct BASIS { char basisname[21];
                       short ord;
                       Shell  *firstshell;
                       Slater *firstslater;
                       struct BASIS *next;
                     } Basis;
typedef struct TEXTURE_STRUCT {
                       GLuint texname;
                       GLubyte *image;
                       GLenum texprop[4];
                       int idepth;
                       int iwidth;
                       int iheight;
                       struct TEXTURE_STRUCT *next;
                     } Texture;

#ifdef WIN32
typedef struct  ATM  { short ord, nbonds, name;
#else
typedef struct ATOM  { short ord, nbonds, name;
#endif
                       unsigned picked      : 1;
                       unsigned het         : 1;
                       unsigned planar      : 1;
                       unsigned main        : 1;
		       unsigned fixed	    : 1;
                       float coord[3], charge, coordination, spin;
                       float force[3];
                       Shell  *firstshell, *actualshell;
                       Valence_shell *valence;
                       Slater *firstslater;
#ifdef WIN32
                       struct ATM *next;
#else
                       struct ATOM *next;
#endif
                     } AtoM;
typedef struct BOND  { AtoM *a, *b;
                       /* float middle[3]; */
                       short type, name, picked;
#ifdef BOND_MATRIX
                       Matrix bondMatrix;
#endif
                       struct BOND *next;
                     } Bond;
typedef struct TER   { AtoM *first, *last;
                       short colindex;
                       struct TER *next;
                     } Ter;
typedef struct RES   { AtoM *first, *last;
                       char type[4];
                       struct RES *next;
                     } Residue;

typedef struct { float v[3]; float n[3];}   Surfdot;
typedef struct { unsigned int p1, p2, p3; } Triangle; 
typedef struct { float x, y, z; } Vector;
typedef struct SURF { Surfdot  *dot;
                      Triangle *tri;
                      Vector   *trinorm;
                      float contour, *val, vmin, vmax;
                      float *val1, vmin1, vmax1; /* second property */
                      int   npts, ntri, matindex;
                      Texture *texture;
                      int textype;
                      GLenum texenv;
                      unsigned char alpha;
                      char  type, *name;
                      float clipplane_tvec[3], clipplane_rvec[4];
                      int surf_transp    ; 
                      int surf_clip      ; 
                      int chickenwire    ;
                      int flatshade      ;
                      int dot_surf       ;
                      struct SURF *next, *second;
                    } Surface;

typedef struct CUTPLANE { float alpha, beta, a[3];
                       int npts, ntri;
                       float dd;
                       float (**plane_point)[7];
//                     x, y, xrot, yrot, zrot, cubeval, visible
                       unsigned int (*tri)[3][2];
                     } Cutplane;

struct MOL;
typedef struct DYNA { Vector **trajectory;
                      AtoM **freeat;
                      struct MOL *molecule;
                      long ntotalsteps, nfreat;
                      long start, end, current;
                      int isrunning, runtype, direction;
                      float stepsize, timestep;
                    } Dynamics;

typedef struct      { float x1, x2, y1, y2, z1, z2, cubesize;
                      short nx, ny, nz;
                      short flag;
                    } Box;
typedef struct MOND { AtoM *a, *b;
                      struct MOND *next;
                    } Mon_dist;
typedef struct MONA { AtoM *a, *b, *c;
                      struct MONA *next;
                    } Mon_ang;
typedef struct MONT { AtoM *a, *b, *c, *d;
                      struct MONT *next;
                    } Mon_tor;

typedef struct { char type[8];
                 int flag;
                 int number;
                 float  occ;
                 double eigenvalue;
                 double *coefficient;
               } MolecularOrbital;
typedef struct { char   type[5];
                 float  frequency;
                 Vector *coord;
               } Vibration;

typedef struct { float start[3];
                 float end[3];
                 float absolute;
               } Dipole;

typedef struct ATOMTYPE {
      short ord;
      int natoms;
      AtoM **atomList;
      struct ATOMTYPE *next;

} AtomType;



typedef struct MOL  { AtoM *firstatom, *lastatom;
                      AtomType *firstAtomType;
                      Bond *firstbond, *lastbond;
                      Surface *firstsurf, *lastsurf;
                      Ter *firstter;
                      Residue *firstresidue;
                      Amoss_basis *firstamoss;
                      Basis *basisset;
                      Mon_dist *firstdist;
                      Mon_ang  *firstang;
                      Mon_tor  *firsttor;
                      Box box;
                      int plane_iz, plane_iy, plane_ix;
                      MolecularOrbital *alphaOrbital, *betaOrbital;
                      int nMolecularOrbitals, nBasisFunctions;
                      int firstOrbital, lastOrbital, n_frequencies;
                      int multiplicity;
                      float **alphaDensity, **betaDensity;
                      float charge, mass;
                      int nAlpha, nBeta, nElectrons, alphaBeta;
                      Vibration *vibration, *freq_arrow;
                      Dipole *dipole;
                      float sc_freq_ar;
                      float sc_dipole_ar;
                      Cutplane *plane;

                      char *filename;
                      void *dvs;
                      void *anim_data;
                      short natoms, nbonds, n_h_bonds, nsurfs, nters, nresidues;
                      float tvec[3], rvec[4], centervec[3];
                      float ***cube_value, cubemin, cubemax, cutoff;
                      int wire_model     ;
                      int sticks         ;
                      int ball_and_stick ;
                      int spacefill      ;
                      int suppress_H     ;
                      int labels         ;
                      int atm_char       ;
                      int atm_spin       ;
                      int transparent    ;
                      int box_on         ;
                      int multiple_bonds ;
                      int got_macufile   ;
                      int cubeplanes     ;
                      int draw_contour   ;
                      int charges        ;
                      int outline        ;
                      int ribbon         ;
                      int main_chain     ;
                      int residues       ;
                      int bond_col       ;
                      int h_bond         ;
                      int show_freq_arrow;
                      int show_dipole    ;
                      struct MOL *next;
                    } Mol;

/* global attributes */

typedef struct INTD { AtoM *a, *b;
                      Mol *ma, *mb;
                      struct INTD *next;
                    } Int_dist;


typedef struct { int persp          ; /* global */
                 int background     ;
                 int both_signs     ;
                 int fullscreen     ;
                 int addsurf        ;
                 int depthcue       ;
                 int smoothline     ;
                 int manip_all      ;
                 int independent    ;
                 int rc_stereo      ;
                 int texture        ;
                 int use_track      ;
                 int high_qual      ;
                 int singlebuf      ;
                 int acbuf          ;
                 int use_suffix     ;
                 int activ_logo     ;
                 int lbl_on_top     ;
                 int no_measure_lbl ;
                 int solid_arrow    ;
                 int keepbonds      ;
                 int superimpose    ;
                 int showfixed      ;
                 int tgl_idle       ;
                 int render         ;
                 int coord_sys      ;
                      int wire_model     ;/* templates for attributes*/
                      int sticks         ;/* per molecule */
                      int ball_and_stick ;
                      int spacefill      ;
                      int suppress_H     ;
                      int labels         ;
                      int atm_char       ;
                      int atm_spin       ;
                      int transparent    ;
                      int box_on         ;
                      int multiple_bonds ;
                      int got_macufile   ;
                      int cubeplanes     ;
                      int draw_contour   ;
                      int outline        ;
                      int ribbon         ;
                      int main_chain     ;
                      int residues       ;
                      int bond_col       ;
                      int h_bond         ;
                      int show_freq_arrow;
                      int show_dipole         ;
                            int surf_transp    ; /* templates for attributes*/
                            int surf_clip      ; /* per surface */
                            int chickenwire    ;
                            int flatshade      ;
                            int dot_surf       ;
               } Bit;


typedef struct { int bury;
                 int  lon;
                 int  bin;
               } Conflag;
typedef struct { float mass, rvdw, coval;
		 short col;
                 Texture *texture;
                 int textype;
                 GLenum texenv;
                 short nbond;
                 char  symbol[3];
               } Element;
typedef float Matrix[4][4];

typedef struct { AtoM **a;
                 short na;
               } Atomcube;


/* definitions */

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define ANGLE_PI 3.14159265358979323846
#define BOHR         0.529177249
#define _1_BOHR      1.88972599

#ifdef WIN32
#define NO_F_FUNC
#define M_PI            3.14159265358979323846
#define M_1_PI          0.31830988618379067154
#define M_2_PI          0.63661977236758134308
#endif

#ifdef SUN
#define NO_F_FUNC
#endif 

#ifdef LINUX
#define NO_F_FUNC
#endif

#ifdef NO_F_FUNC
#define fsin sin
#define fcos cos
#define facos acos
#define fsqrt sqrt
#define fabsf fabs
#define floorf floor
#define ffloor floor
#define fceil ceil
#endif
