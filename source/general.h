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


#include <string.h>
#include <ctype.h>
/*#include <malloc.h>*/

/*----------------------------------------*/
extern Bit bit;
extern FILE *fhelp;
extern char starting_directory[], actual_directory[];
extern Element element[];
extern Mol *firstmol, *actualmol;
extern Surface *actualsurf;
extern Int_dist *first_intermolecular_distance;
extern float sc_freq_ar;
extern float sc_dipole_ar;
extern float freq_pos;
extern Vector *coord_backup;
extern Vibration *actualvib;
extern int quality;
extern int alphablending;
extern int got_dots;
extern int label_size;
extern int cubeplane_axis;
extern int xyzformat;
extern int maplegend;
extern float max_h_bond, min_h_angle;
extern char t41cont[];
extern Dynamics dynamics;

#define NLIST         105
/*----------------------------------------*/

void read_table(void);
Mol *add_mol(char *file);
AtoM *add_atom(int ord, float x, float y, float z);
int updateAtomList(AtoM *ap);
AtomType *addAtomType(int ordinal);
int addAtomToList(AtoM *ap, AtomType *at);
void freeAtomList(AtomType *ap);
void add_bond(register AtoM *i, register AtoM *j);
Ter *add_ter(AtoM *a, AtoM *b);
Residue *add_residue(AtoM *a, AtoM *b, char *s);
Mon_dist *add_distance(AtoM *a, AtoM *b);
Int_dist *add_intermolecular_distance(AtoM *a, Mol *ma, AtoM *b, Mol *mb);
Mon_ang *add_angle(AtoM *a, AtoM *b, AtoM *c);
Mon_tor *add_torsion(AtoM *a, AtoM *b, AtoM *c, AtoM *d);
Dipole *add_dipole(float x, float y, float z);
void removeMolFromList(Mol *x);
void appendMolToList(Mol *x);
void free_mol(Mol *x);
void free_surfaces(Surface *x);
void free_amoss(Amoss_basis *x);
void free_basisset(Basis *x);
void free_slater(Slater *x);
void free_atoms(AtoM *x);
void free_shell(Shell *shell);
void free_gauss(Gauss *gauss);
void free_bonds(Bond *x);
void free_ters(Ter *x);
void free_residues(Residue *x);
void free_distances(Mon_dist *x);
void free_angles(Mon_ang *x);
void free_torsions(Mon_tor *x);
void free__intermolecular_distances(Int_dist *x);
void free_mo(Mol *mp);
void free_dyna(void);
int remIntermolecularDistances(Mol *mp);
void remIntermolecularDistance(Int_dist *x);
float dist(float *a, float *b);
void new_mole(char *file);
void restore_coords(Mol *mp);
int save_coords(Mol *mp);
void remove_values(void);
void rem_bond_from_list(Bond *bp);
void modify_bond_type(short type, Bond *b);
MolecularOrbital *allocOrbital(int nOrbitals, int nBasis, int flag);
MolecularOrbital *reallocOrbital(MolecularOrbital *orbital, int newer, int old);
void *alloc_trimat(int n, size_t size);
Amoss_basis *add_amoss(Mol *mp);
Shell *add_shell_to_amoss(Amoss_basis *ap);
Gauss *add_gauss(Shell *sp);
Shell *add_shell(AtoM *ap);
Slater *add_slater(AtoM *ap);
void help(void);
