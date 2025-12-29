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


/* nongraphic routines : file-I/O, data-structuring ...
 * was main in old molekel
*/

#include "constant.h"
#include "main.h"
#include "molekel.h"
#include "general.h"
#include "drawing.h"
#include "macu.h"
#include "manip.h"
#include "maininterf.h"
#include "chooseinterf.h"

/*----------------------------------------*/
#define NLIST         105

char starting_directory[128], actual_directory[128];
Element element[NLIST];
Mol *firstmol = NULL, *actualmol = NULL;
Surface *actualsurf = NULL;
Int_dist *first_intermolecular_distance = NULL;
float sc_freq_ar = 0.5;
float sc_dipole_ar = 0.5;
float freq_pos;
Vector *coord_backup = NULL;
Vibration *actualvib = NULL;
int quality = 0;
int alphablending = 0;
int got_dots = 0;
int label_size = 0;
int cubeplane_axis = 2;
int xyzformat = 0;
int maplegend = 1;
float max_h_bond = 2.3, min_h_angle = 90.0;
char t41cont[20] = {"Density"};
Dynamics dynamics = {(Vector **)NULL,
                     (AtoM **)NULL,
                     (Mol *)NULL,
                     0,    /* ntotalsteps */
                     0,    /* nfreat */
                     0,    /* start */
                     0,    /* end */
                     0,    /* current */
                     0,    /* isrunning */
                     2,    /* runtype */
                     1,    /* direction */
                     1.0,  /* stepsize */
                     1.0}; /* timestep */

/*----------------------------------------*/
Bit bit = {1,  /* persp          */
           1,  /* background     */
           0,  /* both_signs     */
           0,  /* fullscreen     */
           0,  /* addsurf        */
           0,  /* depthcue       */
           0,  /* smoothline     */
           0,  /* manip_all      */
           1,  /* independent    */
           0,  /* rc_stereo      */
           1,  /* texture        */
           1,  /* use_track      */
           1,  /* high_qual      */
           0,  /* singlebuf      */
           0,  /* acbuf          */
           1,  /* use_suffix     */
           1,  /* activ_logo     */
           1,  /* lbl_on_top     */
           0,  /* no_measure_lbl */
           1,  /* solid_arrow    */
           1,  /* keepbonds      */
           0,  /* superimpose    */
           1,  /* showfixed      */
           0,  /* tgl_idle       */
           0,  /* render         */
           0,  /* coord_sys      */
           0,  /* wire_model     */
           1,  /* sticks         */
           0,  /* ball_and_stick */
           0,  /* spacefill      */
           0,  /* suppress_H     */
           0,  /* labels         */
           0,  /* atm_char       */
           0,  /* atm_spin       */
           0,  /* transparent    */
           0,  /* box_on         */
           1,  /* multiple_bonds */
           0,  /* got_macufile   */
           0,  /* cubeplanes     */
           0,  /* draw_contour   */
           0,  /* outline        */
           0,  /* ribbon         */
           1,  /* main_chain     */
           0,  /* residues       */
           0,  /* bond_col       */
           0,  /* h_bond         */
           0,  /* freq_arrow     */
           0,  /* show dipole    */
           0,  /* surf_transp    */
           0,  /* surf_clip      */
           0,  /* chickenwire    */
           0,  /* flatshade      */
           0,  /* dot_surf       */
          };
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
int remIntermolecularDistances(Mol *mp);
void remIntermolecularDistance(Int_dist *x);
float dist(float *a, float *b);

/*----------------------------------------*/

/* read atomic data from "atoms.data" */
void read_table(void)
{
   char name[100], *s = NULL;
   Element *ep;
   int ord;
   FILE *fp;

   (void)getcwd(starting_directory, 128);
   (void)strcpy(actual_directory, starting_directory);
#ifdef WIN32
   (void)strcat(actual_directory, "\\");
#else
   (void)strcat(actual_directory, "/");
#endif
//printf("%s\n", actual_directory);

/*   define an environment variable for the correct path */
   if((s = getenv("MOLEKEL")) != NULL){
      (void)strcpy(name, s);
   }
   else {
      /* (void)fprintf(stderr, " Can't find environment variable molekel\n"); */
      /* yeah, I know : junior-programming... (greetings to MWT) */
      (void)strcpy(name, "/usr/local/lib/molekel");
   }
#ifdef WIN32
   (void)strcat(name, "\\atoms.data");
#else
   (void)strcat(name, "/atoms.data");
#endif

   if((fp = fopen("atoms.data", "r")) == NULL){
      if((fp = fopen(name, "r")) == NULL){
         (void)fprintf(stderr, " Can't open %s\n", name);
         exit(-1);
      }
   }

   while(fgets(name, 99, fp)){
      if(name[0] == '#') continue;
      (void)sscanf(name, "%d", &ord);
      ep = element + ord;
      (void)sscanf(name, "%*d %s %f %f %hd %hd %f", ep->symbol,
         &ep->coval, &ep->rvdw, &ep->nbond, &ep->col, &ep->mass);
   }
   (void)fclose(fp);

/*
 * moved to beginning of function for debugging purpose
   (void)getcwd(starting_directory, 128);
   (void)strcpy(actual_directory, starting_directory);
   (void)strcat(actual_directory, "/");
printf("%s\n", actual_directory);
*/
}


/* add a molecule to the linked list of molecules, initialize the molecule */
Mol *add_mol(char *file)
{
   Mol *last;

   if((last = (Mol*) malloc(sizeof(Mol))) == NULL){
      (void)printf(" can't allocate enough memory\n");
      exit(-1);
   }

//   if(actualmol) free_3D(actualmol);

   appendMolToList(last);

   last->firstatom = NULL;
   last->firstAtomType = NULL;
   last->firstbond = NULL;
   last->firstsurf = NULL;
   last->firstter = NULL;
   last->firstresidue = NULL;
   last->firstamoss = NULL;
   last->basisset = NULL;
   last->firstdist = NULL;
   last->firstang = NULL;
   last->firsttor = NULL;
   last->alphaOrbital = last->betaOrbital = NULL;
   last->alphaDensity = last->betaDensity = NULL;
   last->vibration = NULL;
   last->freq_arrow = NULL;
   last->sc_freq_ar = 5 * sc_freq_ar;
   last->nMolecularOrbitals = last->nBasisFunctions = last->multiplicity = 0;
   last->n_frequencies = 0;
   last->plane_iz = last->plane_iy = last->plane_ix = 0;
   last->firstOrbital = last->lastOrbital = 1;
   last->nAlpha = last->nBeta = last->nElectrons = last->alphaBeta = 0;
   if(file) last->filename = strdup(file);
   else     last->filename = NULL;
   last->tvec[0] = last->tvec[1] = last->tvec[2] = 0.0;
   last->rvec[0] = last->rvec[1] = last->rvec[2] = 0.0;
   last->rvec[3] = 1.0;
   last->centervec[0] = last->centervec[1] = last->centervec[2] = 0.0;
   last->dipole = NULL;
   last->sc_dipole_ar = 2 * sc_dipole_ar;
   last->plane = NULL;
   last->mass = 0;
   last->charge = 0;
   last->dvs = 0;
   last->natoms = last->nbonds = last->n_h_bonds = last->nsurfs = 0;
   last->nters = last->nresidues = 0;
   last->next = NULL;
   last->wire_model     = bit.wire_model;
   last->sticks         = bit.sticks;
   last->ball_and_stick = bit.ball_and_stick;
   last->spacefill      = bit.spacefill;
   last->suppress_H     = bit.suppress_H;
   last->labels         = bit.labels;
   last->atm_char       = bit.atm_char;
   last->atm_spin       = bit.atm_spin;
   last->transparent    = bit.transparent;
   last->box_on         = bit.box_on;
   last->multiple_bonds = bit.multiple_bonds;
   last->got_macufile   = 0;
   last->cubeplanes     = 0;
   last->draw_contour   = bit.draw_contour;
   last->outline        = bit.outline;
   last->ribbon         = bit.ribbon;
   last->main_chain     = bit.main_chain;
   last->residues       = bit.residues;
   last->bond_col       = bit.bond_col;
   last->h_bond         = bit.h_bond;
   last->show_freq_arrow= 0;
   last->show_dipole    = bit.show_dipole;

   last->cube_value = NULL;
   last->cubemin = last->cubemax = last->cutoff = 0;

   actualmol = last;
   got_dots = 0;
   actualvib = NULL;

   return last;
}



/* add an atom to the linked list of atoms, read the data into structure */
AtoM *add_atom(int ord, float x, float y, float z)
{
   AtoM *temp;

   if((temp = (AtoM*) malloc(sizeof(AtoM))) == NULL){
      (void)printf(" can't allocate enough memory\n");
      exit(-1);
   }

   temp->ord = ord;
   temp->nbonds = 0;
   temp->name = ++(actualmol->natoms);
   temp->picked = 0;
   temp->het = 0;
   temp->planar = 0;
   temp->main = 0;
   temp->fixed = 1;
#ifdef ATOM_NEIGHBORS
   temp->nneighbors = 0;
#endif
   temp->coord[0] = x;
   temp->coord[1] = y;
   temp->coord[2] = z;
   temp->coordination = 0.0;
   temp->charge = 0.0;
   temp->spin = 0.0;
   temp->force[0] = temp->force[1] = temp->force[2] = 0.0;

   temp->firstshell = temp->actualshell = NULL;
   temp->firstslater = NULL;
   temp->valence = NULL;
   temp->next = NULL;

   if(!actualmol->firstatom) actualmol->firstatom = temp;
   else actualmol->lastatom->next = temp;
   actualmol->lastatom = temp;

   return temp;
}




int updateAtomList(AtoM *ap)
{
   AtomType *temp, *at, *previous;

   if(!actualmol->firstAtomType) {
      if((temp = addAtomType(ap->ord)) == NULL) return 0;
      if(!addAtomToList(ap, temp)) return 0;
      actualmol->firstAtomType = temp;
   }
   else {
      for(at=actualmol->firstAtomType; at; at=at->next) {
         if(ap->ord == at->ord) {
            if(!addAtomToList(ap, at)) return 0;
            return 1;
         }
         previous = at;
      }
      
      if((temp = addAtomType(ap->ord)) == NULL) return 0;
      if(!addAtomToList(ap, temp)) return 0;
      previous->next = temp;
   }

   return 1;
}



AtomType *addAtomType(int ordinal)
{
   AtomType *at;

   if((at = (AtomType*) malloc(sizeof(AtomType))) == NULL) {
      (void)printf("can't allocate memory for the atom type %d\n", ordinal);
      return 0;
   }

   at->ord = ordinal;
   at->natoms = 0;
   at->atomList = NULL;
   at->next = NULL;

   return at;
}



int addAtomToList(AtoM *ap, AtomType *at)
{
   if(at->natoms) {
      if((at->atomList = (AtoM**) realloc(at->atomList, 
                         (at->natoms+1)*sizeof(AtoM *))) == NULL) {
         return 0;
      }
   }
   else {
      if((at->atomList = (AtoM**) malloc(sizeof(AtoM *))) == NULL) {
         return 0;
      }
   }

   at->atomList[at->natoms++] = ap;
   return 1;
}


void freeAtomList(AtomType *ap)
{
   if(ap->next) freeAtomList(ap->next);
   free(ap->atomList);
   free(ap);
}


/* add a bond to the linked list of bonds */
void add_bond(register AtoM *i, register AtoM *j)
{
   Bond *temp;

   if((temp = (Bond*) malloc(sizeof(Bond))) == NULL){
      (void)printf(" can't allocate enough memory\n");
      exit(-1);
   }

   temp->a = i;
   temp->b = j;
/*
   temp->middle[0] = (i->coord[0] + j->coord[0])*0.5;
   temp->middle[1] = (i->coord[1] + j->coord[1])*0.5;
   temp->middle[2] = (i->coord[2] + j->coord[2])*0.5;
*/
   temp->type = SINGLEBOND;
   temp->picked = 0;
   temp->next = NULL;
   temp->name = actualmol->natoms + actualmol->nbonds;

   if(!actualmol->firstbond) actualmol->firstbond = temp;
   else actualmol->lastbond->next = temp;
   actualmol->lastbond = temp;

   i->nbonds++;
   j->nbonds++;
}




/* add a ter to the linked list of ters, read the data into structure */
Ter *add_ter(AtoM *a, AtoM *b)
{
   static short col = 0;
   Ter *temp, *tp;

   if((temp = (Ter*) malloc(sizeof(Ter))) == NULL){
      (void)fprintf(stderr, " can't allocate enough memory for ter\n");
      exit(-1);
   }

   temp->first = a;
   temp->last = b;
   temp->colindex = col++;
   temp->next = NULL;

   if(!actualmol->firstter) actualmol->firstter = temp;
   else {
      for(tp=actualmol->firstter; tp->next; tp=tp->next);
      tp->next = temp;
   }

   return temp;
}



/* add a residue to the linked list, read the data into structure */
Residue *add_residue(AtoM *a, AtoM *b, char *s)
{
   Residue *temp, *rp;

   if((temp = (Residue*) malloc(sizeof(Residue))) == NULL){
      (void)fprintf(stderr, " can't allocate enough memory for residue\n");
      exit(-1);
   }

   temp->first = a;
   temp->last = b;
   temp->type[0] = s[0];
   temp->type[1] = s[1];
   temp->type[2] = s[2];
   temp->type[3] = 0;
   temp->next = NULL;

   if(!actualmol->firstresidue) actualmol->firstresidue = temp;
   else {
      for(rp=actualmol->firstresidue; rp->next; rp=rp->next);
      rp->next = temp;
   }

   return temp;
}



/* add a distance monitor to the linked list, read the data into structure */
Mon_dist *add_distance(AtoM *a, AtoM *b)
{
   Mon_dist *temp, *mp;

   if((temp = (Mon_dist*) malloc(sizeof(Mon_dist))) == NULL){
      (void)fprintf(stderr, " can't allocate enough memory for distance\n");
      exit(-1);
   }

   temp->a = a;
   temp->b = b;
   temp->next = NULL;

   if(!actualmol->firstdist) actualmol->firstdist = temp;
   else { 
      for(mp=actualmol->firstdist; mp->next; mp=mp->next);
      mp->next = temp;
   }

   return temp;
}



Int_dist *add_intermolecular_distance(AtoM *a, Mol *ma, AtoM *b, Mol *mb)
{
   Int_dist *temp, *mp;

   if((temp = (Int_dist*) malloc(sizeof(Int_dist))) == NULL){
      (void)fprintf(stderr, " can't allocate enough memory for distance\n");
      exit(-1);
   }

   temp->a  = a;
   temp->ma = ma;
   temp->b  = b;
   temp->mb = mb;
   temp->next = NULL;

   if(!first_intermolecular_distance) first_intermolecular_distance = temp;
   else { 
      for(mp=first_intermolecular_distance; mp->next; mp=mp->next);
      mp->next = temp;
   }

   return temp;
}



/* add an angle monitor to the linked list, read the data into structure */
Mon_ang *add_angle(AtoM *a, AtoM *b, AtoM *c)
{
   Mon_ang *temp, *mp;

   if((temp = (Mon_ang*) malloc(sizeof(Mon_ang))) == NULL){
      (void)fprintf(stderr, " can't allocate enough memory for valence\n");
      exit(-1);
   }

   temp->a = a;
   temp->b = b;
   temp->c = c;
   temp->next = NULL;

   if(!actualmol->firstang) actualmol->firstang = temp;
   else {
      for(mp=actualmol->firstang; mp->next; mp=mp->next);
      mp->next = temp;
   }

   return temp;
}



/* add a torsion monitor to the linked list, read the data into structure */
Mon_tor *add_torsion(AtoM *a, AtoM *b, AtoM *c, AtoM *d)
{
   Mon_tor *temp, *mp;

   if((temp = (Mon_tor*) malloc(sizeof(Mon_tor))) == NULL){
      (void)fprintf(stderr, " can't allocate enough memory for torsion\n");
      exit(-1);
   }

   temp->a = a;
   temp->b = b;
   temp->c = c;
   temp->d = d;
   temp->next = NULL;

   if(!actualmol->firsttor) actualmol->firsttor = temp;
   else {
      for(mp=actualmol->firsttor; mp->next; mp=mp->next);
      mp->next = temp;
   }

   return temp;
}

Dipole *add_dipole(float x, float y, float z)
{
   Dipole *dp;
   char str[30];

   if((dp = (Dipole*) malloc(sizeof(Dipole))) == NULL) {
      fprintf(stderr, " can't allocate memory for the dipole moment\n");
      return NULL;
   }

   dp->start[0] = -actualmol->centervec[0] - x/2;
   dp->start[1] = -actualmol->centervec[1] - y/2;
   dp->start[2] = -actualmol->centervec[2] - z/2;
   dp->end[0] = -actualmol->centervec[0] + x/2;
   dp->end[1] = -actualmol->centervec[1] + y/2;
   dp->end[2] = -actualmol->centervec[2] + z/2;

   dp->absolute = dist(dp->start, dp->end);

/*
printf("Center Vec   : %f %f %f\n", actualmol->centervec[0], actualmol->centervec[1], actualmol->centervec[2]);
printf("Dipole Moment: %f %f %f\n", dp->start[0], dp->start[1], dp->start[2]);
printf("               %f %f %f\n", dp->end[0], dp->end[1], dp->end[2]);
printf("              |%f|\n", dist(dp->start, dp->end));
*/

   return dp;
}


void removeMolFromList(Mol *x)
{
   Mol *i;

   if(!firstmol) return;

   if(x == firstmol) firstmol = x->next;
   else {
      for(i=firstmol; i; i=i->next) {
         if(x == i->next) {
            i->next = x->next;
            break;
         }
      }
   }
   x->next = NULL;
}




void appendMolToList(Mol *x)
{
   Mol *i;

   if(!firstmol) firstmol = x;
   else { /* find the last molecule */
      for(i = firstmol; i->next; i = i->next);
      i->next = x;
   }
   x->next = NULL;
}



/* free one molecule from the linked list */
void free_mol(Mol *x)
{
   AtoM *ap;

   if(!x) return;

   if(x->firstamoss) free_amoss(x->firstamoss);
   for(ap=x->firstatom; ap; ap=ap->next) ap->firstshell = NULL;
   if(x->basisset) free_basisset(x->basisset);
   for(ap=x->firstatom; ap; ap=ap->next) {
      ap->firstshell = NULL;
      ap->firstslater = NULL;
   }
   free_mo(x);

   if(x->firstatom) free_atoms(x->firstatom);
   if(x->firstbond) free_bonds(x->firstbond);
   if(x->firstsurf) free_surfaces(x->firstsurf);
   if(x->firstter)  free_ters(x->firstter);
   if(x->firstresidue) free_residues(x->firstresidue);
   if(x->firstdist) free_distances(x->firstdist);
   if(x->firstang)  free_angles(x->firstang);
   if(x->firsttor)  free_torsions(x->firsttor);
   if(x->cube_value) free_3D(x);
   if(x->filename) free(x->filename);

   removeMolFromList(x);

   actualmol = firstmol;

   if ( actualmol ) { 
      if(actualmol->cube_value){
         cubemin = actualmol->cubemin;
         cubemax = actualmol->cubemax;
         cutoff = actualmol->cutoff = 0;
/* 
         if(actual_menu == &macu_menu ||
            actual_menu == &box_menu) drawfields();
 * ueberpruefen!!!
*/
      }
      actualsurf = actualmol->firstsurf;
   }
   if(x == dynamics.molecule) free_dyna();
   free(x);

   remIntermolecularDistances(x);
}


/* recursively free the linked list of surfaces */
void free_surfaces(Surface *x)
{
   if(x->next) free_surfaces(x->next);

   if(x->dot)  free(x->dot);
   if(x->tri)  free(x->tri);
   if(x->val)  free(x->val);
   if(x->name) free(x->name);
   if(x->trinorm) free(x->trinorm);

   free(x);
}



/* recursively free the linked list of amoss */
void free_amoss(Amoss_basis *x)
{
   if(x->next) free_amoss(x->next);
   if(x->firstshell) free_shell(x->firstshell);
   free(x);
}



/* recursively free the linked list of amoss */
void free_basisset(Basis *x)
{
   if(x->next) free_basisset(x->next);
   if(x->firstshell) free_shell(x->firstshell);
   if(x->firstslater) free_slater(x->firstslater);
   free(x);
}



/* recursively free the linked list of slater shells */
void free_slater(Slater *x)
{
   if(x->next) free_slater(x->next);
   free(x);
}



/* recursively free the linked list of atoms */
void free_atoms(AtoM *x)
{
   if(x->next) free_atoms(x->next);
   if(x->firstshell) free_shell(x->firstshell);
   if(x->valence) free(x->valence);
   if(x->firstslater) free_slater(x->firstslater);
   free(x);
}


void free_shell(Shell *shell)
{
   if(shell->next) free_shell(shell->next);
   if(shell->firstgauss) free_gauss(shell->firstgauss);
   free(shell);
}



void free_gauss(Gauss *gauss)
{
   if(gauss->next) free_gauss(gauss->next);
   free(gauss);
}



/* recursively free the linked list of bonds */
void free_bonds(Bond *x)
{
   if(x->next) free_bonds(x->next);
   free(x);
}



/* recursively free the linked list of ters */
void free_ters(Ter *x)
{
   if(x->next) free_ters(x->next);
   free(x);
}



/* recursively free the linked list of residues */
void free_residues(Residue *x)
{
   if(x->next) free_residues(x->next);
   free(x);
}


void free_distances(Mon_dist *x)
{
   if(x->next) free_distances(x->next);
   free(x);
}


void free_angles(Mon_ang *x)
{
   if(x->next) free_angles(x->next);
   free(x);
}


void free_torsions(Mon_tor *x)
{
   if(x->next) free_torsions(x->next);
   free(x);
}


void free__intermolecular_distances(Int_dist *x)
{
   if(x->next) free__intermolecular_distances(x->next);
   free(x);
}

void free_mo(Mol *mp)
{
   register int i;
   int norbs;

   norbs = mp->nMolecularOrbitals;

   if(mp->alphaOrbital) {
      for(i=0; i<norbs; i++){
         free(mp->alphaOrbital[i].coefficient);
      }
      free(mp->alphaOrbital);
      mp->alphaOrbital = NULL;
   }
   if(mp->betaOrbital) {
      for(i=0; i<norbs; i++){
         free(mp->betaOrbital[i].coefficient);
      }
      free(mp->betaOrbital);
      mp->betaOrbital = NULL;
   }
}


int remIntermolecularDistances(Mol *mp)
{
   Int_dist *i;
   int removed;

   removed = 0;
   for(i=first_intermolecular_distance; i; i=i->next) {
      if(i->ma == mp || i->mb == mp) {
         remIntermolecularDistance(i);
         removed = 1;
      }
   }

   return removed;
}


void remIntermolecularDistance(Int_dist *x)
{
   Int_dist *i;

   if(!first_intermolecular_distance) return;

   if(x == first_intermolecular_distance) {
      first_intermolecular_distance = x->next;
   }
   else {
      for(i=first_intermolecular_distance; i; i=i->next) {
         if(x == i->next) {
            i->next = x->next;
            break;
         }
      }
   }

   free(x);
}


/* return the distance between two points */
float dist(float *a, float *b)
{
   return sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) +
               (a[2]-b[2])*(a[2]-b[2]));
}


void new_mole(char *file)
{
/* taken from grafix.c */
   char str2[30], *str1, str3[30];

   reset();

   if((str1 = strrchr(file, '/')) == NULL) str1 = file;
   else str1++;

   (void)sprintf(str2, "%d atoms and %d bonds",
                     actualmol->natoms, actualmol->nbonds - actualmol->n_h_bonds);
   (void)sprintf(str3, "total mass = %.2f", actualmol->mass);
   logprint(str1);
   logprint(str2);
   logprint(str3);
}

void free_dyna(void)
{
   register int i;

   if(!dynamics.trajectory) return;
   for(i=0; i<dynamics.ntotalsteps; i++) {
      if(dynamics.trajectory[i]) free(dynamics.trajectory[i]);
   }
   if(dynamics.freeat) {
      for(i=0; i<dynamics.nfreat; i++) dynamics.freeat[i]->fixed = 1;
      free(dynamics.freeat);
   }
   dynamics.nfreat = 0;
   dynamics.freeat = NULL;
   dynamics.ntotalsteps = 0;
   dynamics.trajectory = NULL;
   dynamics.start = dynamics.end = dynamics.current = 0;

   restore_coords(dynamics.molecule);
   dynamics.molecule = NULL;
//   recalculate_bondmat(0);
//   drawit();
}

void restore_coords(Mol *mp)
{
   register AtoM *ap;
   register Vector *vw;

   if(!mp || !mp->firstatom || !coord_backup) return;

   for(ap = mp->firstatom, vw=coord_backup; ap; ap=ap->next, vw++){
      ap->coord[0] = vw->x; ap->coord[1] = vw->y; ap->coord[2] = vw->z;
   }
   free(coord_backup);
   coord_backup = NULL;
}

int save_coords(Mol *mp)
{
   register AtoM *ap;
   register Vector *vw;

   if(!coord_backup){
      if((coord_backup = (Vector *)malloc(mp->natoms * sizeof(Vector))) == NULL){
         showinfobox("Can't allocate vector-backup");
         return 0;
      }

      for(ap = mp->firstatom, vw=coord_backup; ap; ap=ap->next, vw++){
         vw->x = ap->coord[0]; vw->y = ap->coord[1]; vw->z = ap->coord[2];
      }
   }
   return 1;
}


void remove_values(void)
{
   int doRedraw;

   if(!actualmol) return;

   doRedraw = 0;

   if(actualmol->firstdist){
      free_distances(actualmol->firstdist);
      actualmol->firstdist = NULL;

      doRedraw = 1;
   }

   if(actualmol->firstang){
      free_angles(actualmol->firstang);
      actualmol->firstang = NULL;

      doRedraw = 1;
   }

   if(actualmol->firsttor){
      free_torsions(actualmol->firsttor);
      actualmol->firsttor = NULL;

      doRedraw = 1;
   }

   if(remIntermolecularDistances(actualmol)) doRedraw = 1;

   if(doRedraw) glutPostRedisplay();
}


void rem_bond_from_list(Bond *bp)
{
   register Bond *bd;

   if(bp == actualmol->firstbond)
      actualmol->firstbond = actualmol->firstbond->next;
   else {
      for(bd=actualmol->firstbond; bd->next!=bp; bd=bd->next);
         bd->next = bd->next->next;
         if(!bd->next) actualmol->lastbond = bd;
   }

   switch(bp->type){
      case TRIPLEBOND : bp->a->nbonds--;
                        bp->b->nbonds--;
      case DOUBLEBOND : bp->a->nbonds--;
                        bp->b->nbonds--;
      case SINGLEBOND : bp->a->nbonds--;
                        bp->b->nbonds--;
   }

   free(bp);
   actualmol->nbonds--;
   if(bp->type == H_BOND) actualmol->n_h_bonds--;
}

void modify_bond_type(short type, Bond *b)
{
   b->a->nbonds += type - b->type;
   b->b->nbonds += type - b->type;
   b->type = type;
}

MolecularOrbital *allocOrbital(int nOrbitals, int nBasis, int flag)
{
   MolecularOrbital *orbital;
   register int i;

   if((orbital = (MolecularOrbital*) malloc(nOrbitals*sizeof(MolecularOrbital))) == NULL){
      return 0;
   }
   for(i=0; i<nOrbitals; i++) {
      if(nBasis) {
         if((orbital[i].coefficient = (double*) malloc(nBasis*sizeof(double))) == NULL) {
            while(i--) free(orbital[i].coefficient);
            free(orbital);
            return 0;
         }
      }
      else orbital[i].coefficient = NULL;

      orbital[i].eigenvalue = 0;
      orbital[i].number = 0;
      orbital[i].occ = 0;
      orbital[i].type[0] = 0;
      orbital[i].flag = flag;
   }

   return orbital;
}

MolecularOrbital *reallocOrbital(MolecularOrbital *orbital, int newer, int old)
{
   MolecularOrbital *neworb;
   register int i;

   if(newer <= 0) return  NULL;

   if(newer == old) return orbital;

   if(newer < old) {
      for(i=newer; i<old; i++) {
         free(orbital[i].coefficient);
      }
      if((neworb = (MolecularOrbital*) realloc(orbital, newer*sizeof(MolecularOrbital))) == NULL){
         for(i=0; i<newer; i++) {
            free(orbital[i].coefficient);
         }
         return 0;
      }
      return neworb;
   }

   if(newer > old) {

      fprintf(stderr, "reallocation for more memory is not yet implemented\n");
      return NULL;

   }
   return NULL;
}


/* allocate dynamically a triangular matrix */
void *alloc_trimat(int n, size_t size)
{
   void **pointerarray;
   char *array;
   register short i;

   if((array = (char*) malloc((n*(n+1))/2*size)) == NULL) return NULL;
      /* array will hold the data */
   if((pointerarray = (void**) malloc(n*sizeof(char *))) == NULL) return NULL;
      /* pointerarray will hold the pointers to the rows of data */
   for(i=0; i<n; i++) pointerarray[i] = array + (i*(i+1))/2*size;

   return pointerarray;
}


Amoss_basis *add_amoss(Mol *mp)
/* add shell to linked list of atom's shells */
{
   Amoss_basis *amoss, *ap;

   if((amoss = (Amoss_basis *)malloc(sizeof(Amoss_basis))) == NULL){
      showinfobox("can't allocate AMOSS basis-set");
      return NULL;
   }

   if(mp->firstamoss == NULL) mp->firstamoss = amoss;
   else {
      for(ap = mp->firstamoss; ap->next; ap = ap->next);
      ap->next = amoss;
   }

   amoss->firstshell = NULL;
   amoss->next = NULL;

   return amoss;
}



Shell *add_shell_to_amoss(Amoss_basis *ap)
/* add shell to linked list of AMOSS basis-set shells */
{
   Shell *shell, *sp;

   if((shell = (Shell *)malloc(sizeof(Shell))) == NULL){
      showinfobox("can't allocate shell");
      return NULL;
   }

   if(ap->firstshell == NULL) ap->firstshell = shell;
   else {
      for(sp = ap->firstshell; sp->next; sp = sp->next);
      sp->next = shell;
   }

   shell->next = NULL;
   shell->firstgauss = NULL;

   return shell;
}

Gauss *add_gauss(Shell *sp)
/* add gaussian primitive to linked list of shell's contracted gaussians */
{
   Gauss *gauss;

   if((gauss = (Gauss*) malloc(sizeof(Gauss))) == NULL){
      showinfobox("can't allocate gaussian primitive function");
      return NULL;
   }

   if(sp->firstgauss == NULL) sp->firstgauss = gauss;
   else sp->actualgauss->next = gauss;

   sp->actualgauss = gauss;
   gauss->coeff = 0;
   gauss->coeff2 = 0;
   gauss->next = NULL;

   return gauss;
}


Shell *add_shell(AtoM *ap)
/* add shell to linked list of atom's shells */
{
   Shell *shell;

   if((shell = (Shell*) malloc(sizeof(Shell))) == NULL){
      showinfobox("can't allocate shell");
      return NULL;
   }

   if(ap->firstshell == NULL) ap->firstshell = shell;
   else ap->actualshell->next = shell;

   ap->actualshell = shell;
   shell->next = NULL;
   shell->firstgauss = NULL;

   return shell;
}


Slater *add_slater(AtoM *ap)
{
   Slater *temp, *sp;

   if((temp = (Slater *)malloc(sizeof(Slater))) == NULL){
      showinfobox("can't allocate slater");
      return NULL;
   }

   if(ap->firstslater == NULL) ap->firstslater = temp;
   else {
      for(sp=ap->firstslater; sp->next; sp=sp->next);
      sp->next = temp;
   }

   temp->n = 0;
   temp->type[0] = ' '; temp->type[1] = 0;
   temp->exponent = 0;
   temp->norm[0] = temp->norm[1] = temp->norm[2] = 1;
   temp->next = NULL;

   return temp;
}


void help(void)
{
#ifndef WIN32
   char command[500];

   sprintf(command, "if { which netscape ;}\nthen\nexit 1\nfi\n");
   if (!system(command)) {
      showinfobox("help:\nCan't find netscape!");
      return;
   }
   sprintf(command, "if { netscape -remote \"openURL(http://www.cscs.ch/molekel/manual_frame.html)\" 2>>/dev/null ;}\n");
   sprintf(command, "%sthen\n   echo\nelse\n   ", command);
   sprintf(command, "%snetscape http://www.cscs.ch/molekel/manual_frame.html & \nfi", command);
   system(command);
#else
   showinfobox("help:\nNot implemented yet!");
#endif
}
