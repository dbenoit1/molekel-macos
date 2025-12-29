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


#ifndef WIN32
#include <signal.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <sys/param.h>
#include <sys/stat.h>
#endif
#include "main.h"
#include "molekel.h"
#include "constant.h"
#include "general.h"
#include "calcdens.h"
#include "box.h"
#include "maininterf.h"
#include "chooseinterf.h"
#include "browser.h"
#include "readgauss.h"
#include "macu.h"
#include "drawing.h"
#include "menu.h"
#include "texture.h"
#include "connolly.h"
#include "manip.h"
#include "readadf.h"

#ifdef LINUX
#define NOSOCK
#endif

#ifdef WIN32
#define NOSOCK
#endif

#define POW_5(x,y) ((y)>0?(x):1)
#define POW_4(x,y) ((y)>0?(x)*POW_5((x),(y)-1):1)
#define POW_3(x,y) ((y)>0?(x)*POW_4((x),(y)-1):1)
#define POW_2(x,y) ((y)>0?(x)*POW_3((x),(y)-1):1)
#define POW_1(x,y) ((y)>0?(x)*POW_2((x),(y)-1):1)
#define POW(x,y)   ((y)>0?(x)*POW_1((x),(y)-1):1)

void sighandler(int sig);

int did_sigset = 0, pipein = 0, pipeout = 0;
int child_action_key = 0;
float **density;
double *chi;
MolecularOrbital *molOrb;
char timestring[300];

void calcdens(int key)
{
   float dim[6];
   int ncub[3], child_id;
   int pipeline[2];
   char command[500];

   child_action_key = action_key;
   define_box(dim, ncub);

   green_note("calculating");

#ifndef NOSOCK
   if(!pipein){
      if(pipe(pipeline) < 0){
         fprintf(stderr, "can't open pipe...\n");
         perror("opening communication pipe...");
         return;
      }
      pipein = pipeline[0];
      pipeout = pipeline[1];
   }

   if((child_id = fork()) == -1){
      fprintf(stderr, "can't create child...\n");
      return;
   }
   if(child_id){                   /* parent */
#ifdef LINUX
         signal(SIGCHLD, sighandler);
#else
      if(!did_sigset){
         sigset(SIGCHLD, sighandler);
         did_sigset = 1;
      }
#endif
   }
   else {                          /* child */
      close(pipein);
      setpriority(PRIO_PROCESS, 0, job_priority_live_var);
      sprintf(command, "%s", filename);
      process_calc(command, dim, ncub, key);
      sprintf(command, "%s\n%s\nMolecule Id %p",
		    filename, timestring, actualmol);
      if(write(pipeout, command, strlen(command) + 1) < 0)
               perror("writing message...");
      close(pipeout);
      exit(0);
   }
#else

   sprintf(command, "%s", filename);
   process_calc(command, dim, ncub, key);

   if (child_action_key != MEP) {
      defaultSurfaces(filename);
   }
   else {
      logprint("");
      logprint("MEP: 3D-grid generation done");
      update_logs();
   }
  
#endif
}
 
#ifndef NOSOCK
void sighandler(int sig)
{
   char message[512], strcopy[512];
   char *str[10];
   int  ntokens;
   struct stat buf;

   putc('\7', stdout);
   fflush(stdout);

   if(read(pipein, message, 511) < 0) perror("reading message..");
   strcpy(strcopy, message);
   showinfobox(strcopy);
   if(message[0]) {
     str[0] = strtok(message, "\n");
     for(ntokens = 1; ntokens<10;) {
        if((str[ntokens] = strtok(NULL, "\n")) == NULL) break;
        ntokens++;
     }
   }

    if(strstr(str[0], ".macu")) {
	int i;

	for(i=1; i<ntokens; i++) {
	    char actual[32];
	    if(strstr(str[i], "Molecule")) {
		sprintf(actual, "%p", actualmol);
		if(strstr(str[i], actual)) {
                    if (child_action_key != MEP) {
                       if(stat(str[0], &buf) == 0) defaultSurfaces(str[0]);
                    }
                    else {
                       logprint("");
                       logprint("MEP: 3D-grid generation done");
                       update_logs();
                    }
		}
		
	    }
	}
    }
    else if(child_action_key == CONNOLLY) {
       default_read_dots();
    }
}
#endif

void process_calc(char *s, float *dim, int *ncubes, int key)
{
   float x, y, z, dx, dy, dz;
   float **slice, *array;
   short i, j, k;
   int len, ncub[3];
   FILE *fp;
#ifndef NOSOCK
   struct tms starttime, endtime;
   float systime, cputime;
#endif
   double (*funct)(float x, float y, float z);

   ncub[0] = *ncubes++;
   ncub[1] = *ncubes++;
   ncub[2] = *ncubes++;

    if((actualmol->alphaOrbital[0].flag == ADF_ORB_A ||
        actualmol->alphaOrbital[0].flag == ADF_ORB_B ) && 
        actualmol->alphaOrbital[0].coefficient == NULL) {
	   executeAdfUtilities(s, dim, ncub, key);
   	   return;
    }

   if((fp = fopen(s, "wb")) == NULL){
      fprintf(stderr, "can't open .macu file\n");
      return;
   }

   if((array = (float *)malloc(ncub[0]*ncub[1]*sizeof(float))) == NULL){
      fprintf(stderr, "can't allocate slice\n");
      exit(-1);
   }
  if((slice = (float **)malloc(ncub[1]*sizeof(float *))) == NULL){
      fprintf(stderr, "can't allocate slice\n");
      exit(-1);
   }
   for(i=0; i<ncub[1]; i++) slice[i] = array + (i * ncub[0]);

    chi = (double *)calloc(actualmol->nBasisFunctions, sizeof(double));
    if(!chi && !(actualmol->alphaOrbital[0].flag == EHT_ORB)){
	fprintf(stderr, "can't allocate chi\n");
	return;
    }


    switch(key){
	case CALC_ORB  : 
	    switch(actualmol->alphaOrbital[0].flag) {
		case GAMESS_ORB :
		case HONDO_ORB  :
		case GAUSS_ORB  : funct = calc_point; break;
/* to be fixed
                case ADF_ORB_A  :
		case ADF_ORB_B  : funct = calc_adf_point; break;

		case EHT_ORB    : funct = calc_eht_point; break;
*/
		case MOS_ORB    :
		case ZINDO_ORB  :
		case PRDDO_ORB  : funct = calc_prddo_point; break;
		case MLD_SLATER_ORB  : funct = calc_sltr_point; break;
	    }
	    break;

	case EL_DENS :
	    switch(actualmol->alphaOrbital[0].flag) {
		case GAMESS_ORB :
		case HONDO_ORB  :
		case GAUSS_ORB  : 
		    if(!generate_density_matrix(key)) {
			fprintf(stderr, "Can't generate the density matrix!\n");
			strcpy(timestring, "Can't generate the density matrix!");
			return;
		    }
		    else {
			printf("density matrix generated...\n");
			funct = calculate_density;
		    }
		    break;

/* to be fixed
                case ADF_ORB_A  :
		case ADF_ORB_B  : funct = calc_adf_density; break;

		case EHT_ORB    : funct = calc_eht_density; break;
*/
		case MOS_ORB    :
		case ZINDO_ORB  :
		case PRDDO_ORB  : funct = calc_prddo_density; break;
		case MLD_SLATER_ORB  : funct = calc_sltr_density; break;
	    }
	    break;

	case SPIN_DENS :
	    switch(actualmol->alphaOrbital[0].flag) {
		case GAMESS_ORB :
		case HONDO_ORB  :
		case GAUSS_ORB  : 
		    if(!actualmol->alphaBeta) {
			funct = calculateSomo;
		    }
		    else if(!generate_density_matrix(key)) {
			fprintf(stderr, "Can't generate the density matrix!\n");
			strcpy(timestring, "Can't generate the density matrix!");
			return;
		    }
		    else {
			printf("density matrix generated...\n");
			funct = calculate_density;
		    }
		    break;

/* to be fixed
		case ADF_ORB_A  :
                case ADF_ORB_B  :
*/
		case MOS_ORB    :
		case ZINDO_ORB  :
		case PRDDO_ORB  : funct = calc_prddo_spindensity; break;
		case MLD_SLATER_ORB  : funct = calc_sltr_spindensity; break;
	    }
	    break;

	case MEP       :
	    funct = calc_mep; break;
   }

   len = 36;
   fwrite(&len, sizeof(int), 1, fp);
   fwrite(dim, sizeof(float), 6, fp);
   fwrite(ncub, sizeof(int), 3, fp);
   fwrite(&len, sizeof(int), 1, fp);

   dx = (dim[1]-dim[0])/(ncub[0]-1);
   dy = (dim[3]-dim[2])/(ncub[1]-1);
   dz = (dim[5]-dim[4])/(ncub[2]-1);
   len = 4*ncub[0]*ncub[1];
#ifndef NOSOCK
   times(&starttime);
#endif

   for(i=0, z=dim[4]; i<ncub[2]; i++, z += dz){
      for(j=0, y=dim[2]; j<ncub[1]; j++, y += dy){
         for(k=0, x=dim[0]; k<ncub[0]; k++, x += dx){
            slice[j][k] = (*funct)(x, y, z);
         }
      }
      fwrite(&len, sizeof(int), 1, fp);
      fwrite(slice[0], sizeof(float), ncub[0]*ncub[1], fp);
      fwrite(&len, sizeof(int), 1, fp);
      fflush(fp);
      printf(".");
      fflush(stdout);
   }
   printf("\n");

#ifndef NOSOCK
   putc('\7', stdout);            /* the machine that goes "beep" */
   fflush(stdout);

   times(&endtime);
   cputime = (endtime.tms_utime - starttime.tms_utime)/(float)HZ;
   systime = (endtime.tms_stime - starttime.tms_stime)/(float)HZ;
   printf("  time : cpu %.2fs, sys %.2fs\n", cputime, systime);
   sprintf(timestring, " : cpu %.2fs, sys %.2fs\0", cputime, systime);
#endif

   fclose(fp);
   free(slice);
   free(array);
   if(density){
      free(density[0]);
      free(density);
   }
   density = NULL;
   if(chi) free(chi);
}




double calc_point(float x, float y, float z)
/* calculate the MO-value at given point */
/* no functions for speed */
{
    register int i;
    double value, *ao_coeff;

    ao_coeff = molOrb->coefficient;
    calc_chi(x, y, z);

    value = 0;
    for(i=0; i<actualmol->nBasisFunctions; i++) {
	value += ao_coeff[i]*chi[i];
    }

    return value;
}






double calculateSomo(float x, float y, float z)
/* calculate the spin density of the singly occ. orbitals */
{
   register short i;
   double value, point;

   value = 0;

   for(i=actualmol->nBeta; i<actualmol->nAlpha; i++){
      molOrb = actualmol->alphaOrbital + i - actualmol->firstOrbital + 1;
      point = calc_point(x, y, z);
      value += point*point;;
   }

   return value;
}






double calculate_density(float x, float y, float z)
/* calculate the electron or spin density at given point */
{
   register short i, j;
   double value;

   value = 0;
   calc_chi(x, y, z);

   for(i=0; i<actualmol->nBasisFunctions; i++){
      value += density[i][i] * chi[i] * chi[i];
      for(j=0; j<i; j++)
         value += density[i][j] * chi[i] * chi[j] * 2.0;
   }

   return value;
}




void calc_chi(float x, float y, float z)
/* calculate sum of AO-contributions chi for each MO at given point */
/* no functions for speed */
{
   AtoM *ap;
   register Shell     *sp;
   register Gauss     *gp;
   double radial_part, *cp;
   float xa, ya, za, ra2;  /* atomic units !! */

   memset(chi, 0, actualmol->nBasisFunctions * sizeof(double));
   cp = chi;
   for(ap = actualmol->firstatom; ap; ap = ap->next){
      xa = (x - ap->coord[0]) * _1_BOHR;
      ya = (y - ap->coord[1]) * _1_BOHR;
      za = (z - ap->coord[2]) * _1_BOHR;
               
      ra2 = xa*xa + ya*ya + za*za;      /* cutoff-distance ? */

      for(sp = ap->firstshell; sp; sp = sp->next){
         switch(sp->n_base){
            case 1  :           /*** S-orbital ***/
               for(gp = sp->firstgauss; gp; gp = gp->next){
                  radial_part = exp(-ra2*gp->exponent);
                  *cp += gp->coeff * radial_part;
               }
               cp++;
               break;

            case 4 :           /*** SP-orbital ***/
               for(gp = sp->firstgauss; gp; gp = gp->next){
                  radial_part = exp(-ra2*gp->exponent);
                  *cp     += gp->coeff * radial_part;
                  *(cp+1) += gp->coeff2 * xa * radial_part;
                  *(cp+2) += gp->coeff2 * ya * radial_part;
                  *(cp+3) += gp->coeff2 * za * radial_part;
               }
               cp += 4;
               break;

            case 3  :           /*** P-orbital ***/
               for(gp = sp->firstgauss; gp; gp = gp->next){
                  radial_part = gp->coeff * exp(-ra2*gp->exponent);
                  *cp     += xa * radial_part;
                  *(cp+1) += ya * radial_part;
                  *(cp+2) += za * radial_part;
               }
               cp += 3;
               break;

            case 5  :           /*** D-orbital (5) ***/
               for(gp = sp->firstgauss; gp; gp = gp->next){
                  radial_part = gp->coeff * exp(-ra2*gp->exponent);
                  *cp     += 0.288675135 *
                            (2*za*za - xa*xa - ya*ya) * radial_part;
                  *(cp+3) += 0.5 * (xa*xa - ya*ya) * radial_part;
                  *(cp+4) += xa * ya * radial_part;
                  *(cp+1) += xa * za * radial_part;
                  *(cp+2) += ya * za * radial_part;
               }
               cp += 5;
               break;

            case 6  :           /*** D-orbital (6) ***/
               for(gp = sp->firstgauss; gp; gp = gp->next){
                  radial_part = gp->coeff * exp(-ra2*gp->exponent);
                  *cp     += radial_part * xa * xa * 0.57735027;
                  *(cp+1) += radial_part * ya * ya * 0.57735027;
                  *(cp+2) += radial_part * za * za * 0.57735027;
                  *(cp+3) += radial_part * xa * ya;
                  *(cp+4) += radial_part * xa * za;
                  *(cp+5) += radial_part * ya * za;
               }
               cp += 6;
               break;

            case 7  :           /*** F-orbital (7) ***/
               for(gp = sp->firstgauss; gp; gp = gp->next){
                  radial_part = exp(-ra2*gp->exponent) * gp->coeff;
                  *cp     += radial_part * za * (5. * za * za - 3. * ra2)/* * k */;
                  *(cp+1) += radial_part * xa * (5. * za * za - ra2)/* * k */;
                  *(cp+2) += radial_part * ya * (5. * za * za - ra2)/* * k */;
                  *(cp+3) += radial_part * za * (xa * xa - ya * ya)/* * k */;
                  *(cp+4) += radial_part * xa * ya * za;
                  *(cp+5) += radial_part * (xa * xa * xa - 3. * xa * ya * ya)/* * k */;
                  *(cp+6) += radial_part * (3. * xa * xa * ya - ya * ya * ya)/* * k */;
               }
               cp += 7;
               break;

            case 10 :           /*** F-orbital (10) ***/
                      /* correct order ??? */
               for(gp = sp->firstgauss; gp; gp = gp->next){
                  radial_part = gp->coeff * exp(-ra2*gp->exponent);
                  *cp     += radial_part * xa * xa * xa * .25819889;
                  *(cp+1) += radial_part * ya * ya * ya * .25819889;
                  *(cp+2) += radial_part * za * za * za * .25819889;
                  *(cp+3) += radial_part * xa * xa * ya * .57735027;
                  *(cp+4) += radial_part * xa * xa * za * .57735027;
                  *(cp+5) += radial_part * xa * ya * ya * .57735027;
                  *(cp+6) += radial_part * ya * ya * za * .57735027;
                  *(cp+7) += radial_part * xa * za * za * .57735027;
                  *(cp+8) += radial_part * ya * za * za * .57735027;
                  *(cp+9) += radial_part * xa * ya * za;
               }
               cp += 10;
               break;

         } /* end of switch */
      } /* end of loop over the shells (for(sp...) */
   } /* end of loop over the atoms (for(ap...)*/

   return;
}




int generate_density_matrix(int key)
{
   register short i, j, k;
   float adder;

   if(!actualmol->alphaDensity && datasource == USE_MATRICES) return 0;

   if(density){
      free(density[0]);
      free(density);
      density = NULL;
   }

   if((density = (float **)alloc_trimat(actualmol->nBasisFunctions, sizeof(float))) == NULL){
      fprintf(stderr, "can't allocate density-matrix\n");
      return 0;
   }

/* use the alpha (and beta) density matrices from the gaussian output file */
   if(datasource == USE_MATRICES){
      if(key == EL_DENS){ 
         if(actualmol->betaDensity){
            for(i=0; i<actualmol->nBasisFunctions; i++){
               for(j=0; j<=i; j++)
                  density[i][j] = actualmol->alphaDensity[i][j] + actualmol->betaDensity[i][j];
            }
         }
         else {
            for(i=0; i<actualmol->nBasisFunctions; i++){
               for(j=0; j<=i; j++) density[i][j] = actualmol->alphaDensity[i][j];
            }
         }
      }
      if(key == SPIN_DENS){
         if(!actualmol->betaDensity) return 0; 
         for(i=0; i<actualmol->nBasisFunctions; i++){
            for(j=0; j<=i; j++)
               density[i][j]   = actualmol->alphaDensity[i][j] - actualmol->betaDensity[i][j];
         }
      }
   }

/* generate the density matrices with the coefficients */
   else if(datasource == USE_COEFFS) {
      if(key == EL_DENS){ 
         for(i=0; i<actualmol->nBasisFunctions; i++){
            for(j=0; j<=i; j++){
               density[i][j] = 0;
               for(k=0; k<actualmol->nBeta; k++){
                  adder = actualmol->alphaOrbital[k].coefficient[i] *
                          actualmol->alphaOrbital[k].coefficient[j];
                  if(actualmol->alphaBeta)
                     adder += actualmol->betaOrbital[k].coefficient[i] *
                              actualmol->betaOrbital[k].coefficient[j];
                  else adder *= 2.0;
                  density[i][j] += adder;
               }
               for(; k<actualmol->nAlpha; k++)
                  density[i][j] += actualmol->alphaOrbital[k].coefficient[i] *
                                   actualmol->alphaOrbital[k].coefficient[j];
            }
         }
      }
      if(key == SPIN_DENS){
         if(!actualmol->alphaBeta) {
            for(i=0; i<actualmol->nBasisFunctions; i++){
               for(j=0; j<=i; j++){
                  density[i][j] = 0;
                  for(k=actualmol->nBeta; k<actualmol->nAlpha; k++)
                     density[i][j] += actualmol->alphaOrbital[k].coefficient[i] *
                                      actualmol->alphaOrbital[k].coefficient[j];
               }
            }
         }
         else {  
            for(i=0; i<actualmol->nBasisFunctions; i++){
               for(j=0; j<=i; j++){
                  density[i][j] = 0;
                  for(k=0; k<actualmol->nAlpha; k++)
                     density[i][j] += actualmol->alphaOrbital[k].coefficient[i] *
                                      actualmol->alphaOrbital[k].coefficient[j];
                  for(k=0; k<actualmol->nBeta; k++)
                     density[i][j] -= actualmol->betaOrbital[k].coefficient[i] *
                                      actualmol->betaOrbital[k].coefficient[j];
               }
            }
         }
      }
   }
   else return 0;

   return 1;
}

void mep_dot_surface(Surface *surf)
{

   if(!surf->val){
      if((surf->val = (float *)malloc(surf->npts * sizeof(float))) == NULL){
         showinfobox("can't allocate memory for the quality");
         return;
      }
   }
   {
      Surfdot *sp;
      float *vp;
      register int i;

      sp = surf->dot;
      vp = surf->val;
      for(i=0; i<surf->npts; i++, sp++, vp++){
          *vp = calc_mep(sp->v[0], sp->v[1], sp->v[2]);
      }

      surf->vmin = surf->vmax = *surf->val;
      for(i=1, vp=surf->val+1; i<surf->npts; i++, vp++){
         if(*vp < surf->vmin) surf->vmin = *vp;
         if(*vp > surf->vmax) surf->vmax = *vp;
      }
      if(!firsttexture) init_texture();

      surf->texture = firsttexture;
      surf->texenv = GL_MODULATE;
      surf->textype = TEX_MAP;
      set_vmin_vmax(surf);
   }

}


void mep_dot_surface_no_pick(void)
{
   Surface *surf;
   int count;
   char str[100];

   surf = actualmol->firstsurf;
   for(count=1; count < surface_live_var; count++) surf = surf->next;

   if(!surf->val){
      if((surf->val = (float *)malloc(surf->npts * sizeof(float))) == NULL){
         showinfobox("can't allocate memory for the quality");
         return;
      }
   }
   {
      Surfdot *sp;
      float *vp;
      register int i;

      sp = surf->dot;
      vp = surf->val;
      for(i=0; i<surf->npts; i++, sp++, vp++){
          *vp = calc_mep(sp->v[0], sp->v[1], sp->v[2]);
      }

      surf->vmin = surf->vmax = *surf->val;
      for(i=1, vp=surf->val+1; i<surf->npts; i++, vp++){
         if(*vp < surf->vmin) surf->vmin = *vp;
         if(*vp > surf->vmax) surf->vmax = *vp;
      }
      if(!firsttexture) init_texture();

      surf->texture = firsttexture;
      surf->texenv = GL_MODULATE;
      surf->textype = TEX_MAP;
   }

   glutPostRedisplay();
   sprintf(str, "vmin = %f   vmax = %f", surf->vmin, surf->vmax);
   logprint("");
   logprint("MEP:");
   logprint(str);
   update_logs();
}


double calc_mep(float x, float y, float z)
/* calculate the MEP at given point based on the point charges 
 * of the atoms
 */
{
   AtoM *ap;
   double nuclear;
   double xa, ya, za, ra2;

   nuclear = 0;
   for(ap = actualmol->firstatom; ap; ap = ap->next){
      xa = (x - ap->coord[0]) * _1_BOHR; /* atomic units */
      ya = (y - ap->coord[1]) * _1_BOHR;
      za = (z - ap->coord[2]) * _1_BOHR;
      ra2 = xa*xa + ya*ya + za*za;
      if(!ra2) {
         if(ap->charge > 0) return 10000;
         else if(ap->charge < 0) return -10000;
      }
      else     nuclear += ap->charge/sqrt(ra2);
   }

   return nuclear;
}

void spin_on_atoms(void)
{
   register AtoM *ap;
   register short i;
   char str[40];

   if(!generate_density_matrix(SPIN_DENS)){
      showinfobox("Can't generate the density matrix!");
      return;
   }

   chi = (double *)calloc(actualmol->nBasisFunctions, sizeof(double));
   if(!chi) showinfobox("Can't allocate chi");

   printf("Spin-densities on the atoms :\n");
   for(ap=actualmol->firstatom, i=0; ap; ap=ap->next, i++){
         ap->spin = calculate_density(ap->coord[0], ap->coord[1], ap->coord[2]);
         sprintf(str, "    %2s%-3d   %12.6f", element[ap->ord].symbol, i+1,
                 ap->spin);
         printf("%s\n", str);
   }
   actualmol->atm_spin = 1;
   update_interface_flags();

   free(chi);
}


double calc_prddo_point(float x, float y, float z)
/* calculate the MO-value at given point */
{
   AtoM *ap;
   Slater *vp;
   double value, *ao_coeff, angular_part;
   float xa, ya, za, ra2, ra;  /* atomic units !! */

   value = 0;
   ao_coeff = molOrb->coefficient;

   for(ap = actualmol->firstatom; ap; ap = ap->next){
      xa = (x - ap->coord[0]) * _1_BOHR;
      ya = (y - ap->coord[1]) * _1_BOHR;
      za = (z - ap->coord[2]) * _1_BOHR;

      ra2 = xa*xa + ya*ya + za*za;
      ra = sqrt(ra2);

      for(vp = ap->firstslater; vp; vp=vp->next){

         switch(vp->type[0]) {
            case 'S' :
               if(*ao_coeff) value += *ao_coeff * POW(ra, vp->n - 1) *
                             exp(-vp->exponent * ra) * vp->norm[0];
               ao_coeff++;
               break;
            case 'P' :
               angular_part = ao_coeff[0] * za + ao_coeff[1] * xa +
                              ao_coeff[2] * ya;
               value += angular_part * POW(ra, vp->n - 2) * 
                        exp(-vp->exponent * ra) * vp->norm[1];
               ao_coeff += 3;
               break;
            case 'D' :
               angular_part =  ao_coeff[0] * (3.*za*za - ra2) * vp->norm[3] +
                               ao_coeff[2] * (xa*xa - ya*ya) * vp->norm[2] +
                              (ao_coeff[1] * xa*za +
                               ao_coeff[3] * ya*za +
                               ao_coeff[4] * xa*ya) * vp->norm[4];

               value += angular_part * POW(ra, vp->n - 3) *
                        exp(-vp->exponent * ra);
               ao_coeff += 5;
               break;
         }
      }
   } /* end of loop over the atoms (for(ap...)*/

   return value;
}




double calc_prddo_density(float x, float y, float z)
{
   int i, norbs;
   double value, contr;

   value = 0;
   for(i=0; i<actualmol->nMolecularOrbitals; i++){
      if(actualmol->alphaOrbital[i].occ > 0){
         molOrb = actualmol->alphaOrbital + i;
         contr = calc_prddo_point(x, y, z);
         value += actualmol->alphaOrbital[i].occ*contr*contr;
      }
      if(actualmol->betaOrbital && actualmol->betaOrbital[i].occ > 0){
         molOrb = actualmol->betaOrbital + i;
         contr = calc_prddo_point(x, y, z);
         value += actualmol->betaOrbital[i].occ*contr*contr;
      }
   }

   return value;
}




double calc_prddo_spindensity(float x, float y, float z)
{
   int i, norbs;
   double value, contr;

   value = 0;
   if(actualmol->betaOrbital) {
      for(i=0; i<actualmol->nMolecularOrbitals; i++){
         if(actualmol->alphaOrbital[i].occ > 0){
            molOrb = actualmol->alphaOrbital + i;
            contr = calc_prddo_point(x, y, z);
            value += actualmol->alphaOrbital[i].occ*contr*contr;
         }
         if(actualmol->betaOrbital[i].occ > 0){
            molOrb = actualmol->betaOrbital + i;
            contr = calc_prddo_point(x, y, z);
            value -= actualmol->betaOrbital[i].occ*contr*contr;
         }
      }
   }
   else if (actualmol->nAlpha != actualmol->nBeta) {
      for(i=0; i<actualmol->nMolecularOrbitals; i++){
         if(actualmol->alphaOrbital[i].occ == 1){
            molOrb = actualmol->alphaOrbital + i;
            contr = calc_prddo_point(x, y, z);
            value += contr*contr;
         }
      }
   }

   return value;
}


double calc_sltr_point(float x, float y, float z)
{
   AtoM *ap;
   Slater *slp;
   float xa, ya, za, ra2, ra;
   double value = 0, exponent = 0, *ao_coeff;

   ao_coeff = molOrb->coefficient;

   for(ap = actualmol->firstatom; ap; ap = ap->next){
      xa = (x - ap->coord[0]) * _1_BOHR;
      ya = (y - ap->coord[1]) * _1_BOHR;
      za = (z - ap->coord[2]) * _1_BOHR;
      ra2 = xa*xa + ya*ya + za*za;
      ra = sqrt(ra2);

      for(slp = ap->firstslater; slp; slp=slp->next){
         if(*ao_coeff) {
            exponent = exp(-(slp->exponent*ra));
            value += *ao_coeff * slp->norm[0] * POW(xa,slp->a) * 
              POW(ya,slp->b) * POW(za,slp->c) * POW(ra,slp->d) * exponent;
         }
         ao_coeff++;
      }

   }
   return value;
}


double calc_sltr_density(float x, float y, float z)
{
   int i, norbs;
   double value, contr;

   value = 0;
   for(i=0; i<actualmol->nMolecularOrbitals; i++){
      if(actualmol->alphaOrbital[i].occ > 0){
         molOrb = actualmol->alphaOrbital + i;
         contr = calc_sltr_point(x, y, z);
         value += actualmol->alphaOrbital[i].occ*contr*contr;
      }
      if(actualmol->betaOrbital && actualmol->betaOrbital[i].occ > 0){
         molOrb = actualmol->betaOrbital + i;
         contr = calc_sltr_point(x, y, z);
         value += actualmol->betaOrbital[i].occ*contr*contr;
      }
   }

   return value;
}

double calc_sltr_spindensity(float x, float y, float z)
{
   int i, norbs;
   double value, contr;

   value = 0;
   if(actualmol->betaOrbital) {
      for(i=0; i<actualmol->nMolecularOrbitals; i++){
         if(actualmol->alphaOrbital[i].occ > 0){
            molOrb = actualmol->alphaOrbital + i;
            contr = calc_sltr_point(x, y, z);
            value += actualmol->alphaOrbital[i].occ*contr*contr;
         }
         if(actualmol->betaOrbital[i].occ > 0){
            molOrb = actualmol->betaOrbital + i;
            contr = calc_sltr_point(x, y, z);
            value -= actualmol->betaOrbital[i].occ*contr*contr;
         }
      }
   }
   else if (actualmol->nAlpha != actualmol->nBeta) {
      for(i=0; i<actualmol->nMolecularOrbitals; i++){
         if(actualmol->alphaOrbital[i].occ == 1){
            molOrb = actualmol->alphaOrbital + i;
            contr = calc_sltr_point(x, y, z);
            value += contr*contr;
         }
      }
   }

   return value;
}
