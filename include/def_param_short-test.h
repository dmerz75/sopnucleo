#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "structure.h"

#define stepsim1 /*1000000*/  100000 /* 100 */ /* 1000 number of simulation steps for heating under zero applied force */
#define stepsim2 /*100000000*/ 1000000 /* 1000 */ /* 10000 LJY  can change criteria in main_*.c */
#define nav      /*100*/  1000   /*10000*/    /* frequency to get information about chain conformation */
#define nav1     /*nav*/  (10*nav)  /*(100*nav)*/  /* frequency to save chain conformation: 2nd for nav=1000;1st and 3rd for nav=10000 */
#define nav_dcd (100*nav) /* frequency to write the chain conformation in dcd format */ 
#define temp     0.6     /* temperature (in units of (eh/kB)); 0.6 kcal/mol = 300 K */
#define zeta     50.0    /* (zeta*(1/tau_L)*(eh/temp) sets the natural time for the overdamped motion */
#define h        /*0.16*/ 0.08 /*0.02*/   /* integration step */
#define R_limit  2.0  /*8.0*/     /* this is the R0 parameter in FENE model */
#define R_limit_bond 8.0 /*9.0*/ /* cut-off distance for 2 beads to form a native non-bonded contact */
#define R_sigma  16.0    /* square of # of stdv tolerated in distance between positions that are taken into account in the LJ pot. */
#define kspring_cov /*40.0*/ 20.0 /* spring constant in FENE model (for the covalently bound beads) */
#define k_trans  /*0.1*/ 0.05 /*0.05*/ /*5.00*/ /*1.50*/     /* spring constant for externally applied force */
#define a        3.8     /* average distance between two consecutive beads along the chain */
#define eps      1.E-6
#define deltax   /*1.00*/ /*0.01*/ /*0.001*/ 0.0005 /*0.0001*/  /* increment in displacement of bead under force (to use in const. loading rate conditions, in A) */
#define eh       /*0.9*/ /*1.5*/ /* 2.0 */ 1.25 /*1.0*/ /*0.8*/ /* energy scale for attractive interactions between beads, i.e., in LJ pot. attract. */
#define el       1.0     /* energy scale for purely repulsive LJ pot. */
#define INF      100000.0
#define scale    1.0     /* factor to scale coordinates of beads in Rasmol */
#define LJ_bond  /*0.2*/ 0.0     /* factor to scale the repulsive LJ pot. for the covalently bound beads */
#define scale_R  1.1 /*0.9*/ /*0.7*/ /*0.8*/ /* scale factor to relate the end-to-end distance with the contour-length=(tot_amino*a) */
#define half     /*2.0*/ 1.0
#define sforce   /*2.0*/ 1.5  /*1.0*/ /*0.50*/ /* value of input force for const. force simulations (in kcal/mol A^2: 1 = 70 pN) */

#define min_force 0.10 /*0.20*/   /* value of the min force to which we reduce (linearly) the initial large force that opened the chain */
#define  bead_select 663 /* position along the chain on which we apply the force */


#define fixed_bead   1  /* position fixed */
#define pulled_bead  663  /* position pulled */

/* LJY20100318 begin calculate between two molecules, i.e., a protein and the nucleotide (ATP or ADP) bound to it */
/* parameters for nucleotide */
#define R_limits      2.0  /*8.0*/     /* R0 parameter in FENE model */
#define kspring_covs /*40.0*/ 20.0 /* spring constant in FENE model (for the covalently bound beads) */
#define as            3.8     /* average distance between two consecutive beads along the chain */
#define ehs          /*0.9*/ /*1.5*/ 2.0 /*1.0*/ /*0.8*/ /* energy scale for attract LJ pot */
#define elss          1.0     /* energy scale for purely repulsive LJ pot. */

/* parameters between protein and nucleotide */
#define ai        3.8  /* average distance between two consecutive beads along the chain */
#define ehi       0.2  /*0.9*/ /*1.5*/ /*2.0*/ /*1.0*/ /*0.8*/ /* energy scale for attract LJ pot */
#define eli       0.1  /*1.0*/     /* energy scale for purely repulsive LJ pot. */
#define zetai     80.0 /* (zeta*(1/tau_L)*(eh/temp) sets the natural time for the overdamped motion */

#define numf      383  /* total number of residues in the protein */
#define rigid       2  /* 0=rigid and static nucleotide; 1=normal nucleotide; 2=CM nucleotide (rigid but CM moves) */
/* LJY20100318 end  */

atom *Atom;
atom *Amino, *Amino_new, *Amino_sec, *d_Amino;

 /* LJY added fdatad, f and r components of interest residues */
FILE *fout, *fout1, *fout_3, *fcd, *fdata, *fdatad ;
char filename_b[100];

struct coord{
  double x, y, z;
};

struct coord c;

int    tot_amino;
int    **connect, **link, **connect_old;
/* int    mseed, step, NC0, NC; */
int    mseed, NC0, NC;
unsigned int    step;
double applied_force, R, epot, rg, **R0, **R0_old;

int    *nr_N, *order_N, stepx;
char   **type_N, **mol_N;
double forcer, forcerx, forcery, forcerz, xt, r0;  /* LJY change forcex x0 -> forcer, r0 */

 /* LJY added r component of r vector between COMs, f and r components of interest residues */
double r38, r201, r240;

int    no_contact, heat;
double epot_fene, epot_LJ, tempav, epot_LJ_att, epot_LJ_nei, epot_LJ_rep;

int    SSbond_nr; /* # disulfide bonds in structure */ 
sec_struct     *SSbond;
