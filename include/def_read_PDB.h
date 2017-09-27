 
               /* this is the file def_read_PDB.h */

/* declarations and assignments of values for amino acids in proteins */

#include "def_const_PDB.h"
#include "structures_protein_PDB.h"
#include "def_parameters_aa.h"

char    Exper_proc; /* Exper_proc = experimental procedure (X-ray or NMR)
		       for producing the PDB file */


int     Helix_nr, Strand_nr; /* # helices and strands in the protein */
int     type[N_class+1][N_class+1], hydro_pot[N_pot+1], hydro_res[N_class+1];
/* the index of each 2-body interaction potential b/ the 20 types of a.a. 
 and the hydrophobicity of each such interaction (1 = HH, 2 = HP, 3 = H+, 
 4 = H-, 5 = PP, 6 = P+, 7 = P-, 8 = ++, 9 = +-, 10 = --)*/
int     Atom_nr, N_resid, Chain_nr; /* # of atoms, residues and chains in 
				       given protein */
int     class_prot; /* the class of each protein (based on its length) */
double  Ad_tot, S_tot, trace_T_tot; /* total asphericity and shape and square 
				       of radius of gyration (trace of inertia
				       tensor) (for a protein, i.e. consider 
				       all its chains together) */

double Ad_tot_1, S_tot_1, Ad_tot_2, S_tot_2; /* shape param. for each
					      half of protein */

double  T_tot[4][4]; /* total tensor of inertia of protein */


extern double Hydro_Eisen[N_class+1]; /* hydrophobicities of a.a. (Eisenberg)*/
extern double ref_area_sc[N_class+1]; /* area of a.a. type X in extended 
					 Gly-X-Gly conformation (Creighton) */
extern double ref_vdw_volume[N_class+1];/* average vdw volume of a.a. inside 
					   a protein (Creighton) */

extern char   symbol[N_class+1][4], symbol1[N_class+1]; 
/* symbol (symbol1) = correspondence between the 3-letter (1-letter) symbol
   and the corresponding order number for each of the 20 types of a.a. */

extern int    number_atom[N_class+1];
/* standard nr. of heavy atoms in each type of residue */

extern int    hydro_res2[N_class+1], hydro_res4[N_class+1];
/* classification of a.a. in 2 or 4 classes based on hydrophobicities */

extern double en[N_class+1][N_class+1]; /* pairwise interaction energies b/w amino acids */

extern atom   *Atom;
extern Amino  *AminoA;
extern sec_struct *Helix;
extern sec_struct *Strand;
extern chain_seq Chain[Chain_max];
extern hydrophob_segm  *H_segm;
