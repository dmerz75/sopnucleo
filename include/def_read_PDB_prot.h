 
               /* this is the file def_read_PDB_prot.h */

/* declarations and assignments of values for amino acids in proteins */

#include "def_const_PDB.h"
#include "structures_protein_PDB.h"
/*#include "parameters_BT.h"*/

/*#include "parameters_KGS.h"*/

extern double Hydro_Eisen[N_class+1];/* hydrophobicities of a.a. (Eisenberg) */
extern double ref_area_sc[N_class+1]; /* area of a.a. type X in extended 
					 Gly-X-Gly conformation (Creighton) */
extern double ref_vdw_volume[N_class+1];/* average vdw volume of a.a. inside 
					   a protein (Creighton) */

extern char   symbol[N_class+1][4], Exper_proc; /* symbol = correspondence 
						   between the 3-letter symbol
						   and the corresponding order 
						   number for each of the 20 
						   types of a.a. ; 
						   Exper_proc = experimental
						   procedure (X-ray or NMR)
						   for producing the PDB file*/
extern char   symbol1[N_class+1];

extern int    hydro_res2[N_class+1], hydro_res4[N_class+1];
/* classification of a.a. in 2 or 4 classes (based on hydrophobicities) */

extern int    hydro_res[N_class+1];

/* charge assignment of residues (when considering Asp,Glu,Lys and Arg 
separately) */
extern int charge_res[N_class+1];

extern  int   classes9_res[N_class+1];
/* classification of a.a. in 9 classes (for obtaining 3-body contacts) */

extern int    Helix_nr, Strand_nr; /* # helices and strands in the protein */

extern int    type[N_class+1][N_class+1], hydro_pot[N_pot+1];
/* the index of each 2-body interaction potential b/ the 20 types of a.a. 
 and the hydrophobicity of each such interaction (1 = HH, 2 = HP, 3 = H+, 
 4 = H-, 5 = PP, 6 = P+, 7 = P-, 8 = ++, 9 = +-, 10 = --)*/

extern int    Atom_nr, N_resid, Chain_nr; 
               /* # of atoms, residues and chains in protein */

extern int    class_prot;/* the class of each protein (based on its length) */
extern double Ad_tot, S_tot, T_tot[4][4]; /* structure factors for a protein */

extern double Ad_tot_1, S_tot_1, Ad_tot_2, S_tot_2;

atom           *Atom;
Amino          *AminoA;
sec_struct     *Helix;
sec_struct     *Strand;
chain_seq      Chain[Chain_max];
hydrophob_segm *H_segm;
