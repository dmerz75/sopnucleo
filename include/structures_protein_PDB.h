
       /* this is the library file "structures_protein_PDB.h" */

/* it contains the definitions of the structures used to store the data about 
a protein from the PDB. */

typedef struct Atom_tot{

  char   Name[4], aa[4], SS, chain;
  int    posn, ind; /* position and index of the whole amino acid along the 
		       sequence */
  int    coord_nr;  /* coordination number for the given atom */
  double x, y, z, mass, radius, radius0, area, volume;

} atom;

typedef struct Aminoacid{

  char   Name[4], SS; /* amino acid name and secondary structure (H = helix,
			 E = strand, T = turn/loop) */
  
  char   chain, **atom; /* label of chain to which this a.a. belongs;
			  and name of each of its atoms */
  
  char   sheet_name[4]; /* name of sheet to which this a.a. belongs to */

  int    ind, posn, ind_H; /* ind = type of a.a. (out of 20) ; 
			      posn = position of aminoacid along sequence;
			      ind_H = 1 if a.a. in segm. of hydro; 
			      0 otherwise*/
  
  int    sec_nr_H, sec_nr_H_new, sec_nr_E, sec_nr_E_new;  
               /* says to which helix or strand it belongs (0 for loops) */
  
  int    start, end, length, length_sc; /* start/end = starting/ending atom of
					   a.a.; whole length of a.a. and
					   # of atoms in its side-chain */

  int    flag_pair; /* 0 means the a.a. has no pair within a given dist.;
		       1 has at least one pair */

  int    coord_nr;  /* coordination number for the given a.a. */

  int   *neighbor_list; /* stores the list of neighbors for the given residue
			 (i.e the coord_nr residues around this one) */

  double *x, *y, *z, *xs, *ys, *zs, xb, yb, zb; 
  /* 3D coord. of all and just side-chain heavy atoms in the amino acid */
  
  double area, ratio; /* the solvent accessible area and ratio of
			 accessibility of the amino acid */
  double volume;

} Amino;


typedef struct Sec_struct{

  char    chain; /* chain = chain to which this sec. struct. elem.
		    belongs to */

  char    Name[4]; /* name of each helix and each sheet in the given protein */
  int     start, end;

  double  energy; /* energy of given element based on summing up the pairwise interactions that each 
		     of its positions have */

} sec_struct;

typedef struct Chain_seq{

  char    name; /* chain identifier */
  int     start, end, sum_H; /* start/end = nr. of 1st and last residue in the
				chain as they appear in the PDB entry
				(e.g.: if the PDB starts from residue 27,
				then the the start for the 1st chain is 27); 
				sum_H = # hydrophobs in each chain */
  int     start_Ind, end_Ind; /* start and end posn. of chain which tell us
				 the actual length of the chain (e.g.: in the
				 above example start_Ind = 1 for the 
				 1st chain of the protein) */ 
  double  T[4][4], trace_T, T_inert[4][4]; /* inertia tensor, its trace and its eigenvectors*/
  double  Ad, lambda[4]; /* asphericity acc. to Dave and eigenvalues of inertia
			    tensor */
  double  S; /* shape parameter acc. to Dave */
  double  xCM, yCM, zCM; /* position of the center of mass of the chain */
  double  T1[4][4], trace_T1, T2[4][4], trace_T2, Ad1, Ad2, S1, S2;
  double  lambda1[4], lambda2[4];
  /* shape parameters for each half of a chain */

} chain_seq;


typedef struct Hydrophob_segm{

  char    chain;
  int     ind, start, end, length, SS, nr_SS, helix, strand, turn;
  
  double  hydro, ratio, ratio_min, ratio_max;
  /* hydro = hydrophobicity of segm (sum of Eisenberg hydrophobicities of its 
     residues); ratio = average ratio of access. of its residues; 
     ratio_min = minimum ratio of access. among its residues; 
     ratio_max = maximum ratio of access. among its residues */
  
} hydrophob_segm;






