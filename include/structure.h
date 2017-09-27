typedef struct Atom{

  char   Name[4], chain, spec;
  int    num, ind, posn; /* position and index of the whole aminoacid along the
			    sequence */
  int    ssbond_nr; /* number of SSbond to which the residue belongs to */
  int    flag_ssbond; /* if 0, site not yet investigated for disulfide bond; 
		       if 1, site has already been considered for disulfide
		       bond */
  int    linker; /* if 0, the position is not in a linker, but in an actual monomer;
		  if 1, the position is in a linker */
  int    monomer; /* it is used for tandems to denote to which monomer does a given position belong to */
  int    p_init;  /* it is used for tandems to store the info about the position of each amino acid from a monomer in the original PDB file (starting the numbering from 1) */
  double x, y, z, BF, mass;
  double fcx, fcy, fcz;

} atom;


typedef struct Nucl{
   char base, **atom;
   int ind, num, size;
   double mass_tot,CMx,CMy,CMz,CPx,CPy,CPz, 
          CMPx,CMPy,CMPz,mass_P,
          CMSx,CMSy,CMSz,mass_S,
          CMBx,CMBy,CMBz,mass_B,
          CBF,BF_P,BF_S,BF_B,
          fcx,fcy,fcz;
   double *x,*y,*z,*mass,*BF;
} nucl;

typedef struct Sec_struct{

  char    Name[4]; /* name of each helix and each sheet in the given protein */
  int     start, end;

} sec_struct;
