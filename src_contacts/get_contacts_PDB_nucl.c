
             /* this is the file get_contacts_PDB_nucl.c */

/* The command line is: ./run_contacts <datafile> (PDB file of the protein) */

/* compile w/ gcc -O2 get_contacts_PDB_nucl.c read_protein_PDB.c identify.c -lm -Wall -o run_contacts */

/* this program first reads the coordinates of the atoms of a given
protein. Then it finds its intra- and inter-chain contacts b/w
positions separated at least by 3 along the chain. The contacts are
determined based on a min cut-off distance b/w pairs of heavy atoms in
the side-chains of the two residues and the distance b/w the C-alpha
atoms of the 2 residues. The output of this file is used to define the
Go-model part of the SOP model of the protein of interest. The
difference with respect to "get_contacts_PDB.c" is that here we
consider a nucleotide (ATP or ADP) bound to the protein and we account
for the nucleotide as another protein chain with all the heavy atoms
listed as CA in the PDB file and starting with the entry ATOM rather
than HETATM. */
/* LJY20100318 : previous code cannot calculate contant j and k 
 if j-k > 2, even if j and k are intercontact */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>

#include "def_read_PDB.h"


extern void    read_protein_PDB(), dist();

//#define DEBUG 

 #define Only_Calpha  /* if this option is turned on, then the
		       contact map is built based only on dist. b/w
		       Calpha atoms; else, we consider both sc-sc
		       contacts and contact b/w Calphas */

FILE *fp, *out, *to, *to2;

int main(int argc, char *argv[]) 
{

  struct tms before, after;
  char   filename_write[120], filename_write2[120], *s, protein_name[10];
  int    i, j, k, p, q, pos;
  int    chain1, chain2, **contact;
  double d1, distance, d2;

  times(&before);

  /* get the name of the protein from the input file */
  s = argv[1];

  s += 3;
  pos = 0;
  while(*s != '.'){
    protein_name[pos] = *s;
    s ++;
    pos ++;
  }
  protein_name[pos] = '\0';


  read_protein_PDB(argv[1]);

#ifdef DEBUG

  sprintf(filename_write,"Output_%4s", protein_name);
  
  fp = fopen(filename_write,"a");
  fprintf(fp,"Atom_nr = %6d N_resid = %6d Chain_nr = %2d\n", Atom_nr, N_resid,
  Chain_nr);
  fprintf(fp,"Helix_nr = %3d Strand_nr = %3d\n\n", Helix_nr, Strand_nr);
  for(i = 1 ; i <= N_resid ; i ++){
    fprintf(fp, "%5d: %4s %c %c %3d :\n", i, AminoA[i].Name, AminoA[i].chain,
	   AminoA[i].SS, AminoA[i].posn);
    fprintf(fp,"start = %3d end = %3d\n", AminoA[i].start, AminoA[i].end);
    /*for(j = 1 ; j <= AminoA[i].length ; j ++)
      fprintf(fp, "%7.3f %7.3f %7.3f\n", AminoA[i].xs[j], AminoA[i].ys[j],
	      AminoA[i].zs[j]);*/
  }
  fclose(fp);

  /*fp = fopen(filename_write,"a");
  for(i = 1 ; i <= Atom_nr ; i ++){
    fprintf(fp, "%5d %7s %7s %10d %7.3f %7.3f %7.3f\n", i, Atom[i].Name, 
	    Atom[i].aa, Atom[i].posn, Atom[i].x, Atom[i].y, Atom[i].z);
  }
  fclose(fp);*/


  exit(100);

#endif

  printf("for %4s: Atom_nr = %6d N_resid = %6d Chain_nr = %2d\n", 
	 protein_name, Atom_nr, N_resid, Chain_nr);


  /* allocate space for the matrix that tells us if 2 positions are in non-bonded contact */
  contact = (int **) malloc((N_resid+1) * sizeof(int *));

  for(i = 1 ; i <= N_resid ; i ++){
  
    contact[i] = (int *) malloc((N_resid+1) * sizeof(int));

    for(j = 1 ; j <= N_resid ; j ++){
      contact[i][j] = 0; /* initially no 2 positions are connected */
    }

  }


    /* ---------------------------------------------------------------------

     go over all chains in the protein and get the intra- and
     inter-chain sc-sc and b-b contacts (that is, any 2 heavy atoms
     within a distance of Dss or any Calpha atoms within a distance DCA)

     --------------------------------------------------------------------- */
  

  /* get the contact map of the protein (for intra- and inter-chain contacts) */
  /* LJY20100318:  calculate all the j, k pairs */
  for(j = 1; j <= (N_resid-1) ; j ++){
    
    chain1 = 0;

    for(k = (j+1) ; k <= N_resid ; k ++){
      
      chain2 = 0;

      /* determine the chain to which these residues belong to */
      for(p = 1 ; p <= Chain_nr ; p ++){
	if(Chain[p].name == AminoA[j].chain) chain1 = p;
	if(Chain[p].name == AminoA[k].chain) chain2 = p;
      }

#ifdef DEBUG
      printf("for:j = %5d (%4s %5d) chain = %d and k = %5d (%4s %5d) chain = %d\n",
	     j, AminoA[j].Name, AminoA[j].posn, chain1, k, AminoA[k].Name, AminoA[k].posn, chain2);
      
#endif
      
      /* check for errors in chain assignments of the residues */
      if((chain1 == 0)){
	printf("for %5d (%5d) %4s %c %c : chain1 = 0 => Error: exit\n",
	       j, AminoA[j].posn, AminoA[j].Name, AminoA[j].SS, 
	       AminoA[j].chain);
	exit(100);
      }

      if((chain2 == 0)){
	printf("for %5d (%5d) %4s %c %c : chain2 = 0 => Error: exit\n",
	       k, AminoA[k].posn, AminoA[k].Name, AminoA[k].SS, 
	       AminoA[k].chain);
	exit(100);
      }

      /* get the min. dist. b/w resid. at posn. j and at posn. k based on side-chain to side-chain dist. */
      d1 = 100000.;
	
      for(p = 1 ; p <= AminoA[j].length_sc ; p ++){
	for(q = 1; q <= AminoA[k].length_sc ; q ++){
	  distance = 0.;
	  dist(AminoA[j].xs[p],AminoA[k].xs[q],AminoA[j].ys[p],AminoA[k].ys[q],AminoA[j].zs[p],AminoA[k].zs[q],&distance);
	  
	  if(distance < d1){
	    d1 = distance;
	  }
	    
	}/* end cycle q (all sc atoms of residue k) */
      }/* end cycle p (all sc atoms of residue j) */

      /* get the min. dist. b/w resid. at posn. j and at posn. k based on C-alpha to C-alpha dist. */
      d2 = 0.;
      dist(AminoA[j].xb,AminoA[k].xb,AminoA[j].yb,AminoA[k].yb,AminoA[j].zb,AminoA[k].zb,&d2);
	      
	
#ifdef DEBUG

      /* if we want info about some specific contacts in the protein */
      if( ( AminoA[j].posn == 170 ) || ( AminoA[k].posn == 170 ) ){
	printf("for %3d(%4s) %3d(%4s) : d_min = %5.3f\n", AminoA[j].posn,
	       AminoA[j].Name, AminoA[k].posn, AminoA[k].Name, d1);
	
	//exit(100);
      }

#endif


      /* get contacts based on cut-off in dist. b/w heavy atoms in the side-chains or b/w Calpha atoms */
      if( ( d1 <= Dss ) || ( d2 <= DCA ) ){
	
	/* say that positions j and k are in contact */
	contact[j][k] = 1;
	contact[k][j] = 1; /* symmetry */

        /* LJY20100318: select intra and inter */
        if( (AminoA[j].chain == AminoA[k].chain) && (k-j) > 2 ){

          sprintf(filename_write2,"Contact_map_intra_sc_and_b_%s", protein_name);
        
	  out = fopen(filename_write2,"a");
	  fprintf(out,"%5d %5d %7.3f ", j, k, d2);
	  fprintf(out,"\n");
	  fclose(out);
        }

        if( AminoA[j].chain != AminoA[k].chain ){
         
	  sprintf(filename_write2,"Contact_map_inter_sc_and_b_%s", protein_name);
        
	  out = fopen(filename_write2,"a");
	  fprintf(out,"%5d %5d %7.3f ", j, k, d2);
	  fprintf(out,"\n");
	  fclose(out);
        }

	/*out = fopen(filename_write2,"a");
	fprintf(out,"%5d (%4s %d %c) %5d (%4s %d %c) %7.3f %7.3f ", AminoA[j].posn, 
		AminoA[j].Name, hydro_res4[AminoA[j].ind], AminoA[j].chain, 
		AminoA[k].posn, AminoA[k].Name, hydro_res4[AminoA[k].ind],
		AminoA[k].chain, d1, d2);
	if( ( d1 <= Dss ) && ( d2 > DCA ) ){
	  fprintf(out,"sc");
	}

	else if( ( d1 > Dss ) && ( d2 <= DCA ) ) {
	  fprintf(out,"b");
	  }*/

      }


#ifdef Only_Calpha

      /* get contacts only based on cut-off in dist. b/w Calpha atoms */
      if( d2 <= DCA ){
	
	/* say that positions j and k are in contact */
	contact[j][k] = 1;
	contact[k][j] = 1; /* symmetry */

        /* LJY20100318: select intra and inter */
        if( (AminoA[j].chain == AminoA[k].chain) && (k-j) > 2 ){
        
	  sprintf(filename_write2,"Contact_map_intra_b_%s", protein_name);

	  out = fopen(filename_write2,"a");
	  fprintf(out,"%5d %5d %7.3f ", j, k, d2);
	  fprintf(out,"\n");
	  fclose(out);

        }

        if( AminoA[j].chain != AminoA[k].chain ){

          sprintf(filename_write2,"Contact_map_inter_b_%s", protein_name);

	  out = fopen(filename_write2,"a");
	  fprintf(out,"%5d %5d %7.3f ", j, k, d2);
	  fprintf(out,"\n");
	  fclose(out);

        }

	/*out = fopen(filename_write2,"a");
	fprintf(out,"%5d (%4s %d %c) %5d (%4s %d %c) %7.3f %7.3f ", AminoA[j].posn, 
		AminoA[j].Name, hydro_res4[AminoA[j].ind], AminoA[j].chain, 
		AminoA[k].posn, AminoA[k].Name, hydro_res4[AminoA[k].ind],
		AminoA[k].chain, d1, d2);
	if( ( d1 <= Dss ) && ( d2 > DCA ) ){
	  fprintf(out,"sc");
	}

	else if( ( d1 > Dss ) && ( d2 <= DCA ) ) {
	  fprintf(out,"b");
	  }*/

      }

#endif
	
    }/* end cycle k (all residues in protein from (j+3) onward)*/

  }/* end cycle j (all residues in the protein) */


  times(&after);

  printf("User time: %d s\n", (after.tms_utime-before.tms_utime));
  printf("System time: %d s\n", (after.tms_stime-before.tms_stime));



  exit(0);

}
