
             /* this is the file read_protein_PDB.c */

/* ------------------------------------------------------------------------
   this program reads the PDB file of a given protein
   that contains 1 or more chains. It assigns to each type of residue a given 
   integer (between 1 and 20). It stores the info about its position along the
   sequence, its sec. struct. assignment, the 3D coordinates of its atoms,
   its number of atoms, its hydrophobicity etc. 

   ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>

#include "def_read_PDB_prot.h"


/*#define  DEBUG*/

#define Superpose_max  100  /* max # superpositions in a protein */
#define line_max       300  /* max # characters on a line */

typedef struct Superpose_seq{

  char    chain;
  int     identity, start, end, index;
  /* identity = original seq. posn. of superposition; index = # of 
     superpositions sharing the same "identity" */

} superpose_seq;

superpose_seq  superpose[Superpose_max];

extern int     string_to_integer();

FILE *fp, *out, *to;

void read_protein_PDB(char *filename_prot) 
{

  char    buffer[line_max], *s, Var_atom1, Var_atom2, Var_aa1, Var_aa2;
  char    file_new[120], filename_Info[120];
  char    filename_write[120], new_indic[5];

  static  double  *D1, *D2, *D3;  /* coord of each atom in the protein */
  static  char    **Type, **Code; /* Type=type of atom (e.g. CA or OXT) */
                                  /* Code=3 letter code of the amino acid */
  static  char    **Number1_h, **Number2_h, **Number1_s, **Number2_s;
  static  char    **Indic, *chain_h, *chain_s, *chain_atom, *chain_pos;
  static char     **sheet_name;
  static  int     *atom_posn, s_elem, *seq_pos;
  
  int     i, j, k, l, p, count, count_h, count_s, pos, tot_atoms;
  int     Ind_pos, tot_chain, posn[L_max];
    
  /* Assign each aminoacid to a hydrophobicity class */
  for(i = 1 ; i <= N_class ; i ++){
    hydro_res[i] = hydro_res2[i];
  }

  /* Assign each a.a. - a.a. interaction a type */
  /*fp = fopen("Index_inter","w");*/
  k = 1;
  for(i = 1 ; i <= N_class ; i ++){
    for(j = i ; j <= N_class ; j ++){
      type[i][j] = k;
      type[j][i] = k;

      if(i <= 7 || i == 9 || i == 10){

	if(j <= 7 || j == 9 || j == 10)
	  hydro_pot[k] = 1; /*hydrophob-hydrophob*/

	else if(j == 16 || j == 19)
	  hydro_pot[k] = 3; /*hydrophob-plus*/

	else if(j == 18 || j == 20)
	  hydro_pot[k] = 4; /*hydrophob-minus*/
      
	else
	  hydro_pot[k] = 2; /*hydrophob-polar*/
      }

      else if(i == 16 || i == 19){
	if(j <= 7 || j == 9 || j == 10)
	  hydro_pot[k] = 3; /*plus-hydrophob*/

	else if(j == 16 || j == 19)
	  hydro_pot[k] = 8; /*plus-plus*/

	else if(j == 18 || j == 20)
	  hydro_pot[k] = 9; /*plus-minus*/
      
	else
	  hydro_pot[k] = 6; /*plus-polar*/
      }

      else if(i == 18 || i == 20){
	if(j <= 7 || j == 9 || j == 10)
	  hydro_pot[k] = 4; /*minus-hydrophob*/

	else if(j == 16 || j == 19)
	  hydro_pot[k] = 9; /*minus-plus*/

	else if(j == 18 || j == 20)
	  hydro_pot[k] = 10; /*minus-minus*/
      
	else
	  hydro_pot[k] = 7; /*minus-polar*/
      }

      else{
	if(j <= 7 || j == 9 || j == 10)
	  hydro_pot[k] = 2; /*polar-hydrophob*/

	else if(j == 16 || j == 19)
	  hydro_pot[k] = 6; /*polar-plus*/

	else if(j == 18 || j == 20)
	  hydro_pot[k] = 7; /*polar-minus*/
      
	else
	  hydro_pot[k] = 5; /*polar-polar*/
      }

      /*fprintf(fp,"%4d %3d %3d %3d\n", k, i, j, hydro_pot[k]);*/
      k ++;

    } 
  }
  /*fclose(fp);*/

  /* every protein has at least 1 chain : we denote it by A */
  Chain_nr = 1;
  Chain[Chain_nr].name = 'A';

#ifdef DEBUG
  sprintf(filename_Info,"Info_protein_%4s", filename_prot);
  
  to = fopen(filename_Info,"a"); 
#endif

  /* ------------- read the structure of the test protein ------------ */

  /* read the helices */
  count_h = 0;
  fp = fopen(filename_prot,"r");  
  while(fgets(buffer,L_max,fp) && (buffer[0] != '\n')){
    if(!strncmp(buffer,"HELIX",5)){
      count_h ++;
    }
  }
  
  if(count_h != 0){

    /* alocate space for the starting/ending points (residues) of each helix */
    Number1_h = (char **) malloc(count_h * sizeof(char *));
    Number2_h = (char **) malloc(count_h * sizeof(char *));

    /* alocate space for the chain identifier of these residues (a helix 
       cannot belong to more than 1 chain) */
    chain_h = (char *) malloc(count_h * sizeof(char));

  }

  for(i = 0 ; i < count_h ; i ++){
    chain_h[i] = 'A';
  }
  
  rewind(fp);
  
  if(count_h != 0){
    
    i = 0;
    while(fgets(buffer,L_max,fp) && (buffer[0] != '\n')){  
      if(!strncmp(buffer,"HELIX",5)){
	
	s = buffer;
	
	s += 19;

	/* if it encounters a chain indicator for i = 0, set Chain_nr to 0 */
	if(*s != ' '){
	  chain_h[i] = *s;

	  if(i == 0){
	    Chain_nr = 0;
	  }

	}

	s += 2;

	Number1_h[i] = (char *) malloc(5 * sizeof(char));
	for(pos = 0 ; pos <= 4 ; Number1_h[i][pos++] = *s++);
	Number1_h[i][pos] = '\0';

#ifdef DEBUG
	for(pos = 0 ; pos <= 5 ; pos ++)
	  printf("Number1[%d] = %c\n", pos, Number1_h[i][pos]);
#endif

	s += 7;

	Number2_h[i] = (char *) malloc( 5 * sizeof(char));
	for(pos = 0 ; pos <= 4 ; Number2_h[i][pos++] = *s++);
	Number2_h[i][pos] = '\0';
	
	i ++;

      }
      
    }
  }
  fclose(fp);

  printf("\n\nfor protein %s\n", filename_prot);

#ifdef DEBUG
     
  printf("\nfor Helices:\n\n");
     
  for(i = 0 ; i < count_h ; i ++){
    printf("%5d: (%c %5s %5s) \n", i, chain_h[i], Number1_h[i], Number2_h[i]);
  }

#endif
  
  /* read the sheets */
  count_s = 0;
  fp = fopen(filename_prot,"r");
  
  while(fgets(buffer,L_max,fp) && (buffer[0] != '\n')){
    if(!strncmp(buffer,"SHEET",5)){
      count_s++;
    }
  }

  if(count_s != 0){

    Number1_s = (char **) malloc(count_s * sizeof(char *));
    Number2_s = (char **) malloc(count_s * sizeof(char *));

    chain_s = (char *) malloc(count_s * sizeof(char));
    sheet_name = (char **) malloc(count_s * sizeof(char *));
      
  }

  for(i = 0 ; i < count_s ; i ++){
    chain_s[i] = 'A';
  }

  rewind(fp);

  if(count_s != 0){

    i = 0;
    while(fgets(buffer,L_max,fp) && (buffer[0] != '\n')){  
      if(!strncmp(buffer,"SHEET",5)){
	
	s = buffer;

	/* read the name of the sheet (based on at most 3 characters) */
	s += 11;

	sheet_name[i] = (char *) malloc(4 * sizeof(char));
	for(pos = 0 ; pos < 3 ; sheet_name[i][pos++] = *s++);
	sheet_name[i][pos] = '\0';

	/*s += 21;*/
	s += 7;

	/* read the name of the chain to which this strand belongs to and,
	 if there are no helices in this protein and we are at i = 0,
	 initialize the names and number of chains in the protein */
	if(*s != ' '){
	  chain_s[i] = *s;
	  
	  if((count_h == 0) && (i == 0)){
	    Chain_nr = 0;
	  }

	}

	s ++;
	
	Number1_s[i] = (char *) malloc(5 * sizeof(char));
	for(pos = 0 ; pos <= 4 ; Number1_s[i][pos++] = *s++);
	Number1_s[i][pos] = '\0';
	
	s += 6;

	Number2_s[i] = (char *) malloc(5 * sizeof(char));
	for(pos = 0 ; pos <= 4 ; Number2_s[i][pos++] = *s++);
	Number2_s[i][pos] = '\0';
	
	i ++;

      }

    }
  }
  fclose(fp);

#ifdef DEBUG

  printf("\nfor Sheets:\n\n");

  for(i = 0 ; i < count_s ; i ++){
    printf("%5d: (%c %5s %5s) \n", i, chain_s[i], Number1_s[i], Number2_s[i]);
  }

#endif

  /* if there is no secondary structure assignment in the header of the 
     PDB file of this protein, exit */
  if((count_h == 0) && (count_s == 0)){

    printf("No sec. struct. assignment in the PDB header of %4s: Exit\n", 
	   filename_prot);

    //exit(100);

    }

  /* -------------- end read the sec. struct. assignment -------------- */

  /* read the positions of the amino acids in the sequence: take only 
   MODEL 1 if NMR, only 1st possible position of atom (from column 13)
   if there are alternatives and only the 1st possible position of amino
   acid (from column 17) if there are alternatives */

  Exper_proc = '\0';

  count = 0;
  fp = fopen(filename_prot,"r");
  
  while(fgets(buffer,L_max,fp) && (buffer[0] != '\n')){
    
    if(!strncmp(buffer,"MODEL        1",14)){
        
      Exper_proc = 'N';

      while(fgets(buffer,L_max,fp) && strncmp(buffer,"ENDMDL",6)){
	if(!strncmp(buffer,"ATOM",4))
	  count ++;
      }

      goto end_read;

    }/* end if we have several NMR models */

    Exper_proc = 'X';

    if(!strncmp(buffer,"ATOM",4)){
      count ++;
    }

  }
   
end_read:
  
  /* position a.a. (to which each atom belongs to) along the sequence */
  Indic      = (char **) malloc((count+1)*sizeof(char *));
  Type       = (char **) malloc((count+1)*sizeof(char *));
  Code       = (char **) malloc((count+1)*sizeof(char *));
  chain_atom = (char *)  malloc((count+1)*sizeof(char *));
  D1         = (double *)malloc((count+1)*sizeof(double));
  D2         = (double *)malloc((count+1)*sizeof(double));
  D3         = (double *)malloc((count+1)*sizeof(double));
    
  /* chain identifier for each atom */
  for(i = 0 ; i <= count ; i ++){
    chain_atom[i] = 'A';
  }

  rewind(fp);

  printf("Exper_proc = %c and count = %7d\n", Exper_proc, count);

#ifdef DEBUG
  fprintf(to,"Exper_proc = %c and count = %7d\n", Exper_proc, count);
#endif

  /* if NMR multi-structures */
  if(Exper_proc == 'N'){

    tot_atoms = 0;

    while(fgets(buffer,L_max,fp) && strncmp(buffer,"ENDMDL",6)){

      if(!strncmp(buffer,"ATOM",4)){

#ifdef DEBUG      
	fprintf(to,"buffer = %s\n", buffer);
#endif

	s = buffer;
	
	s += 12;

	/* if there is just 1 possibility for this atom */
	if(*s == ' ' || *s == '1'){
	  Var_atom1 = 'T';
	}
	else Var_atom1 = 'F';

#ifdef DEBUG
	fprintf(to,"for tot = %7d Var_atom1 = %c\n", tot_atoms, Var_atom1);
#endif

	if(Var_atom1 == 'T'){
	  
	  s ++;

	  /* read the type of atom (N, CA, CB etc.) */
	  pos = 0;
	  Type[tot_atoms+1] = (char *)malloc(4 * sizeof(char));
	  /* get the Type based on non-empty spaces or at most 3 characters */
	  while(*s != ' '){
	    if(pos < 3) Type[tot_atoms+1][pos++] = *s++;
	    else break;
	  }

	  Type[tot_atoms+1][pos] = '\0';

#ifdef DEBUG
	  fprintf(to,"for Type[%d] = %s\n", tot_atoms+1, Type[tot_atoms+1]);
#endif

	  /* check if this is a valid atom */
	  if(!strcmp(Type[tot_atoms+1],"N")||!strcmp(Type[tot_atoms+1],"CA")||
	     !strcmp(Type[tot_atoms+1],"C")||!strcmp(Type[tot_atoms+1],"O") ||
	     !strcmp(Type[tot_atoms+1],"NZ")  ||
	     !strcmp(Type[tot_atoms+1],"ND1") || 
	     !strcmp(Type[tot_atoms+1],"ND2") || 
	     !strcmp(Type[tot_atoms+1],"NE")  || 
	     !strcmp(Type[tot_atoms+1],"NE1") ||
	     !strcmp(Type[tot_atoms+1],"NE2") || 
	     !strcmp(Type[tot_atoms+1],"NH1") ||
	     !strcmp(Type[tot_atoms+1],"NH2") || 
	     !strcmp(Type[tot_atoms+1],"CB")  ||
	     !strcmp(Type[tot_atoms+1],"CD")  || 
	     !strcmp(Type[tot_atoms+1],"CD1") || 
	     !strcmp(Type[tot_atoms+1],"CD2") || 
	     !strcmp(Type[tot_atoms+1],"CE")  || 
	     !strcmp(Type[tot_atoms+1],"CE1") || 
	     !strcmp(Type[tot_atoms+1],"CE2") ||
	     !strcmp(Type[tot_atoms+1],"CE3") || 
	     !strcmp(Type[tot_atoms+1],"CG")  || 
	     !strcmp(Type[tot_atoms+1],"CG1") || 
	     !strcmp(Type[tot_atoms+1],"CG2") ||
	     !strcmp(Type[tot_atoms+1],"CH2") || 
	     !strcmp(Type[tot_atoms+1],"CZ")  || 
	     !strcmp(Type[tot_atoms+1],"CZ2") || 
	     !strcmp(Type[tot_atoms+1],"CZ3") ||
	     !strcmp(Type[tot_atoms+1],"OG")  || 
	     !strcmp(Type[tot_atoms+1],"OH")  ||
	     !strcmp(Type[tot_atoms+1],"OD1") || 
	     !strcmp(Type[tot_atoms+1],"OD2") ||
	     !strcmp(Type[tot_atoms+1],"OG1") || 
	     !strcmp(Type[tot_atoms+1],"OE1") ||
	     !strcmp(Type[tot_atoms+1],"OE2") || 
	     !strcmp(Type[tot_atoms+1],"OXT") ||
	     !strcmp(Type[tot_atoms+1],"SD")  || 
	     !strcmp(Type[tot_atoms+1],"SG")){
	  
	    Var_atom2 = 'T';

	  }

	  else{
	    Var_atom2 = 'F';
	    free(Type[tot_atoms+1]);
	  }

	  /* if this is a standard atom, continue */
	  if(Var_atom2 == 'T'){

	    if(pos == 1) s += 2; /* atom of 1-letter type, e.g. N */
	    else if(pos == 2) s += 1; /* atom of 2-letters type, e.g. CA */
	  
	    /* read the name of the amino acid to which this atom belongs to */

	    /* see if the current amino acid has 1 or more possibilities */
	    if(*s == ' ' || *s == 'A' || *s == '1'){
	      Var_aa1 = 'T';
	    }

	    else{
	      Var_aa1 = 'F';
	      free(Type[tot_atoms+1]);
	    }

#ifdef DEBUG
	    fprintf(to,"for tot_atoms = %7d Var_aa1 = %c\n", tot_atoms, 
		    Var_aa1);
#endif
	  
	    if(Var_aa1 == 'T'){
	    
	      s ++;
	      Code[tot_atoms+1] = (char *) malloc(4 * sizeof(char));
	      for(pos = 0 ; pos < 3 ; Code[tot_atoms+1][pos++] = *s++);
	      Code[tot_atoms+1][3] = '\0';
	      
	      /* check if the type of a.a. is OK, that is if it is one of the
		 20 standard aminoacids */
	      if(strncmp(Code[tot_atoms+1],"CYS",3) && 
		 strncmp(Code[tot_atoms+1],"PHE",3) && 
		 strncmp(Code[tot_atoms+1],"LEU",3) && 
		 strncmp(Code[tot_atoms+1],"TRP",3) && 
		 strncmp(Code[tot_atoms+1],"VAL",3) && 
		 strncmp(Code[tot_atoms+1],"ILE",3) && 
		 strncmp(Code[tot_atoms+1],"MET",3) &&
		 strncmp(Code[tot_atoms+1],"HIS",3) && 
		 strncmp(Code[tot_atoms+1],"TYR",3) &&
		 strncmp(Code[tot_atoms+1],"ALA",3) && 
		 strncmp(Code[tot_atoms+1],"GLY",3) &&
		 strncmp(Code[tot_atoms+1],"PRO",3) && 
		 strncmp(Code[tot_atoms+1],"ASN",3) &&
		 strncmp(Code[tot_atoms+1],"THR",3) && 
		 strncmp(Code[tot_atoms+1],"SER",3) &&
		 strncmp(Code[tot_atoms+1],"ARG",3) && 
		 strncmp(Code[tot_atoms+1],"GLN",3) && 
		 strncmp(Code[tot_atoms+1],"ASP",3) && 
		 strncmp(Code[tot_atoms+1],"LYS",3) &&
		 strncmp(Code[tot_atoms+1],"GLU",3)){
		Var_aa2 = 'F';
		free(Type[tot_atoms+1]);
		free(Code[tot_atoms+1]);
	      }
	    
	      else Var_aa2 = 'T';
	    
	      /* if this is a valid a.a., get the coord. of its atoms */
	      if(Var_aa2 == 'T'){
		
		/* get the chain identifier of the residue */
		s ++;
		if(*s != ' ')
		  chain_atom[tot_atoms+1] = *s;
		
		s ++;

		/* get the position of the residue along the sequence */
		Indic[tot_atoms+1] = (char *) malloc(5 * sizeof(char));
		for(pos = 0 ; pos <= 4 ; Indic[tot_atoms+1][pos++] = *s++);
		Indic[tot_atoms+1][pos] = '\0';
		
#ifdef DEBUG	      
     /*fprintf(to,"for Code[%d] = %s Indic[%d] = %s\n", (tot_atoms+1),
       Code[tot_atoms+1], (tot_atoms+1), Indic[tot_atoms+1]);*/

     printf("for Code[%d] = %s Indic[%d] = %s\n", (tot_atoms+1),
	    Code[tot_atoms+1], (tot_atoms+1), Indic[tot_atoms+1]);
#endif
	      
		s += 3;
	    
		sscanf(s, "%lf %lf %lf", &D1[tot_atoms+1], &D2[tot_atoms+1], &D3[tot_atoms+1]);
	    
		tot_atoms ++;
	  	  
	      }/* end this is a valid amino acid */

	    }/* end just 1 variant of this amino acid */

	  }/* end this is a standard atom */

	}/* end just 1 variant of this atom */

      }/* end look at all atoms of the protein */ 

    }/* end we look up to the end of the 1st model */

  }/* end multiple NMR structures */

  else{

    tot_atoms = 0;
    while(fgets(buffer,L_max,fp) && /*(strncmp(buffer,"HETATM",6))*/
	   (strncmp(buffer,"END",3))){
      
      if(!strncmp(buffer,"ATOM",4)){

#ifdef DEBUG      
	fprintf(to,"\n%s\n", buffer);
#endif

	s = buffer;

	s += 12;

	/* if there is just 1 possibility for this atom */
	if(*s == ' ' || *s == '1'){
	  Var_atom1 = 'T';
	}
	else Var_atom1 = 'F';

#ifdef DEBUG
	fprintf(to,"for tot_atoms = %7d Var_atom = %c\n", tot_atoms, Var_atom);
#endif

	if(Var_atom1 == 'T'){

	  s ++;

	  pos = 0;
	  Type[tot_atoms+1] = (char *)malloc(4 * sizeof(char));
	  /* get the Type based on non-empty spaces or at most 3 characters */

	  while(*s != ' '){
	    if(pos < 3) Type[tot_atoms+1][pos++] = *s++;
	    else break;
	  }

	  Type[tot_atoms+1][pos] = '\0';

#ifdef DEBUG
	  fprintf(to,"for Type[%d] = %s pos = %2d\n", tot_atoms+1,Type[tot_atoms+1], pos);
#endif

	  /* check if this is a standard atom */
	  if(!strcmp(Type[tot_atoms+1],"N")||!strcmp(Type[tot_atoms+1],"CA")||
	     !strcmp(Type[tot_atoms+1],"C")||!strcmp(Type[tot_atoms+1],"O") ||
	     !strcmp(Type[tot_atoms+1],"NZ")  ||
	     !strcmp(Type[tot_atoms+1],"ND1") || 
	     !strcmp(Type[tot_atoms+1],"ND2") || 
	     !strcmp(Type[tot_atoms+1],"NE")  || 
	     !strcmp(Type[tot_atoms+1],"NE1") ||
	     !strcmp(Type[tot_atoms+1],"NE2") || 
	     !strcmp(Type[tot_atoms+1],"NH1") ||
	     !strcmp(Type[tot_atoms+1],"NH2") || 
	     !strcmp(Type[tot_atoms+1],"CB")  ||
	     !strcmp(Type[tot_atoms+1],"CD")  || 
	     !strcmp(Type[tot_atoms+1],"CD1") || 
	     !strcmp(Type[tot_atoms+1],"CD2") || 
	     !strcmp(Type[tot_atoms+1],"CE")  || 
	     !strcmp(Type[tot_atoms+1],"CE1") || 
	     !strcmp(Type[tot_atoms+1],"CE2") ||
	     !strcmp(Type[tot_atoms+1],"CE3") || 
	     !strcmp(Type[tot_atoms+1],"CG")  || 
	     !strcmp(Type[tot_atoms+1],"CG1") || 
	     !strcmp(Type[tot_atoms+1],"CG2") ||
	     !strcmp(Type[tot_atoms+1],"CH2") || 
	     !strcmp(Type[tot_atoms+1],"CZ")  || 
	     !strcmp(Type[tot_atoms+1],"CZ2") || 
	     !strcmp(Type[tot_atoms+1],"CZ3") ||
	     !strcmp(Type[tot_atoms+1],"OG")  || 
	     !strcmp(Type[tot_atoms+1],"OH")  ||
	     !strcmp(Type[tot_atoms+1],"OD1") || 
	     !strcmp(Type[tot_atoms+1],"OD2") ||
	     !strcmp(Type[tot_atoms+1],"OG1") || 
	     !strcmp(Type[tot_atoms+1],"OE1") ||
	     !strcmp(Type[tot_atoms+1],"OE2") || 
	     !strcmp(Type[tot_atoms+1],"OXT") ||
	     !strcmp(Type[tot_atoms+1],"SD")  || 
	     !strcmp(Type[tot_atoms+1],"SG")){

	    Var_atom2 = 'T';

	  }

	  else{
	    Var_atom2 = 'F';
	    free(Type[tot_atoms+1]);
	  }

	  /* if this is a standard atom, continue */
	  if(Var_atom2 == 'T'){

	    if(pos == 1) s += 2;
	    else if(pos == 2) s += 1;
	  	  
	    /* see if the current amino acid has 1 or more possibilities */
	    if(*s == ' ' || *s == 'A' || *s == '1'){

#ifdef DEBUG
	      fprintf(to,"for tot_atoms = %7d s = %c\n", tot_atoms, *s);
#endif

	      Var_aa1 = 'T';
	    }

	    else{
	      Var_aa1 = 'F';
	      free(Type[tot_atoms+1]);
	    }

#ifdef DEBUG
	    fprintf(to,"for tot_atoms = %7d Var_aa1 = %c\n", tot_atoms, Var_aa1);
#endif

	    if(Var_aa1 == 'T'){

	      s ++;
	      Code[tot_atoms+1] = (char *)malloc(4 * sizeof(char));
	      for(pos = 0 ; pos < 3 ; Code[tot_atoms+1][pos++] = *s++);
	      Code[tot_atoms+1][3] = '\0';

	      /* check if the type of a.a. is OK, that is if it is one of
		 the 20 standard aminoacids */
	      if(strncmp(Code[tot_atoms+1],"CYS",3) && 
		 strncmp(Code[tot_atoms+1],"PHE",3) && 
		 strncmp(Code[tot_atoms+1],"LEU",3) && 
		 strncmp(Code[tot_atoms+1],"TRP",3) && 
		 strncmp(Code[tot_atoms+1],"VAL",3) && 
		 strncmp(Code[tot_atoms+1],"ILE",3) && 
		 strncmp(Code[tot_atoms+1],"MET",3) &&
		 strncmp(Code[tot_atoms+1],"HIS",3) && 
		 strncmp(Code[tot_atoms+1],"TYR",3) &&
		 strncmp(Code[tot_atoms+1],"ALA",3) && 
		 strncmp(Code[tot_atoms+1],"GLY",3) &&
		 strncmp(Code[tot_atoms+1],"PRO",3) && 
		 strncmp(Code[tot_atoms+1],"ASN",3) &&
		 strncmp(Code[tot_atoms+1],"THR",3) && 
		 strncmp(Code[tot_atoms+1],"SER",3) &&
		 strncmp(Code[tot_atoms+1],"ARG",3) && 
		 strncmp(Code[tot_atoms+1],"GLN",3) && 
		 strncmp(Code[tot_atoms+1],"ASP",3) && 
		 strncmp(Code[tot_atoms+1],"LYS",3) &&
		 strncmp(Code[tot_atoms+1],"GLU",3)){
		Var_aa2 = 'F';
		free(Type[tot_atoms+1]);
		free(Code[tot_atoms+1]);
	      }
	    
	      else Var_aa2 = 'T';
	    
	      /* if this is a valid a.a., get the coord. of its atoms */
	      if(Var_aa2 == 'T'){

		/* get the chain identifier of the residue */
		s ++;
		
		if(*s != ' ')
		  chain_atom[tot_atoms+1] = *s;
		
		s ++;
	    
		/* get the position of the residue along the sequence */
		Indic[tot_atoms+1] = (char *) malloc(5 * sizeof(char));
		for(pos = 0 ; pos <= 4 ; Indic[tot_atoms+1][pos++] = *s++);
		Indic[tot_atoms+1][pos] = '\0';


#ifdef DEBUG	      
		fprintf(to,"for Code[%d] = %s Indic[%d] = %s\n", tot_atoms+1,
			Code[tot_atoms+1], tot_atoms+1,Indic[tot_atoms+1]);
#endif
	      
		s += 3;
	    
		sscanf(s, "%lf %lf %lf", &D1[tot_atoms+1], &D2[tot_atoms+1], &D3[tot_atoms+1]);
	    
		tot_atoms ++;
	  	  
	      }/* end this is a valid amino acid */

	    }/* end just 1 variant of this amino acid */
	  
	  }/* end this is a standard atom */

	}/* end just 1 variant of this atom */

      }/* end look at all atoms of the protein */ 

    }/* end cycle through the protein */

  }/* end this protein is a valid one */

  fclose(fp);


/*------------end read the structure of the protein---------------------*/

  //printf("\n\nfor protein %s\n", filename_prot);

  printf("nr atoms = %6d nr helix = %3d nr sheet = %3d\n", tot_atoms, count_h, 
	 count_s);

  if(tot_atoms == 0){
    printf("tot_atoms = 0 => Exit\n");
    exit(100);
  }

#ifdef DEBUG

  for(i = 1 ; i <= tot_atoms ; i ++){
    if(Indic[i] != NULL){
      printf("%3d: %4s %4s %c %s => %5.3f %5.3f %5.3f\n", i, Type[i], Code[i],
	     chain_atom[i], Indic[i], D1[i], D2[i], D3[i]);
    }
  }

  exit(0);

#endif 



  /* ---- check what was read and assign each posn. a sec. struct. -------*/
 
#ifdef DEBUG

  fp = fopen("check","a");
  
  fprintf(fp,"\n\nfor protein %s\n", filename_prot);

#endif  

  /* count the number of helices in the protein and determine their 
     starting/ending points */
  Helix_nr = 0;
  for(i = 0 ; i < count_h ; i ++){
    if(Number1_h[i] != NULL){   
      Helix_nr ++;
    }
  }

  /* check for overflow in vector */
  if(Helix_nr > Helix_max){

    printf("Helix_nr (%3d) larger than Helix_max (%3d) =>increase Helix_max\n"
	   , Helix_nr, Helix_max);

    exit(100);

  }

  /* allocate space for these helices */
  Helix     = malloc((Helix_nr+1) * sizeof(sec_struct));

  if(Helix == NULL){
    printf("malloc failed in Helix => Exit\n");
    exit(100);
  }

  for(i = 1 ; i <= Helix_nr ; i ++){
    Helix[i].chain = 'A';
  }

  Helix_nr = 0;
  for(i = 0 ; i < count_h ; i ++){
    if(Number1_h[i] != NULL){
      
      Helix_nr ++;

      strncpy(new_indic,Number1_h[i],4);
      new_indic[4] = '\0';

      Helix[Helix_nr].start = string_to_integer(new_indic);
      Helix[Helix_nr].chain = chain_h[i];
      
#ifdef DEBUG
      fprintf(fp,"%5d %s\n", i, Number1_h[i]);
      fprintf(fp,"%5d start_h = %5d\n", Helix_nr, Helix[Helix_nr].start);
#endif
    
    }

    if(Number2_h[i] != NULL){
    
      strncpy(new_indic,Number2_h[i],4);
      new_indic[4] = '\0';

      Helix[Helix_nr].end = string_to_integer(new_indic);
      
#ifdef DEBUG
      fprintf(fp,"%5d %s\n", i, Number2_h[i]);
      fprintf(fp,"%5d end_h = %5d\n", Helix_nr, Helix[Helix_nr].end);
#endif

    }
  }

  /* count how many strands are in the protein  and determine their 
     starting/ending points */
  Strand_nr = 0;
  for(i = 0 ; i < count_s ; i ++){
    if(Number1_s[i] != NULL){   
      Strand_nr ++;
    }
  }

  /* check for overflow in vector */
  if(Strand_nr > Strand_max){

    printf("Strand_nr(%3d) larger than Strand_max(%3d)=>increase Strand_max\n"
	   , Strand_nr, Strand_max);

    exit(100);

  }

  /* allocate space for these strands */
  Strand     = malloc((Strand_nr+1) * sizeof(sec_struct));

  if(Strand == NULL){
    printf("malloc failed in Strand => Exit\n");
    exit(100);
  }

  for(i = 1 ; i <= Strand_nr ; i ++){
    Strand[i].chain = 'A';
  }
  
  Strand_nr = 0;
  for(i = 0 ; i < count_s ; i ++){
    if(Number1_s[i] != NULL){
      
      Strand_nr ++;

      strncpy(new_indic,Number1_s[i],4);
      new_indic[4] = '\0';

      Strand[Strand_nr].start = string_to_integer(new_indic);
      Strand[Strand_nr].chain = chain_s[i];
      strcpy(Strand[Strand_nr].Name,sheet_name[i]);

#ifdef DEBUG
      fprintf(fp,"%5d %s\n", i, Number1_s[i]);
      fprintf(fp,"%5d start_s = %5d\n", Strand_nr, Strand[Strand_nr].start);
#endif

    }
    
    if(Number2_s[i] != NULL){
      
      strncpy(new_indic,Number2_s[i],4);
      new_indic[4] = '\0';

      Strand[Strand_nr].end = string_to_integer(new_indic);
      
#ifdef DEBUG
      fprintf(fp,"%5d %s\n", i, Number2_s[i]);
      fprintf(fp,"%5d end_s = %5d\n", Strand_nr, Strand[Ind_S].end);
#endif

    }
  }

  /*for(i = 1 ; i <= Helix_nr ; i ++){
    printf("Helix[%d].start = %5d Helix[%d].end = %5d\n", i, Helix[i].start,
	   i, Helix[i].end);
  }*/

  
  seq_pos     = (int *)  malloc((tot_atoms+1) * sizeof(int));
  chain_pos   = (char *) malloc((tot_atoms+1) * sizeof(char));  
  atom_posn   = (int *)  malloc((tot_atoms+1) * sizeof(int));

  /* identify any superpositions along the seq., i.e. both 21 and 21A, 21B */
  s_elem = 0;
  for(i = 1 ; i <= tot_atoms ; i ++){
    
    if(Indic[i] != NULL){

      if(Indic[i][4] != ' '){

	printf("for i = %5d Indic[%d][4] = %c\n", i, i, Indic[i][4]);

	if(strncmp(Indic[i],Indic[i-1],5)){

	  strncpy(new_indic,Indic[i],4);
	  new_indic[4] = '\0';
	    
	  s_elem ++;
	  superpose[s_elem].identity = string_to_integer(new_indic);
       	  superpose[s_elem].chain = chain_atom[i];
	  superpose[s_elem].start = i;

	  if(i > 1){
	    if((Indic[i][0] != Indic[i-1][0]) || (Indic[i][1] != Indic[i-1][1])
	       ||(Indic[i][2] != Indic[i-1][2])||(Indic[i][3]!=Indic[i-1][3]))
	      superpose[s_elem].identity --;

	  }

	  if(i == 1){
	    superpose[s_elem].identity --;
	  }

	}
	
	else if(strncmp(Indic[i],Indic[i+1],5)){
	  superpose[s_elem].end = i;
	}
	
      }

    }

  }

  /* check for overflow in vectors/matrices */
  if(s_elem > Superpose_max){
    printf("s_elem = %4d larger Superpose_max=%4d => increase Superpose_max\n",
	   s_elem, Superpose_max);
    exit(100);
  }

  /* check if we have repeats in superposed, i.e. 21A and 21B etc. */
  for(i = 1 ; i <= s_elem ; i ++){
    superpose[i].index = 1; /* tells us how many identical superposed posn.
			      are before it (e.g. if we have 21A and 21B and 
			      21C the index for 21A is 1, for 21B is 2,
			      for 21C is 3) */
    for(k = 1 ; k <= s_elem ; k ++){
      if((i - k) > 0){
	if(superpose[i].identity == superpose[i-k].identity){
	  printf("We have repeats in superpositions at %3d %3d %c\n",
		 superpose[i].start,superpose[i].identity, superpose[i].chain);
	  superpose[i].index += 0; 
	}
      }
    }
  }

  if(s_elem > 0){
    printf("We have %3d superpositions\n", s_elem);
    for(k = 1 ; k <= s_elem ; k ++)
      printf("%3d %3d %3d (%3d) %c\n", superpose[k].start, superpose[k].end,
	     superpose[k].identity, superpose[k].index, superpose[k].chain);
  }

  if(s_elem > 0){

    printf("s_elem = %3d\n", s_elem);

    for(k = 1 ; k <= s_elem ; k ++){
      
      if(k == 1){
    
	for(i = 1 ; i <= tot_atoms ; i ++){
    
	  if((Indic[i] != NULL) && (chain_atom[i] == superpose[k].chain)){

	    if(i < superpose[k].start){
	      strncpy(new_indic,Indic[i],4);
	      new_indic[4] = '\0';
	      atom_posn[i] = string_to_integer(new_indic);	
	    }

	    else if((i >= superpose[k].start) && (i <= superpose[k].end)){
	      atom_posn[i] = superpose[k].identity + 1;	      
	    }

	    else if(i > superpose[k].end){

	      if((Indic[i][4] == ' ') || ((Indic[i][4] != ' ') && 
					  ((Indic[i][0] == Indic[i-1][0]) &&
					   (Indic[i][1] == Indic[i-1][1]) &&
					   (Indic[i][2] == Indic[i-1][2]) &&
					   (Indic[i][3] == Indic[i-1][3])))){

		strncpy(new_indic,Indic[i],4);
		new_indic[4] = '\0';
		atom_posn[i] = string_to_integer(new_indic) + 1;
		
	      }

	      else if(((Indic[i][4]!=' ')&&((Indic[i][0] != Indic[i-1][0])
					    ||(Indic[i][1] != Indic[i-1][1])|| 
					    (Indic[i][2] != Indic[i-1][2]) || 
					    (Indic[i][3] != Indic[i-1][3])))){
	
		strncpy(new_indic,Indic[i],4);
		new_indic[4] = '\0';
		atom_posn[i] = string_to_integer(new_indic);

	      }

	    }/* end i follows the superposition k */

	  }/* end i is a valid atom of the chain w/ superposition k */

	}/* end cycle i, i.e. going over all atoms */

#ifdef DEBUG
	printf("\n\nafter k = %d : %3d %3d %c :\n", k, superpose[k].start,
	       superpose[k].end, superpose[k].chain);
	for(j = 1; j <= tot_atoms ; j ++){
	  if(Indic[j] != NULL && chain_atom[j] == superpose[k].chain)
	    printf("%5d: %3d %c\n", j, atom_posn[j], chain_atom[j]);
	}
#endif
	
      }/* end k == 1 */
  
      else{

	if(superpose[k].chain == superpose[k-1].chain){

#ifdef DEBUG
	  printf("\n\nfor k = %d (%3d %3d) %c : \n", k, superpose[k].start, 
		 superpose[k].end, superpose[k].chain);
#endif	

	  for(i = superpose[k].start ; i <= tot_atoms ; i ++){
    
	    if((Indic[i] != NULL) && (chain_atom[i] == superpose[k].chain)){

	      atom_posn[i] += 1;
	      
	    }/* end this is a valid atom */

	  }/* end cycle i, i.e. going over all atoms */

#ifdef DEBUG
	  printf("\n\nafter k = %d : %3d %3d %c :\n", k, superpose[k].start,
		 superpose[k].end, superpose[k].chain);
	  for(j = 1; j <= tot_atoms ; j ++){
	    if(Indic[j] != NULL && chain_atom[j] == superpose[k].chain)
	      printf("%5d: %3d %c\n", j, atom_posn[j], chain_atom[j]);
	  }
#endif

	}/* end k is in same chain as k-1 */

	else{

	  for(i = 1 ; i <= tot_atoms ; i ++){
    
	    if((Indic[i] != NULL) && (chain_atom[i] == superpose[k].chain)){

	      if(i < superpose[k].start){
		strncpy(new_indic,Indic[i],4);
		new_indic[4] = '\0';
		atom_posn[i] = string_to_integer(new_indic);	
	      }

	      else if((i >= superpose[k].start) && (i <= superpose[k].end)){
		atom_posn[i] = superpose[k].identity + 1;
	      }
	      
	      else if(i > superpose[k].end){

		strncpy(new_indic,Indic[i],4);
		new_indic[4] = '\0';
		atom_posn[i] = string_to_integer(new_indic) + 1;
				
	      }/* end i follows the superposition k */
	      
	    }/* end i is a valid atom of the chain w/ superposition k */
	    
	  }/* end cycle i, i.e. going over all atoms */

	}/* end k is in different chain from k-1 */

      }/* end k != 1 */

    }/* end going over all superpositions */

    /* if there are no superpositions in a chain */
    for(i = 1 ; i <= tot_atoms ; i ++){
      
      if(Indic[i] != NULL){

	tot_chain = 0;
 
	for(k = 1 ; k <= s_elem ; k ++){
	  if(superpose[k].chain != chain_atom[i]) tot_chain ++;
	}

	if(tot_chain == s_elem){
	  strncpy(new_indic,Indic[i],4);
	  new_indic[4] = '\0';
	  atom_posn[i] = string_to_integer(new_indic); 
	}

      }
    }


#ifdef DEBUG   
    for(i = 1; i <= tot_atoms ; i ++){
      if(Indic[i] != NULL){
	strncpy(new_indic,Indic[i],4);
	new_indic[4] = '\0';
	printf("%5d: (%5d) %5d %c\n", i, string_to_integer(new_indic), 
	       atom_posn[i], chain_atom[i]);
      }
    }
#endif


    /* see if the position of the sec. struct. elem. changes, due to 
     s_elem != 0 */

    /* go over all helices */
    for(p = 1 ; p <= Helix_nr ; p ++){

      for(i = 1 ; i <= tot_atoms ; i ++){
	if(Indic[i] != NULL){
	  strncpy(new_indic,Indic[i],4);
	  new_indic[4] = '\0';
	  if((string_to_integer(new_indic) == Helix[p].start) && 
	     (chain_atom[i] == Helix[p].chain)){
	    Helix[p].start = atom_posn[i];
	    break;
	  }
	}
      }


      for(i = 1 ; i <= tot_atoms ; i ++){
	if(Indic[i] != NULL){
	  strncpy(new_indic,Indic[i],4);
	  new_indic[4] = '\0';
	  if((string_to_integer(new_indic) == Helix[p].end) && 
	     (chain_atom[i] == Helix[p].chain)){
	    Helix[p].end = atom_posn[i];
	    break;
	  }
	}
      }

#ifdef DEBUG
      printf("Helix %3d after : %3d %3d %c\n", p, Helix[p].start, 
	     Helix[p].end, Helix[p].chain);
#endif

    }

    /* go over all strands */
    for(p = 1 ; p <= Strand_nr ; p ++){

      for(i = 1 ; i <= tot_atoms ; i ++){
	if(Indic[i] != NULL){
	  strncpy(new_indic,Indic[i],4);
	  new_indic[4] = '\0';
	  if((string_to_integer(new_indic) == Strand[p].start) && 
	     (chain_atom[i] == Strand[p].chain)){
	    Strand[p].start = atom_posn[i];
	    break;
	  }
	}
      }

      for(i = 1 ; i <= tot_atoms ; i ++){
	if(Indic[i] != NULL){
	  strncpy(new_indic,Indic[i],4);
	  new_indic[4] = '\0';
	  if((string_to_integer(new_indic) == Strand[p].end) && 
	     (chain_atom[i] == Strand[p].chain)){
	    Strand[p].end = atom_posn[i];
	    break;
	  }
	}
      }
      
#ifdef DEBUG
      printf("Strand %3d after : %3d %3d %c\n", p, Strand[p].start, 
	     Strand[p].end, Strand[p].chain);
#endif

    }
  
  }/* end s_elem > 0 , i.e. there are superpositions in this protein */
  
  
  else{

    for(i = 1 ; i <= tot_atoms ; i ++){
      if(Indic[i] != NULL){
	strncpy(new_indic,Indic[i],4);
	new_indic[4] = '\0';
	atom_posn[i] = string_to_integer(new_indic);
      }
    }

  }/* if there are no superpositions in this protein */

  
  Atom = malloc((tot_atoms+1) * sizeof(atom));

  if(Atom == NULL){
    printf("Error in malloc for Atom => Exit\n");
    exit(100);
  }

  /* get the info about each atom in the sequence, the position of each a.a. 
     along the sequence, its chain identifier 
     and secondary structure assignment; also, get the starting and ending
     point for each chain and its identifier */
  Ind_pos = 0;
  Chain_nr = 0;
  Atom_nr = 0;
  for(i = 1 ; i <= tot_atoms ; i ++){
   
    if(Indic[i] != NULL){

#ifdef DEBUG
      fprintf(fp,"%5d %s\n", i, Indic[i]);
#endif

      Atom_nr ++;
      strcpy(Atom[Atom_nr].Name,Type[i]);
      strcpy(Atom[Atom_nr].aa,Code[i]);
      Atom[Atom_nr].posn = atom_posn[i];
      Atom[Atom_nr].chain = chain_atom[i];
      Atom[Atom_nr].x = D1[i];
      Atom[Atom_nr].y = D2[i];
      Atom[Atom_nr].z = D3[i];

      if(i == 1){

	Ind_pos ++;
	seq_pos[Ind_pos] = atom_posn[i];	
	chain_pos[Ind_pos] = chain_atom[i];

	Chain_nr ++;
	Chain[Chain_nr].start = seq_pos[Ind_pos];
	Chain[Chain_nr].start_Ind = Ind_pos;
	Chain[Chain_nr].name  = chain_atom[i];

#ifdef DEBUG
	fprintf(fp,"%5d seq_pos = %5d\n", Ind_pos, seq_pos[Ind_pos]);
#endif

      }
      
      /*monitor the position only for the whole a.a.*/
      else{

	if(atom_posn[i] == atom_posn[i-1]);

	else{

	  Ind_pos ++;
	  seq_pos[Ind_pos] = atom_posn[i];
	  chain_pos[Ind_pos] = chain_atom[i];

#ifdef DEBUG
	  fprintf(fp,"%5d seq_pos = %5d\n", Ind_pos, seq_pos[Ind_pos]);
#endif

	}

	if(i < tot_atoms){
	    
	  if((chain_atom[i] == chain_atom[i-1]) && 
	     (chain_atom[i] != chain_atom[i+1])){
	    Chain[Chain_nr].end = seq_pos[Ind_pos];
	    Chain[Chain_nr].end_Ind = Ind_pos;
	  }

	  else if((chain_atom[i] == chain_atom[i+1]) && 
		  (chain_atom[i] != chain_atom[i-1])){
	    Chain_nr ++;
	    Chain[Chain_nr].start = seq_pos[Ind_pos];
	    Chain[Chain_nr].start_Ind = Ind_pos;
	    Chain[Chain_nr].name  = chain_atom[i];
	  }

	}

	else{
	  Chain[Chain_nr].end = seq_pos[Ind_pos];
	  Chain[Chain_nr].end_Ind = Ind_pos;
	}

      }

    }/* Indic for this atom exists */

  }/* end cycle i, i.e. going over all atoms in protein */

  /* check for overflow in vectors/matrices */
  if(Atom_nr > Atom_max){
    printf("Atom_nr = %6d larger than Atom_max = %6d => increase Atom_max\n",
	   Atom_nr, Atom_max);
    exit(100);
  }

  if(Ind_pos > L_max){
    printf("Ind_pos = %5d larger than L_max = %5d => increase L_max\n",
	   Ind_pos, L_max);
    exit(100);
  }

  if(Chain_nr > Chain_max){
    printf("Chain_nr = %4d larger than Chain_max = %4d =>increase Chain_max\n",
	   Chain_nr, Chain_max);
    exit(100);
  }

#ifdef DEBUG

  for(i = 1 ; i <= Chain_nr ; i ++){
    printf("%2d : %c (%3d %3d)\n", i, Chain[i].name, Chain[i].start, 
	   Chain[i].end);
  }

#endif

  AminoA = malloc((Ind_pos+1) * sizeof(Amino));

  if(AminoA == NULL){
    printf("Error in malloc for AminoA => Exit\n");
    exit(100);
  }

  /* initialize the sec. struct. assignment and chain identifier for each 
     a.a. of the sequence */
  for(i = 1 ; i <= Ind_pos ; i ++){
    AminoA[i].SS = 'T';
    AminoA[i].sec_nr_H = 0;
    AminoA[i].sec_nr_E = 0;
    AminoA[i].chain = chain_pos[i];
    AminoA[i].posn = seq_pos[i];
    strcpy(AminoA[i].sheet_name,"NON");
  }

  for(k = 1 ; k <= Ind_pos ; k ++){
	
    /* identify the a.a. in helices */

    if(Helix_nr != 0){
      for(j = 1 ; j <= Helix_nr ; j ++){
	
	if(AminoA[k].chain == Helix[j].chain){

	  /* !!!! do not neglect the N and C-caps of the helix !!!!!! */
	  if((AminoA[k].posn >= Helix[j].start) && 
	     (AminoA[k].posn <= Helix[j].end)){
	
	    /* !!!!!! if this helix is longer than L_hel residues !!!! */
	    if((Helix[j].end - Helix[j].start + 1) >= L_hel){
	      AminoA[k].SS = 'H';
	      AminoA[k].sec_nr_H = j;
	    }
	  }
	}
      }
    }


    /* identify the a.a. in strands */
    
    if(Strand_nr != 0){
	
      for(j = 1 ; j <= Strand_nr ; j ++){
	  
	if(Strand[j].chain == AminoA[k].chain){

	  if((AminoA[k].posn >= Strand[j].start) && 
	     (AminoA[k].posn <= Strand[j].end)){
	    AminoA[k].SS = 'E';
	    AminoA[k].sec_nr_E = j;
	    strcpy(AminoA[k].sheet_name,Strand[j].Name); 
	  }
	}
      }
      
    }

#ifdef DEBUG    
    fprintf(fp,"%5d %5d %c %s\n", k, AminoA[k].posn, AminoA[k].SS, 
	    AminoA[i].sheet_name);
#endif  

  }/* end cycle k, i.e. going over all residues */

#ifdef DEBUG
  fclose(fp);
#endif  

  printf("length from S.S. is %5d\n", Ind_pos);

  /*------------ assign index to each type of a.a. -----------------------*/

  N_resid= 0;
  /*fp = fopen("test", "a");*/
  for(i = 1 ; i <= tot_atoms ; i ++){

    if(Indic[i] != NULL){

      /* if we are looking at the CA atom of an amino acid */
      if(!strcmp(Type[i],"CA")){

	N_resid ++;
	posn[N_resid] = i /* i-1 */; /* says from where the a.a. #N_resid starts */
                                     /* changed from i-1 to i, i-1 produce error
 *                                      from PDB with only CA atoms    LJY    */
	if(!strcmp(Code[i],"CYS")){
	  AminoA[N_resid].ind = 1;
	}
	else if(!strcmp(Code[i],"PHE")){
	  AminoA[N_resid].ind = 2;
	}
	else if(!strcmp(Code[i],"LEU")){
	  AminoA[N_resid].ind = 3;
	}
	else if(!strcmp(Code[i],"TRP")){ 
	  AminoA[N_resid].ind = 4;
	}
	else if(!strcmp(Code[i],"VAL")){ 
	  AminoA[N_resid].ind = 5;
	}
	else if(!strcmp(Code[i],"ILE")){ 
	  AminoA[N_resid].ind = 6;
	}
	else if(!strcmp(Code[i],"MET")){ 
	  AminoA[N_resid].ind = 7;
	}
	else if(!strcmp(Code[i],"HIS")){ 
	  AminoA[N_resid].ind = 8;
	}
	else if(!strcmp(Code[i],"TYR")){ 
	  AminoA[N_resid].ind = 9;
	}
	else if(!strcmp(Code[i],"ALA")){ 
	  AminoA[N_resid].ind = 10;
	} 
	else if(!strcmp(Code[i],"GLY")){ 
	  AminoA[N_resid].ind = 11;
	}
	else if(!strcmp(Code[i],"PRO")){ 
	  AminoA[N_resid].ind = 12;
	}
	else if(!strcmp(Code[i],"ASN")){ 
	  AminoA[N_resid].ind = 13;
	}
	else if(!strcmp(Code[i],"THR")){ 
	  AminoA[N_resid].ind = 14;
	}
	else if(!strcmp(Code[i],"SER")){ 
	  AminoA[N_resid].ind = 15;
	}
	else if(!strcmp(Code[i],"ARG")){ 
	  AminoA[N_resid].ind = 16;
	}
	else if(!strcmp(Code[i],"GLN")){ 
	  AminoA[N_resid].ind = 17;
	}
	else if(!strcmp(Code[i],"ASP")){ 
	  AminoA[N_resid].ind = 18;
	}
	else if(!strcmp(Code[i],"LYS")){ 
	  AminoA[N_resid].ind = 19;
	}
	else if(!strcmp(Code[i],"GLU")){ 
	  AminoA[N_resid].ind = 20;
	}

	strcpy(AminoA[N_resid].Name,Code[i]);

	/*fprintf(fp, "N_resid= %d %-4s ind[N] = %d\n",N_resid,Code[i],
	  N_resid,ind[N_resid]);*/

      }/* end this atom is an alpha carbon */

    }/* end this atom is valid */

  }/* end going over all atoms */

  /*fclose(fp);*/

/*--------------------end assignment of indices----------------------------*/

#ifdef DEBUG1

  fp = fopen("test","a");
  fprintf(fp,"\n\nfor protein %s\n", filename_prot);
  
  for(i = 1 ; i <= Ind_pos ; i ++){
    fprintf(fp,"%4s %3d %c %c\n", AminoA[i].Name, seq_pos[i], 
	    AminoA[i].chain, AminoA[i].SS);
  }
  
  fclose(fp);
  

  for(j = 1; j <= N_resid; j ++){
    printf("%5d: %4s %c %c %3d\n", j, AminoA[j].Name, AminoA[j].chain,
	   AminoA[j].SS, AminoA[j].posn);
  }
  
#endif

  printf("length from seq. (N_resid) = %5d\n", N_resid);

  /* checking point (see if the info about length of sequence obtained from 
     sec struct and from CA atoms is the same) */
  if(Ind_pos != N_resid){
    
    printf("Ind_pos (%3d) different from N_resid(%3d): Exit\n", Ind_pos, 
	   N_resid);
    
#ifdef DEBUG

    fprintf(to,"Ind_pos (%3d) larger than N_resid(%3d): exiting\n", Ind_pos, 
	    N_resid);
    fclose(to);

#endif

    exit(100);
    
  }

  /* find the start and end position for each a.a. (that is, from the N to the
     C terminal), its total number
     of atoms and init. the 3D coord. for each atom in each a.a. */
  for(i = 1 ; i < N_resid; i ++){

    AminoA[i].start  = posn[i];
    AminoA[i].end    = posn[i+1] - 1;
    AminoA[i].length = (AminoA[i].end - AminoA[i].start) + 1;

    AminoA[i].x  = (double *)malloc((AminoA[i].length+1)*sizeof(double));
    AminoA[i].y  = (double *)malloc((AminoA[i].length+1)*sizeof(double));
    AminoA[i].z  = (double *)malloc((AminoA[i].length+1)*sizeof(double));
    AminoA[i].xs = (double *)malloc((AminoA[i].length+1)*sizeof(double));
    AminoA[i].ys = (double *)malloc((AminoA[i].length+1)*sizeof(double));
    AminoA[i].zs = (double *)malloc((AminoA[i].length+1)*sizeof(double));

    for(j = 1 ; j <= AminoA[i].length ; j ++){
      AminoA[i].x[j] = AminoA[i].y[j] = AminoA[i].z[j] = 1000.;
      AminoA[i].xs[j] = AminoA[i].ys[j] = AminoA[i].zs[j] = 1000.;
    }

    AminoA[i].atom = (char **)malloc((AminoA[i].length+1)*sizeof(char *));

    for(j = 1 ; j <= AminoA[i].length ; j ++){
      AminoA[i].atom[j] = (char *)malloc((3+1)*sizeof(char));
    }

  }
  
  AminoA[N_resid].start  = posn[N_resid];
  AminoA[N_resid].end    = Atom_nr;
  AminoA[N_resid].length = (AminoA[N_resid].end - AminoA[N_resid].start) + 1;
  
  AminoA[N_resid].x = 
    (double *)malloc((AminoA[N_resid].length+1)*sizeof(double));
  AminoA[N_resid].y = 
    (double *)malloc((AminoA[N_resid].length+1)*sizeof(double));
  AminoA[N_resid].z = 
    (double *)malloc((AminoA[N_resid].length+1)*sizeof(double));
  AminoA[N_resid].xs = 
    (double *)malloc((AminoA[N_resid].length+1)*sizeof(double));
  AminoA[N_resid].ys = 
    (double *)malloc((AminoA[N_resid].length+1)*sizeof(double));
  AminoA[N_resid].zs = 
    (double *)malloc((AminoA[N_resid].length+1)*sizeof(double));
  for(j = 1 ; j <= AminoA[N_resid].length ; j ++){
    
    AminoA[N_resid].x[j] = AminoA[N_resid].y[j] = AminoA[N_resid].z[j] = 
      1000.;
    AminoA[N_resid].xs[j] = AminoA[N_resid].ys[j] = AminoA[N_resid].zs[j] = 
      1000.;
  }

  AminoA[N_resid].atom = 
    (char **)malloc((AminoA[N_resid].length+1)*sizeof(char *));

  for(j = 1 ; j <= AminoA[N_resid].length ; j ++){
    AminoA[N_resid].atom[j] = (char *)malloc((3+1)*sizeof(char));
  }
    
#ifdef DEBUG

  fprintf(to,"Intially\n\n");
  for(i = 1 ; i <= N_resid; i ++){
    fprintf(to,"%3d  %5d  %5d  %3d\n", AminoA[i].ind, AminoA[i].start,
	    AminoA[i].end, AminoA[i].length);
  }

#endif

  /* assign to each atom the sec. struct. of its a.a. and its index (ind) */
  for(i = 1 ; i <= N_resid; i ++){    
    for(j = AminoA[i].start ; j <= AminoA[i].end ; j ++){
      Atom[j].SS  = AminoA[i].SS;
      Atom[j].ind = AminoA[i].ind;
    }
  }


  /* find the 3D coord. of all heavy atoms in the side-chain of each a.a.
     and of the CA atom */
  for(i = 1 ; i <= N_resid; i ++){
    
    /* select the CA atom for the backbone */
    for(j = AminoA[i].start ; j <= AminoA[i].end ; j ++){
      
      if(!strcmp(Atom[j].Name,"CA")){
	AminoA[i].xb = Atom[j].x;
	AminoA[i].yb = Atom[j].y;
	AminoA[i].zb = Atom[j].z;
      }/*end select only the CA atom*/

    }/*end look at all atoms in this a.a.*/
    
    k = 0;
    
    /*if we have a GLY select only its CA atom*/
    if(AminoA[i].ind == 11){

      l = 0;
      for(j = AminoA[i].start ; j <= AminoA[i].end ; j ++){
	
	l ++;
	AminoA[i].x[l] = Atom[j].x;
	AminoA[i].y[l] = Atom[j].y;
	AminoA[i].z[l] = Atom[j].z;
	strcpy(AminoA[i].atom[l],Atom[j].Name);

	if(!strcmp(Atom[j].Name,"CA")){
	  k ++;
	  AminoA[i].xs[k] = Atom[j].x;
	  AminoA[i].ys[k] = Atom[j].y;
	  AminoA[i].zs[k] = Atom[j].z;
	  
	}/* end select only the CA atom */
      }/* end look at all atoms in this a.a. */
    }/* end if the a.a. is GLY */
    
    /* for all other a.a. select all heavy atoms in the side-chain */
    else{

      l = 0;
      for(j = AminoA[i].start ; j <= AminoA[i].end ; j ++){

	l ++;

	/* here we store the info about each heavy atom in the a.a. */
	AminoA[i].x[l] = Atom[j].x;
	AminoA[i].y[l] = Atom[j].y;
	AminoA[i].z[l] = Atom[j].z;
	strcpy(AminoA[i].atom[l],Atom[j].Name);


	/* here we only store the info about heavy atoms in the side-chain */
	if(!strcmp(Atom[j].Name,"ND1") || !strcmp(Atom[j].Name,"ND2") || 
	   !strcmp(Atom[j].Name,"NE")  || !strcmp(Atom[j].Name,"NE1") ||
	   !strcmp(Atom[j].Name,"NE2") || !strcmp(Atom[j].Name,"NH1") ||
	   !strcmp(Atom[j].Name,"NH2") || !strcmp(Atom[j].Name,"CB")  ||
	   !strcmp(Atom[j].Name,"NZ")  ||
	   !strcmp(Atom[j].Name,"CD")  || !strcmp(Atom[j].Name,"CD1") || 
	   !strcmp(Atom[j].Name,"CD2") || !strcmp(Atom[j].Name,"CE")  || 
	   !strcmp(Atom[j].Name,"CE1") || !strcmp(Atom[j].Name,"CE2") ||
	   !strcmp(Atom[j].Name,"CE3") || !strcmp(Atom[j].Name,"CG")  || 
	   !strcmp(Atom[j].Name,"CG1") || !strcmp(Atom[j].Name,"CG2") ||
	   !strcmp(Atom[j].Name,"CH2") || !strcmp(Atom[j].Name,"CZ")  || 
	   !strcmp(Atom[j].Name,"CZ2") || !strcmp(Atom[j].Name,"CZ3") ||
	   !strcmp(Atom[j].Name,"OG")  || !strcmp(Atom[j].Name,"OH")  ||
	   !strcmp(Atom[j].Name,"OD1") || !strcmp(Atom[j].Name,"OD2") ||
	   !strcmp(Atom[j].Name,"OG1") || !strcmp(Atom[j].Name,"OE1") ||
	   !strcmp(Atom[j].Name,"OE2") || !strcmp(Atom[j].Name,"OXT") ||
	   !strcmp(Atom[j].Name,"SD")  || !strcmp(Atom[j].Name,"SG")){
	  k ++;
	  AminoA[i].xs[k] = Atom[j].x;
	  AminoA[i].ys[k] = Atom[j].y;
	  AminoA[i].zs[k] = Atom[j].z;
	  
	}/* end select only the heavy atoms in the side-chain of this a.a.*/
	
      }/* end go over all atoms in this a.a. */
    
    }/* end else, i.e. we deal with other a.a. than GLY */
    
    AminoA[i].length_sc = k;
    
  }/* end cycle i, i.e. go over all a.a. in the chain */
  

#ifdef DEBUG

  for(j = 1; j <= N_resid; j ++){
    fprintf(to, "%5d: %4s %c %c %3d\n", j, AminoA[j].Name, AminoA[j].chain,
	   AminoA[j].SS, AminoA[j].posn);
  }
  fclose(to);

#endif

  /* assign to each amino acid its reference vdw volume */
  for(i = 1 ; i <= N_resid ; i ++){
    AminoA[i].volume = ref_vdw_volume[AminoA[i].ind];
  }

  /* init. the index of segment of hydrophobs for each a.a. (0 = not in a 
     segment) */
  for(i = 1 ; i <= N_resid ; i ++){
    AminoA[i].ind_H = 0;
  }

  /* free the assigned vectors and matrices */
  free(D1); free(D2); free(D3);
  for (i = 1 ; i <= tot_atoms; i++) {
    if(Indic[i] != NULL){
      free(Indic[i]);
      free(Type[i]);
      free(Code[i]);
    }
  }
  free(Code); free(Type);
  free(seq_pos); free(chain_atom); free(atom_posn);

  for (i = 0 ; i < count_h; i++) {
    free(Number1_h[i]);
    free(Number2_h[i]);
  }

  for (i = 0 ; i < count_s ; i++) {
    free(Number1_s[i]);
    free(Number2_s[i]);
  }

  free(Number1_h); free(Number2_h);
  free(chain_h);

  free(Number1_s); free(Number2_s);
  free(chain_s);


}



