#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>

#include "def_param.h"
#define max 500

extern int     string_to_integer();

int read_Protein(char *filename_PRO){

  int    i,j, k, pos;
  char   buffer[max], *s;
  FILE   *fp;
  double rg, rgsq;
  static char **Number1_disulf, **Number2_disulf;
  int    count_disulf;
  char   **Indic, new_indic[5];


  /* ---------------- read info about any disulfide bonds --------------------- */

  count_disulf = 0;
  fp = fopen(filename_PRO,"r");  
  while(fgets(buffer,max,fp) && (buffer[0] != '\n')){
    if(!strncmp(buffer,"SSBOND",6)){
      count_disulf ++;
    }
  }
  
  if(count_disulf != 0){

    /* alocate space for the members of each disulfide bond */
    Number1_disulf = (char **) malloc(count_disulf * sizeof(char *));
    Number2_disulf = (char **) malloc(count_disulf * sizeof(char *));

  }

  rewind(fp);
  
  if(count_disulf != 0){
    
    i = 0;
    while(fgets(buffer,max,fp) && (buffer[0] != '\n')){  
      if(!strncmp(buffer,"SSBOND",6)){
	
	s = buffer;
	
	s += 18;

	Number1_disulf[i] = (char *) malloc(5 * sizeof(char));
	for(pos = 0 ; pos <= 4 ; Number1_disulf[i][pos++] = *s++);
	Number1_disulf[i][pos] = '\0';

#ifdef DEBUG
	for(pos = 0 ; pos <= 5 ; pos ++)
	  printf("Number1[%d] = %c\n", pos, Number1_disulf[i][pos]);
#endif

	s += 9;

	Number2_disulf[i] = (char *) malloc( 5 * sizeof(char));
	for(pos = 0 ; pos <= 4 ; Number2_disulf[i][pos++] = *s++);
	Number2_disulf[i][pos] = '\0';
	
	i ++;

      }
      
    }
  }
  fclose(fp);

  //printf("\n\nfor protein %s\n", filename_PRO);

#ifdef DEBUG
     
  printf("\nfor SSbonds:\n\n");
     
  for(i = 0 ; i < count_disulf ; i ++){
    printf("%5d: (%5s %5s) \n", i, Number1_disulf[i], Number2_disulf[i]);
  }

#endif
  

  /* ---------------- end read the disulfide bonds ----------------- */

  fp = fopen(filename_PRO,"r");  

  i = 0;

  while(fgets(buffer,max,fp) && (buffer[0]!='\n')){

    s = buffer;

    if(!strncmp(s,"ATOM",4)){
    
      s += 13;

      if(!strncmp(s,"CA",2)){
        i ++;
      }
    }

    if(!strncmp(s,"TER",3)){
      break;
    }
  }

  tot_amino = i;
  printf("%d\n",tot_amino);
  rewind(fp);
  
  Amino = malloc((tot_amino+1) * sizeof(atom));
  Indic = (char **) malloc((tot_amino+1)*sizeof(char *));
  
  for(i = 1 ; i <= tot_amino; i++){
   
    Amino[i].x = 0.0;
    Amino[i].y = 0.0;
    Amino[i].z = 0.0;
    Amino[i].posn = 0;
    Amino[i].ssbond_nr = 0;

  }

  printf("begin\n");
  i = 0;
  while(fgets(buffer,max,fp)&& (buffer[0]!='\n')){
   
    s = buffer; 
    
    if(!strncmp(s,"ATOM",4)){
      
      s += 13;

      if(!strncmp(s,"CA",2)){ 

        //printf("%s\n",s);
        s += 10;

	/*while(*s == ' ') s ++;

	printf("print s:%s\n", s);*/

        i ++;

	/* get the position of the residue along the sequence */
	Indic[i] = (char *) malloc(5 * sizeof(char));
	for(pos = 0 ; pos <= 4 ; Indic[i][pos++] = *s++);
	Indic[i][pos] = '\0';
      
	s += 2;

        sscanf(s,"%lf %lf %lf",&Amino[i].x,&Amino[i].y,&Amino[i].z);

        s += 29; 
        sscanf(s,"%lf",&Amino[i].BF);

        //printf("%f %f %f %f\n", Amino[i].x,Amino[i].y,Amino[i].z,Amino[i].BF);

      }  
    }
  }

  /* get the square of the radius of gyration of the protein */
  rgsq = 0.0; 
  for(i = 1 ; i <= (tot_amino-1) ; i ++){
    for(j = (i+1) ; j <= tot_amino ; j ++){

      rgsq = rgsq + (Amino[i].x-Amino[j].x)*(Amino[i].x-Amino[j].x) + 
	(Amino[i].y-Amino[j].y)*(Amino[i].y-Amino[j].y) + (Amino[i].z-Amino[j].z)*(Amino[i].z-Amino[j].z);
    }
  }

  rgsq = rgsq/((double)(tot_amino*tot_amino));
  rg = sqrt(rgsq);

  printf("rg = %f\n", rg);

  /* identify the numbering of the residues as it appears in the PDB file */
  for(i = 1 ; i <= tot_amino ; i ++){
    if(Indic[i] != NULL){
      strncpy(new_indic,Indic[i],4);
      new_indic[4] = '\0';
      Amino[i].posn = string_to_integer(new_indic);
    }
  }

  /* count the number of SSbonds in the protein and determine their 
     starting/ending points */
  SSbond_nr = 0;
  for(i = 0 ; i < count_disulf ; i ++){
    if(Number1_disulf[i] != NULL){   
      SSbond_nr ++;
    }
  }

  /* allocate space for these SSbonds */
  SSbond     = malloc((SSbond_nr+1) * sizeof(sec_struct));

  if(SSbond == NULL){
    printf("malloc failed in SSbond => Exit\n");
    exit(100);
  }

  SSbond_nr = 0;
  for(i = 0 ; i < count_disulf ; i ++){
    if(Number1_disulf[i] != NULL){
      
      SSbond_nr ++;

      strncpy(new_indic,Number1_disulf[i],4);
      new_indic[4] = '\0';

      SSbond[SSbond_nr].start = string_to_integer(new_indic);
        
#ifdef DEBUG
      fprintf(fp,"%5d %s\n", i, Number1_disulf[i]);
      fprintf(fp,"%5d start_disulf = %5d\n", SSbond_nr, SSbond[SSbond_nr].start);
#endif
    
    }

    if(Number2_disulf[i] != NULL){
    
      strncpy(new_indic,Number2_disulf[i],4);
      new_indic[4] = '\0';

      SSbond[SSbond_nr].end = string_to_integer(new_indic);
      
#ifdef DEBUG
      fprintf(fp,"%5d %s\n", i, Number2_disulf[i]);
      fprintf(fp,"%5d end_h = %5d\n", SSbond_nr, SSbond[SSbond_nr].end);
#endif

    }
  }

  /* identify the residues in SSbonds */
  if(SSbond_nr != 0){
  
    for(k = 1 ; k <= tot_amino ; k ++){

      for(j = 1 ; j <= SSbond_nr ; j ++){
	
	if((Amino[k].posn == SSbond[j].start) || 
	   (Amino[k].posn == SSbond[j].end)){
	
	  Amino[k].ssbond_nr = j;

	  //printf("k = %5d Amino.posn = %5d Amino.ssbond = %5d\n", k, 
		// Amino[k].posn, Amino[k].ssbond_nr);
	    
	}
      }
        
    }/* end cycle k, i.e., going over all positions in the chain */

  }/* end there are disulfide bonds in this protein structure */

  /* ------ end identify positions in disulfide bonds ------ */

  return(0);

}
  
  
  

  
  
  




















