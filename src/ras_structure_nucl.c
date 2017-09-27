
/* this is the program "ras_structure_nucl.c" */

/* it writes the conformation of the chain in Rasmol format. This
   contains a nucleotide (ATP or ADP) bound to the protein. */

#include "def_param.h"


extern char Protein_name[5];

void ras_structure(){

  int    i;
  double min_x, min_y, min_z;
  char   filenamestr[30];

  min_x = INF;
  min_y = INF;
  min_z = INF;

  /* obtain the minimum x,y,z coordinates of the beads in the chain */
  for(i = 1 ; i <= tot_amino ; i++){

    if(Amino[i].x <= min_x) min_x = Amino[i].x;
    if(Amino[i].y <= min_y) min_y = Amino[i].y;
    if(Amino[i].z <= min_z) min_z = Amino[i].z;

  }

  /* set the (0,0,0) to the minimum (x,y,z) coordinates of the chain */
  for(i = 1 ; i <= tot_amino ; i++){

    Amino[i].x = scale * (Amino[i].x - min_x);
    Amino[i].y = scale * (Amino[i].y - min_y);
    Amino[i].z = scale * (Amino[i].z - min_z);

  }

  for(i = 1 ; i <= tot_amino ; i ++){

    nr_N[i]    = i;
    order_N[i] = i;
    strcpy(mol_N[i],"CA");
    strcpy(type_N[i],"VAL");
  }

  /* ---------------------------------------------------------
     changed step to unsigned int
     --------------------------------------------------------- */
  sprintf(filenamestr,"Struct_data/%s.structure_%u",Protein_name,step);
  fout1 = fopen(filenamestr,"w");

  /* LJY20100318 this is the protein part (treated as chain A) */
  for(i = 1 ; i <= numf ; i++){
    fprintf(fout1,"ATOM   %4d %2s %5s A %3d %11.3f %7.3f %7.3f\n", nr_N[i],
            mol_N[i], type_N[i], order_N[i], Amino[i].x, Amino[i].y, Amino[i].z);
  }

  /* this is the nucleotide part (treated as chain D) */
  for(i = (numf+1) ; i <= tot_amino ; i++){
    fprintf(fout1,"ATOM   %4d %2s %5s D %3d %11.3f %7.3f %7.3f\n", nr_N[i],
            mol_N[i], type_N[i], order_N[i], Amino[i].x, Amino[i].y, Amino[i].z);
  }

  for(i = 1 ; i <= (numf-1) ; i++){
    fprintf(fout1,"CONECT %4d %4d\n", nr_N[i], nr_N[i+1]);
  }

  for(i = (numf+1) ; i <= (tot_amino-1) ; i++){
    fprintf(fout1,"CONECT %4d %4d\n", nr_N[i], nr_N[i+1]);
  }

  fclose(fout1);

}



