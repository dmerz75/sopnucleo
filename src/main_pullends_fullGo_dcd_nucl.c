/* this is the program main_pullends_fullGo_dcd_nucl.c */

/* compile with gcc -O2 main_pullends_fullGo_dcd_nucl.c read_protein.c rforce.c force_pullends_nucl.c iteration_nucl.c update_dcd_fullGo_nucl.c ras_structure_nucl.c identify.c dcdio.c -lm -Wall -o run_protein_Fene */

/* run command is ./run_protein_Fene file1 <PDB file of input protein> file2 <file with native contacts based on full atomic description of the protein> number <number of trajectory> */

/* run Langevin simulation for a coarse-grained model of a protein
   chain with the chain connectivity modeled using FENE-LJ model. Also,
   an external force is applied on the last bead of the chain under
   constant loading rate conditions. Each aminoacid is represented by its
   CA atom. Here the external force is applied on the C-term end of the
   chain in the x-direction with respect to the 1st bead.  The difference
   with respect to "main_pullends.c" is that here we are printing the
   output data trajectory in dcd format as well. The difference with
   respect to "main_pullends_dcd.c" is that here we are reading the
   random seed for each trajectory from the library file
   "def_random_seeds.h" and therefore as input we need to give both the
   PDB file of the protein and the number of the trajectory (integer
   between 1 and 300). The difference with respect to
   "main_pullends_dcd.v3.c" is that now we base our Go-model description
   on native contacts between both Ca-alpha atoms within a cut-off
   dist. of 8 A and between pairs of heavy atoms in the side-chains of
   the residues within a cut-off distance of 5.2 A (these cut-off are in
   "def_const_PDB.h" and the list of these contacts, together with the
   value of the distance between the Calpha atoms of the residues in the
   starting PDB file are in file2 which is produced by first running the
   program "get_contacts_PDB.c" on the PDB file of the protein). The
   difference with respect to "main_pullends_fullGo_dcd.c" is that here
   we also include nucleotides (ATP or ADP) in the simulation. This
   affects the force part of the program and the output format. Here the
   nucleotide is described at atomistic level (through all its heavy
   atoms), rigid inside, and its movement being the movement of its
   center of mass. The nucleotide needs to be part of the input PDB file
   but listed as an amino acid with ATOM instead of HETATM in the first
   column, each of its heavy atoms listed as CA and assigned a different
   chain name from the protein chain (for example chain D when the
   protein chain is A). Obs: the second file must be the results of the
   concatenation of the "Contact_map_inter_sc_and_b_Protname" and
   "Contact_map_intra_sc_and_b_Protname" (if we use the Full Go
   description) or the "Contact_map_inter_b_Protname" and
   "Contact_map_intra_b_Protname" (if we use the simple Go description)
   obtained using the "get_contacts_PDB_nucl.c" program (with the
   "Only_Calpha" option turned on for the second case). */
/* LJY change into that external force is applied between two residue
   groups in the direction of r vector from fixed center of mass to
   pulled center of mass */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "def_param.h"
#include "def_random_seeds.h"

#define line_max       300  /* max # characters on a line */

extern void read_Protein();
extern double distance();
/* LJY
   void ras_structure(), transform(), rforce(), force(), iteration();
*/
void ras_structure(), rforce(), force(), iteration();

extern int r_seed[301]; /* the value of the random seed for ran2() for the given trajectory */

/*
 * Procedures from dcdio.c to write data into dcd
 */
extern FILE* dcd_open_write(char *dcd_filename);
extern void dcd_write_header(FILE* dcd_file, char *dcd_filename, int N, int NFILE, int NPRIV, int NSAVC, double DELTA);
extern void dcd_write_frame(FILE* dcd_file, int N, float *X, float *Y, float *Z);
extern void dcd_close(FILE* dcd_file);

/*
 * Variables needed to write dcd
 */
FILE* dcd_file;
float* X;
float* Y;
float* Z;

FILE *fp;

char Protein_name[20];

int main(int argc, char *argv[]){

    /* LJY added fdatad, f and r components of interest residues */

    char   buffer[line_max], *str, *s, filename[100], filename_data[100], filename_datad[100];
    int    i, j, k, pos, count, nr_pairs, size_Euclid, traj_nr;
    int    *aa1, *aa2;
    /* LJY  added rpx, rpy, rpz, rfx, rfy, rfz */
    double dx, dy, dz, dr, dist_aa, rpx, rpy, rpz, rfx, rfy, rfz;

    double eh_temp; // temporary eh read from topology;
    int midtop; // middle integer 1 in the gsop topology;


    str = argv[1];

    /* read the name of the protein (i.e., the PDB id) */
    pos = 0;
    str += 3;
    //printf("%s\n",str1);
    while(*str != '.'){
        Protein_name[pos] = *str;
        str ++;
        pos ++;
    }
    Protein_name[pos] = '\0';


    str = argv[1];

    /* read the coordinates of the residues in the given protein structure */
    read_Protein(str);

    tempav = 0.0;


    nr_N    = (int*)malloc((tot_amino+1) * sizeof(int));
    order_N = (int*)malloc((tot_amino+1) * sizeof(int));

    mol_N  = (char **)malloc((tot_amino+1) * sizeof(char*));
    type_N = (char **)malloc((tot_amino+1) * sizeof(char*));

    for(i = 1 ; i <= tot_amino ; i ++)
    {
        mol_N[i]  = (char *)malloc((tot_amino+1) * sizeof(char));
        type_N[i] = (char *)malloc((tot_amino+1) * sizeof(char));
    }

    printf("total amino acid check+1: %d\n",tot_amino+1);
    connect = malloc((tot_amino+1) * sizeof(int *));
    connect_eh = malloc((tot_amino+1) * sizeof(double *));
    link    = malloc((tot_amino+1) * sizeof(int *));
    R0      = malloc((tot_amino+1) * sizeof(double *));

    for(i = 1 ; i <= tot_amino ; i ++)
    {
        connect[i] = malloc((tot_amino+1) * sizeof(int));
        connect_eh[i] = malloc((tot_amino+1) * sizeof(double));
        link[i]    = malloc((tot_amino+1) * sizeof(int));
        R0[i]      = malloc((tot_amino+1) * sizeof(double));
    }


    /* connect_eh added here to account for variable eh; */
    /* --Dale, September 27, 2017 */
    for(i = 1 ; i <= tot_amino ; i ++)
    {
        for(j = 1 ; j <= tot_amino ; j ++)
        {
            connect[i][j] = 0;
            connect_eh[i][j] = eh;
            link[i][j] = 0;
            R0[i][j] = 0.;
        }
    }

    d_Amino = malloc((tot_amino+1) * sizeof(atom));

    /* read the file with the information about the positions in contact
       and the initial distance between their C-alpha atoms, i.e., the
       contact map file */
    fp = fopen(argv[2],"r");


    /* Open topology, do line counting, rewind */
    count = 0;
    while(fgets(buffer,line_max,fp) && (buffer[0] != '\n')){
        count ++;
    }
    rewind(fp);
    printf("number of pairs in contact = %5d\n", count);


    /* assign space and initalize the entries that are in contact in the original PDB file */
    aa1 = malloc((count+1) * sizeof(int));
    aa2 = malloc((count+1) * sizeof(int));

    for(i = 0; i <= count; i ++){
        aa1[i] = 0;
        aa2[i] = 0;
    }

    nr_pairs = 0;

    while(fgets(buffer,line_max,fp) &&
          (buffer[0] != '\n'))
        /* && */
          /* (buffer[0] != '#')) */
    {
        if(buffer[0] == '#')
        {
            continue;
        }

        s = buffer;
        dist_aa = 0.;
        /* sscanf(s,"%d %d %lf", &aa1[nr_pairs+1], &aa2[nr_pairs+1], &dist_aa); */
        // emolsion based topology. -- github:dmerz75/emolsion.git
        /* 8    11     1    5.13938    1.25000 */
        /* 8    12     1    4.90770    1.25000 */
        /* 8    13     1    7.75240    1.25000 */

        sscanf(s,"%d %d %d %lf %lf", &aa1[nr_pairs+1],
               &aa2[nr_pairs+1], &midtop, &dist_aa, &eh_temp);

        /* Crucial eh assignment! */
        /* printf(">>  %d %d %d %lf %lf\n",aa1[nr_pairs+1],aa2[nr_pairs+1], */
        /*        midtop,dist_aa,eh_temp); */

        /* printf("eh(1)> %f\n",connect_eh[aa1[nr_pairs+1]+1][aa2[nr_pairs+1]+1]); */
        connect_eh[aa1[nr_pairs+1]+1][aa2[nr_pairs+1]+1] = eh_temp;
        /* printf("eh(2)> %f\n",connect_eh[aa1[nr_pairs+1]+1][aa2[nr_pairs+1]+1]); */

        /* connect_eh[aa1[nr_pairs]][aa2[nr_pairs]] = eh_temp; */

        nr_pairs ++;
    }


    /* The total number of contacts in the topology */
    printf("nr_pairs = %5d\n", nr_pairs);
    fclose(fp);


/* #ifdef EH_REDEFINE */

/*     double array_eh[nr_pairs]; */


/*     for (i = 0; i < nr_pairs; i++) */
/*     { */
/*         printf("nr_pairs: (i):%d    aa1-aa2| %d %d\n",i,aa1[i],aa2[i]); */
/*         array_eh[i] = 1.0; */
/*         /\* std::cout << i << nr_pairs[i] << std::endl; *\/ */
/*     } */
/*     /\* exit(0); *\/ */
/* #endif // EH_REDEFINE */



    /* read the trajectory number */

    traj_nr = atoi(argv[3]);

    /* set the seed of the random number generetor for the current trajectory */
    mseed = r_seed[traj_nr];

    printf("mseed = %5d\n", mseed);

    // exit(100);

    xt = 0.0;

    for(i = 1 ; i <= (tot_amino-1) ; i ++){
        for(j = (i+1) ; j <= tot_amino ; j ++){
            connect[i][j] = 0;
        }
    }

    /* obtain the Euclidean distance between any pair of aminoacids */
    for(i = 1 ; i <= (tot_amino-1) ; i ++){

        for(j = (i+1) ; j <= tot_amino ; j ++){

            //printf("inside the two for cycles:\n\n");

            dx = Amino[j].x - Amino[i].x;
            dy = Amino[j].y - Amino[i].y;
            dz = Amino[j].z - Amino[i].z;
            dr = sqrt(dx*dx+dy*dy+dz*dz);

            /* build the connectivity matrix which tells us which pairs of residues are in bonded or non-bonded
               contact */
            /* declare them in contact based on the data from
               the input file, file2 */

            /* LJY20100318 obtain contact from file2 */
            for(k = 1 ; k <= nr_pairs ; k ++)
            {
                if((i == aa1[k]) && (j == aa2[k]))
                {
                    connect[i][j] = 1;
                    //printf("%d %d\n",i,j);
                }
            }

            /* LJY20100318 obtain contact within the protein and within the nucleotide from file1 */

            if( (i <= numf) && (j <= numf) && ((j-i) <= 2) && (dr <= R_limit_bond) )
            {
                connect[i][j] = 1;
            }
            if( (i  > numf) && (j  > numf) && ((j-i) <= 2) && (dr <= R_limit_bond) )
            {
                connect[i][j] = 1;
            }
        }
    }

    for(i = 1 ; i <= (tot_amino-1) ; i ++){
        for(j = (i+1) ; j <= tot_amino ; j ++){

            if(connect[i][j] == 1){

                dx = Amino[i].x - Amino[j].x;
                dy = Amino[i].y - Amino[j].y;
                dz = Amino[i].z - Amino[j].z;

                R0[i][j] = sqrt(dx*dx+dy*dy+dz*dz);

                //printf("%d %d %f\n",i,j,R0[i][j]);

            }

            /* LJY20100318 discriminate chain1, chain2, and chain1-chain2 */
            else{
                if( (i <= numf) && (j <= numf) ){
                    R0[i][j] = a;
                }
                if( (i > numf) && (j > numf) ){
                    R0[i][j] = as;
                }
                if( (i <= numf) && (j > numf) ){
                    R0[i][j] = ai;
                }
            }

            /* this is changed as compared to the runs from July-August 2005 */
            /* LJY20100318 discriminate chain1 and chain2 */

            if((j-i)==2){
                if( (i <= numf) && (j <= numf) ){
                    R0[i][j] = a/half;
                }
                if( (i > numf) && (j > numf) ){
                    R0[i][j] = as/half;
                }
            }
        }
    }

    //exit(100);

    step = 0;

    /* turn on the random force */
    rforce();

    /* calculate the force that acts on the protein */
    force();

    /* heating is turned on: heat from 0 K to "temp" under zero force */
    heat = 1;

    /* LJY change forcex -> forcer */
    forcer = 0.0;
    forcerx = 0.0;
    forcery = 0.0;
    forcerz = 0.0;

    sprintf(filename,"contact_record_%s_%4.3f_%5.4f_%4.3f.dat",Protein_name,k_trans,deltax,h);
    fout_3=fopen(filename,"w");

    /* sprintf(filename_b,"Coord/coord_%s_%4.3f_%5.4f_%4.3f.bd",Protein_name,k_trans,deltax,h); */
    /* fcd=fopen(filename_b,"w"); */

    sprintf(filename_data,"out%s_%4.3f_%5.4f_%4.3f.dat",Protein_name,k_trans,deltax,h);
    fdata=fopen(filename_data,"w");

    /* LJY added r component of r vector between COMs, f and r components of interest residues */

    sprintf(filename_datad,"out%s_%4.3f_%5.4f_%4.3f_d.dat",Protein_name,k_trans,deltax,h);
    fdatad=fopen(filename_datad,"w");

    while(step <= stepsim1){

        iteration();

        step ++;

    }

    /* heating is off; apply force at constant temperature */
    heat = 0;

    /* LJY  transform(); */

    /* obtain the starting value of the end-to-end distance */
    /* LJY
       x0 = distance(Amino[tot_amino].x,Amino[tot_amino].y,Amino[tot_amino].z,Amino[1].x,Amino[1].y,Amino[1].z);

       printf("at begin pulling x0 = %5.3f\n", x0);
    */

    /* LJY obtain the starting value r0 of the distance between pulled and fixed groups */
    rpx = Amino[pulled_bead].x ;
    rpy = Amino[pulled_bead].y ;
    rpz = Amino[pulled_bead].z ;

    rfx = Amino[fixed_bead].x ;
    rfy = Amino[fixed_bead].y ;
    rfz = Amino[fixed_bead].z ;

    r0 = distance(rpx,rpy,rpz,rfx,rfy,rfz);
    printf("at begin pulling r0 = %5.3f\n", r0);


    ras_structure();

    step = 0;  /* total number of steps for pulling */
    stepx = 0; /* index for number of steps in pulling varying between 1 and "nav" (i.e., up to the point were
                  we output the chain conformation) */


    // Open DCD
    /* before the first step in pulling, initialize the dcd file */
    char* dcd_filename = (char*)malloc(sizeof(Protein_name) + 9);
    sprintf(dcd_filename, "Coord/%s.dcd", Protein_name);


    //printf("Coordinates will be saved as '%s'.\n", dcd_filename);
    dcd_file = dcd_open_write(dcd_filename);
    dcd_write_header(dcd_file, dcd_filename, tot_amino, stepsim2, nav, nav, h);
	// CHECK HERE..
	/* printf("fn: %s   total_amino: %d  stepsim2: %d   nav: %d   h:%lf\n", */
	/* 	   dcd_filename, tot_amino, stepsim2, nav, nav, h); */


    // Allocate memory for temporary data
    size_Euclid = (tot_amino+1)*sizeof(float);

    X = (float*) malloc(size_Euclid);
    Y = (float*) malloc(size_Euclid);
    Z = (float*) malloc(size_Euclid);


    /* while( R <= (tot_amino * R0[1][2] * scale_R) ){ */
    while (step <= stepsim2)
    {    /* LJY  */
        iteration();
        step ++;
        stepx ++;
    }

    fclose(fout_3);
    /* fclose(fcd); */
    dcd_close(dcd_file);

    fclose(fdata);

    /* LJY added fdatad, f and r components of interest residues */

    fclose(fdatad);

    exit(0);

}

/* Euler's angles rotation */
/* LJY
   void transform(){

   int i;
   double min_x = 0.0, min_y = 0.0, min_z = 0.0, cosphi, sinphi;


   Amino_new  = malloc((tot_amino+1) * sizeof(atom));

   min_x = Amino[1].x;
   min_y = Amino[1].y;
   min_z = Amino[1].z;

   for(i = 1 ; i <= tot_amino ; i ++){
   Amino[i].x = Amino[i].x - min_x;
   Amino[i].y = Amino[i].y - min_y;
   Amino[i].z = Amino[i].z - min_z;
   }

   cosphi = Amino[tot_amino].z/sqrt(Amino[tot_amino].y*Amino[tot_amino].y+Amino[tot_amino].z*Amino[tot_amino].z);
   sinphi = Amino[tot_amino].y/sqrt(Amino[tot_amino].y*Amino[tot_amino].y+Amino[tot_amino].z*Amino[tot_amino].z);

   for(i = 1 ; i <= tot_amino ; i ++){
   Amino_new[i].x = Amino[i].x;
   Amino_new[i].y = cosphi * Amino[i].y - sinphi * Amino[i].z;
   Amino_new[i].z = sinphi * Amino[i].y + cosphi * Amino[i].z;
   }

   cosphi =
   Amino_new[tot_amino].x/sqrt(Amino_new[tot_amino].x*Amino_new[tot_amino].x+Amino_new[tot_amino].z*Amino_new[tot_amino].z);

   sinphi =
   Amino_new[tot_amino].z/sqrt(Amino_new[tot_amino].x*Amino_new[tot_amino].x+Amino_new[tot_amino].z*Amino_new[tot_amino].z);

   for(i = 1 ; i <= tot_amino ; i ++){
   Amino[i].x = cosphi * Amino_new[i].x + sinphi * Amino_new[i].z;
   Amino[i].y = Amino_new[i].y;
   Amino[i].z = -sinphi * Amino_new[i].x + cosphi * Amino_new[i].z;
   }

   free(Amino_new);

   }
*/
