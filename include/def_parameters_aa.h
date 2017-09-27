
             /* this is the file "def_parameters_aa.h" */

/* it contains various known numbers (accessible area, volume etc.) for the
20 types of a.a. */

#include "def_const_PDB.h"

/* area of aminod acid X in a GLY-X-GLY extended conformation (Creighton) */
double ref_area_sc[N_class+1] = {0.00, 140., 218., 180., 259., 160., 182, 204., 194., 229., 113., 85., 143., 158., 146., 122., 241., 189., 151., 211., 183.};

/* vdw volume of amino acids inside a protein (Creighton) */
double ref_vdw_volume[N_class+1] = {0.00, 86., 135., 124., 163., 105., 124., 124., 118., 141., 67., 48., 90., 96., 93., 73., 148., 114., 91., 135., 109.};

/* hydrophobicity of each type of a.a. acc. to Eisenberg */
double Hydro_Eisen[N_class+1] = {0.00, 0.29, 1.19, 1.06, 0.81, 1.08, 1.38, 0.64, -0.40, 0.26, 0.62, 0.48, 0.12, -0.78, -0.05, -0.18, -2.53, -0.85, -0.90, -1.50, -0.74};

/* 3-letter symbol for a.a. */
char symbol[N_class+1][4] = {"NOA", "CYS", "PHE", "LEU", "TRP", "VAL", "ILE", "MET", "HIS", "TYR", "ALA", "GLY", "PRO", "ASN", "THR", "SER", "ARG", "GLN", "ASP", "LYS", "GLU"};

/* 1-letter symbol for a.a. */
char symbol1[N_class+1] = {'Z', 'C', 'F', 'L', 'W', 'V', 'I', 'M', 'H', 'Y', 'A', 'G', 'P', 'N', 'T', 'S', 'R', 'Q', 'D', 'K', 'E'};

/* standard nr. heavy atoms in each type of residue */
int number_atom[N_class+1] = {0, 6, 11, 8, 14, 7, 8, 8, 10, 12, 5, 4, 7, 8, 7, 6, 11, 9, 8, 9, 9}; 

/* classification of amino acids in 2 classes: hydrophobs (1) and polars (0) */
int hydro_res2[N_class+1] = {-1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0};

/* classification of a.a. in 4 classes: H(1), P(0), +(2) and -(3) */
int hydro_res4[N_class+1] = {-1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,2,0,3,2,3};

/* introduce the numbering for the charged a.a. (i.e. here, compared to hydro_res4, Asp and Glu are 
   treated independently and Arg and Lys are treated independently) */
int charge_res[N_class+1] = {-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,2,3,4};

/* classification of a.a in 9 classes (used in determining 3-body contacts) */
int classes9_res[N_class+1] = {0,3,5,4,5,4,4,4,7,5,4,1,2,7,6,6,8,7,9,8,9};

/* frequency of finding each of the 20 types of a.a. in a typical protein sequence (acc. to Creighton) */
double freq_aa[N_class+1] = {0.00, 0.017, 0.039, 0.090, 0.013, 0.066, 0.052, 0.024, 0.022, 0.032, 0.083, 0.072, 0.051, 0.044, 0.058, 0.069, 0.057, 0.040, 0.053, 0.057, 0.062};

