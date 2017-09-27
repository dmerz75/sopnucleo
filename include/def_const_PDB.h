
               /* this is the file def_const_PDB.h */

/* constants referring to reading a PDB file */

#define L_max         10000  /* max chain length (nr. of residues) */
#define Atom_max      60000  /* max # atoms in a protein */
#define Chain_max        30  /* max # chains in a protein */
#define n_max            30  /* max # heavy atoms in an a.a. (it's # for TRP)*/
#define tot_atom         36  /* # of types of atoms in a protein */
#define Helix_max       600  /* max # helices in a protein */
#define Strand_max      600  /* max # strands in a protein */
#define hydro_segm_max 5000  /* max # segm. of hydrophobs in a seq. */
#define N_class          20  /* # of classes of a.a. */
#define N_charged         4  /* # classes of a.a. that are charged */ 
#define N_class_3         9  /* # classes of a.a. for 3-body inter. */
#define Nparam_3   (int)((N_class_3*(N_class_3+1)*(N_class_3+2))/6)     
                            /* total # 3body param. */
#define eps          1.0e-7  /* tol for determining a zero in doubles */
#define Dss           2.0  /*5.2*/  /* cut-off distance for sc-sc contacts */
#define Dbs          5.5/*5.2*//*cut-off distance for b-s 2-body contacts */
#define Dss_3        5.2/*4.5*//* cut-off distance for 3-body sc-sc contacts */
#define Dbb          6.5/*6.5*/  /* cut-off distance for 2-body b-b contacts */
#define DCA          8.0 /* cut-off distance for 2body Calpha-Calpha contacts */
#define N_pot           210  /* # of 2-body energy parameters */
#define N_hydro          10  /* # of hydrophobic parameters */
#define L_hel             4  /* min length of a helix */
#define min_length       20  /* min length for a test protein */
#define Segm_max         40  /* max length of segm of hydrophobs */
#define r_short        0.20  /* max percentage seq. sep. for short-range */
#define Short_sep        10  /* max seq. sep. for a short-range contact */
#define max_val        0.06  /* max value of access. area for buried residue */
#define Max_helix_length  100 /* max length of helical fragment from protein */
#define Max_strand_length 100 /* max length of strand fragment from protein */
#define min_exposed_ratio 0.40 /* min ratio of access. for an exposed resid. */
#define dist_hbonds     3.2   /* cut-off distance between N and O involved in a hydrogen bond */
#define ratio_HBond_min -1.00 /* min value of cos angle between HN and HO vectors (corresp. to 180) forming a H bond */
#define ratio_HBond_max -0.76 /* max value of the same cosine (corresp. to angle of 140 degrees) */


#define true              1
#define false             0
