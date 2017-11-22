/* this is the program "update_dcd_fullGo_nucl.c" */

/* compared with "update.c" it prints the output data trajectory in
dcd format as well as in the binary file.  Also, it takes into account
the Go-description of the chain based on the set of contacts obtained
from all atoms.  The difference with respect to "update_dcd_fullGo.c"
is that here we considered that the protein has a nucleotide (ATP or
ADP) bound to it. */

#include "def_param.h"


extern void ras_structure();
double distance();
void   gyration();

/*
 * Procedures from dcdio.c to write data into dcd
 */
extern void dcd_write_frame(FILE* dcd_file, int N, float *X, float *Y, float *Z);

/*
 * Variables needed to write dcd
 */
extern FILE* dcd_file;
extern float* X;
extern float* Y;
extern float* Z;


void update()
{

    int    i, j, count;
    /* LJY added rpx, rpy, rpz, rfx, rfy, rfz */
    double dpossum, dx, dy, dz, dr, rpx, rpy, rpz, rfx, rfy, rfz;

    dpossum = 0.0;

    for(i = 1 ; i <= tot_amino ; i ++)
    {
        dpossum = dpossum +
            d_Amino[i].x*d_Amino[i].x +
            d_Amino[i].y*d_Amino[i].y +
            d_Amino[i].z*d_Amino[i].z;
    }

    dpossum = dpossum*zeta/(2.0*h);
    tempav  = tempav + dpossum/(3.0*((double)(tot_amino)));

    /* every "nav" number of steps update the number of native non-bonded contacts left and output
       the energies */
    if(step%nav == 0)
    {
        /* get the value of the current end-to-end distance */
        /* LJY
           R = distance(Amino[tot_amino].x,Amino[tot_amino].y,
           Amino[tot_amino].z,Amino[1].x,Amino[1].y,Amino[1].z);
        */
        /* LJY get the value R of the current distance between pulled and fixed groups */

        rpx = Amino[pulled_bead].x ;
        rpy = Amino[pulled_bead].y ;
        rpz = Amino[pulled_bead].z ;

        rfx = Amino[fixed_bead].x ;
        rfy = Amino[fixed_bead].y ;
        rfz = Amino[fixed_bead].z ;

        R = distance(rpx,rpy,rpz,rfx,rfy,rfz);

        count = 0;
        for(i = 1 ; i <= (tot_amino-1) ; i ++)
        {
            for(j = (i+1) ; j <= tot_amino ; j ++)
            {

                dx = Amino[j].x - Amino[i].x;
                dy = Amino[j].y - Amino[i].y;
                dz = Amino[j].z - Amino[i].z;
                dr = sqrt(dx*dx+dy*dy+dz*dz);

                /* see how many non-bonded contacts are left */
                /* LJY20100318 consider inter contact i and j even if (i-j) <= 2 */

                if(((dr <= R_limit_bond) && (connect[i][j] == 1) &&
                      ( (j-i) > 2 )) ||
                    ((dr <= R_limit_bond) && (connect[i][j] == 1) &&
                     (i <= numf) && (j > numf) ) ||
                    ((fabs(dr - R0[i][j]) <= R_limit) &&
                     (connect[i][j] == 1) &&
                     ( (j-i) > 2 )))
                {
                    link[i][j] = 1;
                    link[j][i] = 1;
                }
                else
                {
                    link[i][j] = 0;
                    link[j][i] = 0;
                }


                if(link[i][j] == 1)
                {
                    count ++;
                }
            }
        }
        count = count * 2;

        gyration();

        tempav = tempav/((double)nav);

        //printf("%f\n",epot_fene);

        /* LJY forcex -> forcer */

        /*fflush(stdout);
          printf("from update: \n");
          printf(" %d %f %f %f %f %f %f %f %f %d %f %f %d\n", step, tempav, R, forcer, epot_LJ_att, epot_LJ_nei,
          epot_LJ_rep, epot_LJ, epot_fene, count, xt, rg, heat);*/

        /* while heating, print out info every 10^4 steps */
        if(heat == 1){

            if(step%10000 == 0){

                fprintf(fdata, "%u %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5d %5.3f %5.3f %d\n",
                        step, tempav, R, forcer, epot_LJ_att, epot_LJ_nei,
                        epot_LJ_rep, epot_LJ, epot_fene, count, xt, rg, heat);
            }

        }

        /* ---------------------------------------------------------
           changed step just above and just below to unsigned int
           --------------------------------------------------------- */

        /* when heating is off, print out info depending on nav1 */
        else{
            if(step%nav1 == 0){
                fprintf(fdata, "%u %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5d %5.3f %5.3f %d\n",
                        step, tempav, R, forcer, epot_LJ_att, epot_LJ_nei,
                        epot_LJ_rep, epot_LJ, epot_fene, count, xt, rg, heat);
                /* LJY added r component of r vector between COMs, f and r components of interest residues */
                fprintf(fdatad, "%u %5.3f %5.3f %5.3f %5.3f %5.3f\n",
                        step, R, forcer, forcerx, forcery, forcerz);

            }

        }

        /* determine the number of native non-bonded contacts for each bead in the chain */
        for(i = 1 ; i <= tot_amino ; i ++){

            no_contact = 0;

            for(j = 1 ; j <= tot_amino ; j ++){

                /* do not take into account the covalently bound neighbors of "i" */
                if( (i == j) || (abs(i-j) == 1)) continue;

                if(link[i][j] == 1)
				{
                    no_contact ++;
                }
            }

            if(step%(1000*nav)==0){
                fprintf(fout_3,"%d ", no_contact);
            }

        }

        if(step%(1000*nav)==0){
            fprintf(fout_3,"\n");
        }

        // if((ncount%nav_record==0)&&(stretch==1)){

        /* write to the output file the current conformation of the chain */

        /* while heating, print out info every 10^4 steps */
        if(heat == 1){

            if(step%10000 == 0){
                for(i = 1 ; i <= tot_amino ; i ++){

                    c.x = Amino[i].x;
                    c.y = Amino[i].y;
                    c.z = Amino[i].z;
                    /* fwrite(&c,sizeof(struct coord),1,fcd); */
                }
            }

        }

        /* when heating is off, print out info depending on nav1 */
        else{

            if(step%nav1 == 0){

                for(i = 1 ; i <= tot_amino ; i ++)
                {

                    c.x = Amino[i].x;
                    c.y = Amino[i].y;
                    c.z = Amino[i].z;
                    /* fwrite(&c,sizeof(struct coord),1,fcd); */
                }
            }/* end step is multiple of nav1 */


            /* print the dcd file */
            if(step%nav_dcd == 0){

                for(i = 0 ; i < tot_amino ; i ++)
                {
                    X[i] = Amino[i+1].x;
                    Y[i] = Amino[i+1].y;
                    Z[i] = Amino[i+1].z;
                }

                dcd_write_frame(dcd_file, tot_amino, X, Y, Z);
            }


        }/* end heating is off */
        tempav = 0.0;

    }/* end step number is a multiple of "nav" */

}

double distance(double x1, double y1, double z1, double x2,double y2,double z2)
{
    double r, dx, dy, dz;

    dx = x1 - x2;
    dy = y1 - y2;
    dz = z1 - z2;
    r = sqrt(dx*dx+dy*dy+dz*dz);
    return(r);
}

void gyration()
{
    int    i,j;
    double rgsq;

    rgsq = 0.0;

    for(i = 1 ; i <= (tot_amino-1) ; i ++)
    {
        for(j = (i+1) ; j <= tot_amino ; j ++)
        {
            rgsq = rgsq + (Amino[i].x-Amino[j].x)*(Amino[i].x-Amino[j].x) +
                (Amino[i].y-Amino[j].y)*(Amino[i].y-Amino[j].y) +
                (Amino[i].z-Amino[j].z)*(Amino[i].z-Amino[j].z);
        }
    }

    rgsq = rgsq/(double)(tot_amino*tot_amino);
    rg = sqrt(rgsq);

}
