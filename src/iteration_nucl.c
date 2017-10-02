/* this is the program "iteration_nucl.c" */

/* the difference with respect to "iteration.c" is that here we
   consider a nucleotide (ATP or ADP) bound to the protein. */

#include "def_param.h"


extern void rforce(), force(), update();

void iteration(){

  int i;

  for(i = 1 ; i <= tot_amino ; i ++){

    /* LJY20100318 make the nucleotide rigid but its CM moves */
    if( (rigid == 2)  && (i >= (numf+1)) ){

      d_Amino[i].x = (h * Amino[i].fcx)/zetai;
      d_Amino[i].y = (h * Amino[i].fcy)/zetai;
      d_Amino[i].z = (h * Amino[i].fcz)/zetai;
    }

    /* all the other positions can move individually */
    else {

      d_Amino[i].x = (h * Amino[i].fcx)/zeta;
      d_Amino[i].y = (h * Amino[i].fcy)/zeta;
      d_Amino[i].z = (h * Amino[i].fcz)/zeta;
    }

    /* change of position for each bead at each time step due to the action of force */
    Amino[i].x = Amino[i].x + d_Amino[i].x;
    Amino[i].y = Amino[i].y + d_Amino[i].y;
    Amino[i].z = Amino[i].z + d_Amino[i].z;

#ifdef DEBUG

    if( (heat == 0) && (i == 31) ){
      printf("\n from iteration at i = %3d: \n", i);
      printf("dx = %f dy = %f dz = %f\n", d_Amino[i].x, d_Amino[i].y, d_Amino[i].z);
      printf("x = %f y = %f z = %f\n", Amino[i].x, Amino[i].y, Amino[i].z);
    }

#endif

  }

  rforce();

  force();


  /* after each "nav" number of steps, increase the value of "xt" (i.e., the external force is ramped up
     each "nav" number of steps) */
  if(stepx == nav){

    xt = xt + deltax;
    stepx = 0;

  }

  update();

}
