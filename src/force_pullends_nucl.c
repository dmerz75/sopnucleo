
/* this is the program "force_pullends_nucl.c": it calculates the force that acts on the chain at each moment in time */

/* the difference with respect to "force_pullends.c" is that here we
   consider a protein and a nucleotide bound to it (ATP or ADP). The
   nucleotide is considered rigid inside and its movement is described
   by the movement of its center of mass. */


#include "def_param.h"

void fene(),tension(),LJ();


void force(){
  int    i;
  double fsx, fsy, fsz;

  epot = 0.0;

  LJ();
  fene();

  if(heat == 0){
   tension();
  }

  /* LJY20100318 declare the nucleotide as rigid and static */
  if(rigid == 0){

    for(i = (numf+1) ; i <= tot_amino ; i ++){
      Amino[i].fcx = 0.0;
      Amino[i].fcy = 0.0;
      Amino[i].fcz = 0.0;
    }

  }

  /* LJY20100318 make the nucleotide rigid but the CM (geometric center) moves */
  if(rigid == 2){

    fsx = 0.0;
    fsy = 0.0;
    fsz = 0.0;
    
    for(i = (numf+1) ; i <= tot_amino ; i ++){
      fsx = fsx + Amino[i].fcx;
      fsy = fsy + Amino[i].fcy;
      fsz = fsz + Amino[i].fcz;
    }

    fsx = fsx/(tot_amino - numf);
    fsy = fsy/(tot_amino - numf);
    fsz = fsz/(tot_amino - numf);
    
    for(i = (numf+1) ; i <= tot_amino ; i ++){
      Amino[i].fcx = fsx;
      Amino[i].fcy = fsy;
      Amino[i].fcz = fsz;
    }
  }

}

/* get the Lennard-Jones contribution */
void LJ(){

  double dx, dy, dz, rsqi, rsqi6, rsqi10, rsqi12, p12=0.0, p6=0.0, p6_rep=0.0, pot, blub, dphi, fx, fy, fz;
  int    i, j;

  epot_LJ = 0.0;
  epot_LJ_att = 0.0;
  epot_LJ_nei = 0.0;
  epot_LJ_rep = 0.0;

  /* LJY20100318 within the protein set up the attractive and repulsive LJ entries */

  for(i = 1 ; i <= (numf-1) ; i ++){
  
    for(j = (i+1) ; j <= numf ; j ++){
      
      dx = Amino[j].x - Amino[i].x;
      dy = Amino[j].y - Amino[i].y;
      dz = Amino[j].z - Amino[i].z;

      rsqi   = (R0[i][j] * R0[i][j])/(dx*dx+dy*dy+dz*dz);
      
      /* if we want to speed up the computation, introduce a step function
      to cut-off the LJ potential if the new distance between positions 
      i and j becomes more than R_sigma bigger than their original 
      distance in the PDB structure */

      if(rsqi < 1./R_sigma) continue;		

      rsqi6  = rsqi*rsqi*rsqi;
      rsqi10 = rsqi6*rsqi*rsqi;
      rsqi12 = rsqi6*rsqi6;
      
      /* for pair of beads that form a native non-bonded contact */
      if( ( connect[i][j] == 1 ) && ( (j-i) > 2 ) ){
        p12 = eh * rsqi12;
        p6  = -2.0 * eh * rsqi6;
        pot = p12 + p6;
        epot_LJ_att = epot_LJ_att + pot;
      }
      /* for pair of beads that occur consecutively along the chain apply repulsive LJ potential 
	 to prevent overlaps (because FENE potential does not take care of these) */
      else if( ( connect[i][j] == 1 ) && ( (j-i) == 1 ) ){
	p6_rep = LJ_bond * el * rsqi6;
        pot    = p6_rep;
        epot_LJ_nei = epot_LJ_nei + pot;
      
      }
      /* for pairs of beads that are not covalently bound or form a native non-bonded contact,
	 consider only the repulsive LJ potential */
      else{
        p6_rep = el * rsqi6;
        pot    = p6_rep;
        epot_LJ_rep = epot_LJ_rep + pot;
      }
      blub = 1.0/rsqi;
      if(blub < eps) continue;

      /* for beads in native non-bonded contact, derivative of full LJ pot. */
      if( ( connect[i][j] == 1 ) && ( (j-i) > 2 ) ){
	dphi = rsqi * (12.0*p12 + 6.0*p6)/(R0[i][j]*R0[i][j]);
      }

      /* for all other pairs of beads, derivative of only the repulsive part of the LJ pot. */
      else{
        dphi = rsqi*(6.0*p6_rep)/(R0[i][j]*R0[i][j]);
      }

      fx = dx * dphi;
      fy = dy * dphi;
      fz = dz * dphi;

      Amino[i].fcx = Amino[i].fcx - fx;
      Amino[i].fcy = Amino[i].fcy - fy;
      Amino[i].fcz = Amino[i].fcz - fz;
      Amino[j].fcx = Amino[j].fcx + fx;
      Amino[j].fcy = Amino[j].fcy + fy;
      Amino[j].fcz = Amino[j].fcz + fz;

      epot_LJ = epot_LJ + pot;

    }
  }

  /* LJY20100318 within the nucleotide if it is flexible and moving */
  if(rigid == 1){
    
    for(i = (numf+1) ; i <= (tot_amino-1) ; i ++){
      
      for(j = (i+1) ; j <= tot_amino ; j ++){

	dx = Amino[j].x - Amino[i].x;
	dy = Amino[j].y - Amino[i].y;
	dz = Amino[j].z - Amino[i].z;

	rsqi   = (R0[i][j] * R0[i][j])/(dx*dx+dy*dy+dz*dz);

	/* if we want to speed up the computation, introduce a step function
	   to cut-off the LJ potential if the new distance between positions 
	   i and j becomes more than R_sigma bigger than their original 
	   distance in the PDB structure */

      if(rsqi < 1./R_sigma) continue;

      rsqi6  = rsqi*rsqi*rsqi;
      rsqi10 = rsqi6*rsqi*rsqi;
      rsqi12 = rsqi6*rsqi6;

      /* for pair of beads that form a native non-bonded contact */
      if( ( connect[i][j] == 1 ) && ( (j-i) > 2 ) ){
        p12 = ehs * rsqi12;
        p6  = -2.0 * ehs * rsqi6;
        pot = p12 + p6;
        epot_LJ_att = epot_LJ_att + pot;
      }

      /* for pair of beads that occur consecutively along the chain apply repulsive LJ potential 
	 to prevent overlaps (because FENE potential does not take care of these) */
      else if( ( connect[i][j] == 1 ) && ( (j-i) == 1 ) ){
        p6_rep = LJ_bond * elss * rsqi6;
        pot    = p6_rep ;
        epot_LJ_nei = epot_LJ_nei + pot;

      }

      /* for pairs of beads that are not covalently bound or form a native non-bonded contact,
	 consider only the repulsive LJ potential */
      else{
        p6_rep = elss * rsqi6;
        pot    = p6_rep;
        epot_LJ_rep = epot_LJ_rep + pot;
      }
      blub = 1.0/rsqi;
      if(blub < eps) continue;

      /* for beads in native non-bonded contact, derivative of full LJ pot. */
      if( ( connect[i][j] == 1 ) && ( (j-i) > 2 ) ){
        dphi = rsqi * (12.0*p12 + 6.0*p6)/(R0[i][j]*R0[i][j]);
      }

      /* for all other pairs of beads, derivative of only the repulsive part of the LJ pot. */
      else{
        dphi = rsqi*(6.0*p6_rep)/(R0[i][j]*R0[i][j]);
      }

      fx = dx * dphi;
      fy = dy * dphi;
      fz = dz * dphi;

      Amino[i].fcx = Amino[i].fcx - fx;
      Amino[i].fcy = Amino[i].fcy - fy;
      Amino[i].fcz = Amino[i].fcz - fz;
      Amino[j].fcx = Amino[j].fcx + fx;
      Amino[j].fcy = Amino[j].fcy + fy;
      Amino[j].fcz = Amino[j].fcz + fz;

      epot_LJ = epot_LJ + pot;

      }
    }
  } /* end the nucleotide is treated as a flexible and moving molecule */


  /* LJY20100318 interactions between the protein and the nucleotide */
  /* i is a position in the protein */
  for(i = 1 ; i <= numf ; i ++){

    /* j is a position in the nucleotide */
    for(j = (numf+1) ; j <= tot_amino ; j ++){

      dx = Amino[j].x - Amino[i].x;
      dy = Amino[j].y - Amino[i].y;
      dz = Amino[j].z - Amino[i].z;

      rsqi   = (R0[i][j] * R0[i][j])/(dx*dx+dy*dy+dz*dz);

      /* if we want to speed up the computation, introduce a step function
	 to cut-off the LJ potential if the new distance between positions 
	 i and j becomes more than R_sigma bigger than their original 
	 distance in the PDB structure */
      if(rsqi < 1./R_sigma) continue;
      
      rsqi6  = rsqi*rsqi*rsqi;
      rsqi10 = rsqi6*rsqi*rsqi;
      rsqi12 = rsqi6*rsqi6;

      /* for pair of beads that form a native non-bonded contact */
      if( connect[i][j] == 1  ){
        p12 = ehi * rsqi12;
        p6  = -2.0 * ehi * rsqi6;
        pot = p12 + p6;
        epot_LJ_att = epot_LJ_att + pot;
      }
      /* for pairs of beads that are not covalently bound or form a native non-bonded contact,
	 consider only the repulsive LJ potential */
      else{
        p6_rep = eli * rsqi6;
        pot    = p6_rep;
        epot_LJ_rep = epot_LJ_rep + pot;
      }
      blub = 1.0/rsqi;
      if(blub < eps) continue;

      /* for beads in native non-bonded contact, derivative of full LJ pot. */
      if( connect[i][j] == 1  ){
        dphi = rsqi * (12.0*p12 + 6.0*p6)/(R0[i][j]*R0[i][j]);
      }

      /* for all other pairs of beads, derivative of only the repulsive part of the LJ pot. */
      else{
        dphi = rsqi*(6.0*p6_rep)/(R0[i][j]*R0[i][j]);
      }

      fx = dx * dphi;
      fy = dy * dphi;
      fz = dz * dphi;

      Amino[i].fcx = Amino[i].fcx - fx;
      Amino[i].fcy = Amino[i].fcy - fy;
      Amino[i].fcz = Amino[i].fcz - fz;
      Amino[j].fcx = Amino[j].fcx + fx;
      Amino[j].fcy = Amino[j].fcy + fy;
      Amino[j].fcz = Amino[j].fcz + fz;

      epot_LJ = epot_LJ + pot;

    }
  }

}/* end the LJ part of the force */


/* get FENE potential */
void fene(){

  int    i;
  double dx,dy,dz,pot,fx,fy,fz,dr,dR;
 
  epot_fene = 0.0;

  /*
    for(i=1;i<=tot_amino-1;i++){
    printf("R0 %d %f\n",i,R0[i][i+1]);
    }
  */

  /* LJY20100318 within the protein */
  for(i = 1 ; i <= (numf-1) ; i ++){

      dx = Amino[i+1].x - Amino[i].x;
      dy = Amino[i+1].y - Amino[i].y;
      dz = Amino[i+1].z - Amino[i].z;

      dr = sqrt(dx*dx+dy*dy+dz*dz);
      dR = (dr - R0[i][i+1]);

      pot = - (kspring_cov/2.0) * (R_limit*R_limit) * log(1. - dR*dR/(R_limit*R_limit));
      
      if(dR > R_limit){
          /* printf("%d %d %d %f %f %f %f %f %f %f %f %f\n", step, i, i+1, dR, R_limit, pot, */
          printf("%u %d %d %f %f %f %f %f %f %f %f %f\n", step, i, i+1, dR, R_limit, pot,
                 d_Amino[i].x, d_Amino[i].y, d_Amino[i].z, d_Amino[i+1].x, d_Amino[i+1].y, d_Amino[i+1].z);
          exit(0);
      }
                 
      epot_fene = epot_fene + pot;

      /*if(step>14200){
        printf("%d dR %d %f %f %f\n",step,i,dR,pot,epot_fene);
      }
      */

      /* the corresponding force */
      fx = -kspring_cov*dx*dR/dr/(1. - dR*dR/(R_limit*R_limit));
      fy = -kspring_cov*dy*dR/dr/(1. - dR*dR/(R_limit*R_limit));
      fz = -kspring_cov*dz*dR/dr/(1. - dR*dR/(R_limit*R_limit));
     
      Amino[i].fcx = Amino[i].fcx - fx;
      Amino[i].fcy = Amino[i].fcy - fy;
      Amino[i].fcz = Amino[i].fcz - fz;
 
      Amino[i+1].fcx = Amino[i+1].fcx + fx;
      Amino[i+1].fcy = Amino[i+1].fcy + fy;
      Amino[i+1].fcz = Amino[i+1].fcz + fz;
  }

  /* LJY20100318 within the nucleotide treated as flexible and moving */
  if(rigid == 1){

    for(i = (numf+1) ; i <= (tot_amino-1) ; i ++){

      dx = Amino[i+1].x - Amino[i].x;
      dy = Amino[i+1].y - Amino[i].y;
      dz = Amino[i+1].z - Amino[i].z;

      dr = sqrt(dx*dx+dy*dy+dz*dz);
      dR = (dr - R0[i][i+1]);

      pot = - (kspring_covs/2.0) * (R_limits*R_limits) * log(1. - dR*dR/(R_limits*R_limits));
      
      if(dR > R_limits){
          /* printf("%d %d %d %f %f %f %f %f %f %f %f %f\n", step, i, i+1, dR, R_limits, pot, */
          printf("%u %d %d %f %f %f %f %f %f %f %f %f\n", step, i, i+1, dR, R_limits, pot,
                 d_Amino[i].x, d_Amino[i].y, d_Amino[i].z, d_Amino[i+1].x, d_Amino[i+1].y, d_Amino[i+1].z);
          exit(0);
      }

      epot_fene = epot_fene + pot;

      /* the corresponding force */
      fx = -kspring_covs*dx*dR/dr/(1. - dR*dR/(R_limits*R_limits));
      fy = -kspring_covs*dy*dR/dr/(1. - dR*dR/(R_limits*R_limits));
      fz = -kspring_covs*dz*dR/dr/(1. - dR*dR/(R_limits*R_limits));

      Amino[i].fcx = Amino[i].fcx - fx;
      Amino[i].fcy = Amino[i].fcy - fy;
      Amino[i].fcz = Amino[i].fcz - fz;

      Amino[i+1].fcx = Amino[i+1].fcx + fx;
      Amino[i+1].fcy = Amino[i+1].fcy + fy;
      Amino[i+1].fcz = Amino[i+1].fcz + fz;
  
    }
  }
}/* end the FENE part of the force */


/* potential due to the applied force */
/* LJY
void tension(){

  double fx, dx, dy, dz, fy, fz;
*/
  /* this sets the direction and point of application of the external force */
  /* here it is applied on the last bead of the chain in the x-direction; 
     x0 = initial value of end-to-end distance (in starting structure);
     xt = increment in displacement of last bead under the action of applied force 
     (because we work under constant loading rate conditions) */
/* LJY
  dx = Amino[tot_amino].x - (Amino[1].x + x0 + xt);
  dy = Amino[tot_amino].y - Amino[1].y;
  dz = Amino[tot_amino].z - Amino[1].z;

  fx = -k_trans*dx;
  fy = -k_trans*dy;
  fz = -k_trans*dz;


  Amino[tot_amino].fcx = Amino[tot_amino].fcx + fx;
  Amino[tot_amino].fcy = Amino[tot_amino].fcy + fy;
  Amino[tot_amino].fcz = Amino[tot_amino].fcz + fz;
  
  forcex = fx;
*/  
  /* the first bead in the chain is kept fixed!! */
/* LJY
  Amino[1].fcx = 0.0;
  Amino[1].fcy = 0.0;
  Amino[1].fcz = 0.0;
 
}
*/

/* LJY potential due to the applied force */
void tension(){

  double rpx, rpy, rpz, rfx, rfy, rfz, rx, ry, rz, dr, fr, frx, fry, frz;

  /* this sets the direction and point of application of the external force */
  /* here it is applied on the particular beads of the chain in the r-direction; 
     r0 = initial distance value between fixed and pulled groups (in starting structure);
     xt = increment in displacement of pulled group under the action of applied force 
     (because we work under constant loading rate conditions) */

  rpx = Amino[pulled_bead].x ;
  rpy = Amino[pulled_bead].y ;
  rpz = Amino[pulled_bead].z ;

  rfx = Amino[fixed_bead].x ;
  rfy = Amino[fixed_bead].y ;
  rfz = Amino[fixed_bead].z ;

  rx = rpx - rfx;
  ry = rpy - rfy;
  rz = rpz - rfz;

  dr = sqrt(rx*rx+ry*ry+rz*rz) - r0 - xt;
  
  fr = -k_trans*dr;
  
  frx = fr*rx/sqrt(rx*rx+ry*ry+rz*rz);
  fry = fr*ry/sqrt(rx*rx+ry*ry+rz*rz);
  frz = fr*rz/sqrt(rx*rx+ry*ry+rz*rz);

  forcer = fr;
  forcerx = frx; 
  forcery = fry;
  forcerz = frz;

  /* the following beads in the protein chain are pulled!! */

  Amino[pulled_bead].fcx  = Amino[pulled_bead].fcx + frx;
  Amino[pulled_bead].fcy  = Amino[pulled_bead].fcy + fry;
  Amino[pulled_bead].fcz  = Amino[pulled_bead].fcz + frz;

  /* the following beads in the protein chain are kept fixed!! */

  Amino[fixed_bead].fcx = 0.0;
  Amino[fixed_bead].fcy = 0.0;
  Amino[fixed_bead].fcz = 0.0;

}

