/*
 *  model.c
 *  LIME, The versatile 3D line modeling tool 
 *
 *  Created by Christian Brinch on 11/05/07.
 *  Copyright 2006-2013, Christian Brinch, 
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 *  modelHD_CO.c
 *  Edited by I. van Leeuwen and Y. Rusticus
 *  Date: 27/03/19
 *  Input model for disk HD100546 CO 2-1 emission
 *
 */

#include <stdio.h>
#include "lime.h"
#include "structureHD.h"

/******************************************************************************/

void
input(inputPars *par, image *img){
    
  //Basic parameters. See cheat sheet (!updated cheat sheet) for details.
  par->radius = 500*AU;
  par->minScale	= 0.2*AU;
  par->pIntensity = 10000;
  par->sinkPoints = 5000;
  par->dust = "jena_thin_e6.tab";
  par->moldatfile[0] = "data_lamda_12cov2.dat";
  par->moldatfile[1] = "data_lamda_13co.dat";
  par->moldatfile[2] = "data_lamda_c17o.dat";
  par->moldatfile[3] = "data_lamda_c18o.dat";
  par->sampling	= 0;
  par->outputfile = "populationsHD.pop";
  par->gridfile	= "gridHD.vtk";
  par->lte_only = 1; //unreliable, much faster

  //Definitions for image #0. Add blocks for additional images.
  img[0].nchan			= 200;		      // Number of channels
  img[0].velres			= 200.;           // Channel resolution in m/s
  img[0].trans			= 1;              // zero-indexed J quantum number
  img[0].pxls			= 400;	          // Pixels per dimension
  img[0].imgres			= 0.02;		      // Resolution in arc seconds
  img[0].theta			= -44 * PI/180.;  // 0: face-on, pi/2: edge-on
  img[0].distance		= 110*PC;	      // source distance in m
  img[0].source_vel		= 9.1e3;          // source velocity in m/s
  img[0].unit			= 1;		      // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[0].filename		= "HD_CO21_10000p_5000s_lte.fits";  // Output filename
}

/******************************************************************************/

void
density(double x, double y, double z, double *density){	
    
  double r;
  double dist;
  double bestdens;
  double mindist = 9999;

  //Find density in DALI output
  r=sqrt(x*x+y*y); //in m 
  for (int i=0; i < model_r_size; i++) {
      for (int j=0; j < model_z_size; j++) {
          dist = sqrt( (r/AU-model_r_mid[i][j]) * (r/AU-model_r_mid[i][j]) 
                      + (sqrt(z*z)/AU-model_z_mid[i][j]) * (sqrt(z*z)/AU-model_z_mid[i][j]) );
          if (dist < mindist) {
              mindist = dist;
              bestdens = model_H2[i][j]; //cm^-3
          }//if
      }//for
  }//for
  density[0] = 0.25 * bestdens*1e6; //collision partners CO: o-H2 and p-H2
  density[1] = 0.75 * bestdens*1e6; //assume 3:1 ratio of H2 density
}

/******************************************************************************/

void
temperature(double x, double y, double z, double *temperature){
    
  double r; 
  double dist;
  double besttemp0; //t_gas
  double besttemp1; //t_dust
  double mindist = 9999;

  //Find density in DALI output
  r=sqrt(x*x+y*y); //in m
  for (int i=0; i < model_r_size; i++) {
      for (int j=0; j < model_z_size; j++) {
          dist = sqrt( (r/AU-model_r_mid[i][j]) * (r/AU-model_r_mid[i][j]) 
                      + (sqrt(z*z)/AU-model_z_mid[i][j]) * (sqrt(z*z)/AU-model_z_mid[i][j]) );
          if (dist < mindist) {
              mindist = dist;
              besttemp0 = model_t_gas[i][j];  //K
              besttemp1 = model_t_dust[i][j]; //K
          }//if
      }//for
  }//for

 temperature[0] = besttemp0;
 temperature[1] = besttemp1;
}

/******************************************************************************/

void
abundance(double x, double y, double z, double *abundance){
/*  
 * Here we use a constant abundance. Could be a 
 * function of (x,y,z).
 */
  abundance[0] = 1e-4;
}

/******************************************************************************/

void
doppler(double x, double y, double z, double *doppler){
/* 
 * 200 m/s as the doppler b-parameter. This
 * can be a function of (x,y,z) as well.
 * Note that *doppler is a pointer, not an array. 
 * Remember the * in front of doppler.
 */
  *doppler = 200.;
}

/******************************************************************************/

void
velocity(double x, double y, double z, double *vel){
    
  double r, v, theta, rToUse;
  const double rMin = 0.1*AU;

  //avoid singularities
  r=sqrt(x*x+y*y);
  if(r>rMin)
    rToUse = r;
  else
    rToUse = rMin;
    
  //avoid singularities
  if (x == 0 && y >= 0)
      theta = 0.5*PI/180;
  else if (x == 0 && y < 0)
      theta = 1.5*PI/180;
  else
      theta=atan(y/(float)x);
    
  //orbital velocity
  v=sqrt((2.3*1.989e30*6.67e-11)/rToUse); //mass=2.3Msol

  //velocities seem to have the wrong direction: quick fix, needs better look
  //assume flat disk
  if (x>=0) {
    vel[0]=v*sin(theta);  
    vel[1]=v*(-cos(theta));
  }
  else {
    vel[0]=v*(-sin(theta));
    vel[1]=v*cos(theta);
  }
  vel[2]=0; //might need to implement the z-velocity 
}

/******************************************************************************/