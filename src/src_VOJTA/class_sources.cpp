/* standard */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstdlib>


/* standard parralelization */
#include <mpi.h>
#include <omp.h>


/* local from EXTERNAL */
#include "randomc.h"

/* local from HAL */
#include "COMPLEX.h"
#include "HAL_indexes.h"

/* local from VOJTA */
#include "class_global_wrapper.h"
#include "class_sources.h"



void class_sources :: run(){

   generate_point_source(2,1,11);

   generate_wall_source();

   generate_noise_source_vector(1);
   
}


/////////////////////////////////////////////////////////////////

void class_sources :: generate_point_source(int X_coor, int Y_coor, int Z_coor){

   // allocate memory for the sources
   point_source = new double[2 * XYZnodeSites];
   memset(point_source,0,sizeof(point_source));

   generate_point_ixyz ( point_source, X_coor,Y_coor,Z_coor );

}


void class_sources :: generate_wall_source(){

   // allocate memory for the sources
   wall_source = new double[2 * XYZnodeSites];
   memset(wall_source,0,sizeof(wall_source));

   generate_wall_ixyz ( wall_source );

}


void class_sources :: generate_noise_source_vector(int N_noises_in){

   N_noises = N_noises_in;
   
   // allocate memory for the sources
   noise_sources_vector = new double[2* N_noises * XYZnodeSites];
   memset(noise_sources_vector,0,sizeof(noise_sources_vector));

   for(int n=0; n<N_noises;n++){
   
     generate_noise_ixyz ( noise_sources_vector+2*n*XYZnodeSites );
   
   }

}


/////////////////////////////////////////////////////////////////


double * class_sources :: get_point_ixyz(){

  return point_source;

}


double * class_sources :: get_wall_ixyz(){

  return wall_source;

}


double * class_sources :: get_noise_ixyz(int noise_num){

  return noise_sources_vector+2*noise_num*XYZnodeSites;

}


/////////////////////////////////////////////////////////////////


void class_sources :: generate_point_ixyz(double *point, int X_coor, int Y_coor, int Z_coor){

  int X_coor_node, Y_coor_node, Z_coor_node;

  MPI_Barrier(MPI_COMM_WORLD);

  printf("MPI %2i, xyztnode %2i,%2i,%2i,%2i .. coors %2i,%2i,%2i,  \n" ,
          MPI_rank,XnodeCoor,YnodeCoor,ZnodeCoor,TnodeCoor,
          X_coor, Y_coor, Z_coor);

  MPI_Barrier(MPI_COMM_WORLD);

  if (X_coor/XnodeSites == XnodeCoor){
    if (Y_coor/YnodeSites == YnodeCoor){
      if (Z_coor/ZnodeSites == ZnodeCoor){
        
        X_coor_node = (X_coor+100*XnodeSites) % XnodeSites;
        Y_coor_node = (Y_coor+100*YnodeSites) % YnodeSites;
        Z_coor_node = (Z_coor+100*ZnodeSites) % ZnodeSites;

  printf("selected: MPI %2i, xyztnode %2i,%2i,%2i,%2i .. coors %2i,%2i,%2i,  \n" ,
          MPI_rank,XnodeCoor,YnodeCoor,ZnodeCoor,TnodeCoor,
          X_coor_node, Y_coor_node, Z_coor_node);
       
        point[2 * index_xyz(X_coor_node,Y_coor_node,Z_coor_node)] =1.00;

      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  

}


void class_sources :: generate_wall_ixyz(double* wall){
   
  for(int i=0;i<XYZnodeSites;i++){

    wall[2*i]  =1.0/XYZsites;
    wall[2*i+1]= 0.0;
  
  }

}


void class_sources :: generate_noise_ixyz(double* noise){
   
  for(int i=0;i<XYZnodeSites;i++){

    double rand_num=2*PI*RanGen.Random();
    noise[2*i]=cos(rand_num);
    noise[2*i+1]= sin(rand_num);
  
  }

}


/////////////////////////////////////////////////////////////////


void class_sources :: print(){

  if(MPI_rank==0){
//  for(int n=0; n<N_noises;n++){
   
    for(int i=0; i<XYZnodeSites;i++){
      
      printf("MPI %2i ... n= %2i ... i %2i   ran_num %1.16e %1.16e I, %1.16e %1.16e I, \tnorm %1.16e\n", 
               MPI_rank, 0 , i , point_source[2*(i)],point_source[2*(i)+1], 
               Real(((COMPLEX*)point_source)[i]),Imag(((COMPLEX*)point_source)[i]),
               Real(((COMPLEX*)point_source)[i]*Conj(((COMPLEX*)point_source)[i]))
               );

//      printf("MPI %2i ... n= %2i ... i %2i   ran_num %1.16e %1.16e I, %1.16e %1.16e I, \tnorm %1.16e\n", 
//               MPI_rank, n , i , noise_sources_vector[2*(n*XYZnodeSites+i)],noise_sources_vector[2*(n*XYZnodeSites+i)+1], 
//               Real(((COMPLEX*)noise_sources_vector)[n*XYZnodeSites+i]),Imag(((COMPLEX*)noise_sources_vector)[n*XYZnodeSites+i]),
//               Real(((COMPLEX*)noise_sources_vector)[n*XYZnodeSites+i]*Conj(((COMPLEX*)noise_sources_vector)[n*XYZnodeSites+i]))
//               );
    }

//    printf("MPI %2i ................................... n= %2i\n\n",
//               MPI_rank, n );
//  }

//  for(int i=0;i< N_noises * XnodeSites;i++){
//  printf("%1.16e %1.16e I\n",noise_vector[2*i], noise_vector[2*i+1]);
//  }
  }
}

void class_sources :: print_full(){

  if(MPI_rank==0){
  for(int i=0; i<N_noises*XYZnodeSites+10;i++){
      
      printf("MPI %2i ... n= %2i ... i %2i   ran_num %1.16e %1.16e I, %1.16e %1.16e I, norm %1.16e\n", 
               MPI_rank, -1 , i , noise_sources_vector[2*i],noise_sources_vector[2*i+1], 
               Real(((COMPLEX*)noise_sources_vector)[i]),Imag(((COMPLEX*)noise_sources_vector)[i]),
               Real(((COMPLEX*)noise_sources_vector)[i]*Conj(((COMPLEX*)noise_sources_vector)[i]))
               );
   }
//  for(int i=0;i< N_noises * XnodeSites;i++){
//  printf("%1.16e %1.16e I\n",noise_vector[2*i], noise_vector[2*i+1]);
//  }
  }
}


