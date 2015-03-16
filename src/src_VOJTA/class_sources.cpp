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

   generate_noise_source_vector(1,"11");
   
}


/////////////////////////////////////////////////////////////////

void class_sources :: generate_point_source(int X_coor, int Y_coor, int Z_coor){

   if(MPI_rank==0){
     printf(" ++++++ generating point source at x = %2d, y = %2d, z = %2d\n",X_coor, Y_coor, Z_coor);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   
   // allocate memory for the sources
   point_source = new double[2 * XYZnodeSites];
   memset(point_source,0,sizeof(point_source));

   generate_point_ixyz ( point_source, X_coor,Y_coor,Z_coor );

}


void class_sources :: generate_wall_source(){

   if(MPI_rank==0){
     printf(" ++++++ generating wall source\n");
   }
   MPI_Barrier(MPI_COMM_WORLD);
   
   // allocate memory for the sources
   wall_source = new double[2 * XYZnodeSites];
   memset(wall_source,0,sizeof(wall_source));

   generate_wall_ixyz ( wall_source );

}


void class_sources :: generate_noise_source_vector(int N_noises_in, string type){

   N_noises = N_noises_in;
   
   // allocate memory for the sources
   noise_sources_vector = new double[2* N_noises * XYZnodeSites];
   memset(noise_sources_vector,0,sizeof(noise_sources_vector));

   for(int n=0; n<N_noises;n++){

     if(MPI_rank==0){
       printf(" ++++++ generating %s noise source # %2d\n",type.c_str(),n);
     }
     MPI_Barrier(MPI_COMM_WORLD);

     if(type=="U(1)"){   
       generate_noise_ixyz_U1 ( noise_sources_vector+2*n*XYZnodeSites );
     }
     else if(type=="Z(4)"){
       generate_noise_ixyz_Z4 ( noise_sources_vector+2*n*XYZnodeSites );     
     }
     else{
       printf("ERROR - unknown noise type");
       abort();
     }
   }

   double re_av=0.0; ////////////////////////////////////////////debug
   double im_av=0.0;
   for(int i=0;i< N_noises * XYZnodeSites;i++){
     re_av+=noise_sources_vector[2*i];
     im_av+=noise_sources_vector[2*i+1];
   }

   printf("\n MPI %2i, average = %1.16e + I %1.16e\n", 
          MPI_rank, re_av/XYZnodeSites/N_noises, im_av/XYZnodeSites/N_noises);
   MPI_Barrier(MPI_COMM_WORLD);
   

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

  //MPI_Barrier(MPI_COMM_WORLD);  ////////////////////////////////////////////debug

  //printf("MPI %2i, xyztnode %2i,%2i,%2i,%2i .. coors %2i,%2i,%2i,  \n" ,
  //        MPI_rank,XnodeCoor,YnodeCoor,ZnodeCoor,TnodeCoor,
  //        X_coor, Y_coor, Z_coor);

  //MPI_Barrier(MPI_COMM_WORLD);

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


void class_sources :: generate_noise_ixyz_U1(double* noise){
  
  if(TnodeCoor==0){ 
    for(int i=0;i<XYZnodeSites;i++){

      double rand_num=2*PI*RanGen.Random();
      noise[2*i]=cos(rand_num);
      noise[2*i+1]= sin(rand_num);
    }
  }


  for(int n_x=0; n_x<Xnodes;n_x++){
    for(int n_y=0; n_y<Ynodes;n_y++){
      for(int n_z=0; n_z<Znodes;n_z++){

        int rank_from;
        rank_from = rank_from_coor(n_x,n_y,n_z, 0 );

        for(int n_t=1;n_t<Tnodes;n_t++){
          
          int rank_to;
          rank_to   = rank_from_coor(n_x,n_y,n_z,n_t);

          if(MPI_rank==rank_from){
            MPI_Send(noise,
                     2*XYZnodeSites,MPI_DOUBLE,
                     rank_to,1*n_t+100*n_z+10000*n_y+1000000*n_x,
                     MPI_COMM_WORLD);
          }
          
          MPI_Status status;
          if(MPI_rank==rank_to){
            MPI_Recv(noise,
                     2*XYZnodeSites,MPI_DOUBLE,
                     rank_from,1*n_t+100*n_z+10000*n_y+1000000*n_x,
                     MPI_COMM_WORLD,&status);
          }
          

        }
  
      }  
    }  
  }
  MPI_Barrier(MPI_COMM_WORLD);

}


void class_sources :: generate_noise_ixyz_Z4(double* noise){
  
  COMPLEX nums[4]={1.0,-1.0,COMPLEX_I,-COMPLEX_I};
  
  if(TnodeCoor==0){ 
    for(int i=0;i<XYZnodeSites;i++){

      int rand_num=RanGen.IRandom(0,3);
      ((COMPLEX*)noise)[i]=nums[rand_num];
    }
  }


  for(int n_x=0; n_x<Xnodes;n_x++){
    for(int n_y=0; n_y<Ynodes;n_y++){
      for(int n_z=0; n_z<Znodes;n_z++){

        int rank_from;
        rank_from = rank_from_coor(n_x,n_y,n_z, 0 );

        for(int n_t=1;n_t<Tnodes;n_t++){
          
          int rank_to;
          rank_to   = rank_from_coor(n_x,n_y,n_z,n_t);

          if(MPI_rank==rank_from){
            MPI_Send(noise,
                     2*XYZnodeSites,MPI_DOUBLE,
                     rank_to,1*n_t+100*n_z+10000*n_y+1000000*n_x,
                     MPI_COMM_WORLD);
          }
          
          MPI_Status status;
          if(MPI_rank==rank_to){
            MPI_Recv(noise,
                     2*XYZnodeSites,MPI_DOUBLE,
                     rank_from,1*n_t+100*n_z+10000*n_y+1000000*n_x,
                     MPI_COMM_WORLD,&status);
          }
          

        }
  
      }  
    }  
  }
  MPI_Barrier(MPI_COMM_WORLD);

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



