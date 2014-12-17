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
#include "class_noise_source.h"



void class_noise_source :: run(){

   for(int n=0; n<N_noises;n++){
   
     generate_noise_ixyz ( noise_vector+2*n*XYZnodeSites );
   
   }

}

void class_noise_source :: print(){

  if(MPI_rank==0){
  for(int n=0; n<N_noises;n++){
   
    for(int i=0; i<XYZnodeSites+10;i++){
      
      printf("MPI %2i ... n= %2i ... i %2i   ran_num %1.16e %1.16e I, %1.16e %1.16e I, \tnorm %1.16e\n", 
               MPI_rank, n , i , noise_vector[2*(n*XYZnodeSites+i)],noise_vector[2*(n*XYZnodeSites+i)+1], 
               Real(((COMPLEX*)noise_vector)[n*XYZnodeSites+i]),Imag(((COMPLEX*)noise_vector)[n*XYZnodeSites+i]),
               ((COMPLEX*)noise_vector)[n*XYZnodeSites+i]*Conj(((COMPLEX*)noise_vector)[n*XYZnodeSites+i])
               );
    }

    printf("MPI %2i ................................... n= %2i\n\n",
               MPI_rank, n );
  }

//  for(int i=0;i< N_noises * XnodeSites;i++){
//  printf("%1.16e %1.16e I\n",noise_vector[2*i], noise_vector[2*i+1]);
//  }
  }
}

void class_noise_source :: print_full(){

  if(MPI_rank==0){
  for(int i=0; i<N_noises*XYZnodeSites+10;i++){
      
      printf("MPI %2i ... n= %2i ... i %2i   ran_num %1.16e %1.16e I, %1.16e %1.16e I, norm %1.16e\n", 
               MPI_rank, -1 , i , noise_vector[2*i],noise_vector[2*i+1], 
               Real(((COMPLEX*)noise_vector)[i]),Imag(((COMPLEX*)noise_vector)[i]),
               ((COMPLEX*)noise_vector)[i]*Conj(((COMPLEX*)noise_vector)[i])
               );
   }
//  for(int i=0;i< N_noises * XnodeSites;i++){
//  printf("%1.16e %1.16e I\n",noise_vector[2*i], noise_vector[2*i+1]);
//  }
  }
}


double * class_noise_source :: get_noise_ixyz(int noise_num){

  return noise_vector+2*noise_num*XYZnodeSites;

}



void class_noise_source :: generate_noise_ixyz(double* noise){
   
  for(int i=0;i<XYZnodeSites;i++){

    double rand_num=2*PI*RanGen.Random();
    noise[2*i]=cos(rand_num);
    noise[2*i+1]= sin(rand_num);
  
  }


}
