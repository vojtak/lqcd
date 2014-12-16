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



void class_noise_source_propagator :: run(){

   for(int n=0; n<N_noises;n++){
   
     generate_noise_ixyz ( noise_vector+2*n*XYZnodeSites );
   
   }

}

void class_noise_source_propagator :: print(){

  if(MPI_rank==0){
  for(int n=0; n<N_noises;n++){
   
    for(int i=0; i<XYZnodeSites;i++){
      
      printf("MPI %2i ... n= %2i ... i %2i   ran_num %1.16e %1.16e I, %1.16e %1.16e I\n", 
               MPI_rank, n , i , noise_vector[2*(n*XnodeSites+i)],noise_vector[2*(n*XnodeSites+i)+1], 
               Real(((COMPLEX*)noise_vector)[n*XnodeSites+i]),Imag(((COMPLEX*)noise_vector)[n*XnodeSites+i]));
    }

    printf("MPI %2i ................................... n= %2i\n\n",
               MPI_rank, n );
  }

  for(int i=0;i< N_noises * XnodeSites;i++){
  printf("%1.16e %1.16e I\n",noise_vector[2*i], noise_vector[2*i+1]);
  }
  }
}



double * class_noise_source_propagator :: get_noise_ixyz(int noise_num){

  return noise_vector+2*noise_num*XYZnodeSites;

}



void class_noise_source_propagator :: generate_noise_ixyz(double* noise){
   
  for(int i=0;i<XYZnodeSites;i++){

    double rand_num=2*PI*RanGen.Random();
    noise[2*i]=cos(rand_num);
    noise[2*i+1]= sin(rand_num);
  
  }


}
