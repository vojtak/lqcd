/* standard */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstdlib>

#include <mpi.h>
#include <omp.h>

/* local from Bridge++ */
#include "communicator.h"
#include "commonParameters.h"

/* local from HAL */
#include "COMPLEX.h"
//#include "HAL_idx.h"
//#include "NBS_idx.h"

/* local from VOJTA */
#include "HAL_indexes.h"
#include "class_hadron.h"


using namespace HAL_idx;


void class_hadron::run_GF(char hadron_name[]){

if(strcmp(hadron_name,"pion")==0){
  run_GF_meson(prop_ud,prop_ud);
}
else{
printf("ERROR");
}

};

void class_hadron::run_GF_meson(double* prop_up, double* prop_down){

  // complexify propagators
  COMPLEX* Prop_up   = (COMPLEX*)prop_up   ;//+ prop_slv_idx(0,0,0,0,ixyz,it); 
  COMPLEX* Prop_down = (COMPLEX*)prop_down ;//+ prop_slv_idx(0,0,0,0,ixyz,it);

#define back_prop(prop,c,a,cp,ap,ixyzt,it)                               \
        ( ZGM(a,5) * ZGM (IGM(ap,5),5) *                                 \
          Conj(prop[ prop_slv_idx(c,IGM(a,5),cp,IGM(ap,5) ,ixyz,it) ]) )

  double correlator_local[2*Tsites];
  memset(correlator_local,0,sizeof(correlator_local));
  
  double correlator_global[2*Tsites];
  memset(correlator_global,0,sizeof(correlator_global));
  

  int MPI_flag, MPI_size, MPI_rank;
  MPI_Initialized(&MPI_flag);
  MPI_Comm_size(MPI_COMM_WORLD,&MPI_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&MPI_rank);

  printf("is %i, size %i, id %i \n",MPI_flag,MPI_size, MPI_rank);
  
    #pragma omp parallel for
    for(int it = 0; it < TnodeSites; it++){
 

      // ixyz summation
      COMPLEX sum_ixyz = COMPLEX_ZERO;
      for(int ixyz = 0; ixyz < XYZnodeSites; ixyz++){


        // source contraction
        COMPLEX sum_src = COMPLEX_ZERO;
        for(      int alpha1P = 0; alpha1P < 4; alpha1P++){
          for(    int c1P     = 0; c1P     < 3; c1P++){
            
            int alpha2P=IGM(alpha1P,5);
           
            // sink contraction
            COMPLEX sum_snk = COMPLEX_ZERO;
            for(  int alpha1  = 0; alpha1  < 4; alpha1++){
              for(int c1      = 0; c1      < 3; c1++){
                
                int alpha2=IGM(alpha1,5);

                sum_snk +=  ZGM(alpha1,5) 
                         *            Prop_up   [prop_slv_idx(c1,alpha2, c1P,alpha1P,ixyz,it) ]
                         *  back_prop(Prop_down,              c1,alpha1, c1P,alpha2P,ixyz,it)  ;             

              }
            } // sink

            sum_src +=  sum_snk 
                     *  ZGM(alpha1P,5);
 
          }
        } // source

        sum_ixyz += sum_src ;
        
      } // ixyz

    printf("MPI %i OMP %i ... it %i sum %1.16e %1.16e I\n", 
           MPI_rank,omp_get_thread_num(), it,Real(sum_ixyz),Imag(sum_ixyz));
 
    ((COMPLEX*)correlator_local)[it + TnodeSites * TnodeCoor]=sum_ixyz;
    
  } // it, end of omp parallel

  Communicator::sync();
  MPI_Reduce(correlator_local, correlator_global, 2*Tsites, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  Communicator::sync();
  if(Communicator::self()==0){
    for(int it=0;it<Tsites;it++){
      printf("local    ... MPI %i it %3i sum %1.16e %1.16e I\n", 0, it,correlator_local[2*it],correlator_local[2*it+1]);
    }
    printf("---------------\n");
    for(int it=0;it<Tsites;it++){
      printf("global   ... MPI %i it %3i sum %1.16e %1.16e I\n", 0, it,correlator_global[2*it],correlator_global[2*it+1]);
    }
    printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
  }

}




