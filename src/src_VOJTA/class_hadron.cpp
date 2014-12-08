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
#include "class_hadron.h"
#include "HAL_indexes.h"


using namespace HAL_idx;


void class_hadron::run_GF(char hadron_name[]){

if(strcmp(hadron_name,"pion")==0){
  run_GF_meson(prop_ud,prop_ud);
}
else{
printf("ERROR");
}

};

void class_hadron::run_GF_meson(int gamma_src, double* prop_up, int gamma_snk, double* prop_down){

  double correlator_local[2*Tsites];
  memset(correlator_local,0,sizeof(correlator_local));
  
  double correlator_global[2*Tsites];
  memset(correlator_global,0,sizeof(correlator_global));
  
  int id=Communicator::self();
  
    #pragma omp parallel for
    for(int it = 0; it < TnodeSites; it++){
 
      COMPLEX sum = COMPLEX_ZERO;

      for(int ixyz = 0; ixyz < XYZnodeSites; ixyz++){
//      printf("%s MPI id %i,%i  OMP id %i      num %i\n", hadron_name, id, id2, omp_get_thread_num(), ixyz);

        COMPLEX* Prop_up = (COMPLEX*)prop_ud + prop_slv_idx(0,0,0,0,ixyz,it); 
        COMPLEX* Prop_dn = (COMPLEX*)prop_ud + prop_slv_idx(0,0,0,0,ixyz,it);

        for(      int alpha1P = 0; alpha1P < 4; alpha1P++){
          for(    int c1P     = 0; c1P     < 3; c1P++){
            for(  int alpha1  = 0; alpha1  < 4; alpha1++){
              for(int c1      = 0; c1      < 3; c1++){
                sum
                  +=      Prop_up[ prop_slv_cs_idx(c1,alpha1, c1P,alpha1P) ]
                  *  Conj(Prop_dn[ prop_slv_cs_idx(c1,alpha1, c1P,alpha1P) ]);
              }
            }
          }
        }

      } // ixyz

    printf("MPI %i OMP %i ... it %i sum %1.16e %1.16e I\n", id,omp_get_thread_num(), it,Real(sum),Imag(sum));
    ((COMPLEX*)correlator_local)[it + TnodeSites * TnodeCoor]=sum;
    
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




