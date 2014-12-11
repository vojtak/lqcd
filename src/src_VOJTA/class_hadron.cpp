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


void class_hadron::run_GF(string hadron_name){

  double correlator[2*Tsites];
  memset(correlator,0,sizeof(correlator));


  if(hadron_name=="pion"){
  
    run_GF_meson(correlator, prop_ud, prop_ud);

  }
  else if(hadron_name=="kaon"){
  
    run_GF_meson(correlator, prop_ud, prop_s);

  }
  else if(hadron_name=="eta"){
  
    run_GF_meson(correlator, prop_s, prop_s);

  }
  else{
  
    printf("ERROR");

  }
  
  corr_print(correlator, hadron_name);
  
};


void class_hadron::run_GF_meson(double* correlator, double* prop_up, double* prop_down){

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

//    printf("MPI %i OMP %i ... it %i sum %1.16e %1.16e I\n", 
//           MPI_rank,omp_get_thread_num(), it,Real(sum_ixyz),Imag(sum_ixyz));
 
    ((COMPLEX*)correlator_local)[it + TnodeSites * TnodeCoor]=sum_ixyz;
    
  } // it, end of omp parallel

  // reduce from all MPI processes
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(correlator_local, correlator_global, 2*Tsites, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  // output correlator
  for(int it = 0; it < Tsites; it++){
    ((COMPLEX*)correlator)[it]=((COMPLEX*)correlator_global)[it];
  }
  
}



void class_hadron::corr_print(double *correlator, string hadron_name)
{

  int MPI_flag, MPI_size, MPI_rank;
  MPI_Initialized(&MPI_flag);
  MPI_Comm_size(MPI_COMM_WORLD,&MPI_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&MPI_rank);

  MPI_Barrier(MPI_COMM_WORLD);
  
  if(MPI_rank==0){
  
  char wfile[256];
  snprintf(wfile,sizeof(wfile), "results/%s/%s.%02d",
             base.c_str(),
             hadron_name.c_str(),
             iT_src);
  printf(wfile);

  string ofname(wfile);
  std::ofstream fout(ofname.c_str());
  fout.setf(std::ios::scientific);
  fout.precision(16);

  for(int it = 0; it < Tsites; it++){
    int iT2 = (it + iT_src + 100*Tsites) % Tsites;
    char line[1000];
    snprintf(line, sizeof(line), "%4d\t%1.16e %1.16e\n",
             it, Real(correlator[2* iT2]),Imag(correlator[2* iT2+1]));
    fout << line;
  }

  fout.close();

  }

  MPI_Barrier(MPI_COMM_WORLD);

}

