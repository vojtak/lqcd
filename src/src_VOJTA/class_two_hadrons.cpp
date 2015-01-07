/* standard */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstdlib>


/* standard parralelization */
#include <mpi.h>
#include <omp.h>


/* local from HAL */
#include "COMPLEX.h"
#include "HAL_indexes.h"

/* local from VOJTA */
#include "class_global_wrapper.h"
#include "class_two_hadrons.h"


using namespace HAL_idx;

void class_two_hadrons::run_all_GF(){

  run_GF("Pion-Sigma");

}

// ================================================================================================
// wrapper to run all correlators
//
void class_two_hadrons::run_GF(string hadrons_names){

  double correlator[2*Tsites];
  memset(correlator,0,sizeof(correlator));

  if(MPI_rank==0){
    printf("  +++ run_GF : calculate %s propagator        %s\n",
         hadrons_names.c_str(), LocalTime().c_str());
  }
  
  if(hadrons_names=="Pion-Sigma"){
  
    run_GF_meson(correlator, prop_ud, prop_ud);
  }

  else{
  
    printf("ERROR - unknown hadron name");
    abort();
  }

  corr_print(correlator, hadrons_names);
  
};

// ================================================================================================
// calculate pseudoscalar meson propagator < O OBAR >
// O     =  antiquarkBAR  \gamma5  quark
// OBAR  =  quarkBAR      \gamma5  antiquark
//
void class_two_hadrons::run_GF_meson(double* correlator, double* prop_quark, double* prop_antiquark){

  // complexify propagators
  COMPLEX* Prop_quark     = (COMPLEX*)prop_quark     ;//+ prop_slv_idx(0,0,0,0,ixyz,it); 
  COMPLEX* Prop_antiquark = (COMPLEX*)prop_antiquark ;//+ prop_slv_idx(0,0,0,0,ixyz,it);

#define back_prop(prop,c,a,cp,ap,ixyzt,it)                               \
        ( ZGM(a,5) * ZGM (IGM(ap,5),5) *                                 \
          Conj(prop[ prop_slv_idx(c,IGM(a,5),cp,IGM(ap,5) ,ixyz,it) ]) )

  double correlator_local[2*Tsites];
  memset(correlator_local,0,sizeof(correlator_local));
  
  double correlator_global[2*Tsites];
  memset(correlator_global,0,sizeof(correlator_global));
  
  // correlator itself
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
                       *            Prop_quark     [prop_slv_idx(c1,alpha2, c1P,alpha1P,ixyz,it)  ]
                       *  back_prop(Prop_antiquark,              c1,alpha1, c1P,alpha2P,ixyz,it)  ;             

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
  
  
#undef back_prop
}





// ================================================================================================
// print the correlator
//
void class_two_hadrons::corr_print(double *correlator, string hadrons_names)
{

  MPI_Barrier(MPI_COMM_WORLD);
  
  if(MPI_rank==0){

    // create directory and file name
    string dir="results/"+base+"/"+prefix+"correlators/";
    //printf("%s %s\n",prefix.c_str(),dir.c_str());
    create_directory(dir);
  
    char wfile[256];
    snprintf(wfile,sizeof(wfile), "%s%s.%02d",
             dir.c_str(),
             hadrons_names.c_str(),
             iT_src);

 
    // open output file
    string ofname(wfile);
    std::ofstream fout(ofname.c_str());
    fout.setf(std::ios::scientific);
    fout.precision(16);


    // print correlator
    for(int it = 0; it < Tsites; it++){
      int iT2 = (it + iT_src + 100*Tsites) % Tsites;
      //int iT2 = it;

      char line[1000];
      snprintf(line, sizeof(line), "%4d\t%1.16e %1.16e\n",
               it, correlator[2* iT2], correlator[2* iT2+1]);
      fout << line;
    }

    fout.close();

  }

  MPI_Barrier(MPI_COMM_WORLD);

}




