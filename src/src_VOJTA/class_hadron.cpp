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
#include "class_hadron.h"


using namespace HAL_idx;

void class_hadron::run_all_GF(){

  run_GF("pion");
  run_GF("kaon");
  run_GF("eta");
  
/*
  run_GF("proton");
  run_GF("sigma");
  run_GF("xi");
  run_GF("lambda");
*/

}

// ================================================================================================
// wrapper to run all correlators
//
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
  
    printf("ERROR - unknown hadron name");
    abort();
  }
  
  corr_print(correlator, hadron_name);
  
};

// ================================================================================================
// calculate pseudoscalar meson propagator < O OBAR >
// O     =  antiquarkBAR  \gamma5  quark
// OBAR  =  quarkBAR      \gamma5  antiquark
//
void class_hadron::run_GF_meson(double* correlator, double* prop_quark, double* prop_antiquark){

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
// calculate octet-baryon (except Lambda) propagator < O OBAR >
// O     =  \eps_abc  quark1_a     ( quark1_b^T   C\gamma5  quark2_c      )
// OBAR  =  \eps_abc  quark1BAR_a  ( quark1BAR_b  C\gamma5  quark2BAR_c^T )
//
void class_hadron::run_GF_baryon_octet(double* correlator, double* prop_quark1, double* prop_quark2){

  // complexify propagators
  COMPLEX* Prop_quark1  =  (COMPLEX*)prop_quark1  ;//+ prop_slv_idx(0,0,0,0,ixyz,it); 
  COMPLEX* Prop_quark2  =  (COMPLEX*)prop_quark2  ;//+ prop_slv_idx(0,0,0,0,ixyz,it);

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
                       *            Prop_quark1     [prop_slv_idx(c1,alpha2, c1P,alpha1P,ixyz,it)  ]
                       *            Prop_quark1     [prop_slv_idx(c1,alpha2, c1P,alpha1P,ixyz,it)]  ;             

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




// ================================================================================================
// calculate octet-baryon Lambda propagator < O OBAR >
// O     =  \eps_abc  2 s_a    ( u_b^T     C\gamma5  d_c ) 
//                   +  d_a    ( u_b^T     C\gamma5  s_c )    
//                   -  u_a    ( d_b^T     C\gamma5  s_c )
// OBAR  =  \eps_abc  2 sBAR_a ( uBAR_b^T  C\gamma5  dBAR_c ) 
//                   +  dBAR_a ( uBAR_b^T  C\gamma5  sBAR_c )
//                   -  uBAR_a ( dBAR_b^T  C\gamma5  sBAR_c )
//
void class_hadron::run_GF_baryon_octet_lambda(double* correlator, 
                   double* prop_up, double* prop_down, double* prop_strange){


}



// ================================================================================================
// print the correlator
//
void class_hadron::corr_print(double *correlator, string hadron_name)
{

  MPI_Barrier(MPI_COMM_WORLD);
  
  if(MPI_rank==0){

    // create directory and file name
    string dir="results/"+base+"/"+prefix+"correlators/";
    printf("%s %s\n",prefix.c_str(),dir.c_str());
    create_directory(dir);
  
    char wfile[256];
    snprintf(wfile,sizeof(wfile), "%s%s.%02d",
             dir.c_str(),
             hadron_name.c_str(),
             iT_src);

 
    // open output file
    string ofname(wfile);
    std::ofstream fout(ofname.c_str());
    fout.setf(std::ios::scientific);
    fout.precision(16);

    // print correlator
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




