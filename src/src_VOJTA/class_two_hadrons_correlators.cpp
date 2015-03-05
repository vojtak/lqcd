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
#include "class_sources.h"
#include "class_two_hadrons.h"


using namespace HAL_idx;


void class_two_hadrons::run_all_GF(){

  run_GF("pion-sigma");

}


// ================================================================================================
// wrapper to run all correlators
//
void class_two_hadrons::run_GF(string hadron_names){

  double correlator[2*Tsites];
  memset(correlator,0,sizeof(correlator));

  if(MPI_rank==0){
    printf(" ++++++ run_GF : calculate %s propagator        %s\n",
         hadron_names.c_str(), LocalTime().c_str());
  }
  MPI_Barrier(MPI_COMM_WORLD);


  // ==================================
  // pion-sigma two baryon system
    
  if(hadron_names=="pion-sigma"){

/*
    int indexes[8] = {7,8,9,10,13,14,15,16};
    
    for(int rel =1;rel<4;rel++){
      for(int jjj=0;jjj<8;jjj++){

    if(MPI_rank==0){
      printf("       ++++++ run_GF : the LOOP part test %i, %04i    begin %s\n",
             rel, indexes[jjj], LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    double correlator_loop_test[2*Tsites];
    memset(correlator_loop_test,0,sizeof(correlator_loop_test));  
    run_GF_pi_sigma_loop_TEST(correlator_loop_test,rel,indexes[jjj]);
      char pr[50];
      snprintf(pr,sizeof(pr),"%s_loop_TEST_%i_%02i",hadron_names.c_str(),rel,indexes[jjj]); 

    corr_print(correlator_loop_test, pr);
    if(MPI_rank==0){
      printf("       ++++++ run_GF :                        end %s\n", 
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    }}
*/    
    if(MPI_rank==0){
      printf("       ++++++ run_GF : the TREE part        begin %s\n",
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    double correlator_tree[2*Tsites];
    memset(correlator_tree,0,sizeof(correlator_tree));  
    run_GF_pi_sigma_tree(correlator_tree);
    corr_print(correlator_tree     , hadron_names+"_tree"     );
    if(MPI_rank==0){
      printf("       ++++++ run_GF :                        end %s\n", 
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);

/*
    if(MPI_rank==0){
      printf("       ++++++ run_GF : the TREE part NN     begin %s\n",
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    double correlator_tree_NONOISE[2*Tsites];
    memset(correlator_tree_NONOISE,0,sizeof(correlator_tree_NONOISE));  
    run_GF_pi_sigma_tree_NONOISE(correlator_tree_NONOISE);
    corr_print(correlator_tree_NONOISE, hadron_names+"_tree_NONOISE");
    if(MPI_rank==0){
      printf("       ++++++ run_GF :                        end %s\n", 
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);


    if(MPI_rank==0){
      printf("       ++++++ run_GF : the LOOP part        begin %s\n",
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);

    double correlator_loop[2*Tsites];
    memset(correlator_loop,0,sizeof(correlator_loop));  
    run_GF_pi_sigma_loop(correlator_loop);
    GF_print(correlator_loop     , hadron_names+"_loop"     );
    if(MPI_rank==0){
      printf("       ++++++ run_GF :                        end %s\n", 
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);

    
    if(MPI_rank==0){
      printf("       ++++++ run_GF : the LOOP part  NN    begin %s\n",
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);

    double correlator_loop_NONOISE[2*Tsites];
    memset(correlator_loop_NONOISE,0,sizeof(correlator_loop_NONOISE));  
    run_GF_pi_sigma_loop_NONOISE(correlator_loop_NONOISE);
    corr_print(correlator_loop_NONOISE     , hadron_names+"_loop_NONOISE"     );
    if(MPI_rank==0){
      printf("       ++++++ run_GF :                        end %s\n", 
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //for(int i=0;i<2*Tsites;i++){
    //  correlator[i]=correlator_tree[i]+correlator_loop[i];    
    //}
    //corr_print(correlator          , hadron_names             );
*/
  }
  else{
  
    printf("ERROR - unknown hadron name");
    abort();
  }


  //corr_print(correlator, hadron_names);
  
};


// ================================================================================================
// calculate pi-sigma propagator --- the tree part 
//       formulas in notes 
//
void class_two_hadrons::run_GF_pi_sigma_tree(double* correlator){

  // complexify propagators 
  COMPLEX* Prop_ud            = (COMPLEX*)prop_ud        ;//+ prop_slv_idx(0,0,0,0,ixyz,it); 
  COMPLEX* Prop_s             = (COMPLEX*)prop_s         ;//+ prop_slv_idx(0,0,0,0,ixyz,it);

#define back_prop(prop,c,a,cp,ap,ixyz,it)                               \
        ( ZGM(a,5) * ZGM (IGM(ap,5),5) *                                 \
          Conj(prop[ prop_slv_idx(c,IGM(a,5),cp,IGM(ap,5) ,ixyz,it) ]) )

  double correlator_local[2*Tsites];
  memset(correlator_local,0,sizeof(correlator_local));
  
  double correlator_global[2*Tsites];
  memset(correlator_global,0,sizeof(correlator_global));


  
  // correlator itself
  #pragma omp parallel for
  for(int it = 0; it < TnodeSites; it++){

    // noise summation
    COMPLEX sum_N = COMPLEX_ZERO;
    for(int i_noise = 0; i_noise < N_noises; i_noise++){
    
      COMPLEX* Noise      = (COMPLEX*)sources->get_noise_ixyz(i_noise);

      // calculating backward propagator contraction first
      COMPLEX back_prop_contraction[3*4*3*4];
      memset(back_prop_contraction,0,sizeof(back_prop_contraction));

      for(int a=0;a<3;a++){
      for(int alf=0;alf<4;alf++){
      for(int aP=0; aP<3;aP++){
      for(int alfP=0; alfP<4;alfP++){
      
        int index=prop_slv_cs_idx(a, alf, aP, alfP);

        for(int Y_ixyz = 0;  Y_ixyz < XYZnodeSites; Y_ixyz++){            
          back_prop_contraction[index] += 
                                Noise[Y_ixyz] * 
                                back_prop(Prop_ud,       a, alf, aP, alfP, Y_ixyz, it);
        }

      }}}}

      // free Dirac index summation
      COMPLEX sum_freeDI = COMPLEX_ZERO;
      for(int ALPHA = 0; ALPHA < 2; ALPHA++){


        // source summation
        COMPLEX sum_source = COMPLEX_ZERO;
        for(int dP      = 0; dP     < 3; dP++){
        for(int alphaP  = 0; alphaP < 4; alphaP++){
          int betaP=IGM(alphaP,5);
        for(int colorP  = 0; colorP < 6; colorP++){
          int aP    =Eps(0,colorP);
          int bP    =Eps(1,colorP);
          int cP    =Eps(2,colorP);          
        for(int gammaP  = 0; gammaP < 4; gammaP++){
          int deltaP = icg5[gammaP]; 
                    
          
          // sink summation
          COMPLEX sum_sink = COMPLEX_ZERO;
          for(int d       = 0; d      < 3; d ++){
          for(int alpha   = 0; alpha  < 4; alpha++){
            int beta=IGM(alpha,5);
          for(int color   = 0; color  < 6; color++){
            int a    =Eps(0,color);
            int b    =Eps(1,color);
            int c    =Eps(2,color);          
          for(int gamma   = 0; gamma  < 4; gamma++){
            int delta = icg5[gamma]; 


            //sum X_ixyz
            COMPLEX sum_X_ixyz = COMPLEX_ZERO;
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              sum_X_ixyz+=Conj(Noise[X_ixyz]) * 
                          Prop_s[ prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ]*
                          ( -0.5*                                                               //-1/2 T1
                            Prop_ud[ prop_slv_idx(b, gamma,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(a, ALPHA,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(d, beta ,  bP, gammaP,  X_ixyz,it) ]
                            -0.5*                                                               //-1/2 T2
                            Prop_ud[ prop_slv_idx(a, ALPHA,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(d, beta ,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(b, gamma,  bP, gammaP,  X_ixyz,it) ]
                            +                                                              //+1 T3
                            Prop_ud[ prop_slv_idx(d, beta ,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(b, gamma,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(a, ALPHA,  bP, gammaP,  X_ixyz,it) ]
                            -                                                              //+1 T4
                            Prop_ud[ prop_slv_idx(d, beta ,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(a, ALPHA,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(b, gamma,  bP, gammaP,  X_ixyz,it) ]
                            +0.5*                                                              //-1/2 T5
                            Prop_ud[ prop_slv_idx(a, ALPHA,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(b, gamma,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(d, beta ,  bP, gammaP,  X_ixyz,it) ]
                            +0.5*                                                              //-1/2 T6
                            Prop_ud[ prop_slv_idx(b, gamma,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(d, beta ,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(a, ALPHA,  bP, gammaP,  X_ixyz,it) ]
                            
                          )
                          ;
            }//X_ixyz
           

            sum_sink += back_prop_contraction[prop_slv_cs_idx(d, alpha, dP, betaP)] *
                        sum_X_ixyz * 
                        ZGM(alpha,5) * Eps(3,color) * zcg5[gamma];
          }}}}//sink


          sum_source += sum_sink *
                        ZGM(alphaP,5) * Eps(3,colorP) * zcg5[gammaP];
                        
        }}}}//source
      
        sum_freeDI += 0.5 * sum_source;
      }//ALPHA
    
      sum_N += sum_freeDI/N_noises;
    }//i_noise

 
  ((COMPLEX*)correlator_local)[it + TnodeSites * TnodeCoor] = sum_N;
    
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
// calculate pi-sigma propagator --- the loop part 
//       formulas in notes 
//
void class_two_hadrons::run_GF_pi_sigma_loop(double* correlator){

  // complexify propagators 
  COMPLEX* Prop_ud            = (COMPLEX*)prop_ud        ;//+ prop_slv_idx(0,0,0,0,ixyz,it); 
  COMPLEX* Prop_s             = (COMPLEX*)prop_s         ;//+ prop_slv_idx(0,0,0,0,ixyz,it);
  COMPLEX* Prop_noise_full    = (COMPLEX*)prop_noise     ;//+ prop_slv_idx(0,0,0,0,ixyz,it);


  double correlator_local[2*Tsites];
  memset(correlator_local,0,sizeof(correlator_local));
  
  double correlator_global[2*Tsites];
  memset(correlator_global,0,sizeof(correlator_global));


  // calculating source contraction first
  COMPLEX src_ctr_local[3*4*3*4];
  memset(src_ctr_local,0,sizeof(src_ctr_local));
  COMPLEX src_ctr_global[3*4*3*4];
  memset(src_ctr_global,0,sizeof(src_ctr_global));

    
  int iT_src_pos=(iT_src+100*Tsites) % Tsites;
  if (iT_src_pos/TnodeSites == TnodeCoor) {

    int it = iT_src_pos % TnodeSites;

    for(int a=0;a<3;a++){
    for(int alf=0;alf<4;alf++){
    for(int aP=0; aP<3;aP++){
    for(int alfP=0; alfP<4;alfP++){
      
      int index=prop_slv_cs_idx(a, alf, aP, alfP);
      //printf("MPI %i OMP %i ... iT_src %i it %i  .... index %3i .... %i-%i-%i-%i \n", 
      //        MPI_rank,omp_get_thread_num(), iT_src_pos, it, index ,a ,alf,aP,alfP);
      COMPLEX sum =COMPLEX_ZERO;

      for(int Z_ixyz = 0;  Z_ixyz < XYZnodeSites; Z_ixyz++){            
        src_ctr_local[index] += Prop_ud[ prop_slv_idx(a, alf,  aP, alfP,  Z_ixyz,it) ];
        
      }
      src_ctr_local[index] /= XYZsites;

    }}}}

  }

  MPI_Allreduce(src_ctr_local, src_ctr_global, 
                288, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				
  // correlator itself
  #pragma omp parallel for
  for(int it = 0; it < TnodeSites; it++){


    // noise summation
    COMPLEX sum_N = COMPLEX_ZERO;
    for(int i_noise = 0; i_noise < N_noises; i_noise++){
    
      COMPLEX* Noise      = (COMPLEX*)sources->get_noise_ixyz(i_noise);
      COMPLEX* Prop_noise = Prop_noise_full + i_noise * XYZTnodeSites * 3*4*3*4;

      // free Dirac index summation
      COMPLEX sum_freeDI = COMPLEX_ZERO;
      for(int ALPHA = 0; ALPHA < 2; ALPHA++){

        // source summation
        COMPLEX sum_source = COMPLEX_ZERO;
        for(int dP      = 0; dP     < 3; dP++){
        for(int alphaP  = 0; alphaP < 4; alphaP++){
          int betaP=IGM(alphaP,5);
        for(int colorP  = 0; colorP < 6; colorP++){
          int aP    =Eps(0,colorP);
          int bP    =Eps(1,colorP);
          int cP    =Eps(2,colorP);          
        for(int gammaP  = 0; gammaP < 4; gammaP++){
          int deltaP = icg5[gammaP]; 


          // sink summation
          COMPLEX sum_sink = COMPLEX_ZERO;
          for(int d       = 0; d      < 3; d ++){
          for(int alpha   = 0; alpha  < 4; alpha++){
            int beta=IGM(alpha,5);
          for(int color   = 0; color  < 6; color++){
            int a    =Eps(0,color);
            int b    =Eps(1,color);
            int c    =Eps(2,color);          
          for(int gamma   = 0; gamma  < 4; gamma++){
            int delta = icg5[gamma]; 


            //sum X_ixyz
            COMPLEX sum_X_ixyz = COMPLEX_ZERO;
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
            
              sum_X_ixyz+=1.5* Conj(Noise[X_ixyz]) *
                          Prop_s[ prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it)  ]*
                          (                                                               
                            src_ctr_global[ prop_slv_cs_idx(dP, betaP, aP, ALPHA) ]*
                            ( -
                              Prop_noise    [ prop_slv_idx   (b, gamma, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (a, ALPHA, dP, alphaP,  X_ixyz,it) ]*  //-3/2 L7
                                Prop_ud       [ prop_slv_idx   (d, beta , bP, gammaP,  X_ixyz,it) ]
                                +
                                Prop_ud       [ prop_slv_idx   (d, beta , dP, alphaP,  X_ixyz,it) ]*  //+3/2 L8
                                Prop_ud       [ prop_slv_idx   (a, ALPHA, bP, gammaP,  X_ixyz,it) ]                                
                              )
                              +
                              Prop_noise    [ prop_slv_idx   (a, ALPHA, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (b, gamma, dP, alphaP,  X_ixyz,it) ]*  //-3/2 L9
                                Prop_ud       [ prop_slv_idx   (d, beta , bP, gammaP,  X_ixyz,it) ]
                                +
                                Prop_ud       [ prop_slv_idx   (d, beta , dP, alphaP,  X_ixyz,it) ]*  //+3/2 L10
                                Prop_ud       [ prop_slv_idx   (b, gamma, bP, gammaP,  X_ixyz,it) ]                                
                              )
                            )
                            +
                            src_ctr_global[ prop_slv_cs_idx(dP, betaP, bP, gammaP) ]*
                            ( +
                              Prop_noise    [ prop_slv_idx   (b, gamma, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (a, ALPHA, dP, alphaP,  X_ixyz,it) ]*  //-3/2 L13
                                Prop_ud       [ prop_slv_idx   (d, beta , aP, ALPHA ,  X_ixyz,it) ]
                                +
                                Prop_ud       [ prop_slv_idx   (d, beta , dP, alphaP,  X_ixyz,it) ]*  //+3/2 L14
                                Prop_ud       [ prop_slv_idx   (a, ALPHA, aP, ALPHA ,  X_ixyz,it) ]                                
                              )
                              -
                              Prop_noise    [ prop_slv_idx   (a, ALPHA, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (b, gamma, dP, alphaP,  X_ixyz,it) ]*  //-3/2 L15
                                Prop_ud       [ prop_slv_idx   (d, beta , aP, ALPHA ,  X_ixyz,it) ]
                                +
                                Prop_ud       [ prop_slv_idx   (d, beta , dP, alphaP,  X_ixyz,it) ]*  //+3/2 L16
                                Prop_ud       [ prop_slv_idx   (b, gamma, aP, ALPHA ,  X_ixyz,it) ]                                
                              )
                            )
                          )
                          ;
            }//X_ixyz
            

            sum_sink += sum_X_ixyz * 
                        ZGM(alpha,5) * Eps(3,color) * zcg5[gamma];
          }}}}//sink


          sum_source += sum_sink *
                        ZGM(alphaP,5) * Eps(3,colorP) * zcg5[gammaP];
                        
        }}}}//source
      
        sum_freeDI += 0.5* sum_source;
      }//ALPHA
    
      sum_N += sum_freeDI/N_noises;
    }//noise_i

 
  ((COMPLEX*)correlator_local)[it + TnodeSites * TnodeCoor] = sum_N;
    
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
// print the correlator
//
void class_two_hadrons::corr_print(double *correlator, string hadron_names){

  MPI_Barrier(MPI_COMM_WORLD);
  
  if(MPI_rank==0){

    // create directory and file name
    string dir="results/"+base+"/"+prefix+"correlators_two_hardons/";
    //printf("%s %s\n",prefix.c_str(),dir.c_str());
    create_directory(dir);
  
    char wfile[256];
    snprintf(wfile,sizeof(wfile), "%s%s.%02d",
             dir.c_str(),
             hadron_names.c_str(),
             iT_src);

    printf(" ++++++ corr_print : print %10s propagator to file \n                     %s\n                     %s\n",
           hadron_names.c_str(), wfile, LocalTime().c_str());
 
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




//****************###############**************#################*****************###############***
//****************###############**************#################*****************###############***
// test stuff

void class_two_hadrons::run_GF_pi_sigma_tree_TEST(double* correlator){

  // complexify propagators 
  COMPLEX* Prop_ud            = (COMPLEX*)prop_ud        ;//+ prop_slv_idx(0,0,0,0,ixyz,it); 
  COMPLEX* Prop_s             = (COMPLEX*)prop_s         ;//+ prop_slv_idx(0,0,0,0,ixyz,it);

#define back_prop(prop,c,a,cp,ap,ixyz,it)                               \
        ( ZGM(a,5) * ZGM (IGM(ap,5),5) *                                 \
          Conj(prop[ prop_slv_idx(c,IGM(a,5),cp,IGM(ap,5) ,ixyz,it) ]) )

  double correlator_local[2*Tsites];
  memset(correlator_local,0,sizeof(correlator_local));
  
  double correlator_global[2*Tsites];
  memset(correlator_global,0,sizeof(correlator_global));
  
  // correlator itself
  #pragma omp parallel for
  for(int it = 0; it < TnodeSites; it++){

    //if(MPI_rank==0){
    //  printf("OMP= %2d   time %2d\t\t\t%s",omp_get_thread_num(),it, LocalTime().c_str());
    //}
    
    // noise summation
    COMPLEX sum_N = COMPLEX_ZERO;
    for(int i_noise = 0; i_noise < N_noises; i_noise++){
    
      //if(MPI_rank==0){
      //  printf("OMP= %2d   time %2d   noise %2d \n", omp_get_thread_num(),it, i_noise);
      //}
      COMPLEX* Noise      = (COMPLEX*)sources->get_noise_ixyz(i_noise);

    //if(MPI_rank==0){
    //  printf("OMP= %d\ttime %d\t noise %d    after\n",omp_get_thread_num(),it, i_noise);
    //}


      // free Dirac index summation
      COMPLEX sum_freeDI = COMPLEX_ZERO;
      for(int ALPHA = 0; ALPHA < 2; ALPHA++){

    //if(MPI_rank==0){
    //  printf("OMP= %d\ttime %d\t noise %d\t dirac %d\n",omp_get_thread_num(),it, i_noise,ALPHA);
    //}


        // source summation
        COMPLEX sum_source = COMPLEX_ZERO;
        for(int dP      = 0; dP     < 3; dP++){
        for(int alphaP  = 0; alphaP < 4; alphaP++){
          int betaP=IGM(alphaP,5);
        for(int colorP  = 0; colorP < 6; colorP++){
          int aP    =Eps(0,colorP);
          int bP    =Eps(1,colorP);
          int cP    =Eps(2,colorP);          
        for(int gammaP  = 0; gammaP < 4; gammaP++){
          int deltaP = icg5[gammaP]; 
          
    //if(MPI_rank==0){
    //  printf("OMP= %d\ttime %d\t noise %d\t dirac %d\t source %d %d %d %d\n",omp_get_thread_num(),it, i_noise,ALPHA, dP,alphaP,colorP,gammaP);
    //}

          
          
          // sink summation
          COMPLEX sum_sink = COMPLEX_ZERO;
          for(int d       = 0; d      < 3; d ++){
          for(int alpha   = 0; alpha  < 4; alpha++){
            int beta=IGM(alpha,5);
          for(int color   = 0; color  < 6; color++){
            int a    =Eps(0,color);
            int b    =Eps(1,color);
            int c    =Eps(2,color);          
          for(int gamma   = 0; gamma  < 4; gamma++){
            int delta = icg5[gamma]; 

    //if(MPI_rank==0){
    //  printf("                     sink %d %d %d %d\n", d,alpha,color,gamma);
    //}

            //sum Y_ixyz
            COMPLEX sum_Y_ixyz = COMPLEX_ZERO;
            for(int Y_ixyz = 0;  Y_ixyz < XYZnodeSites; Y_ixyz++){
              sum_Y_ixyz += Noise[Y_ixyz] * 
                            back_prop(Prop_ud,       d, alpha, dP, betaP, Y_ixyz, it);
            }//Y_ixyz
            
            
            //sum X_ixyz
            COMPLEX sum_X_ixyz = COMPLEX_ZERO;
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              sum_X_ixyz+=Conj(Noise[X_ixyz]) * 
                          Prop_s[ prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ]*
                          ( -0.5*                                                               //-1/2 T1
                            Prop_ud[ prop_slv_idx(b, gamma,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(a, ALPHA,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(d, beta ,  bP, gammaP,  X_ixyz,it) ]
                            -0.5*                                                               //-1/2 T2
                            Prop_ud[ prop_slv_idx(a, ALPHA,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(d, beta ,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(b, gamma,  bP, gammaP,  X_ixyz,it) ]
                            +                                                              //+1 T3
                            Prop_ud[ prop_slv_idx(d, beta ,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(b, gamma,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(a, ALPHA,  bP, gammaP,  X_ixyz,it) ]
                            -                                                              //+1 T4
                            Prop_ud[ prop_slv_idx(d, beta ,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(a, ALPHA,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(b, gamma,  bP, gammaP,  X_ixyz,it) ]
                            +0.5*                                                              //-1/2 T5
                            Prop_ud[ prop_slv_idx(a, ALPHA,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(b, gamma,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(d, beta ,  bP, gammaP,  X_ixyz,it) ]
                            +0.5*                                                              //-1/2 T6
                            Prop_ud[ prop_slv_idx(b, gamma,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(d, beta ,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(a, ALPHA,  bP, gammaP,  X_ixyz,it) ]
                            
                          )
                           ;
            }//X_ixyz
           

            sum_sink += sum_X_ixyz * sum_Y_ixyz*
                        ZGM(alpha,5) * Eps(3,color) * zcg5[gamma];
          }}}}//sink


          sum_source += sum_sink *
                        ZGM(alphaP,5) * Eps(3,colorP) * zcg5[gammaP];
                        
        }}}}//source
      
        sum_freeDI += 0.5 * sum_source;
      }//ALPHA
    
      sum_N += sum_freeDI/N_noises;
    }//i_noise


//    printf("MPI %i OMP %i ... it %i sum %1.16e %1.16e I\n", 
//           MPI_rank,omp_get_thread_num(), it,Real(sum_ixyz),Imag(sum_ixyz));
 
  ((COMPLEX*)correlator_local)[it + TnodeSites * TnodeCoor]=sum_N;
    
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


void class_two_hadrons::run_GF_pi_sigma_loop_TEST(double* correlator, int AA, int diag){

  // complexify propagators 
  COMPLEX* Prop_ud            = (COMPLEX*)prop_ud        ;//+ prop_slv_idx(0,0,0,0,ixyz,it); 
  COMPLEX* Prop_s             = (COMPLEX*)prop_s         ;//+ prop_slv_idx(0,0,0,0,ixyz,it);
  COMPLEX* Prop_noise_full    = (COMPLEX*)prop_noise     ;//+ prop_slv_idx(0,0,0,0,ixyz,it);



  double correlator_local[2*Tsites];
  memset(correlator_local,0,sizeof(correlator_local));
  
  double correlator_global[2*Tsites];
  memset(correlator_global,0,sizeof(correlator_global));


  // calculating source contraction first
  COMPLEX src_ctr_local[3*4*3*4];
  memset(src_ctr_local,0,sizeof(src_ctr_local));
  COMPLEX src_ctr_global[3*4*3*4];
  memset(src_ctr_global,0,sizeof(src_ctr_global));

    
  int iT_src_pos=(iT_src+100*Tsites) % Tsites;
  if (iT_src_pos/TnodeSites == TnodeCoor) {

    int it = iT_src_pos % TnodeSites;

    for(int a=0;a<3;a++){
    for(int alf=0;alf<4;alf++){
    for(int aP=0; aP<3;aP++){
    for(int alfP=0; alfP<4;alfP++){
      
      int index=prop_slv_cs_idx(a, alf, aP, alfP);
      //printf("MPI %i OMP %i ... iT_src %i it %i  .... index %3i .... %i-%i-%i-%i \n", 
      //        MPI_rank,omp_get_thread_num(), iT_src_pos, it, index ,a ,alf,aP,alfP);
      COMPLEX sum =COMPLEX_ZERO;

      for(int Z_ixyz = 0;  Z_ixyz < XYZnodeSites; Z_ixyz++){            
        src_ctr_local[index] += Prop_ud[ prop_slv_idx(a, alf,  aP, alfP,  Z_ixyz,it) ];
        
      }
      src_ctr_local[index] /= XYZsites;

    }}}}

  }

  MPI_Allreduce(src_ctr_local, src_ctr_global, 
                288, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //for(int mpi=0;mpi<2;mpi++){
  //MPI_Barrier(MPI_COMM_WORLD);
  //if(MPI_rank==mpi){

  //for(int a=0;a<3;a++){
  //for(int alf=0;alf<4;alf++){
  //for(int aP=0; aP<3;aP++){
  //for(int alfP=0; alfP<4;alfP++){
  //    
  //  int index=prop_slv_cs_idx(a, alf, aP, alfP);
  //  printf("MPI %i  corr_glob %1.16e %1.16e I \tindex %3i .... %i-%i-%i-%i \n", 
  //        MPI_rank, Real(source_contraction_global[index]),  Imag(source_contraction_global[index]),
  //        index ,a ,alf,aP,alfP);

  //for(int itt=0;itt<12;itt++){
  
  //  COMPLEX summ =COMPLEX_ZERO;

  //  for(int Z_ixyz = 0;  Z_ixyz < XYZnodeSites; Z_ixyz++){            
  //      COMPLEX num=COMPLEX_ZERO;
  //      for(int iii=0;iii<N_noises;iii++){
  //        COMPLEX* src = (COMPLEX*)sources->get_noise_ixyz(iii);
  //        num+=Prop_noise_full [iii* XYZTnodeSites * 3*4*3*4 +prop_slv_idx(a, alf, aP, alfP,  Z_ixyz,itt) ]*
  //             Conj(src[Z_ixyz]);
  //      }
            
  //      summ += num; 
  //    }

  //  printf("nn= %4i MPI %i, it= %2i,  %i-%i-%i-%i  corr_nois sum %1.16e  %1.16e I\n", 
  //          N_noises, MPI_rank, itt+12*MPI_rank, a ,alf,aP,alfP, Real(summ)/(N_noises*XYZsites),Imag(summ)/(N_noises*XYZsites) );
  
  //}
  //}}}}
  //}}
  //MPI_Barrier(MPI_COMM_WORLD);
  
  int Ast, Afin;
  if(AA==2){
   Ast=0;
   Afin=2;
  }
  else if(AA==3){
   Ast=2;
   Afin=4;
  }
  else{
   Ast=0;
   Afin=4;
  }

  
  // correlator itself
  #pragma omp parallel for
  for(int it = 0; it < TnodeSites; it++){

    //if(MPI_rank==0){
    //  printf("OMP= %2d   time %2d\t\t\t%s",omp_get_thread_num(),it, LocalTime().c_str());
    //}

    // noise summation
    COMPLEX sum_N = COMPLEX_ZERO;
    for(int i_noise = 0; i_noise < N_noises; i_noise++){
    
      COMPLEX* Noise      = (COMPLEX*)sources->get_noise_ixyz(i_noise);
      COMPLEX* Prop_noise = Prop_noise_full + i_noise * XYZTnodeSites * 3*4*3*4;

      // free Dirac index summation
      COMPLEX sum_freeDI = COMPLEX_ZERO;
      for(int ALPHA = Ast; ALPHA < Afin; ALPHA++){

        // source summation
        COMPLEX sum_source = COMPLEX_ZERO;
        for(int dP      = 0; dP     < 3; dP++){
        for(int alphaP  = 0; alphaP < 4; alphaP++){
          int betaP=IGM(alphaP,5);
        for(int colorP  = 0; colorP < 6; colorP++){
          int aP    =Eps(0,colorP);
          int bP    =Eps(1,colorP);
          int cP    =Eps(2,colorP);          
        for(int gammaP  = 0; gammaP < 4; gammaP++){
          int deltaP = icg5[gammaP]; 
          
    //if(MPI_rank==0){
    //  printf("OMP= %d\ttime %d\t noise %d\t dirac %d\t source %d %d %d %d\n",omp_get_thread_num(),it, i_noise,ALPHA, dP,alphaP,colorP,gammaP);
    //}

          
          
          // sink summation
          COMPLEX sum_sink = COMPLEX_ZERO;
          for(int d       = 0; d      < 3; d ++){
          for(int alpha   = 0; alpha  < 4; alpha++){
            int beta=IGM(alpha,5);
          for(int color   = 0; color  < 6; color++){
            int a    =Eps(0,color);
            int b    =Eps(1,color);
            int c    =Eps(2,color);          
          for(int gamma   = 0; gamma  < 4; gamma++){
            int delta = icg5[gamma]; 

    //if(MPI_rank==0){
    //  printf("                     sink %d %d %d %d\n", d,alpha,color,gamma);
    //}

            //sum X_ixyz
            COMPLEX sum_X_ixyz = COMPLEX_ZERO;
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
            
              if(diag==7){ // 7 
              sum_X_ixyz+=1.5* Conj(Noise[X_ixyz]) * 
                          Prop_s[ prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it)  ]*
                          (                                                               
                            src_ctr_global[ prop_slv_cs_idx(dP, betaP, aP, ALPHA) ]*
                            ( -
                              Prop_noise    [ prop_slv_idx   (b, gamma, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (a, ALPHA, dP, alphaP,  X_ixyz,it) ]*  //-3/2 L7
                                Prop_ud       [ prop_slv_idx   (d, beta , bP, gammaP,  X_ixyz,it) ]
                             )
                            )
                          )
                          ;
              }
              if(diag==8){ // 7 8
              sum_X_ixyz+=1.5* Conj(Noise[X_ixyz]) * 
                          Prop_s[ prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it)  ]*
                          (                                                               
                            src_ctr_global[ prop_slv_cs_idx(dP, betaP, aP, ALPHA) ]*
                            ( -
                              Prop_noise    [ prop_slv_idx   (b, gamma, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (d, beta , dP, alphaP,  X_ixyz,it) ]*  //+3/2 L8
                                Prop_ud       [ prop_slv_idx   (a, ALPHA, bP, gammaP,  X_ixyz,it) ]                                
                              )
                            )
                          )
                          ;
              }
              if(diag==9){ // 9
              sum_X_ixyz+=1.5* Conj(Noise[X_ixyz]) * 
                          Prop_s[ prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it)  ]*
                          (                                                               
                            src_ctr_global[ prop_slv_cs_idx(dP, betaP, aP, ALPHA) ]*
                            (+
                              Prop_noise    [ prop_slv_idx   (a, ALPHA, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (b, gamma, dP, alphaP,  X_ixyz,it) ]*  //-3/2 L9
                                Prop_ud       [ prop_slv_idx   (d, beta , bP, gammaP,  X_ixyz,it) ]
                              )
                            )
                          )
                          ;
              }
              if(diag==10){ // 10
              sum_X_ixyz+=1.5* Conj(Noise[X_ixyz]) * 
                          Prop_s[ prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it)  ]*
                          (                                                               
                            src_ctr_global[ prop_slv_cs_idx(dP, betaP, aP, ALPHA) ]*
                            (+
                              Prop_noise    [ prop_slv_idx   (a, ALPHA, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (d, beta , dP, alphaP,  X_ixyz,it) ]*  //+3/2 L10
                                Prop_ud       [ prop_slv_idx   (b, gamma, bP, gammaP,  X_ixyz,it) ]                                
                              )
                            )
                          )
                          ;
              }
              if(diag==13){ // 13 
              sum_X_ixyz+=1.5* Conj(Noise[X_ixyz]) * 
                          Prop_s[ prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it)  ]*
                          (                                                               
                            +
                            src_ctr_global[ prop_slv_cs_idx(dP, betaP, bP, gammaP) ]*
                            ( +
                              Prop_noise    [ prop_slv_idx   (b, gamma, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (a, ALPHA, dP, alphaP,  X_ixyz,it) ]*  //-3/2 L13
                                Prop_ud       [ prop_slv_idx   (d, beta , aP, ALPHA ,  X_ixyz,it) ]
                              )
                            )
                          )
                          ;
              }
              if(diag==14){ //  14
              sum_X_ixyz+=1.5* Conj(Noise[X_ixyz]) * 
                          Prop_s[ prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it)  ]*
                          (                                                               
                            +
                            src_ctr_global[ prop_slv_cs_idx(dP, betaP, bP, gammaP) ]*
                            ( +
                              Prop_noise    [ prop_slv_idx   (b, gamma, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (d, beta , dP, alphaP,  X_ixyz,it) ]*  //+3/2 L14
                                Prop_ud       [ prop_slv_idx   (a, ALPHA, aP, ALPHA ,  X_ixyz,it) ]                                
                              )
                            )
                          )
                          ;
              }
              if(diag==15){ // 15
              sum_X_ixyz+=1.5* Conj(Noise[X_ixyz]) * 
                          Prop_s[ prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it)  ]*
                          (                                                               
                            +
                            src_ctr_global[ prop_slv_cs_idx(dP, betaP, bP, gammaP) ]*
                            ( -
                              Prop_noise    [ prop_slv_idx   (a, ALPHA, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (b, gamma, dP, alphaP,  X_ixyz,it) ]*  //-3/2 L15
                                Prop_ud       [ prop_slv_idx   (d, beta , aP, ALPHA ,  X_ixyz,it) ]
                              )
                            )
                          )
                          ;
              }
              if(diag==16){ // 16 
              sum_X_ixyz+=1.5* Conj(Noise[X_ixyz]) * 
                          Prop_s[ prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it)  ]*
                          (                                                               
                            +
                            src_ctr_global[ prop_slv_cs_idx(dP, betaP, bP, gammaP) ]*
                            ( -
                              Prop_noise    [ prop_slv_idx   (a, ALPHA, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (d, beta , dP, alphaP,  X_ixyz,it) ]*  //+3/2 L16
                                Prop_ud       [ prop_slv_idx   (b, gamma, aP, ALPHA ,  X_ixyz,it) ]                                
                              )
                            )
                          )
                          ;
              }
            }//X_ixyz
            

            sum_sink += sum_X_ixyz * 
                        ZGM(alpha,5) * Eps(3,color) * zcg5[gamma];
          }}}}//sink


          sum_source += sum_sink *
                        ZGM(alphaP,5) * Eps(3,colorP) * zcg5[gammaP];
                        
        }}}}//source
      
        sum_freeDI += 0.5* sum_source;
      }//ALPHA
    
      sum_N += sum_freeDI/N_noises;
    }//noise_i


//    printf("MPI %i OMP %i ... it %i sum %1.16e %1.16e I\n", 
//           MPI_rank,omp_get_thread_num(), it,Real(sum_ixyz),Imag(sum_ixyz));
 
  ((COMPLEX*)correlator_local)[it + TnodeSites * TnodeCoor] = sum_N;
    
  } // it, end of omp parallel

  // reduce from all MPI processes
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(correlator_local, correlator_global, 2*Tsites, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  // output correlator
  for(int it = 0; it < Tsites; it++){
    ((COMPLEX*)correlator)[it]=((COMPLEX*)correlator_global)[it];
  }
  

}


void class_two_hadrons::run_GF_pi_sigma_tree_NONOISE(double* correlator){

  // complexify propagators 
  COMPLEX* Prop_ud            = (COMPLEX*)prop_ud        ;//+ prop_slv_idx(0,0,0,0,ixyz,it); 
  COMPLEX* Prop_s             = (COMPLEX*)prop_s         ;//+ prop_slv_idx(0,0,0,0,ixyz,it);

#define back_prop(prop,c,a,cp,ap,ixyz,it)                               \
        ( ZGM(a,5) * ZGM (IGM(ap,5),5) *                                 \
          Conj(prop[ prop_slv_idx(c,IGM(a,5),cp,IGM(ap,5) ,ixyz,it) ]) )

  double correlator_local[2*Tsites];
  memset(correlator_local,0,sizeof(correlator_local));
  
  double correlator_global[2*Tsites];
  memset(correlator_global,0,sizeof(correlator_global));
  
  // correlator itself
  #pragma omp parallel for
  for(int it = 0; it < TnodeSites; it++){

    //if(MPI_rank==0){
    //  printf("OMP= %2d   time %2d\t\t\t%s",omp_get_thread_num(),it, LocalTime().c_str());
    //}
    
    // noise summation
//    COMPLEX sum_N = COMPLEX_ZERO;
//    for(int i_noise = 0; i_noise < N_noises; i_noise++){
    
      //if(MPI_rank==0){
      //  printf("OMP= %2d   time %2d   noise %2d \n", omp_get_thread_num(),it, i_noise);
      //}
//      COMPLEX* Noise      = (COMPLEX*)sources->get_noise_ixyz(i_noise);

    //if(MPI_rank==0){
    //  printf("OMP= %d\ttime %d\t noise %d    after\n",omp_get_thread_num(),it, i_noise);
    //}


      // free Dirac index summation
      COMPLEX sum_freeDI = COMPLEX_ZERO;
      for(int ALPHA = 0; ALPHA < 4; ALPHA++){

    //if(MPI_rank==0){
    //  printf("OMP= %d\ttime %d\t noise %d\t dirac %d\n",omp_get_thread_num(),it, i_noise,ALPHA);
    //}


        // source summation
        COMPLEX sum_source = COMPLEX_ZERO;
        for(int dP      = 0; dP     < 3; dP++){
        for(int alphaP  = 0; alphaP < 4; alphaP++){
          int betaP=IGM(alphaP,5);
        for(int colorP  = 0; colorP < 6; colorP++){
          int aP    =Eps(0,colorP);
          int bP    =Eps(1,colorP);
          int cP    =Eps(2,colorP);          
        for(int gammaP  = 0; gammaP < 4; gammaP++){
          int deltaP = icg5[gammaP]; 
          
    //if(MPI_rank==0){
    //  printf("OMP= %d\ttime %d\t noise %d\t dirac %d\t source %d %d %d %d\n",omp_get_thread_num(),it, i_noise,ALPHA, dP,alphaP,colorP,gammaP);
    //}

          
          
          // sink summation
          COMPLEX sum_sink = COMPLEX_ZERO;
          for(int d       = 0; d      < 3; d ++){
          for(int alpha   = 0; alpha  < 4; alpha++){
            int beta=IGM(alpha,5);
          for(int color   = 0; color  < 6; color++){
            int a    =Eps(0,color);
            int b    =Eps(1,color);
            int c    =Eps(2,color);          
          for(int gamma   = 0; gamma  < 4; gamma++){
            int delta = icg5[gamma]; 

    //if(MPI_rank==0){
    //  printf("                     sink %d %d %d %d\n", d,alpha,color,gamma);
    //}

            
            
            //sum X_ixyz
            COMPLEX sum_X_ixyz = COMPLEX_ZERO;
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              sum_X_ixyz+=back_prop(Prop_ud,       d, alpha, dP, betaP, X_ixyz, it) * 
                          Prop_s[ prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ]*
                          ( -0.5*                                                               //-1/2 T1
                            Prop_ud[ prop_slv_idx(b, gamma,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(a, ALPHA,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(d, beta ,  bP, gammaP,  X_ixyz,it) ]
                            -0.5*                                                               //-1/2 T2
                            Prop_ud[ prop_slv_idx(a, ALPHA,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(d, beta ,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(b, gamma,  bP, gammaP,  X_ixyz,it) ]
                            +                                                              //+1 T3
                            Prop_ud[ prop_slv_idx(d, beta ,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(b, gamma,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(a, ALPHA,  bP, gammaP,  X_ixyz,it) ]
                            -                                                              //+1 T4
                            Prop_ud[ prop_slv_idx(d, beta ,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(a, ALPHA,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(b, gamma,  bP, gammaP,  X_ixyz,it) ]
                            +0.5*                                                              //-1/2 T5
                            Prop_ud[ prop_slv_idx(a, ALPHA,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(b, gamma,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(d, beta ,  bP, gammaP,  X_ixyz,it) ]
                            +0.5*                                                              //-1/2 T6
                            Prop_ud[ prop_slv_idx(b, gamma,  dP, alphaP,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(d, beta ,  aP, ALPHA ,  X_ixyz,it) ]*
                            Prop_ud[ prop_slv_idx(a, ALPHA,  bP, gammaP,  X_ixyz,it) ]
                            
                          )
                           ;
            }//X_ixyz
           

            sum_sink += sum_X_ixyz * 
                        ZGM(alpha,5) * Eps(3,color) * zcg5[gamma];
          }}}}//sink


          sum_source += sum_sink *
                        ZGM(alphaP,5) * Eps(3,colorP) * zcg5[gammaP];
                        
        }}}}//source
      
        sum_freeDI += 0.5 * sum_source;
      }//ALPHA
    
//      sum_N += sum_freeDI/N_noises;
//    }//i_noise


//    printf("MPI %i OMP %i ... it %i sum %1.16e %1.16e I\n", 
//           MPI_rank,omp_get_thread_num(), it,Real(sum_ixyz),Imag(sum_ixyz));
 
  ((COMPLEX*)correlator_local)[it + TnodeSites * TnodeCoor]=sum_freeDI;
    
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


void class_two_hadrons::run_GF_pi_sigma_loop_NONOISE(double* correlator){

  // complexify propagators 
  COMPLEX* Prop_ud            = (COMPLEX*)prop_ud        ;//+ prop_slv_idx(0,0,0,0,ixyz,it); 
  COMPLEX* Prop_s             = (COMPLEX*)prop_s         ;//+ prop_slv_idx(0,0,0,0,ixyz,it);
  COMPLEX* Prop_noise_full    = (COMPLEX*)prop_noise     ;//+ prop_slv_idx(0,0,0,0,ixyz,it);



  double correlator_local[2*Tsites];
  memset(correlator_local,0,sizeof(correlator_local));
  
  double correlator_global[2*Tsites];
  memset(correlator_global,0,sizeof(correlator_global));


  // calculating source contraction first
  COMPLEX src_ctr_local[3*4*3*4];
  memset(src_ctr_local,0,sizeof(src_ctr_local));
  COMPLEX src_ctr_global[3*4*3*4];
  memset(src_ctr_global,0,sizeof(src_ctr_global));

    
  int iT_src_pos=(iT_src+100*Tsites) % Tsites;
  if (iT_src_pos/TnodeSites == TnodeCoor) {

    int it = iT_src_pos % TnodeSites;

    for(int a=0;a<3;a++){
    for(int alf=0;alf<4;alf++){
    for(int aP=0; aP<3;aP++){
    for(int alfP=0; alfP<4;alfP++){
      
      int index=prop_slv_cs_idx(a, alf, aP, alfP);
      //printf("MPI %i OMP %i ... iT_src %i it %i  .... index %3i .... %i-%i-%i-%i \n", 
      //        MPI_rank,omp_get_thread_num(), iT_src_pos, it, index ,a ,alf,aP,alfP);
      COMPLEX sum =COMPLEX_ZERO;

      for(int Z_ixyz = 0;  Z_ixyz < XYZnodeSites; Z_ixyz++){            
        src_ctr_local[index] += Prop_ud[ prop_slv_idx(a, alf,  aP, alfP,  Z_ixyz,it) ];
        
      }
      src_ctr_local[index] /= XYZsites;

    }}}}

  }

  MPI_Allreduce(src_ctr_local, src_ctr_global, 
                288, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //for(int mpi=0;mpi<2;mpi++){
  //MPI_Barrier(MPI_COMM_WORLD);
  //if(MPI_rank==mpi){

  //for(int a=0;a<3;a++){
  //for(int alf=0;alf<4;alf++){
  //for(int aP=0; aP<3;aP++){
  //for(int alfP=0; alfP<4;alfP++){
  //    
  //  int index=prop_slv_cs_idx(a, alf, aP, alfP);
  //  printf("MPI %i  corr_glob %1.16e %1.16e I \tindex %3i .... %i-%i-%i-%i \n", 
  //        MPI_rank, Real(source_contraction_global[index]),  Imag(source_contraction_global[index]),
  //        index ,a ,alf,aP,alfP);

  //for(int itt=0;itt<12;itt++){
  
  //  COMPLEX summ =COMPLEX_ZERO;

  //  for(int Z_ixyz = 0;  Z_ixyz < XYZnodeSites; Z_ixyz++){            
  //      COMPLEX num=COMPLEX_ZERO;
  //      for(int iii=0;iii<N_noises;iii++){
  //        COMPLEX* src = (COMPLEX*)sources->get_noise_ixyz(iii);
  //        num+=Prop_noise_full [iii* XYZTnodeSites * 3*4*3*4 +prop_slv_idx(a, alf, aP, alfP,  Z_ixyz,itt) ]*
  //             Conj(src[Z_ixyz]);
  //      }
            
  //      summ += num; 
  //    }

  //  printf("nn= %4i MPI %i, it= %2i,  %i-%i-%i-%i  corr_nois sum %1.16e  %1.16e I\n", 
  //          N_noises, MPI_rank, itt+12*MPI_rank, a ,alf,aP,alfP, Real(summ)/(N_noises*XYZsites),Imag(summ)/(N_noises*XYZsites) );
  
  //}
  //}}}}
  //}}
  //MPI_Barrier(MPI_COMM_WORLD);
  
  // correlator itself
  #pragma omp parallel for
  for(int it = 0; it < TnodeSites; it++){

    //if(MPI_rank==0){
    //  printf("OMP= %2d   time %2d\t\t\t%s",omp_get_thread_num(),it, LocalTime().c_str());
    //}

    // noise summation
    
//      COMPLEX* Noise      = (COMPLEX*)sources->get_noise_ixyz(i_noise);
      COMPLEX* Point      = (COMPLEX*)sources->get_point_ixyz();
      COMPLEX* Prop_noise = Prop_noise_full;// + i_noise * XYZTnodeSites * 3*4*3*4;

      // free Dirac index summation
      COMPLEX sum_freeDI = COMPLEX_ZERO;
      for(int ALPHA = 0; ALPHA < 4; ALPHA++){

        // source summation
        COMPLEX sum_source = COMPLEX_ZERO;
        for(int dP      = 0; dP     < 3; dP++){
        for(int alphaP  = 0; alphaP < 4; alphaP++){
          int betaP=IGM(alphaP,5);
        for(int colorP  = 0; colorP < 6; colorP++){
          int aP    =Eps(0,colorP);
          int bP    =Eps(1,colorP);
          int cP    =Eps(2,colorP);          
        for(int gammaP  = 0; gammaP < 4; gammaP++){
          int deltaP = icg5[gammaP]; 
          
    //if(MPI_rank==0){
    //  printf("OMP= %d\ttime %d\t noise %d\t dirac %d\t source %d %d %d %d\n",omp_get_thread_num(),it, i_noise,ALPHA, dP,alphaP,colorP,gammaP);
    //}

          
          
          // sink summation
          COMPLEX sum_sink = COMPLEX_ZERO;
          for(int d       = 0; d      < 3; d ++){
          for(int alpha   = 0; alpha  < 4; alpha++){
            int beta=IGM(alpha,5);
          for(int color   = 0; color  < 6; color++){
            int a    =Eps(0,color);
            int b    =Eps(1,color);
            int c    =Eps(2,color);          
          for(int gamma   = 0; gamma  < 4; gamma++){
            int delta = icg5[gamma]; 

    //if(MPI_rank==0){
    //  printf("                     sink %d %d %d %d\n", d,alpha,color,gamma);
    //}

            //sum X_ixyz
            COMPLEX sum_X_ixyz = COMPLEX_ZERO;
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
            
              sum_X_ixyz+=1.5* Point[X_ixyz] * 
                          Prop_s[ prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it)  ]*
                          (                                                               
                            src_ctr_global[ prop_slv_cs_idx(dP, betaP, aP, ALPHA) ]*
                            ( -
                              Prop_noise    [ prop_slv_idx   (b, gamma, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (a, ALPHA, dP, alphaP,  X_ixyz,it) ]*  //-3/2 L7
                                Prop_ud       [ prop_slv_idx   (d, beta , bP, gammaP,  X_ixyz,it) ]
                                +
                                Prop_ud       [ prop_slv_idx   (d, beta , dP, alphaP,  X_ixyz,it) ]*  //+3/2 L8
                                Prop_ud       [ prop_slv_idx   (a, ALPHA, bP, gammaP,  X_ixyz,it) ]                                
                              )
                              +
                              Prop_noise    [ prop_slv_idx   (a, ALPHA, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (b, gamma, dP, alphaP,  X_ixyz,it) ]*  //-3/2 L9
                                Prop_ud       [ prop_slv_idx   (d, beta , bP, gammaP,  X_ixyz,it) ]
                                +
                                Prop_ud       [ prop_slv_idx   (d, beta , dP, alphaP,  X_ixyz,it) ]*  //+3/2 L10
                                Prop_ud       [ prop_slv_idx   (b, gamma, bP, gammaP,  X_ixyz,it) ]                                
                              )
                            )
                            +
                            src_ctr_global[ prop_slv_cs_idx(dP, betaP, bP, gammaP) ]*
                            ( +
                              Prop_noise    [ prop_slv_idx   (b, gamma, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (a, ALPHA, dP, alphaP,  X_ixyz,it) ]*  //-3/2 L13
                                Prop_ud       [ prop_slv_idx   (d, beta , aP, ALPHA ,  X_ixyz,it) ]
                                +
                                Prop_ud       [ prop_slv_idx   (d, beta , dP, alphaP,  X_ixyz,it) ]*  //+3/2 L14
                                Prop_ud       [ prop_slv_idx   (a, ALPHA, aP, ALPHA ,  X_ixyz,it) ]                                
                              )
                              -
                              Prop_noise    [ prop_slv_idx   (a, ALPHA, d, alpha,  X_ixyz,it) ]*
                              (
                                Prop_ud       [ prop_slv_idx   (b, gamma, dP, alphaP,  X_ixyz,it) ]*  //-3/2 L15
                                Prop_ud       [ prop_slv_idx   (d, beta , aP, ALPHA ,  X_ixyz,it) ]
                                +
                                Prop_ud       [ prop_slv_idx   (d, beta , dP, alphaP,  X_ixyz,it) ]*  //+3/2 L16
                                Prop_ud       [ prop_slv_idx   (b, gamma, aP, ALPHA ,  X_ixyz,it) ]                                
                              )
                            )
                          )
                          ;
            }//X_ixyz
           

            sum_sink += sum_X_ixyz * 
                        ZGM(alpha,5) * Eps(3,color) * zcg5[gamma];
          }}}}//sink


          sum_source += sum_sink *
                        ZGM(alphaP,5) * Eps(3,colorP) * zcg5[gammaP];
                        
        }}}}//source
      
        sum_freeDI += 0.5* sum_source;
      }//ALPHA
    


//    printf("MPI %i OMP %i ... it %i sum %1.16e %1.16e I\n", 
//           MPI_rank,omp_get_thread_num(), it,Real(sum_ixyz),Imag(sum_ixyz));
 
  ((COMPLEX*)correlator_local)[it + TnodeSites * TnodeCoor] = sum_freeDI;
    
  } // it, end of omp parallel

  // reduce from all MPI processes
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(correlator_local, correlator_global, 2*Tsites, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  // output correlator
  for(int it = 0; it < Tsites; it++){
    ((COMPLEX*)correlator)[it]=((COMPLEX*)correlator_global)[it];
  }
  

}
