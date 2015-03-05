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

void class_two_hadrons::run_all_NBSwf(){

  run_GF("pion-sigma");

}

// ================================================================================================
// wrapper to run all wafe functions
//
void class_two_hadrons::run_NBSwf(string hadron_names){

  double correlator[2*Tsites];
  memset(correlator,0,sizeof(correlator));

  if(MPI_rank==0){
    printf(" ++++++ run_NBSwf : calculate %s NBS wave function        %s\n",
         hadron_names.c_str(), LocalTime().c_str());
  }
  MPI_Barrier(MPI_COMM_WORLD);


  // ==================================
  // pion-sigma two baryon system
    
  if(hadron_names=="pion-sigma"){

    
    if(MPI_rank==0){
      printf("       ++++++ run_NBSwf : the TREE part        begin %s\n",
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    double correlator_tree[2*Tsites];
    memset(correlator_tree,0,sizeof(correlator_tree));  
    run_GF_pi_sigma_tree(correlator_tree);
    corr_print(correlator_tree     , hadron_names+"_tree"     );
    if(MPI_rank==0){
      printf("       ++++++ run_NBSwf :                        end %s\n", 
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);


    if(MPI_rank==0){
      printf("       ++++++ run_NBSwf : the LOOP part        begin %s\n",
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);

    double correlator_loop[2*Tsites];
    memset(correlator_loop,0,sizeof(correlator_loop));  
    run_GF_pi_sigma_loop(correlator_loop);
    corr_print(correlator_loop     , hadron_names+"_loop"     );
    if(MPI_rank==0){
      printf("       ++++++ run_NBSwf :                        end %s\n", 
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);

  }
  else{
  
    printf("ERROR - unknown hadron name");
    abort();
  }


  //corr_print(correlator, hadron_names);
  
};


// ================================================================================================
// calculate pi-sigma NBS wave function --- the tree part 
//       formulas in notes 
//
void class_two_hadrons::run_NBSwf_pi_sigma_tree(double* correlator){

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

      // calculating back propagator contraction first
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
// calculate pi-sigma NBS wave function --- the loop part 
//       formulas in notes 
//
void class_two_hadrons::run_NBSwf_pi_sigma_loop(double* correlator){

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
// print the NBS wave function
//
void class_two_hadrons::NBSwf_print(double *correlator, string hadron_names)
{
  
  for(int it=0; it<TnodeSites;it++){

	int iT_glob= it + TnodeSites * TnodeCoor;
	
    // create directory and file name
    char dir[256];
    string dir_base="results/"+base+"/"+prefix+"NBSwaveF_two_hardons/"+hadron_names;
    snprintf(dir,sizeof(dir), "%s.t%02d/",
             dir_base.c_str(),
             iT_glob);

    create_directory(dir);
    //printf("%s %s\n",prefix.c_str(),dir.c_str());
  
    char wfile[256];
    snprintf(wfile,sizeof(wfile), "%sx%02dy%02dz%02d",
             dir,XnodeCoor,YnodeCoor,ZnodeCoor);

    printf(" ++++++ NBSff_print : print %10s NBSdw at t=%i to file \n                     %s\n                     %s\n",
           hadron_names.c_str(), iT_glob, wfile, LocalTime().c_str());
 
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




