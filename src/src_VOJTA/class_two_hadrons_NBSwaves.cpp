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

/* local from EXTERNAL */
#include "fft3d.h"


using namespace HAL_idx;

void class_two_hadrons::run_all_NBSwf(){

  for(int t =0; t<TnodeSites;t++){
    run_NBSwf("pion-sigma",t);
  }
  
}

// ================================================================================================
// wrapper to run all wafe functions
//
void class_two_hadrons::run_NBSwf(string hadron_names, int time_WF){


  double wave_function[2*XYZnodeSites];
  memset(wave_function,0,sizeof(wave_function));  
  
  if(MPI_rank==0){
    printf(" ++++++ run_NBSwf : calculate %s NBS wave function for t_src=%2d and t_WF=%2d+i*%2d  %s\n",
           hadron_names.c_str(), iT_src, time_WF, Tnodes, LocalTime().c_str());
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // ==================================
  // pion-sigma two baryon system
    
  if(hadron_names=="pion-sigma"){
    

    // tree part
    if(MPI_rank==0){
      printf("        ++++++ run_NBSwf : the TREE part        begin %s\n",
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
         
    double wave_function_tree[2*XYZnodeSites];
    memset(wave_function_tree,0,sizeof(wave_function_tree));  

    run_NBSwf_pi_sigma_tree(wave_function_tree, time_WF);
    NBSwf_print(wave_function_tree    , hadron_names+"_tree"    , time_WF );

    if(MPI_rank==0){
      printf("        ++++++ run_NBSwf :                        end %s\n", 
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);


    // loop part
    if(MPI_rank==0){
      printf("        ++++++ run_NBSwf : the LOOP part        begin %s\n",
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);

    double wave_function_loop[2*XYZnodeSites];
    memset(wave_function_loop,0,sizeof(wave_function_loop));  

    run_NBSwf_pi_sigma_loop(wave_function_loop, time_WF);
    NBSwf_print(wave_function_loop     , hadron_names+"_loop"   , time_WF  );
 
    if(MPI_rank==0){
      printf("        ++++++ run_NBSwf :                        end %s\n", 
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);

    
    // combine tree and loop
    for(int i=0;i<2*XYZnodeSites;i++){
      wave_function[i]=wave_function_tree[i]+wave_function_loop[i];    
    }
    NBSwf_print(wave_function          , hadron_names , time_WF            );

    if(MPI_rank==0){
      printf(" ++++++ run_NBSwf : end %s\n", 
             LocalTime().c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
 
  }
  else{
  
    printf("ERROR - unknown hadron name");
    abort();
  }

  
};


// ================================================================================================
// calculate pi-sigma NBS wave function --- the tree part 
//       formulas in notes 
//
void class_two_hadrons::run_NBSwf_pi_sigma_tree(double* wave_function, int time_WF){

  // set wf to zero
  memset(wave_function,0,sizeof(wave_function));
  
  // complexify propagators 
  COMPLEX* Prop_ud            = (COMPLEX*)prop_ud        ;//+ prop_slv_idx(0,0,0,0,ixyz,it); 
  COMPLEX* Prop_s             = (COMPLEX*)prop_s         ;//+ prop_slv_idx(0,0,0,0,ixyz,it);

#define back_prop(prop,c,a,cp,ap,ixyz,it)                               \
        ( ZGM(a,5) * ZGM (IGM(ap,5),5) *                                 \
          Conj(prop[ prop_slv_idx(c,IGM(a,5),cp,IGM(ap,5) ,ixyz,it) ]) )

  int it_full=  (iT_src+time_WF + 100*Tsites)%Tsites;      
  int it = it_full % TnodeSites;     

  // noise summation
  for(int i_noise = 0; i_noise < N_noises; i_noise++){
    
    COMPLEX* Noise      = (COMPLEX*)sources->get_noise_ixyz(i_noise);

    // calculating back propagator contraction first
    COMPLEX back_prop_contraction[3*4*3*4];
    memset(back_prop_contraction,0,sizeof(back_prop_contraction));

//    if(it_full/TnodeSites==TnodeCoor){
      for(int a=0;a<3;a++){
      for(int alf=0;alf<4;alf++){
      for(int aP=0; aP<3;aP++){
      for(int alfP=0; alfP<4;alfP++){
      
        int index=prop_slv_cs_idx(a, alf, aP, alfP);

        for(int Y_ixyz = 0;  Y_ixyz < XYZnodeSites; Y_ixyz++){            
          back_prop_contraction[index] += Noise[Y_ixyz] * 
                                          back_prop(Prop_ud,       a, alf, aP, alfP, Y_ixyz, it);
        }
          
      }}}}
//    }

    // free Dirac index summation
    for(int ALPHA = 0; ALPHA < 2; ALPHA++){

        // source summation
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
        for(int d       = 0; d      < 3; d ++){
        for(int alpha   = 0; alpha  < 4; alpha++){
          int beta=IGM(alpha,5);
        for(int color   = 0; color  < 6; color++){
          int a    =Eps(0,color);
          int b    =Eps(1,color);
          int c    =Eps(2,color);          
        for(int gamma   = 0; gamma  < 4; gamma++){
          int delta = icg5[gamma]; 


          double WF_mom_space[2*XYZnodeSites];
          memset(WF_mom_space,0,sizeof(WF_mom_space));
   
          double WF_left[2*XYZnodeSites];
          memset(WF_left,0,sizeof(WF_mom_space));
          double WF_right[2*XYZnodeSites];
          memset(WF_right,0,sizeof(WF_mom_space));
   
          //-1/2 T1
//          if(it_full/TnodeSites==TnodeCoor){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){

              ((COMPLEX*)WF_left)[X_ixyz]  =  Conj(Noise[X_ixyz]) * 
                                              Prop_ud[ prop_slv_idx(d, beta ,  bP, gammaP,  X_ixyz,it) ] ;

              ((COMPLEX*)WF_right)[X_ixyz] =  Prop_ud[ prop_slv_idx(a, ALPHA,  aP, ALPHA ,  X_ixyz,it) ] *
                                              Prop_ud[ prop_slv_idx(b, gamma,  dP,  betaP,  X_ixyz,it) ] *
                                              Prop_s[  prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ] ;
            }
//          } 
            
          FFT3D (WF_left  , FFTW_BACKWARD );
          FFT3D (WF_right , FFTW_FORWARD  );
 
//          if(it_full/TnodeSites==TnodeCoor){
            for(int X_ixyz = 0;  X_ixyz < 2*XYZnodeSites; X_ixyz++){
              WF_mom_space[X_ixyz] += -0.5 * WF_left[X_ixyz] * WF_right[X_ixyz];
            } 
//          } 
                           

          //-1/2 T2
//         if( it_full/TnodeSites==TnodeCoor){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              
              ((COMPLEX*)WF_left)[X_ixyz]  =  Conj(Noise[X_ixyz]) * 
                                              Prop_ud[ prop_slv_idx(d, beta ,  aP, ALPHA ,  X_ixyz,it) ] ;

              ((COMPLEX*)WF_right)[X_ixyz] =  Prop_ud[ prop_slv_idx(a, ALPHA,  dP,  betaP,  X_ixyz,it) ] *
                                              Prop_ud[ prop_slv_idx(b, gamma,  bP, gammaP,  X_ixyz,it) ] *
                                              Prop_s[  prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ] ;
            } 
//          } 
            
          FFT3D (WF_left  , FFTW_BACKWARD );
          FFT3D (WF_right , FFTW_FORWARD  );
 
//          if( it_full/TnodeSites==TnodeCoor){
            for(int X_ixyz = 0;  X_ixyz < 2*XYZnodeSites; X_ixyz++){
              WF_mom_space[X_ixyz] += -0.5 * WF_left[X_ixyz] * WF_right[X_ixyz];
            } 
//          } 

          //+1 T3
//          if( it_full/TnodeSites==TnodeCoor){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              
            ((COMPLEX*)WF_left)[X_ixyz]  =  Conj(Noise[X_ixyz]) * 
                                            Prop_ud[ prop_slv_idx(d, beta ,  dP,  betaP,  X_ixyz,it) ] ;

            ((COMPLEX*)WF_right)[X_ixyz] =  Prop_ud[ prop_slv_idx(a, ALPHA,  bP, gammaP,  X_ixyz,it) ] *
                                            Prop_ud[ prop_slv_idx(b, gamma,  aP, ALPHA ,  X_ixyz,it) ] *
                                            Prop_s[  prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ] ;
            } 
//          } 
            
          FFT3D (WF_left  , FFTW_BACKWARD );
          FFT3D (WF_right , FFTW_FORWARD  );
 
//          if( it_full/TnodeSites==TnodeCoor){
            for(int X_ixyz = 0;  X_ixyz < 2*XYZnodeSites; X_ixyz++){
              WF_mom_space[X_ixyz] +=  WF_left[X_ixyz] * WF_right[X_ixyz];
            } 
//          } 

          //+1 T4
//          if( it_full/TnodeSites==TnodeCoor){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
             
            ((COMPLEX*)WF_left)[X_ixyz]  =  Conj(Noise[X_ixyz]) * 
                                            Prop_ud[ prop_slv_idx(d, beta ,  dP,  betaP,  X_ixyz,it) ] ;

            ((COMPLEX*)WF_right)[X_ixyz] =  Prop_ud[ prop_slv_idx(a, ALPHA,  aP, ALPHA ,  X_ixyz,it) ] *
                                            Prop_ud[ prop_slv_idx(b, gamma,  bP, gammaP,  X_ixyz,it) ] *
                                            Prop_s[  prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ] ;
            } 
//          } 
            
          FFT3D (WF_left  , FFTW_BACKWARD );
          FFT3D (WF_right , FFTW_FORWARD  );
 
//          if( it_full/TnodeSites==TnodeCoor){
            for(int X_ixyz = 0;  X_ixyz < 2*XYZnodeSites; X_ixyz++){
              WF_mom_space[X_ixyz] += - WF_left[X_ixyz] * WF_right[X_ixyz];
            } 
//          } 

          //-1/2 T5
//          if( it_full/TnodeSites==TnodeCoor){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              
              ((COMPLEX*)WF_left)[X_ixyz]  =  Conj(Noise[X_ixyz]) * 
                                              Prop_ud[ prop_slv_idx(d, beta ,  bP, gammaP,  X_ixyz,it) ] ;

              ((COMPLEX*)WF_right)[X_ixyz] =  Prop_ud[ prop_slv_idx(a, ALPHA,  dP,  betaP,  X_ixyz,it) ] *
                                              Prop_ud[ prop_slv_idx(b, gamma,  aP, ALPHA ,  X_ixyz,it) ] *
                                              Prop_s[  prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ] ;
            } 
//          } 
            
          FFT3D (WF_left  , FFTW_BACKWARD );
          FFT3D (WF_right , FFTW_FORWARD  );
 
//          if( it_full/TnodeSites==TnodeCoor){
            for(int X_ixyz = 0;  X_ixyz < 2*XYZnodeSites; X_ixyz++){
              WF_mom_space[X_ixyz] += 0.5 * WF_left[X_ixyz] * WF_right[X_ixyz];
            } 
//          } 

          //-1/2 T6
//          if( it_full/TnodeSites==TnodeCoor){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              
              ((COMPLEX*)WF_left)[X_ixyz]  =  Conj(Noise[X_ixyz]) * 
                                              Prop_ud[ prop_slv_idx(d, beta ,  aP, ALPHA ,  X_ixyz,it) ] ;

              ((COMPLEX*)WF_right)[X_ixyz] =  Prop_ud[ prop_slv_idx(a, ALPHA,  bP, gammaP,  X_ixyz,it) ] *
                                              Prop_ud[ prop_slv_idx(b, gamma,  dP, betaP,  X_ixyz,it) ] *
                                              Prop_s[  prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ] ;
            } 
//          } 
            
          FFT3D (WF_left  , FFTW_BACKWARD );
          FFT3D (WF_right , FFTW_FORWARD  );
 
//          if( it_full/TnodeSites==TnodeCoor){
            for(int X_ixyz = 0;  X_ixyz < 2*XYZnodeSites; X_ixyz++){
              WF_mom_space[X_ixyz] += 0.5 * WF_left[X_ixyz] * WF_right[X_ixyz];
            } 
//          } 
                            

          FFT3D(WF_mom_space, FFTW_BACKWARD);

          COMPLEX factor = back_prop_contraction[prop_slv_cs_idx(d, alpha, dP, alphaP)] *
                           ZGM(alpha,5)  * Eps(3,color)  * zcg5[gamma] * 
                           ZGM(alphaP,5) * Eps(3,colorP) * zcg5[gammaP];
 
          for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
            ((COMPLEX*)wave_function)[X_ixyz] += factor *  ((COMPLEX*)WF_mom_space[X_ixyz]);
          } 
            
        }}}}//sink
                                             
      } }}}//source

    }//ALPHA
    
  }//i_noise

  double norm_factor = 0.5 / (N_noises* XYZsites);
 
  for(int X_ixyz = 0;  X_ixyz < 2*XYZnodeSites; X_ixyz++){
    wave_function[X_ixyz] += norm_factor * wave_function[X_ixyz];
  } 
    
  MPI_Barrier(MPI_COMM_WORLD);
  
#undef back_prop  
}


// ================================================================================================
// calculate pi-sigma NBS wave function --- the loop part 
//       formulas in notes 
//

void class_two_hadrons::run_NBSwf_pi_sigma_loop(double* wave_function, int time_WF){

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

    
  int it_full=  (iT_src+time_WF + 100*Tsites)%Tsites;      
  int it = it_full % TnodeSites;     


  for(int i_noise = 0; i_noise < N_noises; i_noise++){
    
    COMPLEX* Noise      = (COMPLEX*)sources->get_noise_ixyz(i_noise);
    COMPLEX* Prop_noise = Prop_noise_full + i_noise * XYZTnodeSites * 3*4*3*4;

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
        for(int d       = 0; d      < 3; d ++){
        for(int alpha   = 0; alpha  < 4; alpha++){
          int beta=IGM(alpha,5);
        for(int color   = 0; color  < 6; color++){
          int a    =Eps(0,color);
          int b    =Eps(1,color);
          int c    =Eps(2,color);          
        for(int gamma   = 0; gamma  < 4; gamma++){
          int delta = icg5[gamma]; 



          double WF_mom_space[2*XYZnodeSites];
          memset(WF_mom_space,0,sizeof(WF_mom_space));
   
          double WF_left[2*XYZnodeSites];
          memset(WF_left,0,sizeof(WF_mom_space));
          double WF_right[2*XYZnodeSites];
          memset(WF_right,0,sizeof(WF_mom_space));
   
          //-3/2 L7
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              
            ((COMPLEX*)WF_left)[X_ixyz]  =  Conj(Noise[X_ixyz]) * 
                                            Prop_ud   [  prop_slv_idx(d, beta ,  bP, gammaP,  X_ixyz,it) ] ;

            ((COMPLEX*)WF_right)[X_ixyz] =  Prop_ud   [  prop_slv_idx(a, ALPHA,  dP, betaP ,  X_ixyz,it) ] *
                                            Prop_noise[  prop_slv_idx(b, gamma,  d , alpha ,  X_ixyz,it) ] *
                                            Prop_s    [  prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ] ;
            } 
//          }            

          FFT3D (WF_left  , FFTW_BACKWARD );
          FFT3D (WF_right , FFTW_FORWARD  );
 
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              ((COMPLEX*)WF_mom_space)[X_ixyz] += -1.5 *  src_ctr_global[ prop_slv_cs_idx(dP, alphaP, aP, ALPHA) ] *
                                                  ((COMPLEX*)WF_left)[X_ixyz] * ((COMPLEX*)WF_right)[X_ixyz];
            } 
//          }                          

          //+3/2 L8
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              
            ((COMPLEX*)WF_left)[X_ixyz]  =  Conj(Noise[X_ixyz]) * 
                                            Prop_ud   [  prop_slv_idx(d, beta ,  dP, betaP ,  X_ixyz,it) ] ;

            ((COMPLEX*)WF_right)[X_ixyz] =  Prop_ud   [  prop_slv_idx(a, ALPHA,  bP, gammaP ,  X_ixyz,it) ] *
                                            Prop_noise[  prop_slv_idx(b, gamma,  d , alpha ,  X_ixyz,it) ] *
                                            Prop_s    [  prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ] ;
            } 
//          }            
         
          FFT3D (WF_left  , FFTW_BACKWARD );
          FFT3D (WF_right , FFTW_FORWARD  );
 
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              ((COMPLEX*)WF_mom_space)[X_ixyz] += -1.5 *  src_ctr_global[ prop_slv_cs_idx(dP, alphaP, aP, ALPHA) ] *
                                                  ((COMPLEX*)WF_left)[X_ixyz] * ((COMPLEX*)WF_right)[X_ixyz];
            } 
//          }

          //-3/2 L9
//          if( it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              
              ((COMPLEX*)WF_left)[X_ixyz]  =  Conj(Noise[X_ixyz]) * 
                                              Prop_ud   [  prop_slv_idx(d, beta ,  bP, gammaP,  X_ixyz,it) ] ;

              ((COMPLEX*)WF_right)[X_ixyz] =  Prop_ud   [  prop_slv_idx(b, gamma,  dP, betaP ,  X_ixyz,it) ] *
                                              Prop_noise[  prop_slv_idx(a, ALPHA,  d , alpha ,  X_ixyz,it) ] *
                                              Prop_s    [  prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ] ;
            } 
//          }           
  
          FFT3D (WF_left  , FFTW_BACKWARD );
          FFT3D (WF_right , FFTW_FORWARD  );
 
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              ((COMPLEX*)WF_mom_space)[X_ixyz] += +1.5 *  src_ctr_global[ prop_slv_cs_idx(dP, alphaP, aP, ALPHA) ] *
                                                  ((COMPLEX*)WF_left)[X_ixyz] * ((COMPLEX*)WF_right)[X_ixyz];
            } 
//          }

          //+3/2 L10
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              
              ((COMPLEX*)WF_left)[X_ixyz]  =  Conj(Noise[X_ixyz]) * 
                                              Prop_ud   [  prop_slv_idx(d, beta ,  dP, betaP ,  X_ixyz,it) ] ;

              ((COMPLEX*)WF_right)[X_ixyz] =  Prop_ud   [  prop_slv_idx(b, gamma,  bP, gammaP,  X_ixyz,it) ] *
                                              Prop_noise[  prop_slv_idx(a, ALPHA,  d , alpha ,  X_ixyz,it) ] *
                                              Prop_s    [  prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ] ;
            } 
//          }           
   
          FFT3D (WF_left  , FFTW_BACKWARD );
          FFT3D (WF_right , FFTW_FORWARD  );
   
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              ((COMPLEX*)WF_mom_space)[X_ixyz] += +1.5 *  src_ctr_global[ prop_slv_cs_idx(dP, alphaP, aP, ALPHA) ] *
                                                  ((COMPLEX*)WF_left)[X_ixyz] * ((COMPLEX*)WF_right)[X_ixyz];
            } 
//          }
          
          //-3/2 L13
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              
              ((COMPLEX*)WF_left)[X_ixyz]  =  Conj(Noise[X_ixyz]) * 
                                              Prop_ud   [  prop_slv_idx(d, beta ,  aP, ALPHA,  X_ixyz,it) ] ;

              ((COMPLEX*)WF_right)[X_ixyz] =  Prop_ud   [  prop_slv_idx(a, ALPHA,  dP, betaP ,  X_ixyz,it) ] *
                                              Prop_noise[  prop_slv_idx(b, gamma,  d , alpha ,  X_ixyz,it) ] *
                                             Prop_s    [  prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ] ;
            } 
//          }
                     
          FFT3D (WF_left  , FFTW_BACKWARD );
          FFT3D (WF_right , FFTW_FORWARD  );
 
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              ((COMPLEX*)WF_mom_space)[X_ixyz] += +1.5 *  src_ctr_global[ prop_slv_cs_idx(dP, alphaP, bP, gammaP) ] *
                                                  ((COMPLEX*)WF_left)[X_ixyz] * ((COMPLEX*)WF_right)[X_ixyz];
            } 
//          } 
          
          //+3/2 L14
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              
              ((COMPLEX*)WF_left)[X_ixyz]  =  Conj(Noise[X_ixyz]) * 
                                              Prop_ud   [  prop_slv_idx(d, beta ,  dP, betaP ,  X_ixyz,it) ] ;

              ((COMPLEX*)WF_right)[X_ixyz] =  Prop_ud   [  prop_slv_idx(a, ALPHA,  aP, ALPHA ,  X_ixyz,it) ] *
                                              Prop_noise[  prop_slv_idx(b, gamma,  d , alpha ,  X_ixyz,it) ] *
                                              Prop_s    [  prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ] ;
            } 
//          }
                     
          FFT3D (WF_left  , FFTW_BACKWARD );
          FFT3D (WF_right , FFTW_FORWARD  );
        
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              ((COMPLEX*)WF_mom_space)[X_ixyz] += +1.5 *  src_ctr_global[ prop_slv_cs_idx(dP, alphaP, bP, gammaP) ] *
                                               ((COMPLEX*)WF_left)[X_ixyz] * ((COMPLEX*)WF_right)[X_ixyz];
            } 
//          }
          
          //-3/2 L15
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              
              ((COMPLEX*)WF_left)[X_ixyz]  =  Conj(Noise[X_ixyz]) * 
                                              Prop_ud   [  prop_slv_idx(d, beta ,  aP, ALPHA ,  X_ixyz,it) ] ;

              ((COMPLEX*)WF_right)[X_ixyz] =  Prop_ud   [  prop_slv_idx(b, gamma,  bP, betaP ,  X_ixyz,it) ] *
                                              Prop_noise[  prop_slv_idx(a, ALPHA,  d , alpha ,  X_ixyz,it) ] *
                                              Prop_s    [  prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ] ;
            } 
//          }           

          FFT3D (WF_left  , FFTW_BACKWARD );
          FFT3D (WF_right , FFTW_FORWARD  );
 
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              ((COMPLEX*)WF_mom_space)[X_ixyz] += -1.5 *  src_ctr_global[ prop_slv_cs_idx(dP, alphaP, bP, gammaP) ] *
                                                  ((COMPLEX*)WF_left)[X_ixyz] * ((COMPLEX*)WF_right)[X_ixyz];
            } 
//          }
          
          //+3/2 L16
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              
              ((COMPLEX*)WF_left)[X_ixyz]  =  Conj(Noise[X_ixyz]) * 
                                              Prop_ud   [  prop_slv_idx(d, beta ,  bP, betaP ,  X_ixyz,it) ] ;

              ((COMPLEX*)WF_right)[X_ixyz] =  Prop_ud   [  prop_slv_idx(b, gamma,  aP, ALPHA ,  X_ixyz,it) ] *
                                              Prop_noise[  prop_slv_idx(a, ALPHA,  d , alpha ,  X_ixyz,it) ] *
                                              Prop_s    [  prop_slv_idx(c, delta,  cP, deltaP,  X_ixyz,it) ] ;
            } 
//          }           
       
          FFT3D (WF_left  , FFTW_BACKWARD );
          FFT3D (WF_right , FFTW_FORWARD  );
  
//          if(  it_full/TnodeSites == TnodeCoor ){
            for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
              ((COMPLEX*)WF_mom_space)[X_ixyz] += -1.5 *  src_ctr_global[ prop_slv_cs_idx(dP, alphaP, bP, gammaP) ] *
                                                  ((COMPLEX*)WF_left)[X_ixyz] * ((COMPLEX*)WF_right)[X_ixyz];
            } 
//          }

          FFT3D(WF_mom_space, FFTW_BACKWARD);

          COMPLEX factor = ZGM(alpha,5)  * Eps(3,color)  * zcg5[gamma] * 
                           ZGM(alphaP,5) * Eps(3,colorP) * zcg5[gammaP];
   
          for(int X_ixyz = 0;  X_ixyz < XYZnodeSites; X_ixyz++){
            ((COMPLEX*)wave_function)[X_ixyz] += factor *  ((COMPLEX*)WF_mom_space)[X_ixyz];
          } 
            

        }}}}//sink
                     
                        
      } }}}//source

    }//ALPHA
    
  }//i_noise


  double norm_factor = 0.5 / (N_noises* XYZsites);
 
  for(int X_ixyz = 0;  X_ixyz < 2*XYZnodeSites; X_ixyz++){
    wave_function[X_ixyz] += norm_factor * wave_function[X_ixyz];
  } 
    

  // reduce from all MPI processes
  MPI_Barrier(MPI_COMM_WORLD);

}



// ================================================================================================
// print the NBS wave function
//
void class_two_hadrons::NBSwf_print(double *wave_function, string hadron_names, int time_WF)
{
  
  int it_full=  (iT_src + time_WF + 100*Tsites)%Tsites;      
  int it = it_full % TnodeSites;     
  
  int time_WF_full = (it + TnodeSites*TnodeCoor - iT_src + 100*Tsites)%Tsites;

//  if( it_full/TnodeSites==TnodeCoor){


    // create directory and file name
    char dir[256];
    string dir_base="results/"+base+"/"+prefix+"NBSwaveF_two_hardons/"+hadron_names;
    snprintf(dir,sizeof(dir), "%s.tS%02dtWF%02d/",
             dir_base.c_str(),iT_src,
             time_WF_full);

    create_directory(dir);
    //printf("%s %s\n",prefix.c_str(),dir.c_str());
  
    char wfile[256];
    snprintf(wfile,sizeof(wfile), "%sx%02dy%02dz%02d",
             dir,XnodeCoor,YnodeCoor,ZnodeCoor);

    printf("        ++++++ NBSff_print : print %10s NBSwf for t_src=%2d and t_WF=%2d to file     %s\n                     %s\n",
           hadron_names.c_str(), iT_src, time_WF_full, LocalTime().c_str(), wfile);
 
    // open output file
    string ofname(wfile);
    std::ofstream fout(ofname.c_str());
    fout.setf(std::ios::scientific);
    fout.precision(16);


    // print NBS wave function
    for(int index = 0; index < XYZnodeSites; index++){

      int xcoor =  index % XnodeSites                 +  XnodeSites*XnodeCoor;
      int ycoor =  (index / XnodeSites) % YnodeSites  +  YnodeSites*YnodeCoor;
      int zcoor =  index / (XnodeSites * YnodeSites)  +  ZnodeSites*ZnodeCoor;

      char line[1000];
      snprintf(line, sizeof(line), "%7d,%3d,%3d,%3d,  %1.16e %1.16e\n",
               xcoor*xcoor+ycoor*ycoor+zcoor*zcoor,
               xcoor, ycoor, zcoor, 
               wave_function[2* index], wave_function[2* index+1]);
      fout << line;
    }



    fout.close();

//  }

  MPI_Barrier(MPI_COMM_WORLD);

}




