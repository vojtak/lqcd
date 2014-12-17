#ifndef IS_INCLUDED_CLASS_NOISE_SOURCE_H
#define IS_INCLUDED_CLASS_NOISE_SOURCE_H

/* local from EXTERNAL */
#include "randomc.h"


class class_noise_source : public class_global_wrapper {

  public:
  
  // =======================================
  // constructor
  //
  class_noise_source(int N_noises_in){
  
    N_noises=N_noises_in;

    // initialize random number generator
    int seed[2] = {(int)time(0),MPI_rank};
    RanGen.RandomInitByArray(seed,2);

    // allocate memory for the noise_vectors
    noise_vector = new double[2* N_noises * XYZnodeSites];
    memset(noise_vector,0,sizeof(noise_vector));

  };


  // =======================================
  // functions
  //
  
  void run();
  void print();
  void print_full();
  
  double *get_noise_ixyz(int noise_num);

  // =======================================
  // data members
  //
  int N_noises;
 
  double *noise_vector;
 
  CRandomMersenne RanGen;

  private:

  // =======================================
  // data members
  //


  // =======================================
  // key functions
  //

  void generate_noise_ixyz(double *noise);

  void calculate_propagator();

  // =======================================
  // assistant functions
  //
  
 
};

#endif
