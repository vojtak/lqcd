#ifndef IS_INCLUDED_CLASS_NOISE_SOURCE_H
#define IS_INCLUDED_CLASS_NOISE_SOURCE_H

/* local from EXTERNAL */
#include "randomc.h"


class class_sources : public class_global_wrapper {

  public:
  
  // =======================================
  // constructor
  //
  class_sources(){
  
    // initialize random number generator
    int seed[2] = {(int)time(0),MPI_rank};
    RanGen.RandomInitByArray(seed,2);

  };


  // =======================================
  // functions
  //
  
  void run();
   
  void generate_wall_source();
  void generate_noise_source_vector(int N_noises_in);

  void print();       //for debugging
  void print_full();  //for debugging
  

  
  double *get_noise_ixyz(int noise_num);

  double *get_wall_ixyz();

  // =======================================
  // data members
  //

  int N_noises;
 
  double *noise_sources_vector;

  double *wall_source;
 
  CRandomMersenne RanGen;

  private:

  // =======================================
  // data members
  //


  // =======================================
  // key functions
  //

  void generate_noise_ixyz(double *noise);
  void generate_wall_ixyz(double *wall);


  // =======================================
  // assistant functions
  //
  
 
};

#endif
