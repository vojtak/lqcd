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

    // allocate memory for the sources
    wall_source = new double[2 * XYZnodeSites];
    memset(wall_source,0,sizeof(wall_source));

    noise_vector = new double[2* N_noises * XYZnodeSites];
    memset(noise_vector,0,sizeof(noise_vector));

  };


  // =======================================
  // functions
  //
  
  void run();
  void print();
  void print_full();  //for debugging
  
  
  double *get_noise_ixyz(int noise_num);

  double *get_wall_ixyz();

  // =======================================
  // data members
  //

  int N_noises;
 
  double *noise_vector;

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
