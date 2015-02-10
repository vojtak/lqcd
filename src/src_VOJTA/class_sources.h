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
   
  void generate_point_source(int X_coor, int Y_coor, int Z_coor);
  void generate_wall_source();
  void generate_noise_source_vector(int N_noises_in);


  void print();       //for debugging
  void print_full();  //for debugging
  

  double *get_point_ixyz();
  double *get_wall_ixyz();
  double *get_noise_ixyz(int noise_num);

  // =======================================
  // data members
  //

  int N_noises;
 

  double *point_source;

  double *wall_source;

  double *noise_sources_vector;
 
  CRandomMersenne RanGen;

  private:

  // =======================================
  // data members
  //


  // =======================================
  // key functions
  //

  void generate_point_ixyz(double *point, int X_coor, int Y_coor, int Z_coor);
  void generate_wall_ixyz(double *wall);
  void generate_noise_ixyz(double *noise);

int rank_from_coor1(int x,int y, int z, int t){
    int rank;
    int coor[4]={x,y,z,t};
    Communicator::grid_rank(&rank,coor);
    return rank;
}


  // =======================================
  // assistant functions
  //
  
 
};

#endif
