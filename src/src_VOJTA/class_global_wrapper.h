#ifndef IS_INCLUDED_CLASS_GLOBAL_WRAPPER_H
#define IS_INCLUDED_CLASS_GLOBAL_WRAPPER_H

/* local from Bridge++ */
#include "communicator.h"
#include "commonParameters.h"

class class_global_wrapper{

  public:
  
  // =======================================
  // constructor
  //
  class_global_wrapper(){

    MPI_Initialized(&MPI_flag);
    MPI_Comm_size(MPI_COMM_WORLD,&MPI_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&MPI_rank);

  }

  // =======================================
  // data members
  //

  // MPI parameters
  int MPI_flag, MPI_size, MPI_rank;       


  // name of the gauge configuration
  string base; 	        
  void set_base_name(string base_in){
    base = base_in;
  };
  

  // source position in time
  int iT_src;		
  void set_source_position(int iT_src_in){
    iT_src = iT_src_in;
  };


  // =======================================
  // global parameters
  //
  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  int XnodeSites = CommonParameters::Nx();
  int YnodeSites = CommonParameters::Ny();
  int ZnodeSites = CommonParameters::Nz();
  int TnodeSites = CommonParameters::Nt();

  int Xsites     = CommonParameters::Lx();
  int Ysites     = CommonParameters::Ly();
  int Zsites     = CommonParameters::Lz();
  int Tsites     = CommonParameters::Lt();

  int XnodeCoor  = Communicator::ipe(0);
  int YnodeCoor  = Communicator::ipe(1);
  int ZnodeCoor  = Communicator::ipe(2);
  int TnodeCoor  = Communicator::ipe(3);

  int Xnodes     = Communicator::npe(0);
  int Ynodes     = Communicator::npe(1);
  int Znodes     = Communicator::npe(2);
  int Tnodes     = Communicator::npe(3);

  int XYZnodeSites = XnodeSites * YnodeSites * ZnodeSites;
  int XYZsites     = Xsites     * Ysites     * Zsites;
  
  int XYZTnodeSites = XYZnodeSites * TnodeSites;
  int XYZTsites     = XYZsites     * Tsites;

  int XYZnodeCoor = XnodeCoor + Xnodes * (YnodeCoor + Ynodes * ZnodeCoor);


  // =======================================
  // other useful functions
  //
 
  void create_directory(string path);
  string LocalTime();

};

#endif
