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

    prefix="";
    
    Nc   = CommonParameters::Nc();
    Nd   = CommonParameters::Nd();
    Ndim = CommonParameters::Ndim();
    Nvol = CommonParameters::Nvol();

    XnodeSites = CommonParameters::Nx();
    YnodeSites = CommonParameters::Ny();
    ZnodeSites = CommonParameters::Nz();
    TnodeSites = CommonParameters::Nt();

    Xsites     = CommonParameters::Lx();
    Ysites     = CommonParameters::Ly();
    Zsites     = CommonParameters::Lz();
    Tsites     = CommonParameters::Lt();

    XnodeCoor  = Communicator::ipe(0);
    YnodeCoor  = Communicator::ipe(1);
    ZnodeCoor  = Communicator::ipe(2);
    TnodeCoor  = Communicator::ipe(3);

    Xnodes     = Communicator::npe(0);
    Ynodes     = Communicator::npe(1);
    Znodes     = Communicator::npe(2);
    Tnodes     = Communicator::npe(3);

    XYZnodeSites = XnodeSites * YnodeSites * ZnodeSites;
    XYZsites     = Xsites     * Ysites     * Zsites;
  
    XYZTnodeSites = XYZnodeSites * TnodeSites;
    XYZTsites     = XYZsites     * Tsites;

    XYZnodeCoor = XnodeCoor + Xnodes * (YnodeCoor + Ynodes * ZnodeCoor);

  }

  // =======================================
  // data members
  //

  // pi
  #define PI  3.1415926535897932385

  // MPI parameters
  int MPI_flag, MPI_size, MPI_rank;       


  // name of the gauge configuration and prefix for measurements
  string base; 	        
  void set_base_name(string base_in){
    base = base_in;
  };
  string prefix; 	        
  void set_prefix_name(string prefix_in){
    prefix = prefix_in;
  };
  

  // source position in time
  int iT_src;		
  void set_source_position(int iT_src_in){
    iT_src = iT_src_in;
  };


  // =======================================
  // global parameters
  //
  int Nc;
  int Nd;
  int Ndim;
  int Nvol;

  int XnodeSites;
  int YnodeSites;
  int ZnodeSites;
  int TnodeSites;

  int Xsites;
  int Ysites;
  int Zsites;
  int Tsites;

  int XnodeCoor;
  int YnodeCoor;
  int ZnodeCoor;
  int TnodeCoor;

  int Xnodes;
  int Ynodes;
  int Znodes;
  int Tnodes;

  int XYZnodeSites;
  int XYZsites;
  
  int XYZTnodeSites;
  int XYZTsites;

  int XYZnodeCoor;


  // =======================================
  // other useful functions
  //

 
  void create_directory(string path);
  string LocalTime();

};

#endif
