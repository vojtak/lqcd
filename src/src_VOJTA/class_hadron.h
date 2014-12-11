#ifndef IS_INCLUDED_CLASS_HADRON_H
#define IS_INCLUDED_CLASS_HADRON_H


class class_hadron {

  public:
  
  // =======================================
  // constructor
  class_hadron(double *prop_ud_in, double *prop_s_in){
    prop_ud=prop_ud_in;
    prop_s=prop_s_in;
  };

  // =======================================
  // data members

  double *prop_ud; 	// ud propagator
  double *prop_s;  	//  s propagator
  
  string base; 	        // name of the gauge configuration
  
  int iT_src;		// source position in time


  // =======================================
  // functions
  //

  void run_GF(string hadron_name);

  void run_GF_meson(double* correlator, double* prop_up, double* prop_down);

  void corr_print(double *correlator, string hadron_name);

  // =======================================
  // common parameters and assistant functions (move to some parent class eventually)
  //
  
  void set_base_name(string base_in){
    base = base_in;
  };
  void set_source_position(int iT_src_in){
    iT_src = iT_src_in;
  };
  
  private:
  
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
   
  #define Dirac(     alpha, x)  (alpha +  4*(x))
  #define Color(     c,     x)  (c     +  3*(x))

  // format of indexing is (sink color, sink dirac, source color, source dirac, sink xyz, sink t) 
  size_t prop_slv_idx(int c, int d, int cP, int dP, int ixyz, int it)
  {
    return Color(c, Dirac(d, 
                          Color(cP, Dirac(dP, 
                                          ixyz + XYZnodeSites * (it)))));
  }
  int prop_slv_cs_idx(int c, int d, int cP, int dP)
  {
    // the following two implementations are equivalent, but K-computer prefers the latter
    //return prop_slv_idx(c, d, cP, dP, 0, 0);
      return Color(c, Dirac(d, 
                            Color(cP, Dirac(dP, 0))));
  }

};

#endif
