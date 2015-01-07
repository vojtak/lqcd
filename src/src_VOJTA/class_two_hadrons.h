#ifndef IS_INCLUDED_CLASS_TWO_HADRONS_H
#define IS_INCLUDED_CLASS_TWO_HADRONS_H


class class_two_hadrons : public class_global_wrapper {

  public:
  
  // =======================================
  // constructor
  //
  
  class_two_hadrons(double *prop_ud_in, double *prop_s_in){
    prop_ud=prop_ud_in;
    prop_s=prop_s_in;
  };


  // =======================================
  // functions
  //
  
  void run_all_GF();

  void run_GF(string hadrons_names);

  void run_GF_meson(double* correlator, double* prop_quark, double* prop_antiquark);


  private:

  // =======================================
  // data members
  //
  
  double *prop_ud; 	// ud propagator
  double *prop_s;  	//  s propagator


  // =======================================
  // key functions
  //

  void corr_print(double *correlator, string hadrons_names);  


  // =======================================
  // assistant functions
  //
  
  #define Dirac(     alpha, x)  (alpha +  4*(x))
  #define Color(     c,     x)  (c     +  3*(x))

  // format of indexing is (sink color, sink dirac, source color, source dirac, sink xyz, sink t) 
  size_t prop_slv_idx(int c, int d, int cP, int dP, int ixyz, int it)
  {
    return Color(c, Dirac(d, 
                          Color(cP, Dirac(dP, 
                                          ixyz + XYZnodeSites * (it)))));
  }

  int prop_slv_cs_idx(int c, int d, int cP, int dP) // not needed
  {
    // the following two implementations are equivalent, but K-computer prefers the latter
    //return prop_slv_idx(c, d, cP, dP, 0, 0);
      return Color(c, Dirac(d, 
                            Color(cP, Dirac(dP, 0))));
  }


};

#endif
