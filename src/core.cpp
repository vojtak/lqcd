/**
 *        @file    $Id:: main.cpp #$
 *
 *        @brief
 *
 *        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
 *                 $LastChangedBy: sueda $
 *
 *        @date    $LastChangedDate:: 2013-07-12 16:56:41 #$
 *
 *        @version $LastChangedRevision: 930 $
 **/

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef KEKBGQ
#include "bgnet.h"
#endif

//---------------------------------------------------------------------
//---------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <libgen.h>
#include <string.h>

#include <string>
#include <vector>

//#ifdef KEKBGQ
/* Define to 1 if you have the <sys/time.h> header file. */
#//define HAVE_SYS_TIME_H 1
/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
//#define TIME_WITH_SYS_TIME 1
//#endif

//#ifdef TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
//#else
//# ifdef HAVE_SYS_TIME_H
//#  include <sys/time.h>
//# else
//#  include <time.h>
//# endif
//#endif

#include "parameterManager_YAML.h"
#include "parameters_factory.h"

#include "bridgeIO.h"
using Bridge::vout;

#include "gaugeConfig.h"
#include "staples.h"
#include "polyakovLoop.h"
#include "wilsonLoop.h"

#include "gaugeFixing.h"
#include "gaugeFixing_Coulomb.h"
#include "gaugeFixing_Landau.h"

#include "randomNumbers_Mseries.h"

#include "fopr_Clover.h"
#include "fopr_Clover_eo.h"
#include "fprop_Standard_lex.h"
#include "fprop_Standard_eo.h"

#include "director_Smear.h"
#include "fopr_Smeared.h"

#include "gammaMatrixSet.h"
#include "gaugeFixing.h"
#include "projection.h"
#include "smear.h"
#include "solver.h"
#include "source.h"

#include "corr2pt_4spinor.h"


#include "source_Local.h"
#include "source_Wall.h"
//#include "source_Noise.h"


/* -------------- VOJTA code */
#include "X_misc_V.h"

/* local from HAL */
#include "COMPLEX.h"
#include "HAL_indexes.h"

/* local from VOJTA */

#include "class_global_wrapper.h"
#include "class_sources.h"
#include "class_hadron.h"
#include "class_two_hadrons.h"

/* local from EXTERNAL */
#include "fft3d.h"
#include "fftw3.h"



static void propagators_solve(string label,
                              Fprop* fprop, 
                              double * prop, 
                              int it_src, double * source);

static void iT_slice_selection(double *prop_final,
                               double *prop_tmp,
                               int noise_i,
                               int iT_pos);


//---------------------------------------------------------------------
//---------------------------------------------------------------------

static int Nc;
static int Nd;
static int Ndim;
static int Nvol;

static int XnodeSites, YnodeSites, ZnodeSites, TnodeSites;
static int Xnodes,     Ynodes,     Znodes,     Tnodes;
static int Xsites,     Ysites,     Zsites,     Tsites;
static int XnodeCoor,  YnodeCoor,  ZnodeCoor,  TnodeCoor;

static int XYZTnodeSites, XYZTsites;
static int XYZnodeSites,  XYZsites;

static int XYZnodeCoor;

//---------------------------------------------------------------------
//---------------------------------------------------------------------
// T. Doi

static const char* LocalTime();


//---------------------------------------------------------------------
//---------------------------------------------------------------------


static void converter(std::valarray<Field_F>& sq, double prop[]);

static void usage(const char* msg=NULL,...)
{
  Bridge::vout.general(
                       "usage: main [OPTIONS] (gconf1) (gconf2) ...\n"
                       "OPTIONS:\n"
                       "     -quark_ud (fname[quark_ud.yaml])\n"
                       "     -quark_s  (fname[quark_s.yaml])\n"
                    // "     -mapf     (fname[mapf])\n"
                       "     -f        (fname[gfile_list])\n"
                       "     -h/--help\n"
                       );
  if(NULL != msg){
    char str[4096];
    va_list args;
    va_start(args, msg);
    vsprintf(str, msg, args);
    va_end(args);

    Bridge::vout.general("\n%s\n", str);
  }
}

static const char* fname_quark_ud = NULL;
static const char* fname_quark_s  = NULL;
//static const char* fname_mapf     = NULL;
static const char* fname_gfile_list = NULL;

static int read_options(int argc,char* argv[])
{
  if(strcmp(argv[0],"-quark_ud")==0){
    fname_quark_ud = argv[1];
    return 2;
  }
  if(strcmp(argv[0],"-quark_s")==0){
    fname_quark_s = argv[1];
    return 2;
  }
  /*
  if(strcmp(argv[0],"-mapf")==0){
    fname_mapf = argv[1];
    return 2;
  }
  */
  if(strcmp(argv[0],"-f")==0){
    fname_gfile_list = argv[1];
    return 2;
  }
  if(strcmp(argv[0],"-h")==0 || strcmp(argv[0],"--help")==0){
    usage();
    exit(1);
    return 1;
  }

  return 0;
}

//====================================================================
//====================================================================
//------------ CORE
//====================================================================

int core(int argc,char** argv)
{
  argc--; argv++;
  
  // #### read parameters ####
  if(argc == 0){
    usage("ERROR: provide arguments");
    Communicator::sync();
    abort();
  }

  while(argc > 0 && argv[0][0]=='-'){
    int shift;
    if((shift = read_options(argc,argv)) != 0){
      Bridge::vout.general("@@@OPTIONS: %s\n", argv[0]);
      argc -= shift; argv += shift;
      continue;
    }

    usage("ERROR: invalid option '%s'\n", argv[0]);
    exit(1);
  }

  // ####  parameter setup  ####
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

  Bridge::VerboseLevel  vl = vout.set_verbose_level("General");

  vout.general("\n\t@@@ JOB(START):\t%s\n\t@@@ compilation date: %s %s\n\n", LocalTime(), __DATE__, __TIME__);

  //---------------------------------------------------------------------
  // common part ends here
  //---------------------------------------------------------------------

  
  //---------------------------------------------------------------------
  // load Gauge fixing parameters
  //---------------------------------------------------------------------
  Parameters_GaugeFixing_Coulomb params_gauge_fixing;
  {
    Parameters params_all;
    params_all.Register_Parameters("GaugeFixing", &params_gauge_fixing);

    ParameterManager_YAML params_manager;
    params_manager.read_params("gauge.fixing.yaml", &params_all);

    //delete params_all;
  }

  //---------------------------------------------------------------------
  // load parameters ud-quark propagator
  //---------------------------------------------------------------------
  Parameters *params_solver_ud = ParametersFactory::New("Solver");
  Parameters *params_clover_ud = ParametersFactory::New("Fopr.Clover");
  Parameters *params_source_ud = ParametersFactory::New("Source");
  {
    ParameterManager_YAML params_manager;
  
    Parameters *params_all    = new Parameters;

    params_all->Register_Parameters("Solver",      params_solver_ud);
    params_all->Register_Parameters("Fopr_Clover", params_clover_ud);
    params_all->Register_Parameters("Source",      params_source_ud);

    params_manager.read_params(fname_quark_ud, params_all);
    delete params_all;
  }

  //---------------------------------------------------------------------
  // load parameters s-quark propagators
  //---------------------------------------------------------------------
  Parameters *params_solver_s = ParametersFactory::New("Solver");
  Parameters *params_clover_s = ParametersFactory::New("Fopr.Clover");
  Parameters *params_source_s = ParametersFactory::New("Source");
  {
    ParameterManager_YAML params_manager;
  
    Parameters *params_all    = new Parameters;

    params_all->Register_Parameters("Solver",      params_solver_s);
    params_all->Register_Parameters("Fopr_Clover", params_clover_s);
    params_all->Register_Parameters("Source",      params_source_s);

    params_manager.read_params(fname_quark_s, params_all);
    delete params_all;
  }


  //---------------------------------------------------------------------
  // Extract parameters
  //---------------------------------------------------------------------

  // -- gauge field  
  Field_G        *U         = NULL;
  Field_G        *U_fixed   = NULL;

  U          = new Field_G(Nvol, Ndim);
  U_fixed    = new Field_G(Nvol, Ndim);
  
  // -- ud quark fermion field  
  GammaMatrixSet    *gmset_ud     = NULL;
  Fopr_Clover_eo    *fopr_c_ud    = NULL;
  Solver            *solver_ud    = NULL;
  Fprop             *fprop_ud     = NULL;
  Source            *source_ud    = NULL;

  {
    string str_solver_type = params_solver_ud -> get_string("solver_type");
    string str_gmset_type  = params_clover_ud -> get_string("gamma_matrix_type");
    string str_source_type = params_source_ud -> get_string("source_type");

    Bridge::vout.general("solver_type(ud)      : %s\n", str_solver_type.c_str());
    Bridge::vout.general("gamma_matrix_type(ud): %s\n", str_gmset_type.c_str());
    Bridge::vout.general("source_type(ud)      : %s\n", str_source_type.c_str());

    gmset_ud      = GammaMatrixSet::New(str_gmset_type);
    
    fopr_c_ud     = new Fopr_Clover_eo(str_gmset_type);
    solver_ud     = Solver::New(str_solver_type, fopr_c_ud);
    fprop_ud      = new Fprop_Standard_eo(solver_ud);
    source_ud     = Source::New(str_source_type);
  }
  
  // -- s quark fermion field
  GammaMatrixSet    *gmset_s      = NULL;
  Fopr_Clover_eo    *fopr_c_s     = NULL;
  Solver            *solver_s     = NULL;
  Fprop             *fprop_s      = NULL;
  Source            *source_s     = NULL;

  {
    string str_solver_type = params_solver_s -> get_string("solver_type");
    string str_gmset_type  = params_clover_s -> get_string("gamma_matrix_type");
    string str_source_type = params_source_s -> get_string("source_type");

    Bridge::vout.general("solver_type(s)      : %s\n", str_solver_type.c_str());
    Bridge::vout.general("gamma_matrix_type(s): %s\n", str_gmset_type.c_str());
    Bridge::vout.general("source_type(s)      : %s\n", str_source_type.c_str());

    gmset_s      = GammaMatrixSet::New(str_gmset_type);
    
    fopr_c_s     = new Fopr_Clover_eo(str_gmset_type);
    solver_s     = Solver::New(str_solver_type, fopr_c_s);
    fprop_s      = new Fprop_Standard_eo(solver_s);
    source_s     = Source::New(str_source_type);
  }


  //=====================================================================
  // end of initialization part
  //=====================================================================



  //////////// T. Doi ////////////
  // setup gauge config file names
  std::vector<std::string> gfile_list;
  {
    std::string str;
    std::ifstream ifs(fname_gfile_list);
    if ( ! ifs.is_open() ) {
      vout.crucial(vl, "cannot open %s\n",fname_gfile_list);
      abort();
    }

    while( getline(ifs, str) ) {
      gfile_list.push_back(str);
    }
  }


  // ========================================
  // ++++++++++++++++++++++ VOJTA PART BEGINS

  int noises[4] = {1,4,10};
  for(int i=0;i<3;i++){

  vout.general("\n\n\t ======================================");
  vout.general("\n\t ============== NUMBER OF NOISES = %2d ",noises[i]);
  vout.general("\n\t                beginning at           %s\n\n", LocalTime());

  if(Communicator::self()==0){
    int grid_coor[4];
    printf("\n");
    for(int ii=0;ii<Xnodes*Ynodes*Znodes*Tnodes;ii++){
      Communicator::grid_coord(grid_coor, ii);
      printf("MPI = %3i, grid coor = %2i-%2i-%2i-%2i\n", 
             ii, grid_coor[0],grid_coor[1],grid_coor[2],grid_coor[3]);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  int N_noises=noises[i];


  // ========================================
  // initialize propagator arrays

  // global source propagator
  double *prop_ud_wall   = new double[XYZTnodeSites * 3*4*3*4 *2];
  memset(prop_ud_wall,0,sizeof(prop_ud_wall));
  
  double *prop_s_wall    = new double[XYZTnodeSites * 3*4*3*4 *2];
  memset(prop_s_wall,0,sizeof(prop_s_wall));

  // noise vector propagators
  double *prop_noise  = new double[N_noises * XYZTnodeSites * 3*4*3*4 *2];
  memset(prop_noise,0,sizeof(prop_noise));


  // ========================================
  //initialize source class and generate point, wall and noise sources :)

  class_sources *sources = new class_sources();
 
  sources->generate_wall_source();
  sources->generate_noise_source_vector(N_noises,"Z(4)");

  //---------------------------------------------------------------------
  // main loop w.r. to gauge configurations
  //---------------------------------------------------------------------

  for(int iarg = 0; iarg <1; iarg++){
//  for(int iarg = 0; iarg < gfile_list.size(); iarg++){


    string ifname(gfile_list[iarg]);
    string base = ifname.substr(ifname.find_last_of('/')+1);

 
    //---------------------------------------------------------------------
    // ####  Set up a gauge configuration et al ####
    //---------------------------------------------------------------------
    {
      // not sure whether barrier is necessary, but just in case...
      //BGNET_GlobalBarrier();
      //BGNET_SwitchMPI();
      //MPI_Barrier(MPI_COMM_WORLD);

      vout.general("\n\t@@@ reading %s\n\t@@@ read conf(start) %s\n", ifname.c_str(), LocalTime());

      GaugeConfig *gconf_read = new GaugeConfig("ILDG_Parallel");
      gconf_read->read_file((Field *)U, ifname);
      delete gconf_read;

      vout.general("\n\t@@@ reading %s DONE\n\t@@@ read conf(end) %s\n\n", ifname.c_str(), LocalTime());

      //MPI_Barrier(MPI_COMM_WORLD);
      //BGNET_SwitchMPI();
      //BGNET_GlobalBarrier();
    }
      
    // ---- Vojta
    // gauge fixing
      
    {
      RandomNumbers_Mseries *rand = new RandomNumbers_Mseries(1100);
      
      GaugeFixing_Coulomb *gauge_fixing = new GaugeFixing_Coulomb(rand);        
      gauge_fixing->set_parameters(params_gauge_fixing);
      
      //gauge_fixing->fix(*U_fixed,*U);
      *U_fixed=*U;
    }

    // ========================================
    //checking the plaquette  and Polyakov loop before and after gauge fixing
    {
      Staples        *staple     = new Staples;
      double   plaq   = staple -> plaquette(*U);
      vout.general(vl, "plaq (original) = %18.14f\n", plaq);
      plaq   = staple -> plaquette(*U_fixed);
      vout.general(vl, "plaq (fixed)    = %18.14f\n\n", plaq);
      delete staple;

      PolyakovLoop       *pl     = new PolyakovLoop;
      dcomplex   ploop   = pl -> measure_ploop(*U);
      vout.general(vl, "Polyakov Loop (original) = %18.14f\n", ploop);
      ploop   = pl -> measure_ploop(*U_fixed);
      vout.general(vl, "Polyakov Loop (fixed)    = %18.14f\n\n\n", ploop);
      delete pl;
    }

    // ========================================
    // setting solver parameters...
 
        fopr_c_ud    -> set_parameters(*params_clover_ud);
        fopr_c_ud    -> set_config    (U_fixed, U_fixed);
        solver_ud    -> set_parameters(*params_solver_ud);

        fopr_c_s    -> set_parameters(*params_clover_s);
        fopr_c_s    -> set_config    (U_fixed, U_fixed);
        solver_s    -> set_parameters(*params_solver_s);  

    //---------------------------------------------------------------------
    // ####  Gauge configuration et al set ####
    //---------------------------------------------------------------------


    // ========================================
    // calculation of noise propagators 


    for(int i_noise=0; i_noise<N_noises; i_noise++){


      // dummy propagator
      double *prop_ud   = new double[XYZTnodeSites * 3*4*3*4 * 2];

      // loop over nT
      for(int iT_noise_src_pos=0;iT_noise_src_pos<CommonParameters::Lt();iT_noise_src_pos++){
//      for(int iT_noise_src_pos=7;iT_noise_src_pos<8;iT_noise_src_pos+=15){

        vout.general("\n\n ++++++ NOISE num %2d, from time slice %2d\n",
                  i_noise, iT_noise_src_pos);
        
        char label[256];
        snprintf(label,sizeof(label), "NOISE num %2d, from time slice %2d",
                 i_noise, iT_noise_src_pos);

    
        // solve the propagator
        propagators_solve(label,
                          fprop_ud, prop_ud,  
                          iT_noise_src_pos, sources->get_noise_ixyz(i_noise));          
        
       // selecto only iT=iT slice        
        iT_slice_selection(prop_noise,prop_ud, i_noise ,iT_noise_src_pos);  
        MPI_Barrier(MPI_COMM_WORLD);  

   
      }  //iT_noise_src_pos
          
      delete[] prop_ud;
    } //i_noise 




    //---------------------------------------------------------------------
    // loop w.r. to source positions
    //---------------------------------------------------------------------

//    for(int iT_src_pos=0;iT_src_pos<CommonParameters::Lt();iT_src_pos++){
    for(int iT_src_pos=2;iT_src_pos<3;iT_src_pos+=14){

      vout.general("\n\t@@@ calculation for source position at %2d start: \t%s @@@\n\n",
                   iT_src_pos, LocalTime());
 
      // ========================================
      // solve propagators
     {
      propagators_solve("WALL ud",
                        fprop_ud, 
                        prop_ud_wall,  
                        iT_src_pos, sources->get_wall_ixyz());          
      propagators_solve("WALL s",
                        fprop_s, 
                        prop_s_wall,  
                        iT_src_pos, sources->get_wall_ixyz());          
     }

      
      // ========================================
      // set prefix for this particular run
      char pr[50];
      snprintf(pr,sizeof(pr),"nn_%02d_",N_noises); 
      //snprintf(pr,sizeof(pr),""); 

      // ========================================
      //initialize single hadron class and set all the parameters

      class_hadron Hadron_wall(prop_ud_wall,prop_s_wall);
      
      Hadron_wall.set_base_name(base);
      Hadron_wall.set_prefix_name(pr);
      Hadron_wall.set_source_position(iT_src_pos);

      // ========================================
      //initialize two-hadron class and set all the parameters
      class_two_hadrons Two_hadrons(prop_ud_wall,prop_s_wall, prop_noise,sources);

      Two_hadrons.set_base_name(base);
      Two_hadrons.set_prefix_name(pr);
      Two_hadrons.set_source_position(iT_src_pos);
      Two_hadrons.set_noise_number(N_noises);


      // ========================================
      // run all green functions 


      //Hadron_wall.run_all_GF();

      Two_hadrons.run_all_GF();

 
      ///////////////////////////////////////////////////////////
   

    } // iT_src_pos - loop over source positions
  } // iarg - loop over configurations
  
  // ####  tidy up  ####
  // VOJTA
  delete sources;

  delete[] prop_s_wall;
  delete[] prop_ud_wall;

  delete[] prop_noise;
  } //N_noises
  
  // from template
  delete U_fixed;
  delete U;

  delete fopr_c_s;
  delete solver_s;
  delete fprop_s;
  delete source_s;

  delete fopr_c_ud;
  delete solver_ud;
  delete fprop_ud;
  delete source_ud;

  delete params_solver_s;
  delete params_clover_s;
  delete params_source_s;

  delete params_solver_ud;
  delete params_clover_ud;
  delete params_source_ud;

  vout.general("\n\t@@@ JOB(END), %s\n\n", LocalTime());

  return 0;
}

//---------------------------------------------------------------------
/**
 * @brief propagator converter: bridge++ to als
 */
//---------------------------------------------------------------------

static void converter(std::valarray<Field_F>& sq, double prop[])
{
  typedef dcomplex COMPLEX;

  for(  int alphaI = 0; alphaI < 4; alphaI++){
    for(int cI     = 0; cI     < 3; cI++){
      int idxI = cI + 3*alphaI;

      double *tmp = sq[idxI].ptr(0);

      for(    int ixyzt  = 0; ixyzt < XYZTnodeSites; ixyzt++){
        for(  int alphaF = 0; alphaF < 4; alphaF++){
          for(int cF     = 0; cF     < 3; cF++){
            int idxF = cF + 3*alphaF;

            ((COMPLEX*)prop)[cF + 3*(alphaF + 4*(cI + 3*(alphaI + 4*(ixyzt))))]
              = ((COMPLEX*)tmp)[idxF + 12*ixyzt];
          }
        }
      }

    }
  }

}


//---------------------------------------------------------------------
/**
 calculating propagators from a given source --- almost all copied from template
 */
//---------------------------------------------------------------------

static void propagators_solve(string label,
                              Fprop* fprop, 
                              double * prop, 
                              int it_src, double * source){

      int iT_source=(it_src+100*Tsites) % Tsites;

      typedef std::valarray<Field_F> PropagatorSet;

      PropagatorSet sq(Nc * Nd);

      {
        vout.general("\n\t@@@ %s solver (start):\t%s\n", label.c_str(),LocalTime());

        for(int indx = 0; indx < Nc*Nd; indx++){
          sq[indx]  = 0.0;
        } 
      
        Field_F b;
        
        int    Nconv;
        double diff;
      
        vout.general("  color spin   Nconv      diff \n");
        for (  int ispin  = 0; ispin  < Nd; ++ispin){
          for (int icolor = 0; icolor < Nc; ++icolor) {
            int idx = icolor + Nc * ispin;

            b = 0.0;
            if (it_src/TnodeSites == TnodeCoor) {
 
              int it_node = it_src % TnodeSites;

              for (int ixyz = 0; ixyz < XYZnodeSites; ixyz++) {

                int isite = ixyz+it_node*XYZnodeSites;

                b.set(2 * idx + 0, isite, 0, source[2*ixyz]);
                b.set(2 * idx + 1, isite, 0, source[2*ixyz+1]);
              }
            }
          
            vout.general("\n\t +++ solver  color %d spin %d\n\n",icolor, ispin);
            fprop -> invert_D(sq[idx], b, Nconv, diff); 
            vout.general("   %2d   %2d   %6d   %12.4e\n",
                           icolor, ispin, Nconv, diff);


          }
          vout.general("\n");
        }
        vout.general("\n\t@@@ solver(end):  \t%s\n\n", LocalTime());
      }
      //---------------------------------------------------------------------
      // solvers end
      //---------------------------------------------------------------------

  converter(sq, prop);

}      


//---------------------------------------------------------------------
/**
 selecting only t=iT slice of the propagator 
 */
//---------------------------------------------------------------------

static void iT_slice_selection(double *prop_final,
                               double *prop_tmp,
                               int noise_i, int iT_pos){

  vout.general("\n ++++++ time slice selection at time %2d and index %d, \t%s\n",
                iT_pos, noise_i, LocalTime());

  int iT_noise_pos=(iT_pos+100*Tsites) % Tsites;

  if (iT_noise_pos/TnodeSites == TnodeCoor) {

    int iT_node = iT_noise_pos % TnodeSites;

//    printf("MPI %2i, xyztnode %2i,%2i,%2i,%2i .. time coords  %2i,%2i,  \n" ,
//          Communicator::self(),XnodeCoor,YnodeCoor,ZnodeCoor,TnodeCoor,
//          iT_noise_pos, iT_node);

    int noise_index = noise_i*XYZTnodeSites * 3*4*3*4;

    #define Dirac(     alpha, x)  (alpha +  4*(x))
    #define Color(     c,     x)  (c     +  3*(x))

    for(        int c     = 0; c    < 3; c++){
      for(      int d     = 0; d    < 4; d++){
        for(    int cP    = 0; cP   < 3; cP++){
          for(  int dP    = 0; dP   < 4; dP++){
            for(int ixyz  = 0; ixyz < XYZnodeSites; ixyz++){
 
              ((COMPLEX*)prop_final)[noise_index +
                                     Color(c, Dirac(d, 
                                           Color(cP, Dirac(dP, 
                                                 ixyz + XYZnodeSites * (iT_node)))))]
              = ((COMPLEX*)prop_tmp)[Color(c, Dirac(d, 
                                           Color(cP, Dirac(dP, 
                                                 ixyz + XYZnodeSites * (iT_node)))))];
          
            }
          }
        }
      }
    }

    #undef Dirac  
    #undef Color          
  }
                       
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

static const char* LocalTime()
{
  struct timeval tp;
  time_t         ptm;
  static char    str[128];

  gettimeofday(&tp,NULL);
  ptm = tp.tv_sec;
  strcpy(str,asctime(localtime(&ptm)));
  str[strlen(str)-1]='\0';

  return str;
}




