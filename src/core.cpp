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


static void propagators_solve(string label,
                              Fprop* fprop_ud, Fprop* fprop_s, 
                              double * prop_ud, double * prop_s, 
                              int it_src, double * source);

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

static void construct_iT_src_lst(std::vector<int> &iT_src_lst, 
                                 std::vector<int> &Dirichlet_BC_lst,
                                 int iT_src_num, int iT_src_shift);

//------------Vojt
extern void hal_run_V(int *Nodes, int *NodeSites, int *NodeCoor,
                      double *prop_ud,  double *prop_s,
                      std::vector<int> iT_src_lst, int iX_src, int iY_src, int iZ_src,
                      std::string conf_name);

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

  double *prop_ud        = new double[XYZTnodeSites * 3*4*3*4 *2];
  double *prop_s         = new double[XYZTnodeSites * 3*4*3*4 *2];

  double *prop_ud_wall   = new double[XYZTnodeSites * 3*4*3*4 *2];
  double *prop_s_wall    = new double[XYZTnodeSites * 3*4*3*4 *2];

  double *prop_ud_noise  = new double[XYZTnodeSites * 3*4*3*4 *2];
  double *prop_s_noise   = new double[XYZTnodeSites * 3*4*3*4 *2];


  int N_sources=1;

  //initialize source class and generate wall and noise sources :)
  class_sources *sources = new class_sources();
      
  sources->generate_wall_source();
  sources->generate_noise_source_vector(N_sources);
 

  //---------------------------------------------------------------------
  // main loop w.r. to gauge configurations
  //---------------------------------------------------------------------

//  for(int iarg = 0; iarg <1; iarg++){
  for(int iarg = 0; iarg < gfile_list.size(); iarg++){
   printf("size = %d, \n", gfile_list.size());

    string ifname(gfile_list[iarg]);
    string base = ifname.substr(ifname.find_last_of('/')+1);

 
    //---------------------------------------------------------------------
    // ####  Set up a gauge configuration  ####
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

    // setting solver parameters...
 
        fopr_c_ud    -> set_parameters(*params_clover_ud);
        fopr_c_ud    -> set_config    (U_fixed, U_fixed);
        solver_ud    -> set_parameters(*params_solver_ud);

        fopr_c_s    -> set_parameters(*params_clover_s);
        fopr_c_s    -> set_config    (U_fixed, U_fixed);
        solver_s    -> set_parameters(*params_solver_s);  




    //---------------------------------------------------------------------
    // loop w.r. to source positions
    //---------------------------------------------------------------------

    for(int iT_src_pos=0;iT_src_pos<CommonParameters::Lt();iT_src_pos++){
    //for(int iT_src_pos=11;iT_src_pos<13;iT_src_pos++){

      vout.general("\n\t@@@ calculation for source position at %2d start: \t%s @@@\n\n",
                   iT_src_pos, LocalTime());
 

     {
      propagators_solve("WALL",
                        fprop_ud, fprop_s, 
                        prop_ud_wall,  prop_s_wall,  
                        iT_src_pos, sources->get_wall_ixyz());          
      propagators_solve("NOISE",
                        fprop_ud, fprop_s, 
                        prop_ud_noise, prop_s_noise, 
                        iT_src_pos, sources->get_noise_ixyz(0));          
    }


      //initialize hadron class
   //   class_hadron Hadron(prop_ud,prop_s);
      class_hadron Hadron_wall(prop_ud_wall,prop_s_wall);
      class_hadron Hadron_noise(prop_ud_noise,prop_s_noise);

   //   Hadron.set_base_name(base);
      Hadron_wall.set_base_name(base);
        char pr[50];
        snprintf(pr,sizeof(pr),"w_");
      Hadron_wall.set_prefix_name(pr);
      Hadron_noise.set_base_name(base);
        snprintf(pr,sizeof(pr),"n_");
      Hadron_noise.set_prefix_name(pr);
      
   //   Hadron.set_source_position(iT_src_pos);
      Hadron_wall.set_source_position(iT_src_pos);
      Hadron_noise.set_source_position(iT_src_pos);

      // run all green functions
   //   Hadron.run_all_GF();
      Hadron_wall.run_all_GF();
      Hadron_noise.run_all_GF();

      //
      
      // class_NBS_WF NBS_WF(prop_ud, prop_s);
      // NBS_WF.set_base_name(base);
      // NBS_WF.set_source_position(it_src_pos);
      // NBS_WF.set_noise_propagators(prop_noise);
      //
      // NBS_WF.calculate();

   ///////////////////====================================================
   // remnants of other codes
   ///*
    string dir_base="results/"+base+"/";
    //vout.general("dir = %s",dir_base);
    if(Communicator::nodeid() == 0){
      if(mkdir(dir_base.c_str(),0755)!=0){
        if(errno != EEXIST){
          vout.crucial("ERROR: cannot create a directory 'corr'\n");
          abort();
        }
      }
    }

   
      
      // -- Vojta
      // set the source position    
      {
        valarray<int> new_source(0,4);
        //vout.general("\nnewsource %i %i %i %i",new_source[0], new_source[1], new_source[2], new_source[3]);
        new_source[3]=iT_src_pos;

        params_source_ud->set_int_vector("source_position",new_source);
        params_source_s->set_int_vector("source_position",new_source);
        valarray<int> source_position_ud = params_source_ud->get_int_vector("source_position");
        vout.general("\n\nnew source position set:   %i %i %i %i\n\n",source_position_ud[0], source_position_ud[1], source_position_ud[2], source_position_ud[3]);
      }

      //////////// T. Doi ////////////
      // setup the source and boundary condition information
      std::vector<int> iT_src_lst;
      std::vector<int> Dirichlet_BC_lst;
      int iX_src, iY_src, iZ_src;
      {
        valarray<int> source_position_ud = params_source_ud->get_int_vector("source_position");
        valarray<int> source_position_s  = params_source_s ->get_int_vector("source_position");
        for(int mu = 0; mu < Ndim; mu++) {
          if ( source_position_ud[mu] != source_position_s[mu] ) {
            vout.crucial(vl, "up/down and strange have different centers.\n");
            abort();
          }
        }

        //  ishift[]: (t, x, y, z)-order
        int ishift[] = { source_position_ud[3], source_position_ud[0], source_position_ud[1], source_position_ud[2] };
   
        iX_src = (ishift[1] + 100*Xsites) % Xsites;
        iY_src = (ishift[2] + 100*Ysites) % Ysites;
        iZ_src = (ishift[3] + 100*Zsites) % Zsites;
   
        construct_iT_src_lst( iT_src_lst, Dirichlet_BC_lst, 1, ishift[0] );
      }

      ////////////////////////////////
      
      //---------------------------------------------------------------------
      // solvers begin 
      //---------------------------------------------------------------------
 
      typedef std::valarray<Field_F> PropagatorSet;


      PropagatorSet sq_ud(Nc * Nd);
      PropagatorSet sq_s(Nc * Nd);

      // ud solver
      {
        vout.general("\n\t@@@ solver ud(start):\t%s,\tkappa=XXX, Csw=YYY\n", LocalTime());

   //     fopr_c_ud    -> set_parameters(*params_clover_ud);
     //   fopr_c_ud    -> set_config    (U_fixed, U_fixed);
       //   solver_ud    -> set_parameters(*params_solver_ud);
        source_ud    -> set_parameters(*params_source_ud);

        for(int indx = 0; indx < Nc*Nd; indx++) sq_ud[indx] = 0.0;
      
        Field_F b;
        b = 0.0;

        int    Nconv;
        double diff;
      
        vout.general(vl, "  color spin   Nconv      diff \n");
        for (  int ispin  = 0; ispin  < Nd; ++ispin){
          for (int icolor = 0; icolor < Nc; ++icolor) {
            int idx = icolor + Nc * ispin;
            source_ud->set(b, idx);
	  
            fprop_ud -> invert_D(sq_ud[idx], b, Nconv, diff); 

            vout.general(vl, "   %2d   %2d   %6d   %12.4e\n",
                         icolor, ispin, Nconv, diff);
          }
          vout.general(vl, "\n");
        }
        vout.general("\n\t@@@ solver(end):  \t%s,\titer/max_iter=XXX/YYY\n\n", LocalTime());
      }
      // solver s
      {
        vout.general("\n\t@@@ solver s(start):\t%s,\tkappa=XXX, Csw=YYY\n", LocalTime());

  //      fopr_c_s    -> set_parameters(*params_clover_s);
    //    fopr_c_s    -> set_config    (U_fixed, U_fixed);
      //  solver_s    -> set_parameters(*params_solver_s);
        source_s    -> set_parameters(*params_source_s);
 
        for(int indx = 0; indx < Nc*Nd; indx++) sq_s[indx] = 0.0;
      
        Field_F b;
        b = 0.0;
      
        int    Nconv;
        double diff;
      
        vout.general(vl, "  color spin   Nconv      diff\n");
        for (  int ispin  = 0; ispin  < Nd; ++ispin) {
          for (int icolor = 0; icolor < Nc; ++icolor) {
            int idx = icolor + Nc * ispin;
            source_s->set(b, idx);
	  
            fprop_s -> invert_D(sq_s[idx], b, Nconv, diff);
          }
          vout.general(vl, "\n");
        }
        vout.general("\n\t@@@ solver(end):  \t%s,\titer/max_iter=XXX/YYY\n\n", LocalTime());
      }


      converter(sq_ud, prop_ud);
      converter(sq_s, prop_s);


      //---------------------------------------------------------------------
      // solvers end
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // measurement begin
      //---------------------------------------------------------------------

      
      //-- two point correlators
      //---------------------------------------------------------------------
      vout.general(vl, "\n\t@@@ 2-point correlators(begin):  \t%s\n", LocalTime());

      //set directory name
      string dirname=dir_base+"corr/";
      //if(Communicator::nodeid() == 0){
        if(mkdir(dirname.c_str(),0755)!=0){
          if(errno != EEXIST){
            vout.crucial("ERROR: cannot create a directory 'corr'\n");
            abort();
          }
        }
      //}

      Corr2pt_4spinor  corr(gmset_ud);

      {
        vout.general(vl, "correlator:  pion\n");
        std::valarray<dcomplex> pion(CommonParameters::Lt());
        corr.meson_corr(pion,
                        gmset_ud->get_GM(gmset_ud->GAMMA5),
                        gmset_ud->get_GM(gmset_ud->GAMMA5),
                        sq_ud,
                        sq_ud);
	
        //vout.general(vl, "create file name\n");

        char wfile[FILENAME_MAX];
        snprintf(wfile,	sizeof(wfile), "%spion%02d",
                 dirname.c_str(),iT_src_pos);
        if(Communicator::nodeid() == 0){
          corr_print(pion, wfile, iT_src_pos);
        }
      }

      {
        vout.general(vl, "correlator:  kaon\n");
        std::valarray<dcomplex> kaon(CommonParameters::Lt());
        corr.meson_corr(kaon,
                        gmset_ud->get_GM(gmset_ud->GAMMA5),
                        gmset_ud->get_GM(gmset_ud->GAMMA5),
                        sq_s,
                        sq_ud);
	
        //vout.general(vl, "create file name\n");
        char wfile[FILENAME_MAX];
        snprintf(wfile,	sizeof(wfile), "%skaon%02d",
                 dirname.c_str(),iT_src_pos);

        if(Communicator::nodeid() == 0){
          corr_print(kaon, wfile, iT_src_pos);
        }
      }

      {
        vout.general(vl, "correlator:  eta\n");
        std::valarray<dcomplex> eta(CommonParameters::Lt());
        corr.meson_corr(eta,
                        gmset_ud->get_GM(gmset_ud->GAMMA5),
                        gmset_ud->get_GM(gmset_ud->GAMMA5),
                        sq_s,
                        sq_s);
	
        //vout.general(vl, "create file name\n");
        char wfile[FILENAME_MAX];
        snprintf(wfile,	sizeof(wfile), "%seta%02d",
                 dirname.c_str(),iT_src_pos);

        if(Communicator::nodeid() == 0){
          corr_print(eta, wfile, iT_src_pos);
        }
      }

/*
      {
        vout.general(vl, "correlator:  proton\n");
        std::valarray<dcomplex> proton(CommonParameters::Lt());
        corr.proton_corr(proton,
                        gmset_ud->get_GM(gmset_ud->GAMMA5),  // DONT WORK PROPERLY
                        sq_ud,
                        sq_ud);
	
        //vout.general(vl, "create file name\n");
        char wfile[FILENAME_MAX];
        snprintf(wfile,	sizeof(wfile), "%sproton%02d",
                 dirname.c_str(),iT_src_pos);

        if(Communicator::nodeid() == 0){
          corr_print(proton, wfile, iT_src_pos);
        }
      }
*/

      vout.general(vl, "\n\t@@@ 2-point correlators(end):  \t%s\n", LocalTime());
      


      //-- HAL code
      //---------------------------------------------------------------------
      //////////// T. Doi ////////////  HAL RUN PART
      vout.general(vl, "\n\t@@@ HAL_RUN part(begin):  \t%s\n", LocalTime());
      {
        // not sure whether barrier is necessary, but just in case...
        //BGNET_GlobalBarrier();
        //BGNET_SwitchMPI();
        //MPI_Barrier(MPI_COMM_WORLD);

        int Nodes    [] = { Xnodes,     Ynodes,     Znodes,     Tnodes     };
        int NodeSites[] = { XnodeSites, YnodeSites, ZnodeSites, TnodeSites };
        int NodeCoor [] = { XnodeCoor,  YnodeCoor,  ZnodeCoor,  TnodeCoor  };

     
        //hal_run(Nodes, NodeSites, NodeCoor,
        //        prop_ud,  prop_s,
        //        smear_ud, smear_s,
        //        iT_i_BS, iT_f_BS, flg_iT_BS_rev,
        //        iT_src_lst, iX_src, iY_src, iZ_src,
        //        base);

        hal_run_V(Nodes, NodeSites, NodeCoor,
                  prop_ud,  prop_s,
                  iT_src_lst, iX_src, iY_src, iZ_src,
                  base);

        //MPI_Barrier(MPI_COMM_WORLD);
        //BGNET_SwitchMPI();
        //BGNET_GlobalBarrier();
        vout.general(vl, "\n\t@@@ HAL_RUN part(end):  \t%s\n", LocalTime());
      }
      ////////////////////////////////
      //---------------------------------------------------------------------
      // measurement end
      //---------------------------------------------------------------------
 
   //*/
   ///////////////////====================================================
   // remnants of other codes END
   

    } // iT_src_pos - loop over source positions
  } // iarg - loop over configurations

  // ####  tydy up  ####
  delete[] sources;
  
  delete[] prop_s;
  delete[] prop_ud;

  delete[] prop_s_wall;
  delete[] prop_ud_wall;
  delete[] prop_s_noise;
  delete[] prop_ud_noise;

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
 calculating ud and s propagators from a given source
 */
//---------------------------------------------------------------------

static void propagators_solve(string label,
                              Fprop* fprop_ud, Fprop* fprop_s, 
                              double * prop_ud, double * prop_s, 
                              int it_src, double * source){

      int iT_source=(it_src+100*Tsites) % Tsites;

      typedef std::valarray<Field_F> PropagatorSet;

      PropagatorSet sq_ud(Nc * Nd);
      PropagatorSet sq_s(Nc * Nd);

      {
        vout.general("\n\t@@@ with %s solver ud and s(start):\t%s,\tkappa=XXX, Csw=YYY\n", label.c_str(),LocalTime());

        for(int indx = 0; indx < Nc*Nd; indx++){
          sq_ud[indx] = 0.0;
          sq_s[indx]  = 0.0;
        } 
      
        Field_F b;
        
        int    Nconv;
        double diff;
      
        vout.general("  color spin   Nconv      diff \n");
        for (  int ispin  = 0; ispin  < Nd; ++ispin){
          for (int icolor = 0; icolor < Nc; ++icolor) {
            int idx = icolor + Nc * ispin;

            b = 0.0;
            if (it_src/TnodeSites == Communicator::ipe(3)) {
 
              int it_node = it_src % TnodeSites;

              for (int ixyz = 0; ixyz < XYZnodeSites; ixyz++) {

                int isite = ixyz+it_node*XYZnodeSites;

                b.set(2 * idx + 0, isite, 0, source[2*ixyz]);
                b.set(2 * idx + 1, isite, 0, source[2*ixyz+1]);
              }
            }
          
            vout.general("\n\t+++ ud solver  color %d spin %d\n\n",icolor, ispin);
            fprop_ud -> invert_D(sq_ud[idx], b, Nconv, diff); 
            vout.general("   %2d   %2d   %6d   %12.4e\n",
                           icolor, ispin, Nconv, diff);

            vout.general("\n\t+++ s  solver  color %d spin %d\n\n",icolor, ispin);
            fprop_s  -> invert_D(sq_s[idx],  b, Nconv, diff); 
            vout.general("   %2d   %2d   %6d   %12.4e\n",
                           icolor, ispin, Nconv, diff);

          }
          vout.general("\n");
        }
        vout.general("\n\t@@@ solver(end):  \t%s,\titer/max_iter=XXX/YYY\n\n", LocalTime());
      }
      //---------------------------------------------------------------------
      // solvers end
      //---------------------------------------------------------------------

  converter(sq_ud, prop_ud);
  converter(sq_s,  prop_s );


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



///////////////////====================================================
// remnants of other codes END
   

//---------------------------------------------------------------------
//---------------------------------------------------------------------

static void construct_iT_src_lst(std::vector<int> &iT_src_lst,
                                 std::vector<int> &Dirichlet_BC_lst, 
                                 int iT_src_num, int iT_src_shift)
{
  if ( Tsites % iT_src_num != 0 ){
    fprintf(stderr, "Tsites %% iT_src_num = %d %% %d != 0\n", Tsites, iT_src_num);
    exit(1);
  }

  int iT_src_sep = Tsites / iT_src_num;

  int iT_src, iT_dbc;
  for(int i = 0; i < iT_src_num; i++){
    iT_src = (i * iT_src_sep + iT_src_shift + 100*Tsites) % Tsites;

    iT_src_lst.push_back(iT_src);

    if ( iT_src_sep % 2 == 0 ){
      iT_dbc = (iT_src + iT_src_sep / 2 -1 + 100*Tsites) % Tsites;
      Dirichlet_BC_lst.push_back(iT_dbc);

      iT_dbc = (iT_src + iT_src_sep / 2    + 100*Tsites) % Tsites;
      Dirichlet_BC_lst.push_back(iT_dbc);
    }
    else {
      iT_dbc = (iT_src + iT_src_sep / 2    + 100*Tsites) % Tsites;
      Dirichlet_BC_lst.push_back(iT_dbc);
    }
  }
}


