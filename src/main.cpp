/**
*        @file    $Id:: main.cpp #$
*
*        @brief
*
*        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
*                 $LastChangedBy: matufuru $
*
*        @date    $LastChangedDate:: 2013-10-20 19:28:34 #$
*
*        @version $LastChangedRevision: 972 $
**/

#ifdef KEKBGQ
#include <mpi.h>
#endif

#include "main.h"
#include "threadManager_OpenMP.h"

extern int core(int argc,char* argv[]);

//====================================================================
int main(int argc, char *argv[])
{
  // ###  initial setup  ###
  Bridge::vout.init(std::cout);
  Bridge::vout.ildg_init(std::cout);
  Bridge::VerboseLevel vl = Bridge::PARANOIAC;

  // this includes MPI_Init_thread() & BGNET_Init()
  Communicator::init(&argc, &argv);


  // ####  parameter setup  ####
  Parameters_Main params_main;

  Parameters params_all;

  params_all.Register_Parameters("Main", &params_main);

  string filename_input = filename_main_input;
  if (filename_input == "stdin") {
    vout.general(vl, "input filename : ");
    std::cin >> filename_input;
    vout.general(vl, "%s\n", filename_input.c_str());
  } else {
    vout.general(vl, "input filename : %s\n", filename_input.c_str());
  }
  vout.general(vl, "\n");

  ParameterManager_YAML params_manager;
  params_manager.read_params(filename_input, &params_all);

  const valarray<int> lattice_size     = params_main.get_int_vector("lattice_size");
  const valarray<int> grid_size        = params_main.get_int_vector("grid_size");
  const int           number_of_thread = params_main.get_int       ("number_of_thread");
  const int           number_of_color  = params_main.get_int       ("number_of_color");
  const string        str_logfile      = params_main.get_string    ("log_filename");
  const string        str_ildg_logfile = params_main.get_string    ("ildg_log_filename");
  const string        str_vlevel       = params_main.get_string    ("verbose_level");


  //- initializations
  vl = vout.set_verbose_level(str_vlevel);

  std::ofstream logfile(str_logfile.c_str());
  std::ofstream ildg_logfile(str_ildg_logfile.c_str());

  if (str_logfile != "stdout") {
    Bridge::vout.init(logfile);
  }
  if (str_ildg_logfile != "stdout") {
    Bridge::vout.ildg_init(ildg_logfile);
  }

  //  CommonParameters::init(lattice_size, grid_size);
  CommonParameters::init(lattice_size, grid_size, number_of_color);
  Communicator::setup();

  // thread manager initialized.
  ThreadManager_OpenMP::init(number_of_thread);


  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();


  //- print input parameters
  vout.general(vl, "Parameters of Main:\n");
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(vl, "  lattice_size[%d] = %d\n", mu, lattice_size[mu]);
    vout.general(vl, "  grid_size[%d]    = %d\n", mu, grid_size[mu]);
  }
  vout.general(vl, "  number of thread = %d\n", number_of_thread);
  vout.general(vl, "  number of color  = %d\n", number_of_color);
  vout.general(vl, "  logfile          = %s\n", str_logfile.c_str());
  vout.general(vl, "  ildg_logfile     = %s\n", str_ildg_logfile.c_str());
  vout.general(vl, "  vlevel           = %s\n", str_vlevel.c_str());
  vout.general(vl, "\n");

  //- input parameter check
  int err = 0;
  err += ParameterCheck::non_NULL(str_logfile);
  err += ParameterCheck::non_NULL(str_ildg_logfile);

  if (err) {
    vout.crucial(vl, "Main: input parameters have not been set.\n");
    abort();
  }


  //- timestamp (starting time)
  Timer *timer = new Timer;
  timer->timestamp();


// #ifdef USE_TESTMANAGER
//   run_testmanager(argc, argv);
// #else
//   run_test();
// #endif
  core(argc,argv);

  //- timestamp (end time)
  timer->timestamp();


  // ####  tydy up  ####
  delete timer;


  // this includes MPI_Finalize()
  Communicator::finalize();

  return EXIT_SUCCESS;
}


//====================================================================
//============================================================END=====
