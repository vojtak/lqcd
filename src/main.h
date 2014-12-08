/**
*        @file    $Id:: main.h #$
*
*        @brief
*
*        @author  $LastChangedBy: matufuru $
*
*        @date    $LastChangedDate:: 2013-09-30 23:02:19 #$
*
*        @version $LastChangedRevision: 961 $
**/

#ifndef MAIN_INCLUDED
#define MAIN_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
using std::string;
#include <valarray>
using std::valarray;

#include "configure.h"
#include "defs.h"
#include "parameters.h"
#include "parameterManager_YAML.h"

#include "bridgeIO.h"
using Bridge::vout;

#include "timer.h"

const string filename_main_input = "main.yaml";
// const string filename_main_input = "stdin";

//- prototype declaration
#ifdef USE_TESTMANAGER
int run_testmanager(int argc, char **argv);

#else
int run_test();
#endif

class Parameters_Main : public Parameters
{
 public:
  Parameters_Main()
  {
    Register_int_vector("lattice_size", valarray<int>());
    Register_int_vector("grid_size", valarray<int>());

    Register_int("number_of_thread", 1);
    Register_int("number_of_color", 3);
    Register_string("log_filename", "NULL");
    Register_string("ildg_log_filename", "NULL");

    Register_string("verbose_level", "NULL");
  }
};
#endif
