/* standard */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstdlib>


/* standard parralelization */
#include <mpi.h>
#include <omp.h>


/* local from HAL */
#include "COMPLEX.h"
#include "HAL_indexes.h"

/* local from VOJTA */
#include "class_global_wrapper.h"
#include "class_sources.h"
#include "class_two_hadrons.h"


using namespace HAL_idx;

void class_two_hadrons::run_all_NBSwf(){

  run_GF("pion-sigma");

}

// ================================================================================================
// wrapper to run all correlators
//



// ================================================================================================
// print the correlator
//
void class_two_hadrons::NBSwf_print(double *correlator, string hadron_names)
{

  MPI_Barrier(MPI_COMM_WORLD);
  
  if(MPI_rank==0){

    // create directory and file name
    string dir="results/"+base+"/"+prefix+"NBSwaveF_two_hardons/";
    //printf("%s %s\n",prefix.c_str(),dir.c_str());
    create_directory(dir);
  
    char wfile[256];
    snprintf(wfile,sizeof(wfile), "%s%s.%02d",
             dir.c_str(),
             hadron_names.c_str(),
             iT_src);

    printf(" ++++++ corr_print : print %10s propagator to file \n                     %s\n                     %s\n",
           hadron_names.c_str(), wfile, LocalTime().c_str());
 
    // open output file
    string ofname(wfile);
    std::ofstream fout(ofname.c_str());
    fout.setf(std::ios::scientific);
    fout.precision(16);


    // print correlator
    for(int it = 0; it < Tsites; it++){
      int iT2 = (it + iT_src + 100*Tsites) % Tsites;
      //int iT2 = it;

      char line[1000];
      snprintf(line, sizeof(line), "%4d\t%1.16e %1.16e\n",
               it, correlator[2* iT2], correlator[2* iT2+1]);
      fout << line;
    }

    fout.close();

  }

  MPI_Barrier(MPI_COMM_WORLD);

}




