/* standard */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstdlib>

#include <sys/time.h>
#include <time.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

/* standard parralelization */
#include <mpi.h>
#include <omp.h>


/* local from HAL */
#include "COMPLEX.h"

/* local from VOJTA */
#include "HAL_indexes.h"

#include "class_global_wrapper.h"
#include "class_hadron.h"



void class_global_wrapper::create_directory(string path){

  string s_path=path;
  string dir = "";

  int pos = 0;

  while ((pos = s_path.find("/")) != string::npos) {

    dir+=s_path.substr(0, pos)+"/";

    if(mkdir(dir.c_str(),0755)!=0){
      if(errno != EEXIST){
        printf("ERROR: cannot create a directory %s\n",dir.c_str());
        abort();
      }
    }

    s_path.erase(0, pos + 1);
  }

}


string class_global_wrapper::LocalTime()
  {
    struct timeval tp;
    time_t         ptm;
    string         str;

    gettimeofday(&tp,NULL);
    ptm = tp.tv_sec;
    str=asctime(localtime(&ptm));

    return str;
  }

