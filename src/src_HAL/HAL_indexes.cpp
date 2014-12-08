//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup HAL_idx
 * @brief   miscellaneous index implementations
 * @brief   The definition of gamma matrices are Montvay's Dirac rep in Minkowski
 *          (which is same as Matsufuru's Dirac rep in Minkowski except for the sign of \gamma5)
 * @author  HAL QCD Collaboration
 * @author  N. Ishii
 * @author  T. Doi
 */
//--------------------------------------------------------------------------

#include "COMPLEX.h"
#include "HAL_indexes.h"
#include <stdio.h>


int HAL_idx::eps[]={
  0,1,2, 1,
  1,2,0, 1,
  2,0,1, 1,
  
  0,2,1, -1,
  2,1,0, -1,
  1,0,2, -1
};

int   HAL_idx::icg5[] = {1,   0,    3,   2};
double HAL_idx::zcg5[] = {1.0, -1.0, 1.0, -1.0};


//--------------------------------------------------------------------------
// originally written in hadron.C by N. Ishii
// However, the definition is changed by T. Doi !!!!

int HAL_idx::igm[]={
  0, 1, 2, 3, // unit matrix
  3, 2, 1, 0, // \gamma_x
  3, 2, 1, 0, // \gamma_y
  2, 3, 0, 1, // \gamma_z
  0, 1, 2, 3, // \gamma_t
  2, 3, 0, 1, // \gamma_5
};

/* Montvay's Dirac rep */
COMPLEX HAL_idx::zgm[]={
  1.0,         1.0,        1.0,         1.0,
 -COMPLEX_I,  -COMPLEX_I,  COMPLEX_I,   COMPLEX_I,
 -1.0,         1.0,        1.0,        -1.0,
 -COMPLEX_I,   COMPLEX_I,  COMPLEX_I,  -COMPLEX_I,
  1.0,         1.0,       -1.0,        -1.0,
 -1.0,        -1.0,       -1.0,        -1.0,
};

/* Ishii's Dirac rep */
/*
static COMPLEX zgm[]={
  1,0,         1,0,         1,0,         1,0,
  COMPLEX_I,   COMPLEX_I,  -COMPLEX_I,  -COMPLEX_I,
 -1.0,         1.0,         1.0,        -1.0,
  COMPLEX_I,  -COMPLEX_I,  -COMPLEX_I,   COMPLEX_I,
  1.0,         1.0,        -1.0,        -1.0,
 -1.0,        -1.0,        -1.0,        -1.0,
};
*/

//--------------------------------------------------------------------------

