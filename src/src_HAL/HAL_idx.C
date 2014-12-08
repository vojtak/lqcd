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


#include <util/halqcd/HAL_config.h>

//////////////////////////////
// avoid the compile error
#if 1==1
#undef  COMPLEX_TYPE
#define COMPLEX_TYPE COMPLEX_TYPE_CPP
#endif
#include <util/halqcd/COMPLEX.h>
//////////////////////////////

#include <util/halqcd/HAL_idx.h>


HAL_START_NAMESPACE


//--------------------------------------------------------------------------

int HAL_idx::eps[]={
  0,1,2, 1,
  1,2,0, 1,
  2,0,1, 1,
  
  0,2,1, -1,
  2,1,0, -1,
  1,0,2, -1
};

int   HAL_idx::icg5[] = {1,   0,    3,   2};
Float HAL_idx::zcg5[] = {1.0, -1.0, 1.0, -1.0};

int   HAL_idx::icg4[] = {3,   2,    1,    0};
Float HAL_idx::zcg4[] = {1.0, -1.0, -1.0, 1.0};

int   HAL_idx::icg45[] = {1,    0,   3,   2};
Float HAL_idx::zcg45[] = {-1.0, 1.0, 1.0, -1.0};

// gamma0 = unit matrix
// gamma1 = gamma_x
// gamma2 = gamma_y
// gamma3 = gamma_z
// gamma4 = gamma_t

// C = + \gamma4 * \gamma2

// complex version for [C \gamma_mu]
int     HAL_idx::icgm[] = {
  3, 2, 1, 0, // C
  0, 1, 2, 3, // C * \gamma_1
  0, 1, 2, 3, // C * \gamma_2
  1, 0, 3, 2, // C * \gamma_3
  3, 2, 1, 0, // C * \gamma_4
  1, 0, 3, 2  // C * \gamma_5
};

COMPLEX HAL_idx::zcgm[] = {
  -1.0,        1.0,       -1.0,        1.0,        // C
  -COMPLEX_I,  COMPLEX_I,  COMPLEX_I, -COMPLEX_I,  // C * \gamma_1
   1.0,        1.0,       -1.0,       -1.0,        // C * \gamma_2
   COMPLEX_I,  COMPLEX_I, -COMPLEX_I, -COMPLEX_I,  // C * \gamma_3
   1.0,       -1.0,       -1.0,        1.0,        // C * \gamma_4
   1.0,       -1.0,        1.0,       -1.0         // C * \gamma_5
};

// complex version for [\gamma_mu C]
int     HAL_idx::igmc[] = {
  3, 2, 1, 0, //            C
  0, 1, 2, 3, // \gamma_1 * C
  0, 1, 2, 3, // \gamma_2 * C
  1, 0, 3, 2, // \gamma_3 * C
  3, 2, 1, 0, // \gamma_4 * C
  1, 0, 3, 2  // \gamma_5 * C
};

COMPLEX HAL_idx::zgmc[] = {
  -1.0,        1.0,       -1.0,        1.0,        //            C
  -COMPLEX_I,  COMPLEX_I,  COMPLEX_I, -COMPLEX_I,  // \gamma_1 * C
  -1.0,       -1.0,        1.0,        1.0,        // \gamma_2 * C
   COMPLEX_I,  COMPLEX_I, -COMPLEX_I, -COMPLEX_I,  // \gamma_3 * C
  -1.0,        1.0,        1.0,       -1.0,        // \gamma_4 * C
   1.0,       -1.0,        1.0,       -1.0         // \gamma_5 * C
};

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

HAL_END_NAMESPACE
