//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup HAL_idx
 * @brief   miscellaneous index implementations
 * @author  HAL QCD Collaboration
 * @author  N. Ishii
 * @author  T. Doi
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_HAL_IDX_H
#define IS_INCLUDED_HAL_IDX_H


#include <util/halqcd/HAL_config.h>
#include <util/halqcd/COMPLEX.h>

  int XnodeSites = CommonParameters::Nx();
  int YnodeSites = CommonParameters::Ny();
  int ZnodeSites = CommonParameters::Nz();
  int TnodeSites = CommonParameters::Nt();
  int XYZnodeSites = XnodeSites * YnodeSites * ZnodeSites;


HAL_START_NAMESPACE


//--------------------------------------------------------------------------

namespace HAL_idx {
  extern int   eps[];
  extern int   icg5[];
  extern Float zcg5[];
  extern int   icg4[];
  extern Float zcg4[];
  extern int   icg45[];
  extern Float zcg45[];

  extern int     icgm[];
  extern COMPLEX zcgm[];
  extern int     igmc[];
  extern COMPLEX zgmc[];

  extern int     igm[];
  extern COMPLEX zgm[];
}

inline int LessThanN(int n, int i)
{
  if(HAL_DEBUG)
    if(i < 0 || n <= i){
      ERR.General("","","ERROR:[HAL_idx.h] LessThanN: i = %d >= n = %d\n",i, n);
    }
  return i;
}


#define Spinor(nd, alpha, x)  (LessThanN(nd, alpha) + nd*(x))
#define Dirac(     alpha, x)  (LessThanN( 4, alpha) +  4*(x))
#define Pauli(     alpha, x)  (LessThanN( 2, alpha) +  2*(x))
#define Color(     c,     x)  (LessThanN( 3, c    ) +  3*(x))


#define Eps(i,n)  (HAL_idx::eps[LessThanN(4,i) + 4*LessThanN(6,n)])

#define ICG5(i)   (HAL_idx::icg5[LessThanN(4,i)])
#define ZCG5(i)   (HAL_idx::zcg5[LessThanN(4,i)])

#define ICG4(i)   (HAL_idx::icg4[LessThanN(4,i)])
#define ZCG4(i)   (HAL_idx::zcg4[LessThanN(4,i)])

#define ICG45(i)  (HAL_idx::icg45[LessThanN(4,i)])
#define ZCG45(i)  (HAL_idx::zcg45[LessThanN(4,i)])

// complex version for [C * \gamma_mu] or [\gamma_mu * C]
#define ICGM(i,n) (HAL_idx::icgm[LessThanN(4,i) + 4*LessThanN(6,n)])
#define ZCGM(i,n) (HAL_idx::zcgm[LessThanN(4,i) + 4*LessThanN(6,n)])

#define IGMC(i,n) (HAL_idx::igmc[LessThanN(4,i) + 4*LessThanN(6,n)])
#define ZGMC(i,n) (HAL_idx::zgmc[LessThanN(4,i) + 4*LessThanN(6,n)])

// IGM, ZGM are different from those written in hadron.C by N. Ishii !!!
#define IGM(i,n)  (HAL_idx::igm[LessThanN(4,i) + 4*LessThanN(6,n)])
#define ZGM(i,n)  (HAL_idx::zgm[LessThanN(4,i) + 4*LessThanN(6,n)])

//--------------------------------------------------------------------------


HAL_END_NAMESPACE

#endif
