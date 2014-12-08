#ifndef IS_INCLUDED_HAL_IDX_H
#define IS_INCLUDED_HAL_IDX_H




#define Eps(i,n)  (HAL_idx::eps[i + 4*n])
#define IGM(i,n)  (HAL_idx::igm[i + 4*n])
#define ZGM(i,n)  (HAL_idx::zgm[i + 4*n])

namespace HAL_idx {

  extern int   eps[];
  extern int   icg5[];
  extern double zcg5[];

  extern int     igm[];
  extern COMPLEX zgm[];
}



#endif
