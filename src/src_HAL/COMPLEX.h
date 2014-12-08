//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup
 * @brief   definition for COMPLEX
 * @author  HAL QCD Collaboration
 * @author  T. Doi
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_COMPLEX_H
#define IS_INCLUDED_COMPLEX_H


#if COMPLEX_TYPE == COMPLEX_TYPE_C

#include <complex.h>
typedef double _Complex COMPLEX;
#define Real(x) (__real__ (x))
#define Imag(x) (__imag__ (x))
#define Conj(x) (~(x))
#define COMPLEX_I    (    1.0i)
#define COMPLEX_ZERO (0.0+0.0i)

#else // COMPLEX_TYPE == COMPLEX_TYPE_CPP

#include <complex>
typedef std::complex<double> COMPLEX;
#define Real(x) (real(x))
#define Imag(x) (imag(x))
#define Conj(x) (conj(x))
#define COMPLEX_I    (std::complex<double>(0.0,1.0))
#define COMPLEX_ZERO (std::complex<double>(0.0,0.0))

#endif


#endif
