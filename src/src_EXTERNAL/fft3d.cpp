//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup FFT_3D
 * @author  T.Inoue based on N.Ishii
 * @brief   3D FFT (Izubuchi method) for MPI
 * 
 * Implemented via MPISCU
 */
//--------------------------------------------------------------------------


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "COMPLEX.h"

//#include <util/gjp.h>
//#include <util/smalloc.h>
//#include <util/error.h>
//#include <util/verbose.h>

#include "fft3d.h"


/* local from Bridge++ */
#include "communicator.h"
#include "commonParameters.h"


//using namespace cps;


/* ------------------------------------------------ */
/* the following is from util/FFT/mpi/FFT_3D_IZBC.C */
/* ------------------------------------------------ */


#ifndef NDEBUG
#define NDEBUG 0
#endif


//typedef std::complex<double> COMPLEX;

extern char* LocalTime();

static int is_initialized=0;

static int Xnodes,     Ynodes,     Znodes;
static int XnodeCoor,  YnodeCoor,  ZnodeCoor;
static int XnodeSites, YnodeSites, ZnodeSites;
static int Xsites,     Ysites,     Zsites;

static int XYZnodes, XYZnodeSites, XYZsites;

static int TnodeCoor;

static fftw_plan plan_x_forward;
static fftw_plan plan_y_forward;
static fftw_plan plan_z_forward;

static fftw_plan plan_x_backward;
static fftw_plan plan_y_backward;
static fftw_plan plan_z_backward;

//static COMPLEX* tmp = NULL;

static int nblk_x;
static int nblk_y;
static int nblk_z;

static COMPLEX* fftw_x_in;
static COMPLEX* fftw_y_in;
static COMPLEX* fftw_z_in;

static COMPLEX* fftw_x_out;
static COMPLEX* fftw_y_out;
static COMPLEX* fftw_z_out;

static MPI_Comm Comm_x;
static MPI_Comm Comm_y;
static MPI_Comm Comm_z;

static void initialize_static_variables()
{
  //if(!NDEBUG){ sync(); VRB.Func("FFT_3D_IZBC.C","initialize_static_variables"); }

  if (is_initialized!=0) return;

    Xnodes     = Communicator::npe(0);
    Ynodes     = Communicator::npe(1);
    Znodes     = Communicator::npe(2);
  
    XnodeCoor  = Communicator::ipe(0);
    YnodeCoor  = Communicator::ipe(1);
    ZnodeCoor  = Communicator::ipe(2);
    TnodeCoor  = Communicator::ipe(3);

    XnodeSites = CommonParameters::Nx();
    YnodeSites = CommonParameters::Ny();
    ZnodeSites = CommonParameters::Nz();
  
  Xsites     = XnodeSites * Xnodes;
  Ysites     = YnodeSites * Ynodes;
  Zsites     = ZnodeSites * Znodes;
  
  XYZnodes     = Xnodes     * Ynodes     * Znodes;
  XYZnodeSites = XnodeSites * YnodeSites * ZnodeSites;
  XYZsites     = Xsites     * Ysites     * Zsites;
  
  //tmp = (COMPLEX*) smalloc(sizeof(COMPLEX) * XYZnodeSites);

  nblk_x = YnodeSites * ZnodeSites / Xnodes;
  nblk_y = ZnodeSites * XnodeSites / Ynodes;
  nblk_z = XnodeSites * YnodeSites / Znodes;

  if (YnodeSites * ZnodeSites % Xnodes != 0 ) nblk_x = nblk_x + 1 ;
  if (ZnodeSites * XnodeSites % Ynodes != 0 ) nblk_y = nblk_y + 1 ;
  if (XnodeSites * YnodeSites % Znodes != 0 ) nblk_z = nblk_z + 1 ;
       
  fftw_x_in  = (COMPLEX*) fftw_malloc(sizeof(COMPLEX) * Xsites);
  fftw_y_in  = (COMPLEX*) fftw_malloc(sizeof(COMPLEX) * Ysites);
  fftw_z_in  = (COMPLEX*) fftw_malloc(sizeof(COMPLEX) * Zsites);

  fftw_x_out = (COMPLEX*) fftw_malloc(sizeof(COMPLEX) * Xsites);
  fftw_y_out = (COMPLEX*) fftw_malloc(sizeof(COMPLEX) * Ysites);
  fftw_z_out = (COMPLEX*) fftw_malloc(sizeof(COMPLEX) * Zsites);

  plan_x_forward  = fftw_plan_dft_1d(Xsites, (fftw_complex*)fftw_x_in, (fftw_complex*)fftw_x_out, FFTW_FORWARD,  FFTW_ESTIMATE);
  plan_y_forward  = fftw_plan_dft_1d(Ysites, (fftw_complex*)fftw_y_in, (fftw_complex*)fftw_y_out, FFTW_FORWARD,  FFTW_ESTIMATE);
  plan_z_forward  = fftw_plan_dft_1d(Zsites, (fftw_complex*)fftw_z_in, (fftw_complex*)fftw_z_out, FFTW_FORWARD,  FFTW_ESTIMATE);
  
  plan_x_backward = fftw_plan_dft_1d(Xsites, (fftw_complex*)fftw_x_in, (fftw_complex*)fftw_x_out, FFTW_BACKWARD, FFTW_ESTIMATE);
  plan_y_backward = fftw_plan_dft_1d(Ysites, (fftw_complex*)fftw_y_in, (fftw_complex*)fftw_y_out, FFTW_BACKWARD, FFTW_ESTIMATE);
  plan_z_backward = fftw_plan_dft_1d(Zsites, (fftw_complex*)fftw_z_in, (fftw_complex*)fftw_z_out, FFTW_BACKWARD, FFTW_ESTIMATE);
  
  int XYTnodeCoor = XnodeCoor + Xnodes * (YnodeCoor + Ynodes * TnodeCoor);
  int YZTnodeCoor = YnodeCoor + Ynodes * (ZnodeCoor + Znodes * TnodeCoor);
  int ZXTnodeCoor = ZnodeCoor + Znodes * (XnodeCoor + Xnodes * TnodeCoor);

  MPI_Comm_split(MPI_COMM_WORLD, YZTnodeCoor, XnodeCoor, &Comm_x);
  MPI_Comm_split(MPI_COMM_WORLD, ZXTnodeCoor, YnodeCoor, &Comm_y);
  MPI_Comm_split(MPI_COMM_WORLD, XYTnodeCoor, ZnodeCoor, &Comm_z);

  is_initialized = 1;

  //if(!NDEBUG){ sync(); VRB.FuncEnd("FFT_3D_IZBC.C","initialize_static_variables"); }
}

//--------------------------------------------------------------------------
/**
 * @brief
 *
 *
 */
//--------------------------------------------------------------------------

void FFT3D(double* buf,int sign)
{
  //if(!NDEBUG) { sync(); VRB.Func("FFT_3D_IZBC.C","FFT3D(double*,int)"); }
  if (is_initialized==0) initialize_static_variables();

  COMPLEX *tmp_x = new COMPLEX[Xsites*nblk_x];
  COMPLEX *tmp_y = new COMPLEX[Ysites*nblk_y];
  COMPLEX *tmp_z = new COMPLEX[Zsites*nblk_z];

  COMPLEX *send_x = new COMPLEX[Xsites];
  COMPLEX *send_y = new COMPLEX[Ysites];
  COMPLEX *send_z = new COMPLEX[Zsites];

  int pos;

  //--------------------------------------------------------------------------
  // FFT along x direction
  //--------------------------------------------------------------------------
  //if(!NDEBUG) { sync(); VRB.Flow("FFT_3D_IZBC.C","FFT3D","X-direction"); }
  for(int iblk = 0; iblk < nblk_x; iblk++){
    // phase 1
    for( int i = 0; i < Xsites; i++){
      pos = i + iblk*XnodeSites*Xnodes;
      if(pos < XYZnodeSites){ 
        send_x[i] = ((COMPLEX*)buf)[pos];
      }	else{
        send_x[i] = 0.0;   //Inputing as a Complex.
      }
    }
    MPI_Alltoall(send_x,    2*XnodeSites, MPI_DOUBLE,
		 fftw_x_in, 2*XnodeSites, MPI_DOUBLE,
		 Comm_x);
    // phase 2
    if(sign == FFTW_FORWARD){
      fftw_execute(plan_x_forward);
    }
    else{
      fftw_execute(plan_x_backward);
    }
    // phase 3
    MPI_Alltoall(fftw_x_out,                                 2*XnodeSites, MPI_DOUBLE,
		 &((COMPLEX*)tmp_x)[iblk*XnodeSites*Xnodes], 2*XnodeSites, MPI_DOUBLE,
		 Comm_x);
  }

  for(    int iz = 0; iz < ZnodeSites; iz++)
    for(  int iy = 0; iy < YnodeSites; iy++)
      for(int ix = 0; ix < XnodeSites; ix++){
	const int ixyz = ix  +  XnodeSites * (iy  +  YnodeSites * iz);
	const int iyzx = iy  +  YnodeSites * (iz  +  ZnodeSites * ix);
	//((COMPLEX*)buf)[iyzx] = ((COMPLEX*)tmp)[ixyz];
	double* dst = (double*)&((COMPLEX*)buf)[iyzx];
	double* src = (double*)&((COMPLEX*)tmp_x)[ixyz];
	dst[0] = src[0];
	dst[1] = src[1];
      }
  //--------------------------------------------------------------------------
  // FFT along y direction
  //--------------------------------------------------------------------------
  //if(!NDEBUG) { sync(); VRB.Flow("FFT_3D_IZBC.C","FFT3D","Y-direction"); }
  for(int iblk = 0; iblk < nblk_y; iblk++){
    // phase 1
    for( int i = 0; i < Ysites;  i++){
      pos = i + iblk*YnodeSites*Ynodes;
      if(pos < XYZnodeSites){ 
        send_y[i] = ((COMPLEX*)buf)[pos];
      }	else{
        send_y[i] = 0.0;   //Inputing as a Complex. 
      }
    }
    MPI_Alltoall(send_y,    2*YnodeSites, MPI_DOUBLE,
		 fftw_y_in, 2*YnodeSites, MPI_DOUBLE,
		 Comm_y);
    // phase 2
    if(sign == FFTW_FORWARD){
      fftw_execute(plan_y_forward);
    }
    else{
      fftw_execute(plan_y_backward);
    }
    // phase 3
    MPI_Alltoall(fftw_y_out,                                 2*YnodeSites, MPI_DOUBLE,
		 &((COMPLEX*)tmp_y)[iblk*YnodeSites*Ynodes], 2*YnodeSites, MPI_DOUBLE,
		 Comm_y);
  }
  
  for(    int ix = 0; ix < XnodeSites; ix++)
    for(  int iz = 0; iz < ZnodeSites; iz++)
      for(int iy = 0; iy < YnodeSites; iy++){
	const int iyzx = iy  +  YnodeSites * (iz  +  ZnodeSites * ix);
	const int izxy = iz  +  ZnodeSites * (ix  +  XnodeSites * iy);
	//((COMPLEX*)buf)[izxy] = ((COMPLEX*)tmp)[iyzx];
	double* dst = (double*)&((COMPLEX*)buf)[izxy];
	double* src = (double*)&((COMPLEX*)tmp_y)[iyzx];
	dst[0] = src[0];
	dst[1] = src[1];
      }
  //--------------------------------------------------------------------------
  // FFT along z direction
  //--------------------------------------------------------------------------
  //if(!NDEBUG) { sync(); VRB.Flow("FFT_3D_IZBC.C","FFT3D","Z-direction"); }
  for(int iblk = 0; iblk < nblk_z; iblk++){
    // phase 1
    for( int i = 0; i < Zsites;  i++){
      pos = i + iblk*ZnodeSites*Znodes;
      if(pos < XYZnodeSites){ 
        send_z[i] = ((COMPLEX*)buf)[pos];
      }	else{
        send_z[i] = 0.0;   //Inputing as a Complex. 
      }
    }
    MPI_Alltoall(send_z,    2*ZnodeSites, MPI_DOUBLE,
		 fftw_z_in, 2*ZnodeSites, MPI_DOUBLE,
		 Comm_z);
    // phase 2
    if(sign == FFTW_FORWARD){
      fftw_execute(plan_z_forward);
    }
    else{
      fftw_execute(plan_z_backward);
    }
    // phase 3
    MPI_Alltoall(fftw_z_out,                                 2*ZnodeSites, MPI_DOUBLE,
		 &((COMPLEX*)tmp_z)[iblk*ZnodeSites*Znodes], 2*ZnodeSites, MPI_DOUBLE,
		 Comm_z);
  }
  
  for(    int iy = 0; iy < YnodeSites; iy++)
    for(  int ix = 0; ix < XnodeSites; ix++)
      for(int iz = 0; iz < ZnodeSites; iz++){
	const int izxy = iz  +  ZnodeSites * (ix  +  XnodeSites * iy);
	const int ixyz = ix  +  XnodeSites * (iy  +  YnodeSites * iz);
	//((COMPLEX*)buf)[ixyz] = ((COMPLEX*)tmp)[izxy];
	double* dst = (double*)&((COMPLEX*)buf)[ixyz];
	double* src = (double*)&((COMPLEX*)tmp_z)[izxy];
	dst[0] = src[0];
	dst[1] = src[1];
      }
  //if(!NDEBUG) { sync(); VRB.FuncEnd("FFT_3D_IZBC.C","FFT3D"); }

  delete [] tmp_x;
  delete [] tmp_y;
  delete [] tmp_z;

  delete [] send_x;
  delete [] send_y;
  delete [] send_z;

}
