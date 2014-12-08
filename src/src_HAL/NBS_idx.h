//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup NBS_idx
 * @brief   miscellaneous index implementations related to NBS wave functions
 * @author  HAL QCD Collaboration
 * @author  T. Doi
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_NBS_IDX_H
#define IS_INCLUDED_NBS_IDX_H

#include <util/halqcd/HAL_config.h>
#include <util/halqcd/HAL_idx.h>


HAL_START_NAMESPACE


//--------------------------------------------------------------------------
/**
 * @brief collections of index implementations related to NBS wave functions
 */
//--------------------------------------------------------------------------

namespace NBS_idx {

  //! idx for propagators obtained by solver (c,d: sink, cP, dP: source)
  size_t prop_slv_idx   (int c, int d, int cP, int dP, int ixyz, int it);
  int    prop_slv_cs_idx(int c, int d, int cP, int dP);

  //! idx for propagators used in NBS calc (c,d: sink, cP, dP: source)
  size_t prop_idx         (int c, int d, int ixyz, int cP, int dP, int it);
  int    prop_snk_idx_no3D(int c, int d);

  //! idx for (sink) hadron block
  //                    [ snk indices   ]  [            src indices                           ]
  int Boct_blk_idx     (int ixyz, int alp, int c1, int alp1, int c2, int alp2, int c3, int alp3, bool flg_NR_limit);
  int Boct_blk_idx_no3D(          int alp, int c1, int alp1, int c2, int alp2, int c3, int alp3, bool flg_NR_limit);
  int Boct_blk_src_idx (                   int c1, int alp1, int c2, int alp2, int c3, int alp3, bool flg_NR_limit);

  //                                    [ snk indices   ]  [src indices] 
  int Boct_blk_idx     (                int ixyz, int alp, int c_alp);
  int Boct_blk_idx_no3D(                          int alp, int c_alp);
  // (Here, effect of "flg_NR_limit" should be already taken into account in c_alp idx)

  //! idx related to NBS wave functions
  int wave_idx_fin (int dof_snk, int isnk, int ixyz,           int isrc); //!< final        idx for wave
  int wave_idx     (int dof_snk,           int ixyz, int isnk, int isrc); //!< intermediate idx for wave

  //                               (snk spinor)     (src spinor)
  int wave_idx_2Boct     (int ixyz, int*,            int*);
  int wave_idx_2Boct     (int ixyz, int,int,         int,int);
  int wave_idx_2Boct_no3D(          int*,            int*);
  int wave_idx_2Boct_no3D(          int,int,         int,int);

  int wave_idx_3Boct     (int ixyz, int*,            int*);
  int wave_idx_3Boct     (int ixyz, int,int,int,     int,int,int);
  int wave_idx_3Boct_no3D(          int*,            int*);
  int wave_idx_3Boct_no3D(          int,int,int,     int,int,int);

  int wave_idx_4Boct     (int ixyz, int*,            int*);
  int wave_idx_4Boct     (int ixyz, int,int,int,int, int,int,int,int);
  int wave_idx_4Boct_no3D(          int*,            int*);
  int wave_idx_4Boct_no3D(          int,int,int,int, int,int,int,int);

}


//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
inline size_t NBS_idx::prop_slv_idx(int c, int d, int cP, int dP, int ixyz, int it)
{
  return Color(c, Dirac(d, 
                        Color(cP, Dirac(dP, 
                                        ixyz + GJP.XYZnodeSites() * (it)))));
}
//--------------------------------------------------------------------------
inline int NBS_idx::prop_slv_cs_idx(int c, int d, int cP, int dP)
{
  // the following two implementations are equivalent, but K-computer prefers the latter
  //return prop_slv_idx(c, d, cP, dP, 0, 0);
  return Color(c, Dirac(d, 
                        Color(cP, Dirac(dP, 0))));
}
//--------------------------------------------------------------------------
inline size_t NBS_idx::prop_idx(int c, int d, int ixyz, int cP, int dP, int it)
{
  return Color(c, Dirac(d, 
                        ixyz + GJP.XYZnodeSites() * 
                        Color(cP, Dirac(dP, it))));
}
//--------------------------------------------------------------------------
inline int NBS_idx::prop_snk_idx_no3D(int c, int d)
{
  // the following two implementations are equivalent, but K-computer prefers the latter
  //return prop_idx(c, d, 0, 0, 0, 0);
  return Color(c, Dirac(d, 0));
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
inline int NBS_idx::Boct_blk_idx(int ixyz, int alp, // snk indices
                                 int c1, int alp1, int c2, int alp2, int c3, int alp3, // src indices
                                 bool flg_NR_limit)
{
  return ixyz + GJP.XYZnodeSites()*
    Boct_blk_idx_no3D(alp, c1, alp1, c2, alp2, c3, alp3, flg_NR_limit);
}
//--------------------------------------------------------------------------
inline int NBS_idx::Boct_blk_idx_no3D(int alp, // snk indices
                                      int c1, int alp1, int c2, int alp2, int c3, int alp3, // src indices
                                      bool flg_NR_limit)
{
  return Pauli(alp, Boct_blk_src_idx(c1, alp1, c2, alp2, c3, alp3, flg_NR_limit));
}
//--------------------------------------------------------------------------
inline int NBS_idx::Boct_blk_src_idx(int c1, int alp1, int c2, int alp2, int c3, int alp3, bool flg_NR_limit)
{
  int idx;
  if ( ! flg_NR_limit ) { idx = Color(c1, Dirac(alp1, Color(c2, Dirac(alp2, Color(c3, Dirac(alp3, 0)))))); }
  else                  { idx = Color(c1, Pauli(alp1, Color(c2, Pauli(alp2, Color(c3, Pauli(alp3, 0)))))); }
  return idx;
}
//--------------------------------------------------------------------------
inline int NBS_idx::Boct_blk_idx(int ixyz, int alp, // snk indices
                                 int c_alp)         // src indices (all color/spin combined)
{
  return ixyz + GJP.XYZnodeSites()*
    Boct_blk_idx_no3D(alp, c_alp);
}
//--------------------------------------------------------------------------
inline int NBS_idx::Boct_blk_idx_no3D(int alp,   // snk indices
                                      int c_alp) // src indices (all color/spin combined)
{
  return Pauli(alp, c_alp);
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
// final idx for wave
inline int NBS_idx::wave_idx_fin(int dof_snk, int isnk, int ixyz, int isrc)
{ 
  return isnk + dof_snk*(ixyz + GJP.XYZnodeSites() * (isrc)); 
}
//--------------------------------------------------------------------------
// itermediate idx for wave
inline int NBS_idx::wave_idx    (int dof_snk, int ixyz, int isnk, int isrc)
{ 
  return ixyz + GJP.XYZnodeSites() * (isnk + dof_snk * (isrc)); 
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
inline int NBS_idx::wave_idx_2Boct(int ixyz,
                                   int *snk_spin,
                                   int *src_spin)
{
  return wave_idx_2Boct(ixyz,
                        snk_spin[0], snk_spin[1], 
                        src_spin[0], src_spin[1]);
}
//--------------------------------------------------------------------------
inline int NBS_idx::wave_idx_2Boct(int ixyz,
                                   int alpha,  int beta,
                                   int alphaP, int betaP)
{
  return ixyz + GJP.XYZnodeSites() * 
    wave_idx_2Boct_no3D(alpha,  beta, 
                        alphaP, betaP);
}
//--------------------------------------------------------------------------
inline int NBS_idx::wave_idx_2Boct_no3D(int *snk_spin,
                                        int *src_spin)
{
  return wave_idx_2Boct_no3D(snk_spin[0], snk_spin[1], 
                             src_spin[0], src_spin[1]);
}
//--------------------------------------------------------------------------
inline int NBS_idx::wave_idx_2Boct_no3D(int alpha,  int beta,
                                        int alphaP, int betaP)
{
  return Pauli(alpha, Pauli(beta, 
                            Pauli(alphaP, Pauli(betaP, 0))));
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
inline int NBS_idx::wave_idx_3Boct(int ixyz, 
                                   int *snk_spin,
                                   int *src_spin)
{
  return wave_idx_3Boct(ixyz,
                        snk_spin[0], snk_spin[1], snk_spin[2], 
                        src_spin[0], src_spin[1], src_spin[2]);
}
//--------------------------------------------------------------------------
inline int NBS_idx::wave_idx_3Boct(int ixyz, 
                                   int alpha,  int beta,  int gamma,
                                   int alphaP, int betaP, int gammaP)
{
  return ixyz + GJP.XYZnodeSites() * 
    wave_idx_3Boct_no3D(alpha,  beta,  gamma, 
                        alphaP, betaP, gammaP);
}
//--------------------------------------------------------------------------
inline int NBS_idx::wave_idx_3Boct_no3D(int *snk_spin,
                                        int *src_spin)
{
  return wave_idx_3Boct_no3D(snk_spin[0], snk_spin[1], snk_spin[2], 
                             src_spin[0], src_spin[1], src_spin[2]);
}
//--------------------------------------------------------------------------
inline int NBS_idx::wave_idx_3Boct_no3D(int alpha,  int beta,  int gamma,
                                        int alphaP, int betaP, int gammaP)
{
  return Pauli(alpha, Pauli(beta, Pauli(gamma, 
                                        Pauli(alphaP, Pauli(betaP, Pauli(gammaP, 0))))));
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
inline int NBS_idx::wave_idx_4Boct(int ixyz,
                                   int *snk_spin,
                                   int *src_spin)
{
  return wave_idx_4Boct(ixyz,
                        snk_spin[0], snk_spin[1], snk_spin[2], snk_spin[3], 
                        src_spin[0], src_spin[1], src_spin[2], src_spin[3]);
}
//--------------------------------------------------------------------------
inline int NBS_idx::wave_idx_4Boct(int ixyz,
                                   int alpha,  int beta,  int gamma,  int delta,
                                   int alphaP, int betaP, int gammaP, int deltaP)
{
  return ixyz + GJP.XYZnodeSites() * 
    wave_idx_4Boct_no3D(alpha,  beta,  gamma,  delta, 
                        alphaP, betaP, gammaP, deltaP);
}
//--------------------------------------------------------------------------
inline int NBS_idx::wave_idx_4Boct_no3D(int *snk_spin,
                                        int *src_spin)
{
  return wave_idx_4Boct_no3D(snk_spin[0], snk_spin[1], snk_spin[2], snk_spin[3], 
                             src_spin[0], src_spin[1], src_spin[2], src_spin[3]);
}
//--------------------------------------------------------------------------
inline int NBS_idx::wave_idx_4Boct_no3D(int alpha,  int beta,  int gamma,  int delta,
                                        int alphaP, int betaP, int gammaP, int deltaP)
{
  return Pauli(alpha, Pauli(beta, Pauli(gamma, Pauli(delta, 
                                                     Pauli(alphaP, Pauli(betaP, Pauli(gammaP, Pauli(deltaP, 0))))))));
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------


HAL_END_NAMESPACE

#endif

