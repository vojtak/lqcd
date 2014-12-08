#include <ctype.h>

#include <util/halqcd/HAL_config.h>
#include <util/halqcd/Hadron.h>

using namespace hal;


void hal_run_V(int *Nodes, int *NodeSites, int *NodeCoor,
             Float *prop_ud,  Float *prop_s,
             vector<int> iT_src_lst, int iX_src, int iY_src, int iZ_src,
             string conf_name)
{
  //-------------------------------------------------------------------------
  // initialization
  //-------------------------------------------------------------------------

  std::new_handler old_handler = std::set_new_handler(out_of_memory);

  GJP.initialize(Nodes, NodeSites, NodeCoor);
 
  VRB.Level          (VERBOSE_FUNC_LEVEL);
  VRB.ActivateLevel  (VERBOSE_FLOW_LEVEL);

  int Tsites = GJP.Tsites();

  //-------------------------------------------------------------------------
  // convert the indices of propagators
  //-------------------------------------------------------------------------

  //Prop_Hadron::conv_idx(prop_ud, prop_s, From_CPS);

  string corr_wdir;

  corr_wdir = "2p_correlators";  // name of the folder to save 2p correlators

    {
      Hadron hadron(prop_ud, prop_s);
      hadron.set_src_params(iT_src_lst, iX_src, iY_src, iZ_src);
      hadron.set_wfile_params(corr_wdir, conf_name);

      hadron.run_2pt();
    }

 // idx_snk

  //-------------------------------------------------------------------------
  // finalization
  //-------------------------------------------------------------------------

  std::set_new_handler(old_handler);
}


