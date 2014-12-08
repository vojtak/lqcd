

#include "X_misc_V.h"

void corr_print(std::valarray<dcomplex>& corr, char* wfile, int iT_src_pos)
{
        int Tsites=CommonParameters::Lt();
          string ofname(wfile);
          std::ofstream fout(ofname.c_str());
          fout.setf(std::ios::scientific);
          fout.precision(16);
          for(int it = 0; it < CommonParameters::Lt(); it++){
            int iT2 = (it + iT_src_pos + 100*Tsites) % Tsites;
            char line[1000];
            snprintf(line, sizeof(line), "%4d\t%1.16e %1.16e\n",
                 it, real(corr[iT2]),imag(corr[iT2]));
            fout << line;
          }
          fout.close();

}

