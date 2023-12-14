#include "common/common.h"

namespace cpp_nav_filt
{
    Common::Common()
    {
        sv_ephem.setZero();
    }

    Common::~Common()
    {

    }    

    void Common::passSvEphem(vec_1_27& ephem_in,const int& sv_in)
    {
        sv_ephem.block<1,27>(sv_in-1,0) = ephem_in;
    }
}
