#include "lc_ins/lc_ins.h"

namespace cpp_nav_filt
{
    LooselyCoupledIns::LooselyCoupledIns()
    {
        std::cout<<"Starting Loosely Coupled INS!"<<std::endl;

        // zeroing state estimates
        x_hat_.setZero();
        P_hat_.setZero();

        filt_init_ = false;
        pos_init_  = false;
        vel_init_  = false;
    }

    LooselyCoupledIns::~LooselyCoupledIns()
    {
        std::cout<<"Loosely Coupled INS Shuting Down!"<<std::endl;
    }

    void LooselyCoupledIns::setInitialPosState(vec_3_1& pos_init,mat_3_3& pos_P)
    {
        x_hat_.block<3,1>(3,0) = pos_init;
        P_hat_.block<3,3>(3,3) = pos_P;
        pos_init_ = true;
        filt_init_ = pos_init_&&vel_init_; // true if pos and vel have been initialized
    }

    void LooselyCoupledIns::setInitialVelState(vec_3_1& vel_init,mat_3_3& vel_P)
    {
        x_hat_.block<3,1>(0,0) = vel_init;
        P_hat_.block<3,3>(0,0) = vel_P;
        vel_init_ = true;
        filt_init_ = pos_init_&&vel_init_; // true if pos and vel have been initialized
    }

    void LooselyCoupledIns::mechanizeSolution()
    {

    }

} // end of namespace
