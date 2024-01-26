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

    void LooselyCoupledIns::getPositionSoln(vec_3_1& pos)
    {
        pos = x_hat_.block<3,1>(3,0);
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

    void LooselyCoupledIns::somiglianaGravityModel(vec_3_1& gamma_b_n)
    {   
        vec_3_1 pos,inner;
        double term_x,term_y,term_z,common_term;
        getPositionSoln(pos);

        common_term = 5*pow((pos[2]/pos.norm()),2);
        term_x = 1-common_term*pos[0];
        term_y = 1-common_term*pos[1];
        term_z = 1-common_term*pos[2];

        inner<<term_x,term_y,term_z;

        inner = (1.5*cpp_nav_filt::J2*pow(cpp_nav_filt::Ro,2)/(pow(pos.norm(),2)))*inner;
        inner = pos + inner;
        gamma_b_n = -(cpp_nav_filt::mu_g/pow(pos.norm(),3))*inner;
    }

} // end of namespace
