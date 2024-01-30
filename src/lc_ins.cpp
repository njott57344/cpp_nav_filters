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

    // =============== Getters ============ //
    void LooselyCoupledIns::getPositionSoln(vec_3_1& pos)
    {
        pos = x_hat_.block<3,1>(3,0);
    }

    void LooselyCoupledIns::getImuMeasurements(vec_3_1& f,vec_3_1& ar,double& t)
    {
        wb_b_ = ar;
        fb_b_ = f;

        dt_ = t - time_;
        time_ = t;
    }

    // =============== Init =============== //
    void LooselyCoupledIns::setInitialPosState(vec_3_1& pos_init,mat_3_3& pos_P)
    {
        x_hat_.block<3,1>(3,0) = pos_init;
        P_hat_.block<3,3>(3,3) = pos_P;
        pos_init_ = true;
    }

    void LooselyCoupledIns::setInitialVelState(vec_3_1& vel_init,mat_3_3& vel_P)
    {
        x_hat_.block<3,1>(0,0) = vel_init;
        P_hat_.block<3,3>(0,0) = vel_P;
        vel_init_ = true;
    }

    void LooselyCoupledIns::setInitialAttState(vec_3_1& att_init,mat_3_3& att_P)
    {
        x_hat_.block<3,1>(6,0) = att_init;
        P_hat_.block<3,3>(6,6) = att_P;
        att_init_ = true;
    }

    void LooselyCoupledIns::setInitialTime(double& init_time)
    {
        time_ = init_time;
        time_init_ = true;
    }

    void LooselyCoupledIns::setInitialBgState(vec_3_1& bg_init,mat_3_3& bg_P)
    {
        x_hat_.block<3,1>(9,0) = bg_init;
        P_hat_.block<3,3>(9,9) = bg_P;
        bg_init_ = true;
    }

    void LooselyCoupledIns::setInitialBaState(vec_3_1& ba_init,mat_3_3& ba_P)
    {
        x_hat_.block<3,1>(12,0) = ba_init;
        P_hat_.block<3,3>(12,12) = ba_P;
        ba_init_ = true;
    }

    // =============== Loose INS ============== //
    void LooselyCoupledIns::mechanizeSolution()
    {
        // propagate full state
        if(filt_init_)
        {

        }
    }

    void LooselyCoupledIns::checkInitStatus()
    {
        filt_init_ = pos_init_&&vel_init_&&att_init_&&time_init_&&bg_init_&&ba_init_;
    }

} // end of namespace
