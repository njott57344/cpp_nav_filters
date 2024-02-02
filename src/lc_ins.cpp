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

        we_i_<<0,0,cpp_nav_filt::w_e;

        I_3_.setZero();
    }

    LooselyCoupledIns::~LooselyCoupledIns()
    {
        std::cout<<"Loosely Coupled INS Shuting Down!"<<std::endl;
    }

    // =============== Getters ============ //
    void LooselyCoupledIns::setPosSol(vec_3_1& pos)
    {
        pos = x_hat_.block<3,1>(3,0);
    }

    void LooselyCoupledIns::setAttSol(vec_3_1& att)
    {
        att = x_hat_.block<3,1>(6,0);
    }

    void LooselyCoupledIns::setVelSol(vec_3_1& vel)
    {
        vel = x_hat_.block<3,1>(0,0);
    }

    void LooselyCoupledIns::getImuMeasurements(vec_3_1& f,vec_3_1& ar,double& t)
    {
        wb_b_ = ar;
        fb_b_ = f;

        dt_ = t - time_;
        time_ = t;
    }

    // =============== Setters =============== //
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

    void LooselyCoupledIns::setFullStateEstimate(vec_3_1& pos,vec_3_1& vel,vec_3_1& att)
    {
        x_hat_.block<3,1>(0,0) = vel;
        x_hat_.block<3,1>(3,0) = pos;
        x_hat_.block<3,1>(6,0) = att;
    }

    void LooselyCoupledIns::setCommonClass(Common& common)
    {
        common_ = common;
        common_.makeSkewSymmetic(we_i_,Omega_e_);
    }

    // =============== Loose INS ============== //
    void LooselyCoupledIns::mechanizeSolution()
    {
        // propagate full state
        if(filt_init_)
        {
            /*
                See groves p.173-175 for equations (ECEF full state propagation)
            */

            setPosSol(minus_pos_);
            setAttSol(minus_vel_);
            setVelSol(minus_att_);

            common_.eul2Rotm(minus_att_,C_n_b_); // rotation matrix based on current solution of attitude
            common_.makeSkewSymmetic(wb_b_,Omega_b_); // skew symmetric of body frame angular rates
            common_.somiglianaGravityModel(minus_pos_,gamma_b_n_); // gravity from somigliana model

            C_b_n_ = C_n_b_.transpose();

            // state propagation
            C_b_n_ = C_b_n_*(I_3_ + Omega_b_*dt_) - Omega_e_*C_b_n_*dt_; // Attitude Update
            fb_n_ = C_b_n_*fb_b_; // rotating specific force into nav frame
            plus_vel_ = minus_vel_ + ((fb_n_ + gamma_b_n_) - 2*Omega_e_*minus_vel_)*dt_; // velocity update
            plus_pos_ = minus_pos_ + (minus_vel_ + plus_vel_)*0.5*dt_; // position update
                // note: pos update assumes velocity varies linearly over integration period [groves 175]

            C_n_b_ = C_b_n_.transpose();

            common_.rotm2Eul(C_n_b_,plus_att_);

            setFullStateEstimate(plus_pos_,plus_vel_,plus_att_);
        }
    }

    void LooselyCoupledIns::checkInitStatus()
    {
        filt_init_ = pos_init_&&vel_init_&&att_init_&&time_init_&&bg_init_&&ba_init_;
    }

} // end of namespace
