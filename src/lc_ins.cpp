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
        att_init_  = false;
        bg_init_   = false;
        ba_init_   = false;
        time_init_ = false;

        we_i_<<0,0,cpp_nav_filt::w_e;

        I_3_.setIdentity();
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
        
        if(dt_>0)
        {
            mechanizeSolution();
        }
    }

    // =============== Setters =============== //
    void LooselyCoupledIns::setInitialPosState(vec_3_1& pos_init,mat_3_3& pos_P)
    {
        x_hat_.block<3,1>(3,0) = pos_init;
        std::cout<<"pos: "<<x_hat_.block<3,1>(3,0).transpose()<<std::endl;
        P_hat_.block<3,3>(3,3) = pos_P;
        pos_init_ = true;
        checkInitStatus();
    }

    void LooselyCoupledIns::setInitialVelState(vec_3_1& vel_init,mat_3_3& vel_P)
    {
        x_hat_.block<3,1>(0,0) = vel_init;
        std::cout<<"vel: "<<x_hat_.block<3,1>(0,0).transpose()<<std::endl;
        P_hat_.block<3,3>(0,0) = vel_P;
        vel_init_ = true;
        checkInitStatus();
    }

    void LooselyCoupledIns::setInitialAttState(vec_3_1& att_init,mat_3_3& att_P)
    {
        x_hat_.block<3,1>(6,0) = att_init;
        std::cout<<"att: "<<x_hat_.block<3,1>(6,0).transpose()<<std::endl;
        P_hat_.block<3,3>(6,6) = att_P;
        att_init_ = true;
        checkInitStatus();
    }

    void LooselyCoupledIns::setInitialTime(double& init_time)
    {
        time_ = init_time;
        time_init_ = true;
        checkInitStatus();
    }

    void LooselyCoupledIns::setInitialBgState(vec_3_1& bg_init,mat_3_3& bg_P)
    {
        x_hat_.block<3,1>(9,0) = bg_init;
        P_hat_.block<3,3>(9,9) = bg_P;
        bg_init_ = true;
        checkInitStatus();
    }

    void LooselyCoupledIns::setInitialBaState(vec_3_1& ba_init,mat_3_3& ba_P)
    {
        x_hat_.block<3,1>(12,0) = ba_init;
        P_hat_.block<3,3>(12,12) = ba_P;
        ba_init_ = true;
        checkInitStatus();
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

            setPosSol(minus_pos_); // get position state before propagation
            setAttSol(minus_att_); // get attitude state before propagation
            setVelSol(minus_vel_); // get velocity state before propagation

            common_.eul2Rotm(minus_att_,C_n_b_minus_); // rotation matrix based on current solution of attitude                      
            common_.makeSkewSymmetic(wb_b_,Omega_b_); // skew symmetric of body frame angular rates
            common_.somiglianaGravityModel(minus_pos_,gamma_b_n_); // gravity from somigliana model

            C_b_n_minus_ = C_n_b_minus_.transpose();

            // state propagation
            C_b_n_plus_ = C_b_n_minus_*(I_3_ + Omega_b_*dt_); // - Omega_e_*C_b_n_minus_*dt_; // Attitude Update
            fb_n_ = 0.5*(C_b_n_minus_ + C_b_n_plus_)*fb_b_; // rotating specific force into nav frame
            
            gamma_b_n_ << 0,0,9.81;

            plus_vel_ = minus_vel_ + (fb_n_ + gamma_b_n_)*dt_; // velocity update
            plus_pos_ = minus_pos_ + plus_vel_*dt_; // position update
                // note: pos update assumes velocity varies linearly over integration period [groves 175]

            C_n_b_plus_ = C_b_n_plus_.transpose();

            common_.rotm2Eul(C_n_b_plus_,plus_att_);

            setFullStateEstimate(plus_pos_,plus_vel_,plus_att_);
        }
    }

    void LooselyCoupledIns::checkInitStatus()
    {
        filt_init_ = pos_init_&&vel_init_&&att_init_&&time_init_&&bg_init_&&ba_init_;
        if(filt_init_)
        {
            std::cout<<"Loosely Coupled INS has been initialized!"<<std::endl;
        }
    }

} // end of namespace
