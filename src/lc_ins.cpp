#include "lc_ins/lc_ins.h"

namespace cpp_nav_filt
{
    LooselyCoupledIns::LooselyCoupledIns(LooselyCoupledInsSettings& settings_in)
    {
        std::cout<<"Starting Loosely Coupled INS!"<<std::endl;

        ins_settings = settings_in;

        // zeroing state estimates
        x_hat_.setZero();
        P_hat_.setZero();
        dx_hat_.setZero(); // zeroize error state
        Phi_.setZero();
        Q_.setZero();

        filt_init_ = false;
        pos_init_  = false;
        vel_init_  = false;
        att_init_  = false;
        bg_init_   = false;
        ba_init_   = false;
        time_init_ = false;

        we_i_<<0,0,cpp_nav_filt::w_e;

        I_3_.setIdentity();
        I_15_.setIdentity();

        Omega_e_ = cpp_nav_filt::makeSkewSymmetic(we_i_);
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

    void LooselyCoupledIns::setBgSol(vec_3_1& bg_hat)
    {
        bg_hat = x_hat_.block<3,1>(9,0);
    }

    void LooselyCoupledIns::setBaSol(vec_3_1& ba_hat)
    {
        ba_hat = x_hat_.block<3,1>(12,0);
    }

    void LooselyCoupledIns::setErrorCoverEst(mat_15_15& P_hat)
    {
        P_hat = P_hat_;
    }

    void LooselyCoupledIns::getImuMeasurements(vec_3_1& f,vec_3_1& ar,double& t)
    {
        wb_b_ = ar;
        fb_b_ = f;

        dt_ = t - time_;
        
        time_ = t;
        
        if(dt_>0)
        {
            correctImuMeasurements(); // perform bias compensation
            mechanizeSolution(); // mechanize full state PVA solution 

            generateErrorStateSTM(); // generate error state state transition matrix
            generateProcessCovarMat(); // generate current process covariance matrix

            propagateErrorState();

            setFullStateEstimate(plus_pos_,plus_vel_,plus_att_);
        }
    }

    void LooselyCoupledIns::getGnssMeasurements(vec_3_1& gnss_measurement,const int& meas_type,mat_3_3& R)
    {
        Eigen::MatrixXd S; // innovation Covariance
        Eigen::MatrixXd K; // kalman gain

        gnss_meas_ = gnss_measurement;
        meas_type_ = meas_type;

        estimateGnssMeasurement(y_hat_,H_); // gives full state measurement estimate and jacobian
        
        innov_ = gnss_meas_ - y_hat_; // full state innovation
        del_innov_ = innov_ - H_*dx_hat_; // error state innovation

        S = H_*P_hat_*H_.transpose() + R;
        K = ( (S.transpose()).fullPivLu().solve((P_hat_ * H_.transpose()).transpose()) ).transpose();
        dx_hat_ = dx_hat_ + K*del_innov_;

        P_hat_ = (I_15_ - K*H_)*P_hat_*(I_15_ - K*H_).transpose() + K*R*K.transpose(); // joseph form

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

            C_n_b_minus_ = cpp_nav_filt::eul2Rotm(minus_att_);
            Omega_b_ = cpp_nav_filt::makeSkewSymmetic(wb_b_);
            gamma_b_n_ = cpp_nav_filt::somiglianaGravityModel(minus_pos_);

            C_b_n_minus_ = C_n_b_minus_.transpose();

            // state propagation
            C_b_n_plus_ = C_b_n_minus_*(I_3_ + Omega_b_*dt_); // - Omega_e_*C_b_n_minus_*dt_; // Attitude Update
            fb_n_ = 0.5*(C_b_n_minus_ + C_b_n_plus_)*fb_b_; // rotating specific force into nav frame
            
            gamma_b_n_ << 0,0,9.81; // this is wrong look into a better model

            plus_vel_ = minus_vel_ + (fb_n_ + gamma_b_n_)*dt_; // velocity update
            plus_pos_ = minus_pos_ + plus_vel_*dt_; // position update
                // note: pos update assumes velocity varies linearly over integration period [groves 175]

            C_n_b_plus_ = C_b_n_plus_.transpose();
            plus_att_ = cpp_nav_filt::rotm2Eul(C_n_b_plus_);
        }
    }

    void LooselyCoupledIns::generateErrorStateSTM()
    {
        /*
            See Groves p. 582-584 (ECEF error state INS state propogation)
        */

        vec_3_1 current_pos,lla,current_att;
        mat_3_3 Cbn,Cnb;

        setPosSol(current_pos);
        setAttSol(current_att);
        Cnb = cpp_nav_filt::eul2Rotm(current_att);
        Cbn = Cnb.transpose();
        fb_n_ = Cbn*fb_b_; // nav frame specific fornces

        gamma_b_n_ = cpp_nav_filt::somiglianaGravityModel(current_pos);
        lla = cpp_nav_filt::ecef2llaPos(current_pos);
        lla.block<2,1>(0,0) = lla.block<2,1>(0,0)*cpp_nav_filt::D2R;
        geocentric_radius_ = cpp_nav_filt::geocentricRadius(lla[0]);

        F23_ = (-2*gamma_b_n_*(current_pos.transpose()))/(geocentric_radius_*current_pos.norm());
        F21_ = cpp_nav_filt::makeSkewSymmetic(fb_n_);
        F21_ = -F21_;

        // velocity
        Phi_.block<3,3>(0,0) = I_3_ - 2*dt_*Omega_e_; // del V w/ del V
        Phi_.block<3,3>(0,3) = F23_*dt_;
        Phi_.block<3,3>(0,6) = F21_*dt_;
        Phi_.block<3,3>(0,12) = Cbn*dt_;
        
        // Position
        Phi_.block<3,3>(3,0) = I_3_*dt_;
        Phi_.block<3,3>(3,3) = I_3_;
        
        // Attitude
        Phi_.block<3,3>(6,6) = I_3_ - Omega_e_*dt_;
        Phi_.block<3,3>(6,9) = Cbn*dt_;

        //Gyro Bias
        // assume is constant for now (may change to first order GMP in the future)
        Phi_.block<3,3>(9,9) = I_3_;

        //Accel Bias
        // assume is constant for now (may change to first order GMP in the future)
        Phi_.block<3,3>(12,12) = I_3_;
    }

    void LooselyCoupledIns::propagateErrorState()
    {
        dx_hat_ = Phi_*dx_hat_; // state update equation
        P_hat_ = Phi_*P_hat_*(Phi_.transpose()) + Q_;
    }

    void LooselyCoupledIns::generateProcessCovarMat()
    {
        // Groves p. 592

        // this is an approximation, assumes dt less than or equal to 0.2 s
        // I will fix this when Loosely Coupled INS is working

        Q_.block<3,3>(0,0)   = ins_settings.psd_accel_noise*I_3_;
        Q_.block<3,3>(6,6)   = ins_settings.psd_gyro_noise*I_3_;
        Q_.block<3,3>(9,9)   = ins_settings.psd_gyro_bias*I_3_;
        Q_.block<3,3>(12,12) = ins_settings.psd_accel_bias*I_3_;

        Q_ = Q_*dt_;
    }

    void LooselyCoupledIns::correctImuMeasurements()
    {
        setBgSol(bg_hat_);
        setBaSol(ba_hat_);

        wb_b_ = wb_b_ - bg_hat_; // subtract bias estimate from gyroscopes
        fb_b_ = fb_b_ - ba_hat_; // subtract bias estimate from accelerometers
    }

    void LooselyCoupledIns::estimateGnssMeasurement(vec_3_1& y_hat,mat_3_15& H)
    {
        vec_3_1 state,att;
        
        setAttSol(att);

        if(meas_type_ == POS_TYPE)
        {
            setPosSol(state);
        }
        else if(meas_type_ == VEL_TYPE)
        {
            setVelSol(state);
        }
        
        measurementModel(state,att,H);
        rigidBodyTransform(state,att,y_hat); // gets current estimate of measured quantity
    }

    void LooselyCoupledIns::rigidBodyTransform(vec_3_1& X,vec_3_1& att,vec_3_1& Y)
    {
        // see groves p. 598 eq 14.102
        mat_3_3 Cbn,Cnb; // current estimate of DCM
        vec_3_1 lb_b; // current estimate of lever arm (body frame)

        Cnb = cpp_nav_filt::eul2Rotm(att);
        Cbn = Cnb.transpose();

        lb_b = ins_settings.lever_arm;

        if(meas_type_ == POS_TYPE)
        {
            Y = X + Cbn*lb_b;
        }
        else if(meas_type_ == VEL_TYPE)
        {
            mat_3_3 skew_ang_rate;
            skew_ang_rate = cpp_nav_filt::makeSkewSymmetic(wb_b_);
            Y = X + Cbn*(skew_ang_rate*lb_b) - Omega_e_*Cbn*lb_b;
        }
    }

    void LooselyCoupledIns::measurementModel(vec_3_1& X,vec_3_1& att_sol,mat_3_15& H)
    {
        H.setZero();

        // see groves p. 600 eq 14.111-14.112
        mat_3_3 Cnb,Cbn;
        vec_3_1 lb_b = ins_settings.lever_arm;
        vec_3_1 lb_n = Cbn*lb_b;
        Cnb = cpp_nav_filt::eul2Rotm(att_sol);
        Cbn = Cnb.transpose();

        if(meas_type_ == POS_TYPE)
        {
            mat_3_3 Hr1;
            Hr1 = cpp_nav_filt::makeSkewSymmetic(lb_n);
            
            H.block<3,3>(0,6) = Hr1; // attitude 
            H.block<3,3>(0,3) = -1*I_3_; // position 
        }
        else if(meas_type_ == VEL_TYPE)
        {
            mat_3_3 Hv1,Hv2,skew_ang_rate,hv2;
            vec_3_1 hv1;

            skew_ang_rate = cpp_nav_filt::makeSkewSymmetic(wb_b_);
            hv1 = Cbn*(skew_ang_rate*lb_b) - Omega_e_*Cbn*lb_b;
            Hv1 = cpp_nav_filt::makeSkewSymmetic(hv1);
            
            hv2 = cpp_nav_filt::makeSkewSymmetic(lb_b);
            Hv2 = Cbn*hv2;

            H.block<3,3>(0,6) = Hv1; // attitude
            H.block<3,3>(0,0) = -1*I_3_; // velocity
            H.block<3,3>(0,9) = Hv2;
        }
        else
        {
            H = NAN*H;
        }
    }

    void LooselyCoupledIns::errorStateCorrection()
    {
        // see groves p. 564 eq 14.7-14.9
        vec_3_1 pos,vel,att,att_err;
        mat_3_3 Cnb,Cbn,Cbn_plus,skew_att_err;

        setPosSol(pos);
        setVelSol(vel);
        setAttSol(att);

        pos = pos - dx_hat_.block<3,1>(3,0); // position correction
        vel = vel - dx_hat_.block<3,1>(0,0); // velocity correction

        // attitude correction
        Cnb = cpp_nav_filt::eul2Rotm(att);
        Cbn = Cnb.transpose();

        att_err = dx_hat_.block<3,1>(6,0);
        skew_att_err = cpp_nav_filt::makeSkewSymmetic(att_err);

        Cbn_plus = (I_3_ - skew_att_err)*Cbn;
        Cnb = Cbn_plus.transpose();

        att = cpp_nav_filt::rotm2Eul(Cnb);

        dx_hat_.block<9,1>(0,0).setZero();

        x_hat_.block<3,1>(9,0) = dx_hat_.block<3,1>(9,0);
        x_hat_.block<3,1>(12,0) = dx_hat_.block<3,1>(12,0);
        
        setFullStateEstimate(pos,vel,att);
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
