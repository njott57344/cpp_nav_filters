#include "imu_mech.h"

namespace cpp_nav_filt
{
    ImuMechanization::ImuMechanization(ImuMechanizationSettings& settings_in)
    {
        settings_ = settings_in;

        state_is_init_ = false;
        pos_is_init_ = false;
        vel_is_init_ = false;
        att_is_init_ = false;

        I3.setIdentity();
        w_ie_e_ << 0.0,0.0,cpp_nav_filt::w_e;
        Omega_ie_e_ = cpp_nav_filt::makeSkewSymmetic(w_ie_e_);

    }

    ImuMechanization::~ImuMechanization()
    {
        std::cout<<"Imu Mechanization Shutting Down"<<std::endl;
    }

    // ================== Getters ================== //

    void ImuMechanization::getPos(vec_3_1& pos)
    {
        r_ib_e_posterior_ = pos;
        pos_is_init_ = true;
        checkInitStatus();
    }

    void ImuMechanization::getVel(vec_3_1& vel)
    {
        v_ib_e_posterior_ = vel;
        vel_is_init_ = true;
        checkInitStatus();
    }

    void ImuMechanization::getAtt(vec_3_1& att,vec_3_1& init_lla)
    {
        C_be_posterior_ = eul2EcefDCM(att,init_lla);
        att_is_init_ = true;
        checkInitStatus();
    }

    void ImuMechanization::getTime(double& time)
    {
        t_new_ = time;
        time_is_init_ = true;
        checkInitStatus();
    }

    void ImuMechanization::getInertialMeasurements(vec_3_1& f_ib_b,vec_3_1& w_ib_b,vec_3_1& fbias_ib_b,vec_3_1& wbias_ib_b,double& t)
    {
        f_ib_b_ = f_ib_b - fbias_ib_b;
        w_ib_b_ = w_ib_b - wbias_ib_b;
        Omega_ib_b_ = cpp_nav_filt::makeSkewSymmetic(w_ib_b_);

        t_old_ = t_new_;
        t_new_ = t;
        dt_ = t_new_ - t_old_;

        if(state_is_init_)
        {
            mechanizePVA();
        }
        else
        {
            std::cout<<"WARNING: PVA State is Not Initialized"<<std::endl;
        }

    }

    // ===================== Senders ==================== //
    void ImuMechanization::sendPos(vec_3_1& pos)
    {

    }

    void ImuMechanization::sendVel(vec_3_1& vel)
    {

    }

    void ImuMechanization::sendAtt(vec_3_1& att)
    {

    }

    void ImuMechanization::sendCbm(mat_3_3& Cbm)
    {
        
    }

    // =================== Interaction Functions ============ //
    void ImuMechanization::correctPVA(vec_3_1& del_pos,vec_3_1& del_vel,vec_3_1& del_att)
    {

    }

    // =================== Internal Functions =============== //
    void ImuMechanization::checkInitStatus()
    {
        state_is_init_ = pos_is_init_ && vel_is_init_ && att_is_init_;
    }

    void ImuMechanization::mechanizePVA()
    {
        r_ib_e_prior_ = r_ib_e_posterior_;
        v_ib_e_prior_ = v_ib_e_posterior_;
        C_be_prior_ = C_be_posterior_;

        g_ib_e_ = cpp_nav_filt::ecefGravity(r_ib_e_prior_);

        // 5.27
        C_be_posterior_ = C_be_prior_*(I3 + Omega_ib_b_*dt_) - 
                          Omega_ie_e_*C_be_prior_*dt_;

        // 5.28
        f_ib_e_ = 0.5*(C_be_posterior_+C_be_prior_)*f_ib_b_;

        // 5.36
        v_ib_e_posterior_ = v_ib_e_prior_ + (f_ib_e_ + g_ib_e_ - 
                            2*Omega_ie_e_*v_ib_e_prior_)*dt_;

        // 5.38
        r_ib_e_posterior_ = r_ib_e_prior_ + (v_ib_e_prior_+v_ib_e_posterior_)*0.5*dt_;
    }

}// end of namespace