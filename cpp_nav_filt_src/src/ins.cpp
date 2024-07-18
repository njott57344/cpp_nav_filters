#include "ins.h"

namespace cpp_nav_filt
{

    Ins::Ins()
    {
        omega_ie_e_ << 0,0,w_e;
        Omega_ie_e_ = makeSkewSymmetic(omega_ie_e_);
        I3_.setIdentity();
    }

    Ins::~Ins()
    {

    }

    void Ins::getInitialPosition(vec_3_1& pos)
    {
        r_eb_e_ = pos;
        pos_init_ = true;
        initCheck();
    }

    void Ins::getInitialVelocity(vec_3_1& vel)
    {
        v_eb_e_ = vel;
        vel_init_ = true;
        initCheck();
    }

    void Ins::getInitialAttitude(vec_3_1& att)
    {
        if(pos_init_)
        {
            vec_3_1 init_lla;
            init_lla = ecef2llaPos(r_eb_e_);
            C_be_ = eul2EcefDCM(att,init_lla);
            att_init_ = true;
            initCheck();
        }
        else
        {
            std::cout<<"To Initiate Attitude Please Initial Position First"<<std::endl;
        }
    }

    void Ins::getInertialMeasurements(vec_3_1& f_ib_b,vec_3_1& w_ib_b,double& dt)
    {
        f_ib_b_ = f_ib_b;
        w_ib_b_ = w_ib_b;
        dt_ = dt;

        mechanizeFullState();
    }

    void Ins::sendPosition(vec_3_1& pos_out)
    {
        pos_out = r_eb_e_;
    }

    void Ins::sendVelocity(vec_3_1& vel_out)
    {
        vel_out = v_eb_e_;
    }

    void Ins::sendEcefECM(mat_3_3& dcm_out)
    {
        dcm_out = C_be_;
    }

    void Ins::sendEulerAngles(mat_3_3& eul_out)
    {
        
    }

    void Ins::initCheck()
    {
        full_state_init_ = pos_init_&&vel_init_&&att_init_;
    }

    void Ins::mechanizeFullState()
    {
        if(full_state_init_)
        {
            // prior state (before propagation)
            vec_3_1 r_eb_e_min,v_eb_e_min;
            mat_3_3 C_be_min;

            r_eb_e_min = r_eb_e_;
            v_eb_e_min = v_eb_e_;
            C_be_min   = C_be_;

            Omega_ib_b_ = makeSkewSymmetic(w_ib_b_);
            g_be_ = ecefGravity(r_eb_e_min);

            // DCM Propagation
            C_be_ = C_be_min*(I3_ + Omega_ib_b_*dt_) - Omega_ie_e_*C_be_min*dt_;
            
            // Specific Force Rotation
            f_ib_e_ = 0.5*(C_be_ + C_be_min)*f_ib_b_;

            v_eb_e_ = v_eb_e_min + (f_ib_e_ + g_be_ - 2*Omega_ie_e_*v_eb_e_min)*dt_;

            r_eb_e_ = r_eb_e_min + 0.5*dt_*(v_eb_e_ + v_eb_e_min);
        }
    }

}// end of namespace