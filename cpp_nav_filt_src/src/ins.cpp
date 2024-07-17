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
    }

    void Ins::initCheck()
    {
        full_state_init_ = pos_init_&&vel_init_&&att_init_;
    }

    void Ins::mechanizeFullState()
    {
        if(full_state_init_)
        {
            Omega_ib_b_ = makeSkewSymmetic(w_ib_b_);
            g_be_ = ecefGravity(r_eb_e_);

            C_be_ = C_be_*(I3_ + Omega_ib_b_*dt_) - Omega_ie_e_*C_be_*dt_;
            v_
        }
    }

}// end of namespace