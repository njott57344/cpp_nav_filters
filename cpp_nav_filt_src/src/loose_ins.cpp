#include "loose_ins.h"

namespace cpp_nav_filt
{

    LooseIns::LooseIns()
    {
        omega_ie_e_ << 0,0,w_e;
        Omega_ie_e_ = makeSkewSymmetic(omega_ie_e_);
        I3_.setIdentity();
    }

    LooseIns::~LooseIns()
    {

    }

    void LooseIns::getInitialPosition(vec_3_1& pos)
    {
        r_eb_e_ = pos;
        pos_init_ = true;
        initCheck();
    }

    void LooseIns::getInitialVelocity(vec_3_1& vel)
    {
        v_eb_e_ = vel;
        vel_init_ = true;
        initCheck();
    }

    void LooseIns::getInitialAttitude(vec_3_1& att)
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

    void LooseIns::getInertialMeasurements(vec_3_1& f_ib_b,vec_3_1& w_ib_b,double& dt)
    {
        f_ib_b_ = f_ib_b;
        w_ib_b_ = w_ib_b;
        dt_ = dt;

        mechanizeFullState();
    }

    void LooseIns::sendPosition(vec_3_1& pos_out)
    {
        pos_out = r_eb_e_;
    }

    void LooseIns::sendVelocity(vec_3_1& vel_out)
    {
        vel_out = v_eb_e_;
    }

    void LooseIns::sendEcefDCM(mat_3_3& dcm_out)
    {
        dcm_out = C_be_;
    }

    void LooseIns::sendEulerAngles(vec_3_1& eul_out)
    {
        
    }

    void LooseIns::initCheck()
    {
        full_state_init_ = pos_init_&&vel_init_&&att_init_;
    }

    void LooseIns::mechanizeFullState()
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
            
            C_be_ = normalizeDCM(C_be_);

            // Specific Force Rotation
            f_ib_e_ = 0.5*(C_be_ + C_be_min)*f_ib_b_;

            v_eb_e_ = v_eb_e_min + (f_ib_e_ + g_be_ - 2*Omega_ie_e_*v_eb_e_min)*dt_;

            r_eb_e_ = r_eb_e_min + 0.5*dt_*(v_eb_e_ + v_eb_e_min);
        }
    }

    mat_3_3 LooseIns::normalizeDCM(mat_3_3& dcm_in)
    {
        mat_3_3 normalized_dcm;

        vec_3_1 x,y,z;
        vec_3_1 x_ort,y_ort,z_ort;
        vec_3_1 x_new,y_new,z_new;
        double error;

        x = dcm_in.block(0,0,3,1);
        y = dcm_in.block(0,1,3,1);
        z = dcm_in.block(0,2,3,1);

        error = x.dot(y);
    
        x_ort = x - (0.5*error)*y;
        y_ort = y - (0.5*error)*x;
        z_ort = x_ort.cross(y_ort);

        x_new = 0.5*(3-x_ort.dot(x_ort))*x_ort;
        y_new = 0.5*(3-y_ort.dot(y_ort))*y_ort;
        z_new = 0.5*(3-z_ort.dot(z_ort))*z_ort;
        
        normalized_dcm.block(0,0,3,1) = x_new;
        normalized_dcm.block(0,1,3,1) = y_new;
        normalized_dcm.block(0,2,3,1) = z_new;
        
        return normalized_dcm;
    }

    
}// end of namespace