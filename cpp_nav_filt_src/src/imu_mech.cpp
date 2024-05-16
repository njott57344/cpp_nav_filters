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
    }

    ImuMechanization::~ImuMechanization()
    {
        std::cout<<"Imu Mechanization Shutting Down"<<std::endl;
    }

    void ImuMechanization::getPos(vec_3_1& pos)
    {
        pos_posterior_ = pos;
        pos_is_init_ = true;
        checkInitStatus();
    }

    // ================== Getters ================== //
    void ImuMechanization::getVel(vec_3_1& vel)
    {
        vel_posterior_ = vel;
        vel_is_init_ = true;
        checkInitStatus();
    }

    void ImuMechanization::getAtt(vec_3_1& att,vec_3_1& init_lla)
    {
        C_bm_posterior_ = eul2EcefDCM(att,init_lla);
        att_is_init_ = true;
        checkInitStatus();
    }

    void ImuMechanization::getInertialMeasurements(vec_3_1& f_ib_b,vec_3_1& w_ib_b,vec_3_1& fbias_ib_b,vec_3_1& wbias_ib_b)
    {
        f_ib_b_ = f_ib_b - fbias_ib_b;
        w_ib_b_ = w_ib_b - wbias_ib_b;

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

    // =================== Internal Functions =============== //
    void ImuMechanization::checkInitStatus()
    {
        state_is_init_ = pos_is_init_ && vel_is_init_ && att_is_init_;
    }
}// end of namespace