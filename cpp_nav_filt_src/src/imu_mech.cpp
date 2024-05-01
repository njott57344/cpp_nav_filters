#include "imu_mech.h"

namespace cpp_nav_filt
{
    ImuMechanization::ImuMechanization(ImuMechanizationSettings& settings_in)
    {
        settings_ = settings_in;
    }

    ImuMechanization::~ImuMechanization()
    {
        std::cout<<"Imu Mechanization Shutting Down"<<std::endl;
    }

    void ImuMechanization::getPos(vec_3_1& pos)
    {

    }

    void ImuMechanization::getVel(vec_3_1& vel)
    {

    }

    void ImuMechanization::getAtt(vec_3_1& att)
    {
        
    }

    void ImuMechanization::getInertialMeasurements(vec_3_1& f_ib_b,vec_3_1& w_ib_b)
    {

    }

    void ImuMechanization::sendPos(vec_3_1& pos)
    {

    }

    void ImuMechanization::sendVel(vec_3_1& vel)
    {

    }

    void sendAtt(vec_3_1& att)
    {

    }

    void sendCbm(mat_3_3& Cbm)
    {
        
    }
}// end of namespace