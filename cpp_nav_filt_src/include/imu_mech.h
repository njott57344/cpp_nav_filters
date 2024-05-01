#ifndef IMU_MECH_H
#define IMU_MECH_H

#include "cpp_nav_filt_lib.h"

namespace cpp_nav_filt
{
    typedef struct
    {
        std::string mechanization_frame;
    } ImuMechanizationSettings;

    class ImuMechanization
    {
        public:
            
            ImuMechanization(ImuMechanizationSettings& settings_in);
            ~ImuMechanization();

            // Getters (to class) and senders (from class)
            void getPos(vec_3_1& pos); // set position state (for init or hard reset)
            void getVel(vec_3_1& vel); // set velocity state (for init or hard reset)
            void getAtt(vec_3_1& att); // set local tangent frame attitude state (for init or hard reset)
            void getInertialMeasurements(vec_3_1& f_ib_b,vec_3_1& w_ib_b); // receive current inertial measurements

            void sendPos(vec_3_1& pos); // send pos to outside class
            void sendVel(vec_3_1& vel); // send vel to outside class
            void sendAtt(vec_3_1& att); // send att to outside class
            void sendCbm(mat_3_3& Cbm); // send Cbm (rotation from body to mech frame) to outside class


        private:
            
            ImuMechanizationSettings settings_;

            vec_9_1 xhat_prior_; // Estimate of PVA prior to propagation
            vec_9_1 xhat_posterior_; // Estimate of PVA posterior to propagation

            vec_3_1 f_ib_b_; // specific force of body wrt inertial frame in body frame
            vec_3_1 f_ib_m_; // specific force of body wrt inertial frame in mechanization frame 

            vec_3_1 w_ib_b_; // angular rates of body wrt inertial frame in body frame

            bool state_is_init_;
            bool pos_is_init_;
            bool vel_is_init_;
            bool att_is_init_;
            
            void mechanizePVA(); // function that runs mechanization

        protected:
    }; // end of class

} // end of namespace
#endif