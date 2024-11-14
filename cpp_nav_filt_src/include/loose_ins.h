#ifndef LOOSE_INS_H
#define LOOSE_INS_H

#include "cpp_nav_filt_lib.h"

namespace cpp_nav_filt
{

    class LooseIns
    {
        public:
            
            LooseIns();
            ~LooseIns();

            /* convention: 
                a get function is the class GETTING something
                a send function is the class SENDING something
            */

            void getInitialPosition(vec_3_1& pos);
            void getInitialVelocity(vec_3_1& vel);
            void getInitialAttitude(vec_3_1& att);

            // to-do: correct for bias estimate
            void getInertialMeasurements(vec_3_1& f_ib_b,vec_3_1& w_ib_b,double& dt);

            void getPVMeasurements(vec_3_1& r_eb_e_tilde,vec_3_1& v_eb_e_tilde,vec_3_1 w_ib_b,
                                   vec_3_1& r_ab_b,mat_3_3& R_pos,mat_3_3& R_vel);

            void sendPosition(vec_3_1& pos_out);
            void sendVelocity(vec_3_1& vel_out);
            void sendEcefDCM(mat_3_3& dcm_out);
            void sendEulerAngles(vec_3_1& eul_out);

        private:

            void mechanizeFullState();
            void mechanizeErrorState();

            void kalmanUpdate();
            void feedbackErrorState();

            void initCheck();
        
            // full states
            vec_3_1 r_eb_e_;
            vec_3_1 v_eb_e_;
            mat_3_3 C_be_,C_eb_;
            vec_3_1 g_be_;

            // inertial measurements
            vec_3_1 f_ib_b_,f_ib_e_;
            vec_3_1 w_ib_b_;
            mat_3_3 Omega_ib_b_;
            double dt_;

            bool pos_init_,vel_init_,att_init_;
            bool full_state_init_;

            // misc
            mat_3_3 I3_;
            mat_3_3 Omega_ie_e_;
            vec_3_1 omega_ie_e_;

            mat_3_3 normalizeDCM(mat_3_3& dcm_in);
            
        protected:

    };

}

#endif