#ifndef LC_INS_H
#define LC_INS_H

#include "common/common.h"

/*
    Inputs:
    1) Body Frame specific force, angular rates, time associated with imu observables
    2) GPS Position Solution and machine time associated with GPS pos sol
        need to pass lever arm FROM antenna TO imu 
    3) GPS Velocity Solution and machine time associated with GPS vel sol
        need to pass lever arm FROM antenna TO imu

    Outputs:
    1) Output PVA solution in ECEF at any point on body 
        need to pass output function lever arm FROM imu TO position on body
    2) Covariance associated with PVA state (P)

    Notes:
    1) Navigator is in ECEF
    2) Body frame lever arms are:
        x -> forward
        y -> left of vehicle
        z -> completing right hand triad through floor
    
*/


namespace cpp_nav_filt
{

    class LooselyCoupledIns
    {
        public:

            LooselyCoupledIns();
            ~LooselyCoupledIns();

            void getImuMeasurements(vec_3_1& f,vec_3_1& ar,double& t);

            void setInitialPosState(vec_3_1& pos_init,mat_3_3& pos_P);
            void setInitialVelState(vec_3_1& vel_init,mat_3_3& vel_P);
            void setInitialAttState(vec_3_1& att_init,mat_3_3& att_P);
            void setInitialTime(double& init_time);
            void setInitialBgState(vec_3_1& bg_init,mat_3_3& bg_P);
            void setInitialBaState(vec_3_1& ba_init,mat_3_3& ba_P);

            void getPositionSoln(vec_3_1& pos);
            
        private:
            
            bool filt_init_; // boolean to check if filter has initial full state estimate
            bool pos_init_;
            bool vel_init_;
            bool att_init_;
            bool time_init_;
            bool bg_init_;
            bool ba_init_;

            /*
                X => [dx,dy,dx,x,y,z,roll,pitch,yaw,bgx,bgy,bgz,bax,bay,baz]^T
            */

            vec_15_1 x_hat_; // full state mean solution
            vec_15_1 dx_hat_; // estimate of error in full state
            mat_15_15 P_hat_; // full state error covariance
            
            vec_3_1 fb_b_; // body frame specific force in body frame
            vec_3_1 fb_n_; // body frame specific force in nav frame
            vec_3_1 gamma_b_n_; // gravity acting on body in nav frame

            vec_3_1 wb_b_; // body frame angular rate in body frame
            vec_3_1 wb_n_; // body frame angular rate in nav frame
            mat_3_3 Omegab_n; // skew symetric of angular rate in nav frame

            mat_3_3 C_b_n_; // rotation matrix from body to nav frame (based on attitude estimate)

            vec_3_1 pos_; // ECEF GPS position solution
            vec_3_1 vel_; // ECEF GPS velocity solution

            double time_;
            double dt_;

            void mechanizeSolution();
            void checkInitStatus();
            
        protected:


    }; // end of class
} // end of namespace

#endif