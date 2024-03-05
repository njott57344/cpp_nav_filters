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
    3) IMU assumed to be aligned with body frame 
        if above is not true, user needs to apply constant rotation before using ins filter
        x -> forward
        y -> left of vehicle
        z -> completing right hand traid through floor   
*/


namespace cpp_nav_filt
{
    typedef struct
    {
        double psd_gyro_noise;
        double psd_gyro_bias;
        double psd_accel_noise;
        double psd_accel_bias;

        vec_3_1 lever_arm;

    }LooselyCoupledInsSettings;
    
    class LooselyCoupledIns
    {
        public:



            LooselyCoupledIns(LooselyCoupledInsSettings& settings_in);
            ~LooselyCoupledIns();

            void getImuMeasurements(vec_3_1& f,vec_3_1& ar,double& t);

            void setInitialPosState(vec_3_1& pos_init,mat_3_3& pos_P);
            void setInitialVelState(vec_3_1& vel_init,mat_3_3& vel_P);
            void setInitialAttState(vec_3_1& att_init,mat_3_3& att_P);
            void setInitialTime(double& init_time);
            void setInitialBgState(vec_3_1& bg_init,mat_3_3& bg_P);
            void setInitialBaState(vec_3_1& ba_init,mat_3_3& ba_P);

            void setCommonClass(Common& common);
            void setPosSol(vec_3_1& pos);
            void setAttSol(vec_3_1& att);
            void setVelSol(vec_3_1& vel);
            void setBgSol(vec_3_1& bg_hat);
            void setBaSol(vec_3_1& ba_hat);

        private:
            
            LooselyCoupledInsSettings ins_settings;

            Common common_;

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

            mat_15_15 Phi_; // error state STM
            mat_15_15 Q_; // discrete process covariance

            vec_3_1 fb_b_; // body frame specific force in body frame
            vec_3_1 fb_n_; // body frame specific force in nav frame
            vec_3_1 gamma_b_n_; // gravity acting on body in nav frame

            vec_3_1 plus_pos_;
            vec_3_1 plus_vel_;
            vec_3_1 plus_att_;

            vec_3_1 minus_pos_;
            vec_3_1 minus_vel_;
            vec_3_1 minus_att_;

            vec_3_1 bg_hat_; // estimate of bias in gyro
            vec_3_1 ba_hat_; // estimate of bias in accel

            vec_3_1 wb_b_; // body frame angular rate in body frame
            mat_3_3 Omega_b_; // skew symetric of angular rate in body frame

            vec_3_1 we_i_; // rotation rate of earth
            mat_3_3 Omega_e_; // skew symmetric of rotation rate of earth

            mat_3_3 C_b_n_plus_; // rotation matrix from body to nav frame (based on attitude estimate)
            mat_3_3 C_n_b_plus_; // rotation matrix from nav to body frame (based on attitude estimate)

            mat_3_3 C_b_n_minus_;
            mat_3_3 C_n_b_minus_;

            vec_3_1 pos_; // ECEF GPS position solution
            vec_3_1 vel_; // ECEF GPS velocity solution

            // identity matrices and such
            mat_3_3 I_3_;
            
            double time_;
            double dt_;
            double geocentric_radius_;

            void mechanizeSolution();
            void checkInitStatus();
            void generateErrorStateSTM(); // generate error state state state transition matrix
            void generateProcessCovarMat(); 
            void propagateErrorState();

            // private setters
            void setFullStateEstimate(vec_3_1& pos,vec_3_1& vel,vec_3_1& att);

        protected:


    }; // end of class
} // end of namespace

#endif