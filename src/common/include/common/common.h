#ifndef COMMON_H
#define COMMON_H

// Eigen
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"

// CPP Headers
#include <iostream>
#include <cmath>

typedef Eigen::Matrix<double,3,1> vec_3_1;
typedef Eigen::Matrix<double,32,3> mat_32_3;
typedef Eigen::Matrix<double,32,27> mat_32_27;
typedef Eigen::Matrix<double,1,27> vec_1_27;
typedef Eigen::Matrix<double,6,1> vec_6_1;
typedef Eigen::Matrix<double,7,1> vec_7_1;
typedef Eigen::Matrix<double,4,1> vec_4_1;

/*
    To-Do:
    1) Function to Calculate SV States
        [needs testing]
*/

/*
    Order of Ephemeris
    IODE_sf1
    IODE_sf2
    gpsWeek
    t_oe
    A
    deltaN
    M_0
    e
    omega
    C_uc
    C_us
    C_rc
    C_rs
    C_ic
    C_is
    i_0
    iDot
    omega_O
    omegaDot
    IODC
    t_oc
    T_GD
    a_f0
    a_f1
    a_f2

*/



namespace cpp_nav_filt
{
    class Common
    {
        public:

            Common();
            ~Common();

            // Getter and Passer for Common
            void receiveSvEphem(vec_1_27& ephem_in,const int& sv_in);
            void sendSvEphem(vec_1_27& ephem_out,const int& desired_sv);

            void sendUnitVectors(vec_4_1& X_hat,Eigen::MatrixXd& sv_pos,Eigen::MatrixXd& H);

            vec_7_1 sendSvStates(const int& sv_in,const double& transmit_time,const double& transit_time);

            // GPS Constants
            const double gps_pi = M_PI; // PI
            const double omega_e_dot = 7.29211561467*pow(10,-5); // rotation rate of earth
            const double GM = 3.986005*pow(10,14);  
            const double F = -4.442807633*pow(10,-10);
            const int half_week = 302400; // [s]
            const double c = 299792458.0; // [m/s] speed of light
            
        private:


            mat_32_27 sv_ephem; // matrix of satellite ephemeris
            vec_7_1 sv_state_; // is the pos and vel of a satellite we care about given by sv in

            // Internal Variables
            int desired_sv_;
            double T_transmit_;
            double T_transit_;
            double dt;
            double tk;
            double time;

            Eigen::MatrixXd psr_tilde_;
            Eigen::MatrixXd sv_pos_;
            Eigen::MatrixXd H_;
            vec_4_1 x_hat_;

            // Internal Functions
            void calcSvPVStates(vec_7_1& sv_state); // this is adapted from Dr. Bevly's provided class code
            void setCurrentEphem(vec_1_27& ephem,const int& sv);
            double checkT(double time);
            void calcUnitVectors(double num_measurements);
            
        protected:


    }; // end of class

} // end of namespace

#endif