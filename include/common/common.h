#ifndef COMMON_H
#define COMMON_H

// Eigen
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"

// CPP Headers
#include <iostream>
#include <cmath>
#include <vector>
#include <map>

#include "frame_conversions/frame_conversions.h"

typedef Eigen::Matrix<double,3,1> vec_3_1;
typedef Eigen::Matrix<double,32,3> mat_32_3;
typedef Eigen::Matrix<double,32,27> mat_32_27;
typedef Eigen::Matrix<double,1,27> vec_1_27;
typedef Eigen::Matrix<double,6,1> vec_6_1;
typedef Eigen::Matrix<double,7,1> vec_7_1;
typedef Eigen::Matrix<double,4,1> vec_4_1;
typedef Eigen::Matrix<double,8,1> vec_8_1;
typedef Eigen::Matrix<double,32,1> vec_32_1;
typedef Eigen::Matrix<double,4,4> mat_4_4;

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

    // GPS Constants
    const double gps_pi = M_PI; // PI
    const double omega_e_dot = 7.29211561467*pow(10,-5); // rotation rate of earth
    const double GM = 3.986005*pow(10,14);  
    const double F = -4.442807633*pow(10,-10);
    const int half_week = 302400; // [s]
    const double c = 299792458.0; // [m/s] speed of light
    
    // Frequency Constant
    const double f_l1 = 1.57542*pow(10,9);
    const double f_l2 = 1.2276*pow(10,9);
    const double f_l5 = 1.176*pow(10,9);

    class Common
    {
        public:

            Common();
            ~Common();

            // Getter and Passer for Common
            void receiveSvEphem(vec_1_27& ephem_in,const int& sv_in);
            void sendSvEphem(vec_1_27& ephem_out,const int& desired_sv);

            void sendUnitVectors(vec_3_1& pos,double& clk,Eigen::MatrixXd& SvPVT,Eigen::MatrixXd& H);
            void sendMeasEst(vec_3_1& pos,vec_3_1& vel,double& clk,double& clk_drift,Eigen::MatrixXd& SvPVT,Eigen::MatrixXd& Yhat);

            vec_7_1 sendSvStates(const int& sv_in,const double& transmit_time,const double& transit_time);

            // Frame Conversion Functions
            void convertECEF2LLA(vec_3_1& ecef_pos,vec_3_1& lla_pos,WgsConversions& frame_conversions);
            void convertLLA2ECEF(vec_3_1& lla_pos,vec_3_1& ecef_pos,WgsConversions& frame_conversions);
            // void convertECEF2NED(vec_3_1& ecef_pos,vec_3_1& )

        private:

            mat_32_27 sv_ephem; // matrix of satellite ephemeris
            vec_7_1 sv_state_; // is the pos and vel of a satellite we care about given by sv idx

            // Internal Variables
            int desired_sv_;
            double T_transmit_;
            double T_transit_;
            double dt;
            double tk;
            double time;
            double x_comp,y_comp,z_comp,psr_hat;
            double x_comp_vel,y_comp_vel,z_comp_vel;
            double num_sv_;
            double psr_rate_hat;

            Eigen::MatrixXd Y_;
            Eigen::MatrixXd Yhat_;
            Eigen::MatrixXd sv_pvt_;
            Eigen::MatrixXd H_;
            vec_3_1 pos_; // [Position;clk;Velocity;clk drift]
            vec_3_1 vel_;
            double clk_;
            double clk_drift_;

            vec_32_1 ones_32_1;

            // Ephemeris std::map
            std::map<std::string,double> current_ephem_;

            // Temporary Variables for doing the frame conversions
            double ecef_pos_[3];
            double lla_pos_[3];

            // Internal Functions
            void calcSvPVStates(vec_7_1& sv_state); // this is adapted from Dr. Bevly's provided class code
            void setCurrentEphem(const int& sv);
            double checkT(double time);

            void calcUnitVectors();
            void calcPsr(double sv_id);
            void calcPsrRate(double sv_id);
            void calcMeasEst();

            void eigen2array(double array[3],vec_3_1& eigen);
            void array2eigen(vec_3_1& eigen,double array[3]);

        protected:


    }; // end of class

} // end of namespace

#endif