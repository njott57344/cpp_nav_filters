#ifndef CPP_NAV_FILT_LIB_H
#define CPP_NAV_FILT_LIB_H

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
typedef Eigen::Matrix<double,1,3> vec_1_3;
typedef Eigen::Matrix<double,32,3> mat_32_3;
typedef Eigen::Matrix<double,32,27> mat_32_27;
typedef Eigen::Matrix<double,1,27> vec_1_27;
typedef Eigen::Matrix<double,6,1> vec_6_1;
typedef Eigen::Matrix<double,7,1> vec_7_1;
typedef Eigen::Matrix<double,4,1> vec_4_1;
typedef Eigen::Matrix<double,8,1> vec_8_1;
typedef Eigen::Matrix<double,32,1> vec_32_1;
typedef Eigen::Matrix<double,4,4> mat_4_4;
typedef Eigen::Matrix<double,3,3> mat_3_3;
typedef Eigen::Matrix<double,15,1> vec_15_1;
typedef Eigen::Matrix<double,15,15> mat_15_15;

namespace cpp_nav_filt
{
    // ========= Math Constants ============= //

    const double R2D = 180/M_PI; // FROM [rad] TO [deg]
    const double D2R = M_PI/180; // FROM [deg] TO [rad]

    // GPS Constants
    const double gps_pi = M_PI; // PI
    const double omega_e_dot = 7.29211561467*pow(10,-5); // rotation rate of earth
    const double GM = 3.986005*pow(10,14);  
    const double F = -4.442807633*pow(10,-10);
    const int half_week = 302400; // [s]
    const double c = 299792458.0; // [m/s] speed of light
    const double J2 = 1.082627*pow(10,-3); // J2 gravity constant (from Groves p. 72);
    const double mu_g = 3.986004418*pow(10,14); // Earth gravity constant (from Goves p. 71)
    const double Ro = 6378137.0; // WGS84 equatorial radius [m]
    const double w_e = 7.292115*pow(1,-5); // rotation rate of earth [rad/s]
    const double e = 0.0818191908425; // ecentricity of earth

    // Frequency Constant
    const double f_l1 = 1.57542*pow(10,9);
    const double f_l2 = 1.2276*pow(10,9);
    const double f_l5 = 1.176*pow(10,9);

    // ============ Meas Estimate ========== //
    Eigen::MatrixXd calcUnitVectors(Eigen::MatrixXd& SvPVT,vec_3_1& ecef_pos,double& clk_b);
    Eigen::MatrixXd calcPsr(Eigen::MatrixXd& SvPVT,vec_3_1& ecef_pos,double& clk_b,Eigen::MatrixXd& psr_hat);
    Eigen::MatrixXd calcPsrRate(Eigen::MatrixXd& SvPVT,vec_3_1& ecef_pos,vec_3_1& ecef_vel,double& clk_b,double& clk_d);
    Eigen::MatrixXd calcMeasEst(Eigen::MatrixXd& SvPVT,vec_3_1& ecef_pos,vec_3_1& ecef_vel,double& clk_b,double& clk_d);

    void calcElAngle();
  
    // ========= Frame Conversion Functions ============= //
    void setRefLla(vec_3_1& lla_in);

    // ecef to/from lla
    void convertECEF2LLA(vec_3_1& ecef_pos,vec_3_1& lla_pos);
    void convertLLA2ECEF(vec_3_1& lla_pos,vec_3_1& ecef_pos);
    
    // ecef to/from ned
    void convertECEF2NED(vec_3_1& ecef_pos,vec_3_1& ned_pos,vec_3_1& ref_lla);
    void convertNED2ECEF(vec_3_1& ned_pos,vec_3_1& ecef_pos,vec_3_1& ref_lla);
    
    // enu to/from ned
    void convertNED2ENU(vec_3_1& ned_pos,vec_3_1& enu_pos);
    void convertENU2NED(vec_3_1& enu_pos,vec_3_1& ned_pos);

    // ecef to/from enu
    void convertECEF2ENU(vec_3_1& ecef_pos,vec_3_1& enu_pos,vec_3_1& ref_lla);
    void convertENU2ECEF(vec_3_1& enu_pos,vec_3_1& ecef_pos,vec_3_1& ref_lla);

    // ============ Nav Functions ========== //
    
    // dealing with euler angles
    mat_3_3 eul2Rotm(vec_3_1& euler_angles);
    vec_3_1 rotm2Eul(mat_3_3& C);
    mat_3_3 makeSkewSymmetic(vec_3_1& vec_in);

    // gravity model
    vec_3_1 somiglianaGravityModel(vec_3_1& ecef_pos); // see groves p. 72
    double meridianRadius(double& lat,double);
    double geocentricRadius(double& lat);


    // Init Roll Pitch
    void calcSampleMean();
    bool levelInsAccel(vec_3_1& fb_b);
    void initRPfromAccel(vec_3_1& att);
}

#endif