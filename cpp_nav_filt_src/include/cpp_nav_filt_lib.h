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

/*! @brief variable type 3 by 1 matrix */
typedef Eigen::Matrix<double,3,1> vec_3_1;

/*! @brief variable type 1 by 3 matrix */
typedef Eigen::Matrix<double,1,3> vec_1_3;

/*! @brief variable type 32 by 3 matrix */
typedef Eigen::Matrix<double,32,3> mat_32_3;

/*! @brief variable type 32 by 37 matrix*/
typedef Eigen::Matrix<double,32,27> mat_32_27;

/*! @brief variable type 1 by 27 matrix */
typedef Eigen::Matrix<double,1,27> vec_1_27;

/*! @brief variable type 6 by 1 matrix */
typedef Eigen::Matrix<double,6,1> vec_6_1;

/*! @brief variable type 7 by 1 matrix */
typedef Eigen::Matrix<double,7,1> vec_7_1;

/*! @brief variable type 4 by 1 matrix */
typedef Eigen::Matrix<double,4,1> vec_4_1;

/*! @brief variable type 8 by 1 matrix */
typedef Eigen::Matrix<double,8,1> vec_8_1;

/*! @brief variable type 32 by 1 matrix */
typedef Eigen::Matrix<double,32,1> vec_32_1;

/*! @brief variable type 4 by 4 matrix */
typedef Eigen::Matrix<double,4,4> mat_4_4;

/*! @brief variable type 3 by 3 matrix */
typedef Eigen::Matrix<double,3,3> mat_3_3;

/*! @brief variable type 15 by 1 matrix */
typedef Eigen::Matrix<double,15,1> vec_15_1;

/*! @brief variable type 15 by 15 matrix */
typedef Eigen::Matrix<double,15,15> mat_15_15;

/*! @brief variable type 2 by 1 matrix */
typedef Eigen::Matrix<double,2,1> vec_2_1;

/*! @brief variable type 3 by 15 matrix */
typedef Eigen::Matrix<double,3,15> mat_3_15;

/*! @brief variable type 9 by 1 matrix */
typedef Eigen::Matrix<double,9,1> vec_9_1;

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
    const double w_e = 7.292115*pow(10,-5); // rotation rate of earth [rad/s]
    
    // Parameters of WGS84 Ellipsoid
    const double Ro = 6378137.0; // WGS84 equatorial radius [m]
    const double e = 0.0818191908425; // ecentricity of earth
    const double flattening = 1/298.257223563;
    const double Rp = 6356752.31425; // WGS84 polar radius [m]
    const double A = 6378137.0;
    
    // Frequency Constant
    const double f_l1 = 1.57542*pow(10,9);
    const double f_l2 = 1.2276*pow(10,9);
    const double f_l5 = 1.176*pow(10,9);

    // ============ Meas Estimate ========== //
    /*! @brief function to calculate a matrix of unit vectors from an antenna to a set of satellites
        @param[in] SvPVT matrix of SV states ordered [x,y,z,dx,dy,dz,clk_correction]
        @param[in] ecef_pos estimate of Earth Frame Position
        @param[in] clk_b estimate of receiver clock bias
        @return matrix of unit vectors rows = num svs, cols = [ux,uy,uz]
    */
    Eigen::MatrixXd calcUnitVectors(Eigen::MatrixXd& SvPVT,vec_3_1& ecef_pos,double& clk_b);

    /*! @brief function to calculate pseudorange estimates from an antennna to a set of satellites
        @param[in] SvPVT matrix of SV states ordered [x,y,z,dx,dy,dz,clk_correction]
        @param[in] ecef_pos estimate of Earth Frame Position
        @param[in] clk_b estimate of receiver clock bias
        @return matrix of pseudoranges rows = num svs, col = pseudorange
    */
    Eigen::MatrixXd calcPsr(Eigen::MatrixXd& SvPVT,vec_3_1& ecef_pos,double& clk_b);
    
    /*! @brief function to calculate pseudorange rate estimates from an antenna to a set of satellites
        @param[in] SvPVT matrix of SV states ordered [x,y,z,dx,dy,dz,clk_corrections]
        @param[in] ecef_pos estimate of Earth Frame Position
        @param[in] ecef_vel estimate of Earth Frame Velocity
        @param[in] clk_b estimate of receiver clock bias
        @param[in] clk_d estimate of receiver clock drift
        @return matrix of pseudorange rates rows = num svs, col = pseudorange rates
    */
    Eigen::MatrixXd calcPsrRate(Eigen::MatrixXd& SvPVT,vec_3_1& ecef_pos,vec_3_1& ecef_vel,double& clk_b,double& clk_d);

    /*! @brief function to calculate a full set of measurement estimates
        @param[in] SvPVT matrix of SV states ordered [x,y,z,dx,dy,dz,clk_corrections]
        @param[in] ecef_pos estimate of Earth Frame Position
        @param[in] ecef_vel estimate of Earth Frame Velocity
        @param[in] clk_b estimate of receiver clock bias
        @param[in] clk_d estimate of receiver clock drift
        @return stacked vector of [pseudoranges,pseudorange rates]'
    */
    Eigen::MatrixXd calcMeasEst(Eigen::MatrixXd& SvPVT,vec_3_1& ecef_pos,vec_3_1& ecef_vel,double& clk_b,double& clk_d);

    /*! @brief function to calculate elevation angle from a receiver antenna to a set of satellites
        @note this is WIP
    */
    void calcElAngle();
  
    // ========= Frame Conversion Functions ============= //

    /*! @brief function converts Earth Frame Positions to LLA Positions
        @param[in] ecef_pos Earth Frame Position to be converted to LLA position
        @return returns lla output equivalent to Earth Frame input
    */
    vec_3_1 ecef2llaPos(vec_3_1& ecef_pos);
    
    /*! @brief function converts LLA Position to Earth Frame Position
        @param[in] lla_pos lat lon altitude Position to be converted to Earth Frame Position
        @return returns Earth Frame position from lla input
    */
    vec_3_1 lla2ecefPos(vec_3_1& lla_pos);
    
    // ecef to/from ned

    /*! @brief function converts Earth Frame Position to North East Down Position
        @param[in] ecef_pos Earth Frame Position to be converted to North East Down Position
        @param[in] ref_lla reference latitude longitude altitude coordinate defining local tangent frame
        @return returns North East Down position from Earth Frame input
    */
    vec_3_1 ecef2nedPos(vec_3_1& ecef_pos,vec_3_1& ref_lla);

    /*! @brief function converts North East Down to an Earth Frame position
        @param[in] ned_pos North East Down Position to be converted to Earth Frame
        @param[in] ref_lla reference latitude longitude altitude coordinate defining local tangent frame
        @return returns Earth Frame position from North East Down Position
    */
    vec_3_1 ned2ecefPos(vec_3_1& ned_pos,vec_3_1& ref_lla);
    
    // ecef to/from enu
    /*! @brief function converts Earth Frame Position to East North Up Position
        @param[in] ecef_pos Earth Frame Position to be converted to East North Up Position
        @param[in] ref_lla reference latitude longitude altitude coordinate defining local tangent frame
        @return returns East North Up position from Earth Frame Coordinate
    */
    vec_3_1 ecef2enuPos(vec_3_1& ecef_pos,vec_3_1& ref_lla);

    /*! @brief function converts East North Up Position to Earth Frame Position
        @param[in] enu_pos East North Up Position to be converted to Earth Frame Position
        @param[in] ref_lla reference latitude longitude altitude coordinate defining local tangent frame
        @return returns Earth Frame Position from East North Up Position
    */
    vec_3_1 enu2ecefPos(vec_3_1& enu_pos,vec_3_1& ref_lla);

    // enu to/from ned
    /*! @brief function converts East North Up Position to North East Down Position
        @param[in] enu_pos East North Up Position to be converted to North East Down Position
        @return returns North East Down Position from East North Up position
    */
    vec_3_1 enu2nedPos(vec_3_1& enu_pos);
    
    /*! @brief function converts North East Down Position to East North Up Position
        @param[in] ned_pos North East Down Position to be converted to East North Up Position
        @return returns East North Up Position from North East Down Position
    */
    vec_3_1 ned2enuPos(vec_3_1& ned_pos);
    
    // lla to/from ned
    /*! @brief functions converts Latitude Longitude Altitude position to North East Down Position
        @param[in] lla_pos Latitude Longitude Altitude position to be converted to North East Down Position
        @param[in] ref_lla reference latitude longitude altitude coordinate defining local tangent frame
        @return returns North East Down Position from Latitude Longitude Altitude Position
     */
    vec_3_1 lla2nedPos(vec_3_1& lla_pos,vec_3_1& ref_lla);

    /*! @brief function converts North East Down Position to a Latitude Longitude Altitude Position
        @param[in] ned_pos North East Down Position to be converted to Latitude Longitude Altitude Position
        @param[in] ref_lla reference latitude longitude altitude coordinate defining local tangent frame
        @return returns Latitude Longitude Altitude Position from North East Down Position
    */
    vec_3_1 ned2llaPos(vec_3_1& ned_pos,vec_3_1& ref_lla);

    // lla to/from enu
    /*! @brief functions converts Latitude Longitude Altitude position to East North Up Position
        @param[in] lla_pos Latitude Longitude Altitude position to be converted to East North Up Position
        @param[in] ref_lla reference latitude longitude altitude coordinate defining local tangent frame
        @return returns East North Up Position from Latitude Longitude Altitude Position
     */
    vec_3_1 lla2enuPos(vec_3_1& lla_pos,vec_3_1& ref_lla);
    
    /*! @brief function converts East North Up Position to a Latitude Longitude Altitude Position
        @param[in] enu_pos East North Up Position to be converted to Latitude Longitude Altitude Position
        @param[in] ref_lla reference latitude longitude altitude coordinate defining local tangent frame
        @return returns Latitude Longitude Altitude Position from East North Up Position
    */    
    vec_3_1 enu2llaPos(vec_3_1& enu_pos,vec_3_1& ref_lla);
    
    // ---------- Velocities --------------------------- //
    // ecef to/from ned velocities
    /*! @brief function converts a Earth Frame Velocity to a North East Down Velocity
        @param[in] ecef_vel Earth Frame Velocity to be converted to North East Down Velocity
        @param[in] ref_lla reference latitude longitude altitude coordinate defining local tangent frame
        @return returns North East Down Velocity from Earth Frame Velocity
    */
    vec_3_1 ecef2nedVel(vec_3_1& ecef_vel,vec_3_1& ref_lla);

    /*! @brief function converts a North East Down Velocity to an Earth Frame Velocity
        @param[in] ned_vel North East Down Velocity to be converted to Earth Frame Velocity
        @param[in] ref_lla reference latitude longitude altitude coordinate defining local tangent frame
        @return returns Earth Frame Velocity from North East Down Velocity
    */
    vec_3_1 ned2ecefVel(vec_3_1& ned_vel,vec_3_1& ref_lla);

    // ecef to/from enu velocities
    /*! @brief function converts a Earth Frame Velocity to a East North Up Velocity
        @param[in] ecef_vel Earth Frame Velocity to be converted to East North Up Velocity
        @param[in] ref_lla reference latitude longitude altitude coordinate defining local tangent frame
        @return returns East North Up Velocity from Earth Frame Velocity
    */
    vec_3_1 ecef2enuVel(vec_3_1& ecef_vel,vec_3_1& ref_lla);
    
    /*! @brief function converts a East North Up Velocity to an Earth Frame Velocity
        @param[in] enu_vel East North Up Velocity to be converted to Earth Frame Velocity
        @param[in] ref_lla reference latitude longitude altitude coordinate defining local tangent frame
        @return returns Earth Frame Velocity from East North Up Velocity
    */    
    vec_3_1 enu2ecefVel(vec_3_1& enu_vel,vec_3_1& ref_lla);

    // DCM transformations
    /*! @brief function calculates the DCM to go from a local tangent frame in NED to the Earth Frame
        @param[in] lla_pos Latitude Longitude Altitude position to calculate DCM from
        @return returns the DCM to rotate FROM a North East Down reference frame TO Earth Frame
    */
    mat_3_3 ned2ecefDCM(vec_3_1& lla_pos);

    /*! @brief function calculates the DCM to go from Earth Frame to a local tangent frame in NED
        @param[in] lla_pos Latitude Longitude Altitude position to calculate DCM from
        @return returns the DCM to rotate FROM an Earth Frame TO a North East Down frame 
    */
    mat_3_3 ecef2nedDCM(vec_3_1& lla_pos);

    /*! @brief function calculates the DCM to go from an Earth Frame to a local tangent frame in ENU
        @param[in] lla_pos Latitude Longitude Altitude position to calculate DCM from
        @return returns the DCM to rotate FROM an Earth Frame TO a East North Up frame
    */
    mat_3_3 ecef2enuDCM(vec_3_1& lla_pos);

    /*! @brief function to calculate the DCM to go from a local tangent frame in ENU to an Earth Frame
        @param[in] lla_pos Latitude Longitude Altitude position to calculate DCM from
        @return returns the DCM to rotate FROM an East North Up tangent frame TO an Earth Frame
    */
    mat_3_3 enu2ecefDCM(vec_3_1& lla_pos);

    // ============ Nav Functions ========== //
    
    // dealing with euler angles
    /*! @brief function to calculate a DCM from a vector of roll pitch yaw euler angles
        @param[in] euler_angles vector ordered [roll,pitch,yaw]' of euler angles
        @return returns the DCM equivalent of a vector of euler angles
    */
    mat_3_3 eul2Rotm(vec_3_1& euler_angles);

    /*! @brief function to calculate the DCM from a set of roll pitch yaw euler angles and lla position to an Earth Frame
        @param[in] euler_angles vector ordered [roll,pitch,yaw]' of euler angles for body to nav DCM
        @param[in] lla_pos Latitude Longitude Altitude position to calculate DCM for nav to ecef DCM
        @return returns the DCM to go from euler angles to an Earth Frame coordinate
        @note See Paul Groves p. 38 eq 2.22
    */
    mat_3_3 eul2EcefDCM(vec_3_1& euler_angles,vec_3_1& lla_pos);

    /*! @brief function to calculate a set of euler angles from a local tangent frame DCM
        @param[in] C DCM from body frame to local tangent navigation frame 
        @return returns euler angles [roll,pitch,yaw]
        @note see Paul Groves p. 38 eqs 2.24-2.25
    */
    vec_3_1 rotm2Eul(mat_3_3& C);

    /*! @brief function to get the skew symmetric equivalent of a 3x1 vector
        @param[in] vec_in vector to make the skew symmetric of
        @return returns the skew symmetric model
    */
    mat_3_3 makeSkewSymmetic(vec_3_1& vec_in);

    // gravity model
    /*! @brief function to find the gravitational force in the ecef frame at a given position
        @param[in] ecef_pos Earth Frame Postition to find the gravitational force at
        @return returns a 3x1 vector of gravitational force
        @note See Paul Groves p. 72 eq 2.142
    */
    vec_3_1 ecefGravitation(vec_3_1& ecef_pos);

    /*! @brief function to find the acceleration due to gravity in the Earth Frame at an Earth Frame Position
        @param[in] ecef_pos Earth Frame Position to find acceleration due to gravity at
        @return acceleration due to gravity in the Earth Frame
        @note See Paul Groves p. 70 eq 2.133
    */
    vec_3_1 ecefGravity(vec_3_1& ecef_pos);

    // WGS Ellipoidal Earth Stuff
    /*! @brief function to find the Transverse Radius of Curvature 
        @param[in] lat latitude to find Trasverse Radius of Curvature at
        @return Transverse Radius of Curvature
        @note See Paul Groves p. 59 eq 2.106
    */
    double transverseRadiusOfCurvature(double& lat);

    /*! @brief function to calculate Meridian Radius of Curvature
        @param[in] lat latitude to find Meridian Radius of Curvature at
        @return Meridian Radius of Curvature
        @note See Paul Groves p. 59 eq 2.105
    */
    double meridianRadiusOfCurvature(double& lat);

    /*! @brief function to calculate Geocentric Radius
        @param[in] lat latitude to find Geocentric Radius at
        @return Geocentric Radius
        @note See Paul Groves p. 71 eq 2.137
    */
    double geocentricRadius(double& lat);

    // Init Roll Pitch
    /*! @brief function to calculate levelling paramters roll and pitch from specific force measurements
        @param[in] f_ib_b specific force of body wrt inertial frame in body frame
        @return returns roll and pitch from specific forces
        @note See Paul Groves p. 670 eq 16.27
        @note assumes body is static
    */
    vec_2_1 levelInsAccel(Eigen::MatrixXd& f_ib_b); // groves p. 670
}

#endif