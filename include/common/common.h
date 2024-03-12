#ifndef COMMON_H
#define COMMON_H

#include "cpp_nav_filt_lib/cpp_nav_filt_lib.h"

typedef struct
{
    std::map<std::string,double> ephem_map;
    int sv;
} SatEphemeris;

namespace cpp_nav_filt
{

    class Common
    {
        public:

            Common();
            ~Common();

            void sendUnitVectors(vec_3_1& pos,double& clk,Eigen::MatrixXd& SvPVT,Eigen::MatrixXd& H);
            void sendMeasEst(vec_3_1& pos,vec_3_1& vel,double& clk,double& clk_drift,Eigen::MatrixXd& SvPVT,Eigen::MatrixXd& Yhat);
            void sendElAngles(Eigen::MatrixXd& SvPVT,vec_3_1& pos,Eigen::MatrixXd& el_angles);
            vec_7_1 sendSvStates(const int& sv_in,const double& transmit_time,const double& transit_time);

            // Frame Conversion Functions

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

            std::vector<SatEphemeris> ephem_vect;

            void setRefLla(vec_3_1& lla_in);

            void eul2Rotm(vec_3_1& euler_angles,mat_3_3& C);
            void rotm2Eul(mat_3_3& C,vec_3_1& euler_angles);

            // Nav Functions
            void somiglianaGravityModel(vec_3_1& pos,vec_3_1& gamma_b_n); // see groves p. 72
            void makeSkewSymmetic(vec_3_1& vec_in,mat_3_3& skew_out);
            bool levelInsAccel(vec_3_1& fb_b);
            void initRPfromAccel(vec_3_1& att);
            void meridianRadius(double& lat,double& meridian_radius);
            void geocentricRadius(double& lat,double& geocentric_radius);
            
        private:

            WgsConversions fc;

            vec_7_1 sv_state_; // is the pos and vel of a satellite we care about given by sv idx

            // Internal Variables
            int desired_sv_;
            int num_fb_b_meas_;
            
            double T_transmit_;
            double T_transit_;
            double dt;
            double tk;
            double time;
            double x_comp,y_comp,z_comp,psr_hat;
            double x_comp_vel,y_comp_vel,z_comp_vel;
            double num_sv_;
            double psr_rate_hat;

            // accelerometer init
            vec_3_1 d_var_; // delta variances for accelerometer levelling
            vec_3_1 var_;
            vec_3_1 old_var_;
            vec_1_3 samp_mean_;
            Eigen::MatrixXd fb_b_; // specific forces for levelling

            Eigen::MatrixXd Y_;
            Eigen::MatrixXd Yhat_;
            Eigen::MatrixXd sv_pvt_;
            Eigen::MatrixXd H_;
            Eigen::MatrixXd el_angles_;

            vec_3_1 pos_;
            vec_3_1 vel_;
            double clk_;
            double clk_drift_;
            
            mat_3_3 C_ned_enu; // rotation NED to ENU
            mat_3_3 C_enu_ned; // rotation ENU to NED
            mat_3_3 C_x,C_y,C_z;

            vec_32_1 ones_32_1;
            vec_3_1 ones_3_1;

            // Ephemeris std::map
            std::map<std::string,double> current_ephem_;

            // Temporary Variables for doing the frame conversions
            double ecef_pos_[3];
            double lla_pos_[3];
            double ned_pos_[3];
            double enu_pos_[3];

            // Internal Functions
            void calcSvPVStates(vec_7_1& sv_state); // this is adapted from Dr. Bevly's provided class code
            void setCurrentEphem(const int& sv);
            double checkT(double time);

            void calcUnitVectors();
            void calcPsr(double sv_id);
            void calcPsrRate(double sv_id);
            void calcMeasEst();
            void calcElAngle();
            void calcSampleMean();

            void nanEphemerisMap();
            void eigen2array(double array[3],vec_3_1& eigen);
            void array2eigen(vec_3_1& eigen,double array[3]);

        protected:


    }; // end of class

} // end of namespace

#endif