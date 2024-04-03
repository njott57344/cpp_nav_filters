#ifndef SV_MANAGER_H
#define SV_MANAGER_H

#include "cpp_nav_filt_lib.h"

typedef struct
{
    std::map<std::string,double> ephem_map;
    int sv;
} SatEphemeris;

namespace cpp_nav_filt
{

    class SvManager
    {
        public:

            SvManager();
            ~SvManager();

            void sendUnitVectors(vec_3_1& pos,double& clk,Eigen::MatrixXd& SvPVT,Eigen::MatrixXd& H);
            void sendMeasEst(vec_3_1& pos,vec_3_1& vel,double& clk,double& clk_drift,Eigen::MatrixXd& SvPVT,Eigen::MatrixXd& Yhat);
            void sendElAngles(Eigen::MatrixXd& SvPVT,vec_3_1& pos,Eigen::MatrixXd& el_angles);
            vec_7_1 sendSvStates(const int& sv_in,const double& transmit_time,const double& transit_time);

            std::vector<SatEphemeris> ephem_vect;
      
        private:

            vec_7_1 sv_state_; // is the pos and vel of a satellite we care about given by sv idx

            // Internal Variables
            int desired_sv_;
            int num_fb_b_meas_;
            
            double T_transmit_;
            double T_transit_;
            double dt;
            double tk;
            double time;

            vec_32_1 ones_32_1;
            vec_3_1 ones_3_1;

            // Ephemeris std::map
            std::map<std::string,double> current_ephem_;

            // Internal Functions
            void calcSvPVStates(vec_7_1& sv_state); // this is adapted from Dr. Bevly's provided class code
            void setCurrentEphem(const int& sv);
            double checkT(double time);


            void nanEphemerisMap();

        protected:


    }; // end of class

} // end of namespace

#endif