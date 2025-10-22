#ifndef TIGHTLY_COUPLED_INS_H
#define TIGHTLY_COUPLED_INS_H

#include "cpp_nav_filt_lib.h"
#include "sv_manager.h"

namespace cpp_nav_filt{

        typedef struct 
        {
            double lat; // rad
            double lon; // rad
            double alt; // rad
            double vn;  // m/s
            double ve;  // m/s
            double vd;  // m/s
            double r;   // rad
            double p;   // rad
            double y;   // rad
            double clk_b; // m
            double clk_d; // m/s

            double lat_stdev;
            double lon_stdev;
            double alt_stdev;
            double vn_stdev;
            double ve_stdev;
            double vd_stdev;
            double r_stdev;
            double p_stdev;
            double y_stdev;
            double clk_b_stdev;
            double clk_d_stdev;
        }NavState;

    class tightlyCoupledNavigator{
        
        public:

            /*! @brief constructor
            */
            tightlyCoupledNavigator();

            /*! @brief destructor
            */
            ~tightlyCoupledNavigator();

            /*! @brief function to propagate a PVA given IMU measurements
            @param[in] wb angular rates [rad/s]
            @param[in] fb specific forces [m/s/s]
            @param[in] dt time step to propagate over [sec]
            @return bool successful or not
            */
            bool Dynamics(const vec_3_1 & wb, const vec_3_1& fb, const vec_3_1& dt);

            /*! @brief return the current nav frame state of the system
            @return NavState returns out the current nav state (see def of NavState struct)
            */
            NavState returnNavState();

            /*! @brief performs a single antenna tightly coupled correction for a given set of psr, psr rate and SvPVT states
            @param[in] psr pseudo-range measurements from antenna to satellites [m]
            @param[in] psr_rate pseudo-range-rate measurements from antenna to satellites [m/s]
            @param[in] la lever arm from IMU to antenna
            @param[in] SvPVT satellite PVT states (needed for unit vectors)
            */
            bool tightlyCoupledCorrection(const Eigen::MatrixXd & psr, const Eigen::MatrixXd & psr_rate, const vec_3_1 & la, const Eigen::MatrixXd & SvPVT);

            /*! @brief set the initial position of the navigator 
            @param[in] r_nb_n [lat,lon,alt] in [rad,rad,m]
            */
            bool setInitialPosition(const vec_3_1 & r_nb_n);

            /*! @brief set the initial velocity of the navigator 
            @param[in] v_nb_n [vn,ve,vd] in m/s
            */
            bool setInitialVelocity(const vec_3_1 & v_nb_n);
            
            /*! @brief set the initial attitude of the navigator
            @param[in] a_nb [roll,pitch,yaw] euler angles (rad)
             */
            bool setInitialAttitude(const vec_3_1 & a_nb);
            
        private:

            /* Scalar's for Nav*/
            double phi_{0.0};   // Latitude [rad]
            double lamb_{0.0};  // Longitude [rad]
            double h_{0.0};     // Altitude [m]
            double vn_{0.0};    // North velocity [m/s]
            double ve_{0.0};    // East velocity [m/s]
            double vd_{0.0};    // Down velocity [m/s]
            double r_{0.0};     // Roll angle [rad]
            double p_{0.0};     // Pitch angle [rad]
            double y_{0.0};     // Yaw angle [rad]
            double clk_b_{0.0}; // clock bias [m]
            double clk_d_{0.0}; // clock drift [m/s]

            /*
                error state is:
                [r_nb_n cartesian position;
                 v_nb_n cartesian velocity;
                 a_nb r,p,y euler angles
                 bg gyro bias error
                 ba accel bias error
                 clk_b clock bias
                 clk_d clock drift]
            */
            vec_17_1 dx_;

            vec_4_1 q_nb;
            mat_3_3 C_nb_;



            vec_3_1 bg_;
            vec_3_1 ba_;

            mat_17_17 P_;

    }; // end of class
} // end of namespace

#endif