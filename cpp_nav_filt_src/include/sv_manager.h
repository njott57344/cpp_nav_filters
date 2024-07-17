#ifndef SV_MANAGER_H
#define SV_MANAGER_H

#include "cpp_nav_filt_lib.h"

/*! @brief struct for containing all of the ephemeris of a particular sv
    @param ephem_map std::map from an ephermide string name to it's double value
    @param sv integer for which SV this is
*/
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
  
            /*! @brief function to send SV states from Satellite Ephemeris
                @param[in] sv_in index of satellite you want the States of
                @param[in] transmit_time time the satellite transmitted 
                @param[in] transit_time time the satellite message was in transit 
                @return Satellite states [x,y,z,dx,dy,dz,clk_corrections]
            */
            vec_7_1 sendSvStates(const int& sv_in,const double& transmit_time,const double& transit_time);

            /*! @brief vector of Satellite ephemeris structs */
            std::vector<SatEphemeris> ephem_vect;
      
        private:

            vec_7_1 sv_state_; // is the pos and vel of a satellite we care about given by sv idx

            // Internal Variables
            int desired_sv_;
            
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
            /*! @brief function that calculates the SV States
                @param[out] sv_state is the output state of SV pos,vel,clk_corrections
            */
            void calcSvPVStates(vec_7_1& sv_state);

            /*! @brief function that sets the current ephemeris map to the SV under consideration
                @param[in] sv is the integer id for the satellite under current consideration
            */
            void setCurrentEphem(const int& sv);

            /*! @brief function to check half week time
                @param[in] time time to run check on
                @return returns checked time
            */
            double checkT(double time);

            void nanEphemerisMap();

        protected:


    }; // end of class

} // end of namespace

#endif