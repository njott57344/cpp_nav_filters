#include "common/common.h"

#include <fstream>

int main(int argc,char **argv)
{
    cpp_nav_filt::Common common;

    // pointers for file handlers sv ephem and measurements
    std::fstream sv_ephem;
    std::fstream sv_meas;

    // strings of files
    std::string ephem_file = "/home/nicholas/devel/class_data/class_ephem.csv";
    std::string meas_file = "/home/nicholas/devel/class_data/class_meas_data.csv";

    sv_ephem.open(ephem_file,std::ios::in);
    sv_meas.open(meas_file,std::ios::in);

    // ============== Reading Ephemeris ============== //
    std::string ephem_line,ephem_word;
    std::vector<double> ephem_vect;
    double ephemeride;

    std::cout<<"Reading SV Ephemeris and Passing to Common"<<std::endl;

    int i = 0; // ctr

    while(std::getline(sv_ephem,ephem_line))
    {
        std::stringstream s(ephem_line);
        // skip first line
        if(i>0)
        {
            // cast read string to a double vector
            while(std::getline(s,ephem_word,','))
            {                
                ephemeride = std::stod(ephem_word);
                ephem_vect.push_back(ephemeride);
            }

            vec_1_27 ephem_eigen_vect(ephem_vect.data());

            common.receiveSvEphem(ephem_eigen_vect,i);
        }

        i++;
    }
    
    i = 0;

    std::cout<<"Ephemeris is Read"<<std::endl;

    // ============== Reading GPS Psr and Dopplers ===== //

    std::string sv_meas_line,sv_meas_word;
    std::vector<double> sv_measurements;
    std::vector<int> sv_id_vect;
    int sv_id;
    double cur_meas;

    while(std::getline(sv_meas,sv_meas_line))
    {
        std::stringstream ss(sv_meas_line);

        while(std::getline(ss,sv_meas_word,','))
        {
            if(i == 0)
            {
                sv_id = std::stoi(sv_meas_word);
                sv_id_vect.push_back(sv_id);
            }
            else
            {
                cur_meas = std::stod(sv_meas_word);
                sv_measurements.push_back(cur_meas);
            }
        }

        if(i == 1)
        {
            vec_7_1 sv_pvt;
            double transit_time = sv_measurements[1]/common.c;
            double transmit_time = sv_measurements[0] - transit_time;
            sv_pvt = common.sendSvStates(sv_id_vect[1],transmit_time,transit_time);
            std::cout<<transit_time<<" "<<transmit_time<<std::endl;
        }
        i++;
    }

    return 0;
}