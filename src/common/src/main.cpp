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

    std::cout<<"Ephemeris is Read"<<std::endl;

    vec_1_27 ephem_out;

    common.sendSvEphem(ephem_out,2);

    return 0;
}