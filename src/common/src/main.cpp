#include "common/common.h"

#include <fstream>

int main(int argc,char **argv)
{
    cpp_nav_filt::Common common;

    // pointers for file handlers sv ephem and measurements
    std::fstream sv_ephem;
    std::fstream sv_meas;

    // strings of files
    std::string ephem_file = "/home/njo0004/devel/cpp_nav_filters_data/class_ephem.csv";
    std::string meas_file = "/home/njo0004/devel/cpp_nav_filters_data/class_meas_data.csv";

    sv_ephem.open(ephem_file,std::ios::in);
    sv_meas.open(meas_file,std::ios::in);

    // MRI antenna
    vec_8_1 true_x;
    true_x(0) = 422593.629;
    true_x(1) = -5362864.287;
    true_x(2) = 3415493.797;
    true_x(3) = 37.0937;
    true_x(4) = 0.0;
    true_x(5) = 0.0;
    true_x(6) = 0.0;
    true_x(7) = 0.0;
    
    std::cout<<true_x<<std::endl;

    // ============== Reading Ephemeris ============== //
    std::string ephem_line,ephem_word;
    std::vector<double> ephem_vect;
    double ephemeride;

    std::cout<<"Reading SV Ephemeris and Passing to Common"<<std::endl;

    int i = 0; // ctr
    
    vec_1_27 ephem_eigen_vect;

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
            
            Eigen::Map<vec_1_27> temp(ephem_vect.data(),1,27);
            ephem_eigen_vect = temp;
            common.receiveSvEphem(ephem_eigen_vect,i);
            ephem_vect.clear();
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

    Eigen::MatrixXd sv_states;
    Eigen::MatrixXd H;

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

        double num_measurements = sv_measurements.size();
        Eigen::Map<Eigen::MatrixXd> temp(sv_measurements.data(),1,num_measurements);
        Eigen::MatrixXd meas_vect = temp;
        
        if(i == 1)
        {        
            vec_7_1 sv_pvt;

            double num_svs = (sv_measurements.size() - 1)*0.5;

            for(int j = 0;j<num_svs;j++)
            {
                double transit_time = sv_measurements[j+1]/common.c;
                double transmit_time = sv_measurements[0] - transit_time;
                sv_pvt = common.sendSvStates(sv_id_vect[j+1],transmit_time,transit_time);
                sv_states.conservativeResize(j+1,7);
                sv_states.block<1,7>(j,0) = sv_pvt.transpose();
            }

            common.sendUnitVectors(true_x,sv_states,H);
            std::cout<<H<<std::endl;
        
        }

        sv_measurements.clear();

        i++;
    }

    return 0;
}