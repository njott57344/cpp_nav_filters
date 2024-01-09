#include "common/common.h"
#include "gps_least_squares/gps_least_squares.h"

#include <fstream>

int main(int argc,char **argv)
{
    cpp_nav_filt::Common common;
    cpp_nav_filt::GpsLeastSquares gps_least_squares;

    // pointers for file handlers sv ephem and measurements
    std::fstream sv_ephem;
    std::fstream sv_meas;

    // strings of files
    std::string ephem_file = "/home/njo0004/devel/cpp_nav_filters_data/class_ephem.csv";
    std::string meas_file = "/home/njo0004/devel/cpp_nav_filters_data/class_meas_data.csv";

    sv_ephem.open(ephem_file,std::ios::in);
    sv_meas.open(meas_file,std::ios::in);

    // MRI antenna
    vec_3_1 true_x;
    true_x(0) = 422596.629;
    true_x(1) = -5362864.287;
    true_x(2) = 3415493.797;


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

    // ============== GPS Least Squares ===== //

    std::string sv_meas_line,sv_meas_word;
    std::vector<double> sv_measurements;
    std::vector<double> gps_time;
    std::vector<int> sv_id_vect;
    int sv_id;
    double cur_meas,cur_time;

    Eigen::MatrixXd sv_states;
    Eigen::MatrixXd H;

    int j = 0;

    vec_8_1 x_hat;
    vec_7_1 sv_pvt;
    Eigen::MatrixXd meas_vect;
    double num_measurements,num_svs;
    double transit_time,transmit_time;

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
                if(j == 0)
                {
                    cur_time = std::stod(sv_meas_word);
                    gps_time.push_back(cur_time);
                }
                else
                {
                    cur_meas = std::stod(sv_meas_word);
                    sv_measurements.push_back(cur_meas);
                }
                j++;
            }
        }

        j = 0;

        num_measurements = sv_measurements.size();
        Eigen::Map<Eigen::MatrixXd> temp(sv_measurements.data(),1,num_measurements);
        
        meas_vect.resize(num_measurements,1);
        meas_vect = temp.transpose();
                    
        vec_7_1 sv_pvt;

        num_svs = (sv_measurements.size())*0.5;

        for(int j = 0;j<num_svs;j++)
        {
            transit_time = sv_measurements[j]/common.c;
            transmit_time = cur_time - transit_time;
            sv_pvt = common.sendSvStates(sv_id_vect[j+1],transmit_time,transit_time);
            sv_states.conservativeResize(j+1,7);
            sv_states.block<1,7>(j,0) = sv_pvt.transpose();
        }

        gps_least_squares.sendStateEstimate(meas_vect,sv_states,common,x_hat);
        // std::cout<<true_x<<std::endl;

        sv_measurements.clear();

        i++;
    }

    return 0;
}