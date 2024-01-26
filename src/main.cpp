#include "common/common.h"
#include "gps_least_squares/gps_least_squares.h"
#include "matplot/matplot.h"

#include <fstream>

namespace plt = matplot;

int main(int argc,char **argv)
{
    cpp_nav_filt::Common common;
    cpp_nav_filt::GpsLeastSquaresSettings LeastSquaresSettings;
    LeastSquaresSettings.weighted_least_squares = false;
    cpp_nav_filt::GpsLeastSquares gps_least_squares(LeastSquaresSettings);
    
    // pointers for file handlers sv ephem and measurements
    std::fstream sv_ephem;
    std::fstream sv_meas;

    std::ofstream solution_out;

    // strings of files
    std::string ephem_file = "/home/njo0004/devel/cpp_nav_filters_data/dynamic_ephem.csv";
    std::string meas_file = "/home/njo0004/devel/cpp_nav_filters_data/dynamic_measurements.csv";
    std::string output_file = "/home/njo0004/devel/cpp_nav_filters_data/output.csv";

    sv_ephem.open(ephem_file,std::ios::in);
    sv_meas.open(meas_file,std::ios::in);
    solution_out.open(output_file);

    // MRI antenna
    vec_3_1 true_x,true_lla;
    true_x(0) = 422596.629;
    true_x(1) = -5362864.287;
    true_x(2) = 3415493.797;

    // ============== Reading Ephemeris ============== //
    std::string ephem_line,ephem_word;
    std::vector<std::string> ephem_label_vect;
    std::vector<double> ephem_vect;
    double ephemeride;

    std::cout<<"Reading SV Ephemeris and Passing to Common"<<std::endl;

    int i = 0; // ctr
    int j = 0; // ctr 2

    vec_1_27 ephem_eigen_vect;

    while(std::getline(sv_ephem,ephem_line))
    {
        std::stringstream s(ephem_line);

        while(std::getline(s,ephem_word,','))
        {
            if(i == 0)
            {
                ephem_label_vect.push_back(ephem_word);
                //std::cout<<ephem_word;
            }
            else
            {
                if(j == 0)
                {   
                    common.ephem_vect[i-1].sv = std::stoi(ephem_word);
                }
                else
                {
                    common.ephem_vect[i-1].ephem_map[ephem_label_vect[j]] = std::stod(ephem_word);
                }
                j++;
            }
        }
        i++;
        j = 0;
    }
    
    i = 0;

    std::cout<<"Ephemeris is Read"<<std::endl;

    // ============== GPS Least Squares ===== //

    std::string sv_meas_line,sv_meas_word;
    std::vector<double> sv_measurements;
    std::vector<double> gps_time;

    std::vector<double> lat_soln;
    std::vector<double> lon_soln;
    std::vector<double> alt_soln;

    std::vector<int> sv_id_vect;
    std::vector<int> valid_sv_id_vect;
    int sv_id;
    double cur_meas,cur_time;

    Eigen::MatrixXd sv_states;
    Eigen::MatrixXd H;
    Eigen::MatrixXd state_soln;

    j = 0;

    vec_8_1 x_hat;
    mat_4_4 DOP;

    vec_7_1 sv_pvt;
    Eigen::MatrixXd meas_vect;
    double num_measurements,num_svs;
    double transit_time,transmit_time;

    vec_3_1 lla_pos,ecef_pos;

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

                    if(!std::isnan(cur_meas))
                    {
                        sv_measurements.push_back(cur_meas);
                        valid_sv_id_vect.push_back(sv_id_vect[j]);
                    }
                }
                j++;
            }
        }

        j = 0;

        if(i>0)
        {
            num_measurements = sv_measurements.size();
            num_svs = (sv_measurements.size())*0.5;

            Eigen::Map<Eigen::MatrixXd> temp(sv_measurements.data(),1,num_measurements);

            // resizing matrices
            meas_vect.resize(num_measurements,1);
            sv_states.resize(num_svs,7);
            meas_vect = temp.transpose();

            for(int j = 0;j<num_svs;j++)
            {
                transit_time = sv_measurements[j]/cpp_nav_filt::c;
                transmit_time = cur_time - transit_time;
                sv_pvt = common.sendSvStates(valid_sv_id_vect[j],transmit_time,transit_time);
                sv_states.block<1,7>(j,0) = sv_pvt.transpose();
            }

            gps_least_squares.sendStateEstimate(meas_vect,sv_states,common,x_hat);
            gps_least_squares.sendDOPEstimate(x_hat,sv_states,common,DOP);
            ecef_pos = x_hat.block<3,1>(0,0);

            common.convertECEF2LLA(ecef_pos,lla_pos);

            lat_soln.push_back(lla_pos[0]);
            lon_soln.push_back(lla_pos[1]);
            alt_soln.push_back(lla_pos[2]);

            sv_measurements.clear();
            valid_sv_id_vect.clear();
            
            x_hat.setZero();
        }
                
        i++;
    }
        
    solution_out.close();
    
    plt::geoplot(lat_soln,lon_soln,"*");
    plt::geolimits({32.5667,32.6167},{-85.4750,-85.5083}); // geolimits for auburn al
    plt::hold(plt::on);
    plt::show();

    return 0;
}