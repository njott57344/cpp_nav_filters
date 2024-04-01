#include "cpp_nav_filt_lib/cpp_nav_filt_lib.h"
#include "lc_ins/lc_ins.h"

#include "matplot/matplot.h"

#include "matplotlibcpp/matplotlibcpp.h"

#include <fstream>

namespace plt = matplotlibcpp;

int main(int argc,char **argv)
{
    cpp_nav_filt::LooselyCoupledInsSettings lc_ins_settings;

    // fill this in with actual simualted psd's
    lc_ins_settings.psd_accel_bias = 1;
    lc_ins_settings.psd_accel_noise = 1;
    lc_ins_settings.psd_gyro_noise = 1;
    lc_ins_settings.psd_gyro_bias = 1;
    lc_ins_settings.lever_arm.setZero();

    cpp_nav_filt::LooselyCoupledIns ins(lc_ins_settings);

    std::fstream imu_meas;
    std::fstream imu_truth;

    std::string imu_file = "/home/njo0004/devel/cpp_nav_filters_data/sim/wiggly_path_imu.csv";
    std::string imu_truth_file = "/home/njo0004/devel/cpp_nav_filters_data/sim/wiggly_path_truth.csv";

    imu_meas.open(imu_file,std::ios::in);
    imu_truth.open(imu_truth_file,std::ios::in);

// ============== Truth PVA States ========= //

    std::string truth_line,truth_word;
    std::vector<double> truth_states; // current line of truth Time,Pos,Vel,Att
    std::vector<double> lat_true,lon_true,alt_true; // truth LLA
    std::vector<double> n_true,e_true,d_true; // truth NED Pos
    std::vector<double> dn_true,de_true,dd_true; // truth NED vel
    std::vector<double> roll_true,pitch_true,yaw_true; // truth Attitude
    std::vector<double> time;

    Eigen::Matrix<double,3,3> I3;


    while(std::getline(imu_truth,truth_line))
    {
        std::stringstream truth_stream(truth_line);

        while(std::getline(truth_stream,truth_word,','))
        {
            truth_states.push_back(std::stod(truth_word));
        }
        
        time.push_back((truth_states[0]));
        
        lat_true.push_back(truth_states[1]);
        lon_true.push_back(truth_states[2]);
        alt_true.push_back(truth_states[3]);

        dn_true.push_back(truth_states[4]);
        de_true.push_back(truth_states[5]);
        dd_true.push_back(truth_states[6]);

        roll_true.push_back(truth_states[7]);
        pitch_true.push_back(truth_states[8]);
        yaw_true.push_back(truth_states[9]);

        truth_states.clear();
    }

    // ========== IMU Mechanization ========= //
    vec_3_1 init_pos,init_vel,init_att,init_ba,init_bg;
    vec_3_1 ned_att,lla_pos;
    vec_3_1 x_hat,dx_hat,ecef_att_hat,ned_att_hat;
    vec_3_1 ned_hat,ned_truth,lla_0,lla_truth;
    vec_3_1 lla_hat,ned_vel_hat;
    mat_3_3 Cne; // DCM ned to ecef
    double init_time;
    int j = 0;

    I3.setIdentity();
    lla_pos<<lat_true[0],lon_true[0],alt_true[0];
    init_vel<<0,0,0;
    init_ba<<0,0,0;
    init_bg<<0,0,0;
    init_time = time[0];

    lla_0 << lat_true[0],lon_true[0],alt_true[0];

    Cne = cpp_nav_filt::ned2ecefDCM(init_pos);
    ned_att<<roll_true[0],pitch_true[0],yaw_true[0];
    init_att = Cne*ned_att;
    init_pos = cpp_nav_filt::lla2ecefPos(lla_pos);

    ins.setInitialPosState(init_pos,I3);
    ins.setInitialVelState(init_vel,I3);
    ins.setInitialAttState(init_att,I3);
    ins.setInitialBaState(init_ba,I3);
    ins.setInitialBgState(init_bg,I3);
    ins.setInitialTime(init_time);

    std::string imu_line,imu_word;
    vec_3_1 fb_b,wb_b;
    double cur_time;
    std::vector<double> current_imu_meas;

    std::vector<double> lat_hat,lon_hat,alt_hat;
    std::vector<double> dn_hat,de_hat,dd_hat;
    std::vector<double> roll_hat,pitch_hat,yaw_hat;
    std::vector<double> n_hat,e_hat,d_hat;

    std::vector<double> lat_err,lon_err,alt_err;
    std::vector<double> dn_err,de_err,dd_err;
    std::vector<double> roll_err,pitch_err,yaw_err;
    std::vector<double> n_err,e_err,d_err;

    while(std::getline(imu_meas,imu_line))
    {
        std::stringstream imu_meas_stream(imu_line);
        
        while(std::getline(imu_meas_stream,imu_word,','))
        {
            current_imu_meas.push_back(std::stod(imu_word));
        }

        Eigen::Map<Eigen::MatrixXd> temp(current_imu_meas.data(),1,7);
        fb_b = temp.transpose().block<3,1>(4,0);
        wb_b = temp.transpose().block<3,1>(1,0);
        cur_time = temp(0,0);
        current_imu_meas.clear();

        // send imu measurements to class (triggers time update)
        ins.getImuMeasurements(fb_b,wb_b,cur_time); 

        // get current estimate of PVA
        ins.setPosSol(x_hat);
        ins.setVelSol(dx_hat);
        ins.setAttSol(ecef_att_hat);

        lla_hat = cpp_nav_filt::ecef2llaPos(x_hat);
        Cne = cpp_nav_filt::ecef2nedDCM(lla_hat);
        ned_vel_hat = Cne*dx_hat;
        ned_att_hat = Cne*ecef_att_hat;

        lla_truth << lat_true[j],lon_true[j],alt_true[j];
        ned_truth = cpp_nav_filt::lla2nedPos(lla_truth,lla_0);
        ned_hat = cpp_nav_filt::ecef2nedPos(x_hat,lla_0);

        lat_hat.push_back(lla_hat[0]);
        lon_hat.push_back(lla_hat[1]);
        alt_hat.push_back(lla_hat[2]);

        lat_err.push_back(lat_true[j]-lat_hat[j]);
        lon_err.push_back(lon_true[j]-lon_hat[j]);
        alt_err.push_back(alt_true[j]-alt_hat[j]);

        n_hat.push_back(ned_hat[0]);
        e_hat.push_back(ned_hat[1]);
        d_hat.push_back(ned_hat[2]);

        n_true.push_back(ned_truth[0]);
        e_true.push_back(ned_truth[1]);
        d_true.push_back(ned_truth[2]);

        n_err.push_back(n_true[j]-n_hat[j]);
        e_err.push_back(e_true[j]-e_hat[j]);
        d_err.push_back(d_true[j]-d_hat[j]);

        dn_hat.push_back(ned_vel_hat[0]);
        de_hat.push_back(ned_vel_hat[1]);
        dd_hat.push_back(ned_vel_hat[2]);

        dn_err.push_back(dn_true[j]-dn_hat[j]);
        de_err.push_back(de_true[j]-de_hat[j]);
        dd_err.push_back(dd_true[j]-dd_hat[j]);

        roll_hat.push_back(ned_att_hat[0]);
        pitch_hat.push_back(ned_att_hat[1]);
        yaw_hat.push_back(ned_att_hat[2]);

        roll_err.push_back(roll_true[j] - roll_hat[j]);
        pitch_err.push_back(pitch_true[j]-pitch_hat[j]);
        yaw_err.push_back(yaw_true[j]-yaw_hat[j]);

        j++;
    }   

    plt::suptitle("Positioning");
    plt::subplot(3,1,1);
    plt::plot(time,e_err);
    plt::title("East (m)");
    plt::subplot(3,1,2);
    plt::plot(time,n_err);
    plt::title("North (m)");
    plt::subplot(3,1,3);
    plt::plot(time,d_err);
    plt::title("Down (m)");
    plt::show();
    
    plt::named_plot("Truth",e_true,n_true);
    plt::named_plot("Estimate",e_hat,n_hat);
    plt::title("East Vs North Position");
    plt::legend();
    plt::show();

}