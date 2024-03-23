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
    std::fstream imu_init_cdn;

    std::string imu_file = "/home/njo0004/devel/cpp_nav_filters_data/straight_driving.csv";
    std::string imu_init_cdn_file = "/home/njo0004/devel/cpp_nav_filters_data/straight_line_init_cdn.csv";

    imu_meas.open(imu_file,std::ios::in);
    imu_init_cdn.open(imu_init_cdn_file,std::ios::in);

// ============== Initial INS Conditions ========= //

    std::string init_line,init_word;
    std::vector<double> init_cdn;
    Eigen::Matrix<double,10,1> eig_init_cdn;
    Eigen::Matrix<double,3,3> I3;
    vec_3_1 init_pos,init_vel,init_att,init_ba,init_bg;
    vec_3_1 ins_lla_0;

    double init_time;
    I3.setIdentity();

    while(std::getline(imu_init_cdn,init_line))
    {
        std::stringstream init_stream(init_line);

        while(std::getline(init_stream,init_word,','))
        {
            init_cdn.push_back(std::stod(init_word));
        }
    }

    Eigen::Map<Eigen::MatrixXd> temp_init(init_cdn.data(),1,10);
    eig_init_cdn = temp_init.transpose();
        
    // INS initial Conditions
    init_pos = eig_init_cdn.segment<3>(1);
    init_vel = eig_init_cdn.segment<3>(4);
    init_att = eig_init_cdn.segment<3>(7);
    init_ba.setZero();
    init_bg.setZero();
    init_time = eig_init_cdn[0];

    ins_lla_0 = cpp_nav_filt::ecef2llaPos(init_pos);

    ins.setInitialPosState(init_pos,I3);
    ins.setInitialVelState(init_vel,I3);
    ins.setInitialAttState(init_att,I3);
    ins.setInitialBaState(init_ba,I3);
    ins.setInitialBgState(init_bg,I3);
    ins.setInitialTime(init_time);

    // ============== Testing INS Mechanization ====== //

    std::string imu_line,imu_word;
    std::vector<double> imu_vect;
    Eigen::Matrix<double,3,1> att_init;
    Eigen::Matrix<double,1,7> imu_eig_vect;
    bool att_is_init = false;
    vec_3_1 fb_b,wb_b;
    vec_3_1 ins_pos,ins_vel,ins_att; // ecef pos states
    vec_3_1 ins_ned,ins_ned_vel;
    double t_current;

    std::vector<double> ins_x,ins_y,ins_z,ins_dx,ins_dy,ins_dz,ins_r,ins_p,ins_h;

    int imu_ctr = 0;

    while(std::getline(imu_meas,imu_line))
    {        
        std::stringstream imu_stream(imu_line);

        while(std::getline(imu_stream,imu_word,','))
        {
            imu_vect.push_back(std::stod(imu_word));
        }

        Eigen::Map<Eigen::MatrixXd> temp_imu(imu_vect.data(),1,7);
        imu_eig_vect = temp_imu;

        fb_b = imu_eig_vect.segment<3>(1).transpose();
        wb_b = imu_eig_vect.segment<3>(4).transpose();
        t_current = imu_eig_vect[0];

        ins.getImuMeasurements(fb_b,wb_b,t_current);

        ins.setPosSol(ins_pos);
        ins.setVelSol(ins_vel);
        ins.setAttSol(ins_att);

        ins_ned = cpp_nav_filt::ecef2nedPos(ins_pos,ins_lla_0);
        ins_ned_vel = cpp_nav_filt::ecef2nedVel(ins_vel,ins_lla_0);

        ins_x.push_back(ins_ned[0]);
        ins_y.push_back(ins_ned[1]);
        ins_z.push_back(ins_ned[2]);
        
        ins_dx.push_back(ins_ned_vel[0]);
        ins_dy.push_back(ins_ned_vel[1]);
        ins_dz.push_back(ins_ned_vel[2]);

        ins_r.push_back(ins_att[0]);
        ins_p.push_back(ins_att[1]);
        ins_h.push_back(ins_att[2]);

        /*
        if(att_is_init == false)
        {
            if(common.levelInsAccel(fb_b))
            {
                common.initRPfromAccel(att_init);
                att_is_init = true;
                std::cout<<att_init*180/M_PI<<std::endl;
            }
        }
        */
        imu_ctr++;

        imu_vect.clear();
    }    

    plt::plot(ins_y,ins_x);
    plt::title("INS NED E vs N Position (m)");
    plt::show();

    plt::plot(ins_dx);
    plt::title("INS E Vel (m/s)");
    plt::show();

    plt::plot(ins_dy);
    plt::title("INS N Vel (m/s)");
    plt::show();

    plt::plot(ins_r);
    plt::title("INS Roll");
    plt::show();

    plt::plot(ins_p);
    plt::title("INS Pitch");
    plt::show();

    plt::plot(ins_h);
    plt::title("INS Heading");
    plt::show();
}