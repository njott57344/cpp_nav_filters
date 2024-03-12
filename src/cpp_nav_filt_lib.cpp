#include "cpp_nav_filt_lib/cpp_nav_filt_lib.h"

namespace cpp_nav_filt
{
    Eigen::MatrixXd calcPsr(Eigen::MatrixXd& SvPVT,vec_3_1& ecef_pos,double& clk_b)
    {
        int num_sv = SvPVT.rows();

        double x_comp,y_comp,z_comp;

        Eigen::Matrix<double,1,1> psr_hat;
        psr_hat.resize(num_sv,1);

        for(int i = 0;i<num_sv;i++)
        {
            x_comp = SvPVT(i,0) - ecef_pos[0];
            y_comp = SvPVT(i,1) - ecef_pos[1];
            z_comp = SvPVT(i,2) - ecef_pos[2];

            psr_hat(i,0) = sqrt(pow(x_comp,2) + pow(y_comp,2) + pow(z_comp,2)) + clk_b;
        }

        return psr_hat;
    }

    Eigen::MatrixXd calcUnitVectors(Eigen::MatrixXd& SvPVT,vec_3_1& ecef_pos,double& clk_b)
    {
        int num_sv = SvPVT.rows();
        double x_comp,y_comp,z_comp;

        Eigen::Matrix<double,1,1> psr_hat,H;
        psr_hat.resize(num_sv,1); // estimated measurement vector
        H.resize(num_sv,4); // estimated unit vector matrix

        psr_hat = calcPsr(SvPVT,ecef_pos,clk_b);

        for(int i = 0;i<num_sv;i++)
        {
            x_comp = SvPVT(i,0) - ecef_pos[0];
            y_comp = SvPVT(i,1) - ecef_pos[1];
            z_comp = SvPVT(i,2) - ecef_pos[2];

            H(i,0) = -x_comp/psr_hat(i);
            H(i,1) = -y_comp/psr_hat(i);
            H(i,2) = -z_comp/psr_hat(i);
            H(i,3) = 1;
        }

        return H;
    }

    Eigen::MatrixXd calcPsrRate(Eigen::MatrixXd& SvPVT,vec_3_1& ecef_pos,vec_3_1& ecef_vel,double& clk_b,double& clk_d)
    {
        int num_sv = SvPVT.rows();
        double x_comp,y_comp,z_comp;

        Eigen::Matrix<double,1,1> psr_rate_hat,H;
        vec_3_1 relative_vel,u;
        
        psr_rate_hat.resize(num_sv,1);
        H.resize(num_sv,4);

        H = calcUnitVectors(SvPVT,ecef_pos,clk_b);

        for(int i = 0;i<num_sv;i++)
        {
            x_comp = SvPVT(i,3) - ecef_vel(0);
            y_comp = SvPVT(i,4) - ecef_vel(1);
            z_comp = SvPVT(i,5) - ecef_vel(2);

            relative_vel << x_comp,y_comp,z_comp;
            u = H.block<1,3>(i,0);

            psr_rate_hat(i,0) = (-u.dot(relative_vel)) + clk_d;
        }

        return psr_rate_hat;
    }

    Eigen::MatrixXd calcMeasEst(Eigen::MatrixXd& SvPVT,vec_3_1& ecef_pos,vec_3_1& ecef_vel,double& clk_b,double& clk_d)
    {
        int num_sv = SvPVT.rows();

        Eigen::Matrix<double,1,1> psr,psr_rate,Y_hat;
        psr.resize(num_sv,1);
        psr_rate.resize(num_sv,1);
        Y_hat.resize(2*num_sv,1);

        psr = calcPsr(SvPVT,ecef_pos,clk_b);
        psr_rate = calcPsrRate(SvPVT,ecef_pos,ecef_vel,clk_b,clk_d);

        Y_hat.block(0,0,num_sv,1) = psr;
        Y_hat.block(num_sv+1,0,num_sv,1) = psr_rate;

        return Y_hat;
    }

    mat_3_3 eul2Rotm(vec_3_1& euler_angles)
    {
        // assumes euler angles are ordered roll pitch yaw
        // see groves p. 38
        mat_3_3 C_x,C_y,C_z;

        vec_3_1 c;
        vec_3_1 s;

        c << std::cos(euler_angles[0]),std::cos(euler_angles[1]),std::cos(euler_angles[2]);
        s << std::sin(euler_angles[0]),std::sin(euler_angles[1]),std::sin(euler_angles[2]);

        C_x << 1,  0,    0,
                0,  c[0], s[0],
                0, -s[0], c[0];

        C_y << c[1], 0, -s[1],
                0,    1,  0,
                s[1], 0,  c[1];

        C_z << c[2], s[2], 0,
                -s[2], c[2], 0,
                0,    0,    1;

        return C_x*C_y*C_z; // rotation nav to body        
    }

    vec_3_1 rotm2Eul(mat_3_3& C)
    {
        vec_3_1 euler_angles;

        euler_angles[0] = std::atan2(C(1,2),C(2,2)); // roll
        euler_angles[1] = -std::asin(C(0,2));        // pitch
        euler_angles[2] = std::atan2(C(0,1),C(0,0)); // yaw
    
        return euler_angles;
    }

    mat_3_3 makeSkewSymmetic(vec_3_1& vec_in)
    {
        mat_3_3 skew_out;

        skew_out.setZero();

        // Z
        skew_out(0,1) = -vec_in[2];
        skew_out(1,0) = vec_in[2];

        //Y
        skew_out(0,2) = vec_in[1];
        skew_out(2,0) = -vec_in[1];

        //X
        skew_out(2,1) = vec_in[0];
        skew_out(1,2) = -vec_in[0];     

        return skew_out;   
    }

    vec_3_1 somiglianaGravityModel(vec_3_1& ecef_pos)
    {
        vec_3_1 inner;
        double term_x,term_y,term_z,common_term;
        double abs_pos = ecef_pos.norm();

        common_term = 5*pow((ecef_pos[2]/abs_pos),2);
        term_x = (1-common_term)*ecef_pos[0];
        term_y = (1-common_term)*ecef_pos[1];
        term_z = (3-common_term)*ecef_pos[2];

        inner<<term_x,term_y,term_z;

        inner = (1.5*cpp_nav_filt::J2*pow(cpp_nav_filt::Ro,2)/(pow(abs_pos,2)))*inner;
        inner = ecef_pos + inner;
        return -(cpp_nav_filt::mu_g/pow(abs_pos,3))*inner;
    }

    double meridianRadius(double& lat)
    {
        double denom = 1-pow(cpp_nav_filt::e,2)*pow(sin(lat),2);
        return cpp_nav_filt::Ro/sqrt(denom);        
    }

    double geocentricRadius(double& lat)
    {
        // groves p. 71
        double meridian_radius;
        meridian_radius = meridianRadius(lat);

        double sqrt_arg = pow(cos(lat),2) + pow(1-pow(cpp_nav_filt::e,2),2)*pow(sin(lat),2);
        return meridian_radius*sqrt(sqrt_arg);        
    }
} // end of namespace