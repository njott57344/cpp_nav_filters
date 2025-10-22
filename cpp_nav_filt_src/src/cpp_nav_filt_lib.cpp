#include "cpp_nav_filt_lib.h"

namespace cpp_nav_filt
{

    // =========== Measurement Functions ======== //

    Eigen::MatrixXd calcPsr(Eigen::MatrixXd& SvPVT,vec_3_1& ecef_pos,double& clk_b)
    {
        int num_sv = SvPVT.rows();

        double x_comp,y_comp,z_comp;

        Eigen::MatrixXd psr_hat;
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

        Eigen::MatrixXd psr_hat,H;
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

        Eigen::MatrixXd psr_rate_hat,H;
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

        Eigen::MatrixXd psr,psr_rate,Y_hat;
        psr.resize(num_sv,1);
        psr_rate.resize(num_sv,1);
        Y_hat.resize(2*num_sv,1);

        psr = calcPsr(SvPVT,ecef_pos,clk_b);
        psr_rate = calcPsrRate(SvPVT,ecef_pos,ecef_vel,clk_b,clk_d);

        Y_hat.block(0,0,num_sv,1) = psr;
        Y_hat.block(num_sv,0,num_sv,1) = psr_rate;

        return Y_hat;
    }

    // ======== Nav Functions ======== //

    mat_3_3 eul2Rotm(const vec_3_1& euler_angles)
    {
        // assumes euler angles are ordered roll pitch yaw (and in rad)
        // see groves p. 38

        // groves 2.24
        vec_3_1 c;
        vec_3_1 s;

        c << std::cos(euler_angles[0]),std::cos(euler_angles[1]),std::cos(euler_angles[2]);
        s << std::sin(euler_angles[0]),std::sin(euler_angles[1]),std::sin(euler_angles[2]);

        mat_3_3 C;
        C(0,0) =  c(1)*c(2); // cos(pitch)cos(yaw)
        C(1,0) =  c(1)*s(2); // cos(pitch)sin(yaw)
        C(2,0) = -s(1);     // -sin(pitch)

        C(0,1) = -c(0)*s(2)+s(0)*s(1)*c(2); // -cos(roll)sin(yaw)+sin(roll)sin(pitch)cos(yaw)
        C(1,1) =  c(0)*c(2)+s(0)*s(1)*s(2); //  cos(roll)cos(yaw)+sin(roll)sin(pitch)sin(yaw)
        C(2,1) =  s(0)*c(1);                //  sin(roll)cos(pitch)
    
        C(0,2) =  s(0)*s(2)+c(0)*s(1)*c(2); //  sin(roll)sin(yaw)+cos(roll)sin(pitch)cos(yaw)
        C(1,2) = -s(0)*c(2)+c(0)*s(1)*s(2); // -sin(roll)cos(yaw)+cos(roll)sin(pitch)sin(yaw)
        C(2,2) =  c(0)*c(1);                //  cos(roll)cos(pitch)

        return C; // rotation body to nav
    }

    mat_3_3 eul2EcefDCM(vec_3_1& euler_angles,vec_3_1& lla_pos)
    {
        // return rotation matrix from body to ecef
        mat_3_3 C_n_b,C_b_n,C_b_e,C_n_e;

        C_n_b = eul2Rotm(euler_angles); // DCM from local nav to body
        C_b_n = C_n_b.transpose(); // DCM from body to local nav

        C_n_e = ned2ecefDCM(lla_pos); // DCM from local nav to ecef
        C_b_e = C_n_e*C_b_n;
        return C_b_e;
    }

    vec_3_1 rotm2Eul(mat_3_3& C)
    {
        // groves 2.25
        vec_3_1 euler_angles;

        euler_angles[0] =  std::atan2(C(2,1),C(2,2)); // roll
        euler_angles[1] = -std::asin(C(2,0));        // pitch
        euler_angles[2] =  std::atan2(C(1,0),C(0,0)); // yaw
    
        return euler_angles;
    }

    vec_3_1 ecefDCM2EulerAngles(mat_3_3& C_be,vec_3_1& lla_pos)
    {
        vec_3_1 eul_angles;
        mat_3_3 C_ne = ned2ecefDCM(lla_pos);
        mat_3_3 C_bn = C_be*C_ne.transpose();

        eul_angles = rotm2Eul(C_bn);

        return eul_angles;
    }

    vec_4_1 eul2q(const vec_3_1 & eul)
    {
        double ss1 = sin(eul(0) / 2);
        double ss2 = sin(eul(1) / 2);
        double ss3 = sin(eul(2) / 2);
        double cc1 = cos(eul(0) / 2);
        double cc2 = cos(eul(1) / 2);
        double cc3 = cos(eul(2) / 2);

        vec_4_1 q_ret = {cc1 * cc2 * cc3 + ss1 * ss2 * ss3, ss1 * cc2 * cc3 - cc1 * ss2 * ss3,
                                 cc1 * ss2 * cc3 + ss1 * cc2 * ss3, cc1 * cc2 * ss3 - ss1 * ss2 * cc3};

        return q_ret;
    }

    vec_3_1 q2eul(const vec_4_1 & q)
    {
        vec_3_1 eul = {
            std::atan2(2 * (q(0) * q(1) + q(2) * q(3)), 1 - 2 * (q(1) * q(1) + q(2) * q(2))),
            std::asin(2 * (q(0) * q(2) - q(1) * q(3))),
            std::atan2(2 * (q(0) * q(3) + q(1) * q(2)), 1 - 2 * (q(2) * q(2) + q(3) * q(3)))};

        return eul;
    }

    vec_4_1 qMult(const vec_4_1 & p, const vec_4_1 & q)
    {
        Eigen::Vector4d q_ret = {p(0) * q(0) - p(1) * q(1) - p(2) * q(2) - p(3) * q(3),
                                p(0) * q(1) + p(1) * q(0) + p(2) * q(3) - p(3) * q(2),
                                p(0) * q(2) - p(1) * q(3) + p(2) * q(0) + p(3) * q(1),
                                p(0) * q(3) + p(1) * q(2) - p(2) * q(1) + p(3) * q(0)};
        return q_ret;
    }

    vec_4_1 qNormalize(const vec_4_1 & q)
    {
        double quat_norm = qNorm(q);
        vec_4_1 q_ret = q / quat_norm;
        return q_ret;
    }

    double qNorm(const vec_4_1 & q)
    {
        double squared_sum = q(0) * q(0) + q(1) * q(1) + q(2) * q(2) + q(3) * q(3);
        double ret = sqrt(squared_sum);
        return ret;
    }

    mat_3_3 q2DCM(const vec_4_1 & q)
    {
        mat_3_3 DCM_ret;
        DCM_ret << (q(0) * q(0) + q(1) * q(1) - q(2) * q(2) - q(3) * q(3)),
            2 * (q(1) * q(2) - q(0) * q(3)), 2 * (q(1) * q(3) + q(0) * q(2)),
            2 * (q(1) * q(2) + q(0) * q(3)), (q(0) * q(0) - q(1) * q(1) + q(2) * q(2) - q(3) * q(3)),
            2 * (q(2) * q(3) - q(0) * q(1)), 2 * (q(1) * q(3) - q(0) * q(2)),
            2 * (q(2) * q(3) + q(0) * q(1)), (q(0) * q(0) - q(1) * q(1) - q(2) * q(2) + q(3) * q(3));
        return DCM_ret;       
    }

    mat_3_3 makeSkewSymmetic(vec_3_1& vec_in)
    {
        mat_3_3 skew_out;

        skew_out << 0,        -vec_in[2], vec_in[1],
                    vec_in[2], 0,        -vec_in[0],
                   -vec_in[1], vec_in[0], 0;
                   
        return skew_out;   
    }

    vec_3_1 ecefGravitation(vec_3_1& ecef_pos)
    {
        vec_3_1 inner;
        double term_x,term_y,term_z,common_term;
        double abs_pos = ecef_pos.norm();
        double pos_3,R02,pos_2;
        double outer_lump,inner_lump;

        common_term = 5*pow((ecef_pos[2]/abs_pos),2);
        term_x = (1-common_term)*ecef_pos[0];
        term_y = (1-common_term)*ecef_pos[1];
        term_z = (3-common_term)*ecef_pos[2];

        inner<<term_x,term_y,term_z;
    
        pos_3 = pow(abs_pos,3);
        pos_2 = pow(abs_pos,2);
        R02 = pow(Ro,2);

        outer_lump = mu_g/pos_3; // lumped together constants 
        inner_lump = 1.5*J2*R02/pos_2; // lumped together constants

        return -outer_lump*(ecef_pos + inner_lump*inner);
    }

    vec_3_1 ecefGravity(vec_3_1& ecef_pos)
    {
        vec_3_1 ecef_gravitation,ecef_gravity;
        mat_3_3 C;
        
        C.setZero();
        C(0,0) = 1;
        C(1,1) = 1;

        ecef_gravitation = ecefGravitation(ecef_pos);
        
        if(ecef_pos.norm() == 0)
        {
            ecef_gravity.setZero();
        }
        else
        {
            ecef_gravity = ecef_gravitation + w_e*w_e*C*ecef_pos;
        }
        return  ecef_gravity;
    }

    double meridianRadiusOfCurvature(const double& lat)
    {
        double num,den;

        num = cpp_nav_filt::Ro*(1-pow(cpp_nav_filt::e,2));
        den = pow(1-pow(cpp_nav_filt::e,2)*pow(sin(lat),2),1.5);

        return num/den;
    }

    double transverseRadiusOfCurvature(const double& lat)
    {
        double denom = 1-pow(cpp_nav_filt::e,2)*pow(sin(lat),2);
        return cpp_nav_filt::Ro/sqrt(denom);        
    }

    double geocentricRadius(const double& lat)
    {
        double radius_of_curvature;
        radius_of_curvature = transverseRadiusOfCurvature(lat);

        double sqrt_arg = pow(cos(lat),2) + pow(1-pow(cpp_nav_filt::e,2),2)*pow(sin(lat),2);
        return radius_of_curvature*sqrt(sqrt_arg);        
    }

    vec_2_1 levelInsAccel(Eigen::MatrixXd& f_ib_b)
    {
        int num_meas = f_ib_b.rows();
        vec_2_1 roll_pitch;

        vec_1_3 samp_mean = f_ib_b.colwise().mean();

        roll_pitch[0] = std::atan2(-samp_mean[1],-samp_mean[2]);
        roll_pitch[1] = std::atan(samp_mean[0]/sqrt(pow(samp_mean[1],2) + pow(samp_mean[2],2)));
        return roll_pitch;
    }   

    // ========= Frame Transformations ======== //

    vec_3_1 convertEnu2NED(vec_3_1& enu_pos)
    {
        mat_3_3 C;
        C.setZero();

        C(1,0) = 1;
        C(0,1) = 1;
        C(2,2) = -1;
        
        return C*enu_pos;
    }

    vec_3_1 ecef2llaPos(vec_3_1& ecef_pos)
    {
        // groves p 61-62
        vec_3_1 lla;
        lla.setOnes();

        double x,y,z;
        x = ecef_pos[0];
        y = ecef_pos[1];
        z = ecef_pos[2];

        double lat,lon,alt;
        double templat,tempalt;

        double num,den;

        double rho_sqr;
        double dlat = 1,dalt = 1; // init to 1 so while loop runs
        int iter = 0;

        double Re; // transverse radius of curavture
        double e2 = cpp_nav_filt::e*cpp_nav_filt::e;

        // solving for lon (groves 2.113)
        lon = std::atan2(y,x);

        // have to iterate to solve for lat and alt
        // init to spherical earth 
        rho_sqr = x*x + z*z;
        templat = atan2(y,sqrt(rho_sqr));
        tempalt = sqrt(rho_sqr + y*y) - cpp_nav_filt::A;

        while (iter<50 && dlat>1e-12 && dalt>1e-12)
        {
            Re = transverseRadiusOfCurvature(templat);

            // lat
            num = (z*(Re + tempalt));
            den = sqrt(x*x + y*y)*((1-e2)*Re + tempalt);
            lat = atan(num/den);

            // alt
            num = sqrt(x*x + y*y);
            den = std::cos(templat);
            alt = num/den - Re;

            // while loop stuff
            dlat = abs(lat - templat);
            dalt = abs(alt - tempalt);
            templat = lat;
            tempalt = alt;
            iter++;
        }

        if(iter == 49)
        {
            std::cout<<"WARNING: ecef2llaPos could not converge!"<<std::endl;
            lla = NAN*lla;
        }
        else
        {
            lla<<lat*cpp_nav_filt::R2D,lon*cpp_nav_filt::R2D,alt;
        }

        return lla;
    }

    vec_3_1 lla2ecefPos(vec_3_1& lla_pos)
    {
        // see groves 2.112
        double lat,lon,alt;
        double x,y,z;
        vec_3_1 ecef_pos;
        
        lat = lla_pos(0)*cpp_nav_filt::D2R;
        lon = lla_pos(1)*cpp_nav_filt::D2R;
        alt = lla_pos(2);

        double Re = transverseRadiusOfCurvature(lat);
        double e2 = pow(cpp_nav_filt::e,2);

        x = (Re + alt)*std::cos(lat)*std::cos(lon);
        y = (Re + alt)*std::cos(lat)*std::sin(lon);
        z = ((1 - e2)*Re + alt)*sin(lat);

        ecef_pos<<x,y,z;
        return ecef_pos;
    }

    vec_3_1 ecef2nedPos(vec_3_1& ecef_pos,vec_3_1& ref_lla)
    {
        double lat_0,lon_0,alt_0; // origin of NED frame
        lat_0 = ref_lla(0)*cpp_nav_filt::D2R;
        lon_0 = ref_lla(1)*cpp_nav_filt::D2R;
        alt_0 = ref_lla(2);
        
        mat_3_3 Cen; // rotation matrix ecef to ned
        vec_3_1 ref_ecef;
                
        ref_ecef = cpp_nav_filt::lla2ecefPos(ref_lla);

        Cen << -std::sin(lat_0)*std::cos(lon_0), -std::sin(lat_0)*std::sin(lon_0),  std::cos(lat_0),
               -std::sin(lon_0),                  std::cos(lon_0),                  0,
               -std::cos(lat_0)*std::cos(lon_0), -std::cos(lat_0)*std::sin(lon_0), -std::sin(lat_0);


        return Cen*(ecef_pos - ref_ecef);
    }

    vec_3_1 ned2ecefPos(vec_3_1& ned_pos,vec_3_1& ref_lla)
    {
        double lat_0,lon_0,alt_0; // origin of NED frame
        lat_0 = ref_lla(0)*cpp_nav_filt::D2R;
        lon_0 = ref_lla(1)*cpp_nav_filt::D2R;
        alt_0 = ref_lla(2);
        
        mat_3_3 Cen; // rotation matrix ecef to ned
        mat_3_3 Cne; // rotation matrix ned to ecef
        vec_3_1 ref_ecef;
        
        ref_ecef = cpp_nav_filt::lla2ecefPos(ref_lla);

        Cen << -std::sin(lat_0)*std::cos(lon_0), -std::sin(lat_0)*std::sin(lon_0),  std::cos(lat_0),
               -std::sin(lon_0),                  std::cos(lon_0),                  0,
               -std::cos(lat_0)*std::cos(lon_0), -std::cos(lat_0)*std::sin(lon_0), -std::sin(lat_0);
        
        Cne = Cen.transpose();

        return ref_ecef + Cne*ned_pos;
    }

    vec_3_1 ecef2enuPos(vec_3_1& ecef_pos,vec_3_1& ref_lla)
    {
        vec_3_1 ned_pos,enu_pos;
        ned_pos = cpp_nav_filt::ecef2nedPos(ecef_pos,ref_lla);
        enu_pos = cpp_nav_filt::ned2enuPos(ned_pos);
        return enu_pos;
    }

    vec_3_1 enu2ecefPos(vec_3_1& enu_pos,vec_3_1& ref_lla)
    {
        vec_3_1 ned_pos,ecef_pos;
        ned_pos = cpp_nav_filt::enu2nedPos(enu_pos);
        ecef_pos = cpp_nav_filt::ned2ecefPos(ned_pos,ref_lla);
        return ecef_pos;
    }

    vec_3_1 enu2nedPos(vec_3_1& enu_pos)
    {
        mat_3_3 Cen; // rotation matrix enu to ned
        Cen << 0, 1, 0,
               1, 0, 0,
               0, 0,-1;
        return Cen*enu_pos;
    }

    vec_3_1 ned2enuPos(vec_3_1& ned_pos)
    {
        mat_3_3 Cne; // rotation matrix ned to enu
        Cne << 0, 1, 0,
               1, 0, 0,
               0, 0,-1;
        return Cne*ned_pos;
    }

    vec_3_1 lla2nedPos(vec_3_1& lla_pos,vec_3_1& ref_lla)
    {
        vec_3_1 ecef_pos,ned_pos;
        ecef_pos = cpp_nav_filt::lla2ecefPos(lla_pos);
        ned_pos = cpp_nav_filt::ecef2nedPos(ecef_pos,ref_lla);
        return ned_pos;
    }
    
    vec_3_1 ned2llaPos(vec_3_1& ned_pos,vec_3_1& ref_lla)
    {
        vec_3_1 ecef_pos,lla_pos;
        ecef_pos = cpp_nav_filt::ned2ecefPos(ned_pos,ref_lla);
        lla_pos = cpp_nav_filt::ecef2llaPos(ecef_pos);
        return lla_pos;
    }

    vec_3_1 lla2enuPos(vec_3_1 lla_pos,vec_3_1& ref_lla)
    {
        vec_3_1 ecef_pos,enu_pos;
        ecef_pos = cpp_nav_filt::lla2ecefPos(lla_pos);
        enu_pos = cpp_nav_filt::ecef2enuPos(ecef_pos,ref_lla);
        return enu_pos;
    }

    vec_3_1 enu2llaPos(vec_3_1 enu_pos,vec_3_1& ref_lla)
    {
        vec_3_1 ecef_pos,lla_pos;
        ecef_pos = cpp_nav_filt::enu2ecefPos(enu_pos,ref_lla);
        lla_pos = cpp_nav_filt::ecef2llaPos(ecef_pos);
        return lla_pos;
    }
    
    vec_3_1 ecef2nedVel(vec_3_1& ecef_vel,vec_3_1& ref_lla)
    {
        double lat_0,lon_0,alt_0; // origin of NED frame
        lat_0 = ref_lla(0)*cpp_nav_filt::D2R;
        lon_0 = ref_lla(1)*cpp_nav_filt::D2R;
        alt_0 = ref_lla(2);
        
        mat_3_3 Cen; // rotation matrix ecef 2 ned
        vec_3_1 ned_vel;

        Cen << -std::sin(lat_0)*std::cos(lon_0), -std::sin(lat_0)*std::sin(lon_0),  std::cos(lat_0),
               -std::sin(lon_0),                  std::cos(lon_0),                  0,
               -std::cos(lat_0)*std::cos(lon_0), -std::cos(lat_0)*std::sin(lon_0), -std::sin(lat_0);

        ned_vel = Cen*ecef_vel;
        return ned_vel;
    }

    vec_3_1 ned2ecefVel(vec_3_1& ned_vel,vec_3_1& ref_lla)
    {
        double lat_0,lon_0,alt_0; // origin of NED frame
        lat_0 = ref_lla(0)*cpp_nav_filt::D2R;
        lon_0 = ref_lla(1)*cpp_nav_filt::D2R;
        alt_0 = ref_lla(2);
        
        mat_3_3 Cen,Cne; // rotation matrix ecef 2 ned
        vec_3_1 ecef_vel;
        vec_3_1 temp;

        Cen << -std::sin(lat_0)*std::cos(lon_0), -std::sin(lat_0)*std::sin(lon_0),  std::cos(lat_0),
               -std::sin(lon_0),                  std::cos(lon_0),                  0,
               -std::cos(lat_0)*std::cos(lon_0), -std::cos(lat_0)*std::sin(lon_0), -std::sin(lat_0);

        Cne = Cen.transpose();

        ecef_vel = Cne*ned_vel;
        
        return ecef_vel;
    }

    vec_3_1 ecef2enuVel(vec_3_1& ecef_vel,vec_3_1& ref_lla)
    {
        vec_3_1 ned_vel,enu_vel;
        ned_vel = ecef2nedVel(ecef_vel,ref_lla);
        enu_vel = ned2enuPos(ned_vel);
        return enu_vel;
    }

    vec_3_1 enu2ecefVel(vec_3_1& enu_vel,vec_3_1& ref_lla)
    {
        vec_3_1 ned_vel,ecef_vel;
        ned_vel = enu2nedPos(enu_vel);
        ecef_vel = ned2ecefVel(ned_vel,ref_lla);    
        return ecef_vel;
    }

    mat_3_3 ecef2nedDCM(vec_3_1& lla_pos)
    {
        mat_3_3 C;

        double lat,lon,alt;
        double clat,slat,clon,slon;
        lat = lla_pos[0]*D2R;
        lon = lla_pos[1]*D2R;

        clat = std::cos(lat);
        slat = std::sin(lat);
        clon = std::cos(lon);
        slon = std::sin(lon);

        C << -slat*clon, -slat*slon, clat,
             -slon,       clon,      0,
             -clat*clon, -clat*slon, -slat;

        return C;
    }
    
    mat_3_3 ned2ecefDCM(vec_3_1& lla_pos)
    {
        mat_3_3 C;

        double lat,lon;
        double clat,slat,clon,slon; // sin and cos of lat/lon
        lat = lla_pos[0]*D2R;
        lon = lla_pos[1]*D2R;

        clat = std::cos(lat);
        slat = std::sin(lat);
        clon = std::cos(lon);
        slon = std::sin(lon);

        C << -slat*clon, -slat*slon, clat,
             -slon,       clon,      0,
             -clat*clon, -clat*slon, -slat;

        return C.transpose();
    }

    mat_3_3 ecef2enuDCM(vec_3_1& lla_pos)
    {
        mat_3_3 C;

        double lat,lon,alt;
        double clat,slat,clon,slon;

        lat = lla_pos[0]*D2R;
        lon = lla_pos[1]*D2R;
        alt = lla_pos[2];

        clat = std::cos(lat);
        slat = std::sin(lat);
        clon = std::cos(lon);
        slon = std::sin(lon);

        C << -slon,       clon,      0,
             -slat*clon, -slat*slon, clat,
              clat*clon,  clat*slon, slat;

        return C;
    }

    mat_3_3 enu2ecefDCM(vec_3_1& lla_pos)
    {
        mat_3_3 C; 
        
        double lat,lon,alt;
        double clat,slat,clon,slon;
        lat = lla_pos[0]*D2R;
        lon = lla_pos[1]*D2R;
        alt = lla_pos[2];

        clat = std::cos(lat);
        slat = std::sin(lat);
        clon = std::cos(lon);
        slon = std::sin(lon);

        C << -slon,       clon,      0,
             -slat*clon, -slat*slon, clat,
              clat*clon,  clat*slon, slat;

        return C.transpose();
    }

    vec_3_1 transportRate(const double& ve, const double & vn,const double & h, const double & lat)
    {
        vec_3_1 transport_rate;
        double Re = transverseRadiusOfCurvature(lat);
        double Rn = meridianRadiusOfCurvature(lat);
        double tL = std::tan(lat);

        transport_rate(0) = ve / (Re + h);
        transport_rate(1) = -vn / (Rn + h);
        transport_rate(2) = (-ve * tL) / (Re + h);

        return transport_rate;
    }

    vec_3_1 navFrameRotationRate(const double & lat)
    {
        vec_3_1 earth_rotation_rate;

        earth_rotation_rate(0) = cpp_nav_filt::w_e * std::cos(lat);
        earth_rotation_rate(1) = 0;
        earth_rotation_rate(2) = cpp_nav_filt::w_e * std::sin(lat);

        return earth_rotation_rate;
    }

    mat_3_3 normalizeDCM(mat_3_3& dcm_in)
    {
        mat_3_3 normalized_dcm;

        vec_3_1 x,y,z;
        vec_3_1 x_ort,y_ort,z_ort;
        vec_3_1 x_new,y_new,z_new;
        double error;

        x = dcm_in.block(0,0,3,1);
        y = dcm_in.block(0,1,3,1);
        z = dcm_in.block(0,2,3,1);

        error = x.dot(y);

        x_ort = x - (0.5*error)*y;
        y_ort = y - (0.5*error)*x;
        z_ort = x_ort.cross(y_ort);

        x_new = 0.5*(3-x_ort.dot(x_ort))*x_ort;
        y_new = 0.5*(3-y_ort.dot(y_ort))*y_ort;
        z_new = 0.5*(3-z_ort.dot(z_ort))*z_ort;

        normalized_dcm.block(0,0,3,1) = x_new;
        normalized_dcm.block(0,1,3,1) = y_new;
        normalized_dcm.block(0,2,3,1) = z_new;

        return normalized_dcm;
    }

} // end of namespace