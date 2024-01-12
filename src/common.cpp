#include "common/common.h"

namespace cpp_nav_filt
{
    Common::Common()
    {
        ones_32_1.setOnes();

        sv_state_.setZero();

        nanEphemerisMap();
        ephem_vect.resize(32);
    }

    Common::~Common()
    {

    }    

    vec_7_1 Common::sendSvStates(const int& sv_in,const double& transmit_time,const double& transit_time)
    {
        desired_sv_ = sv_in;
        T_transmit_ = transmit_time;
        T_transit_ = transit_time;

        calcSvPVStates(sv_state_);

        return sv_state_;
    }

    void Common::setCurrentEphem(const int& sv)
    {
        for(int i = 0;i<ephem_vect.size();i++)
        {
            if(ephem_vect[i].sv == sv)
            {
                std::cout<<ephem_vect[i].sv<<" "<<sv<<std::endl;
                current_ephem_ = ephem_vect[i].ephem_map;
            }
        }
    }

    void Common::calcSvPVStates(vec_7_1& sv_state)
    {
        setCurrentEphem(desired_sv_);

        // Making current ephemerides their own variables
        double T_GD      = current_ephem_.at("T_GD");
        double a_f2      = current_ephem_.at("a_f2");
        double a_f1      = current_ephem_.at("a_f1");
        double a_f0      = current_ephem_.at("a_f0");
        double C_rc      = current_ephem_.at("C_rc");
        double C_rs      = current_ephem_.at("C_rs");
        double C_uc      = current_ephem_.at("C_uc");
        double C_us      = current_ephem_.at("C_us");
        double C_ic      = current_ephem_.at("C_ic");
        double C_is      = current_ephem_.at("C_is");
        double Delta_n   = current_ephem_.at("deltan");
        double M_0       = current_ephem_.at("M_0");
        double e         = current_ephem_.at("e");
        double A         = pow(current_ephem_.at("A"),2);
        double t_oe      = current_ephem_.at("t_oe");
        double Omega_0   = current_ephem_.at("omega_0");
        double i_0       = current_ephem_.at("i_0");
        double omega     = current_ephem_.at("omega");
        double dot_Omega = current_ephem_.at("omegaDot");
        double t_oc      = current_ephem_.at("t_oc");
        double I_dot     = current_ephem_.at("iDot");

        if(T_GD != NAN)
        {

            dt = checkT(T_transmit_ - t_oc);
            
            sv_state(6) = (a_f2*dt + a_f1)*dt + a_f0 - T_GD;

            time = T_transmit_ - sv_state(6);
            tk = checkT(time - t_oe);
            
            double n0 = sqrt(GM/pow(A,3));
            double n = n0 + Delta_n;
            double M = M_0+n * tk;

            M = remainder(M + 2*gps_pi,2*gps_pi);

            double E = M;
            double E_old,dE;

            for(int i = 0;i<11;i++)
            {
                E_old = E;
                E = M + e*sin(E);
                dE = remainder(E - E_old,2*gps_pi);

                if(abs(dE)<1*pow(10,-12))
                {
                    break;
                }
            }

            E = remainder(E + 2*gps_pi,2*gps_pi);

            double dtr = F*e*sqrt(A)*sin(E);
            double nu  = atan2(sqrt(1-pow(e,2))*sin(E),cos(E)-e);
            double phi = nu + omega;
            phi = remainder(phi,2*gps_pi);

            double u = phi + C_uc*cos(2*phi)+C_us*sin(2*phi);
            double r = A*(1-e*cos(E)) + C_rc*cos(2*phi) + C_rs*sin(2*phi);
            double i = i_0 + I_dot * tk + C_ic * cos(2*phi) + C_is*sin(2*phi);

            double Omega = Omega_0 + (dot_Omega - omega_e_dot)*tk - omega_e_dot*t_oe-omega_e_dot*T_transit_;
            Omega = remainder(Omega + 2*gps_pi,2*gps_pi);

            double X = r*cos(u);
            double Y = r*sin(u);

            double x = X*cos(Omega) - Y*cos(i)*sin(Omega);
            double y = X*sin(Omega) + Y*cos(i)*cos(Omega);
            double z = Y*sin(i);

            sv_state(0) = x;
            sv_state(1) = y;
            sv_state(2) = z;
            
            double Edot = (n0 + Delta_n)/(1-e*cos(E));
            double phi_dot = (sqrt(1-pow(e,2))/(1-e*cos(E)))*Edot;
            double u_dot = (1+2*C_us*cos(2*phi)-2*C_uc*sin(2*phi))*phi_dot;
            double r_dot = 2*(C_rs*cos(2*phi)-C_rc*sin(2*phi))*phi_dot + A*e*sin(E)*Edot;
            double i_dot = 2*(C_is*cos(2*phi)-C_ic*sin(2*phi))*phi_dot + I_dot;
            double Xdot = r_dot*cos(u) - r*sin(u)*u_dot;
            double Ydot = r_dot*sin(u) + r*cos(u)*u_dot;
            double Omega_dot  = dot_Omega - omega_e_dot;

            sv_state(3) = Xdot*cos(Omega) - Ydot*cos(i)*sin(Omega) + Y*sin(i)*sin(Omega)*i_dot - y*Omega_dot;
            sv_state(4) = Xdot*sin(Omega) + Ydot*cos(i)*cos(Omega) - Y*sin(i)*cos(Omega)*i_dot + x*Omega_dot;
            sv_state(5) = Ydot*sin(i) + Y*cos(i)*i_dot;

            sv_state(6) = (a_f2*dt + a_f1)*dt + a_f0 - T_GD + dtr;
        }
        else
        {
            sv_state.setZero();
        }
       
       nanEphemerisMap();
    }

    double Common::checkT(double time)
    {
        /* Kai Borre 04-01-96
        Copyright (c) by Kai Borre

        CVS record:
        Id: check_t.m,v 1.1.1.1.2.4 2006/08/22 13:45:59 dpl Exp
        */

        double dt = time;

        if(dt>half_week)
        {
            dt = time - 2*half_week;
        }
        else if(dt<-1*half_week)
        {
            dt = time + 2*half_week;
        }

        return dt;
    }

    void Common::nanEphemerisMap()
    {
        current_ephem_["T_GD"]      = NAN;
        current_ephem_["t_oc"]      = NAN;
        current_ephem_["a_f2"]      = NAN;
        current_ephem_["a_f1"]      = NAN;
        current_ephem_["a_f0"]      = NAN;
        current_ephem_["C_rc"]      = NAN;
        current_ephem_["C_rs"]      = NAN;
        current_ephem_["C_uc"]      = NAN;
        current_ephem_["C_us"]      = NAN;
        current_ephem_["C_ic"]      = NAN;
        current_ephem_["C_is"]      = NAN;
        current_ephem_["Delta_n"]   = NAN;
        current_ephem_["M_0"]       = NAN;
        current_ephem_["e"]         = NAN;
        current_ephem_["A"]         = NAN;
        current_ephem_["t_oe"]      = NAN;
        current_ephem_["Omega_0"]   = NAN;
        current_ephem_["i_0"]       = NAN;
        current_ephem_["omega"]     = NAN;
        current_ephem_["dot_Omega"] = NAN;
        current_ephem_["I_dot"]     = NAN;
    }

    void Common::sendUnitVectors(vec_3_1& X_hat,double& clk,Eigen::MatrixXd& SvPVT,Eigen::MatrixXd& H)
    {
        // sending position states to state estimate
        pos_ = X_hat;
        clk_ = clk;
        
        num_sv_ = SvPVT.rows();
        
        sv_pvt_.resize(num_sv_,7);
        H_.resize(num_sv_,4); // resize H_ to be the size it needs to be

        sv_pvt_ = SvPVT;
                
        calcUnitVectors();

        H = H_;
    }

    void Common::sendMeasEst(vec_3_1& pos,vec_3_1& vel,double& clk,double& clk_drift,Eigen::MatrixXd& SvPVT,Eigen::MatrixXd& Yhat)
    {
        pos_ = pos;
        vel_ = vel;
        clk_ = clk;
        clk_drift_ = clk_drift;

        num_sv_ = SvPVT.rows();

        H_.resize(num_sv_,4);
        sv_pvt_.resize(num_sv_,7);
        Yhat_.resize(num_sv_*2,1);

        sv_pvt_ = SvPVT;

        calcUnitVectors();
        
        calcMeasEst();

        Yhat = Yhat_;

    }

    void Common::calcUnitVectors()
    {
        for(int i = 0;i<num_sv_;i++)
        {
            calcPsr(i);
            
            H_(i,0) = -x_comp/psr_hat;
            H_(i,1) = -y_comp/psr_hat;
            H_(i,2) = -z_comp/psr_hat;
            H_(i,3) = 1;
        }
    }

    void Common::calcMeasEst()
    {
        for(int i = 0;i<num_sv_;i++)
        {
            calcPsr(i);
            calcPsrRate(i);

            Yhat_(i,0) = psr_hat;
            Yhat_(i+num_sv_,0) = psr_rate_hat;
        }
    }

    void Common::calcPsr(double sv_id)
    {        
        x_comp = sv_pvt_(sv_id,0) - pos_[0];
        y_comp = sv_pvt_(sv_id,1) - pos_[1];
        z_comp = sv_pvt_(sv_id,2) - pos_[2];

        psr_hat = sqrt(pow(x_comp,2) + pow(y_comp,2) + pow(z_comp,2)) + clk_;
    }

    void Common::calcPsrRate(double sv_id)
    {
        vec_3_1 u; // unit vector to the sv_id satellite
        vec_3_1 relative_velocity;

        x_comp_vel = sv_pvt_(sv_id,3) - vel_[0];
        y_comp_vel = sv_pvt_(sv_id,4) - vel_[1];
        z_comp_vel = sv_pvt_(sv_id,5) - vel_[2];
        
        relative_velocity << x_comp_vel,y_comp_vel,z_comp_vel;

        u = H_.block<1,3>(sv_id,0);
        
        psr_rate_hat = (-u.dot(relative_velocity)) + clk_drift_;
    }

    //========= Frame Conversions ==========//

    void Common::convertECEF2LLA(vec_3_1& ecef_pos,vec_3_1& lla_pos,WgsConversions& frame_conversions)
    {
        eigen2array(ecef_pos_,ecef_pos);
        frame_conversions.xyz2lla(lla_pos_,ecef_pos_);
        array2eigen(lla_pos,lla_pos_);
    }

    void Common::convertLLA2ECEF(vec_3_1& lla_pos,vec_3_1& ecef_pos,WgsConversions& frame_conversions)
    {
        eigen2array(lla_pos_,lla_pos);
        frame_conversions.lla2xyz(ecef_pos_,lla_pos_);
        array2eigen(ecef_pos,ecef_pos_);
    }

    void Common::eigen2array(double array[3],vec_3_1& eigen)
    {
        array[0] = eigen[0];
        array[1] = eigen[1];
        array[2] = eigen[2];
    }

    void Common::array2eigen(vec_3_1& eigen,double array[3])
    {
        eigen[0] = array[0];
        eigen[1] = array[1];
        eigen[2] = array[2];
    }
}// end of namespace
