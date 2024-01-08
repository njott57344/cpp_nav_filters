#include "common/common.h"

namespace cpp_nav_filt
{
    Common::Common()
    {
        sv_ephem.setZero();
        sv_state_.setZero();
    }

    Common::~Common()
    {

    }    

    void Common::receiveSvEphem(vec_1_27& ephem_in,const int& sv_in)
    {
        // std::cout<<ephem_in<<std::endl;
        sv_ephem.block<1,27>(sv_in-1,0) = ephem_in;
    }

    void Common::sendSvEphem(vec_1_27& ephem_out,const int& desired_sv)
    {
        ephem_out = sv_ephem.block<1,27>(desired_sv-1,0);
    }

    vec_7_1 Common::sendSvStates(const int& sv_in,const double& transmit_time,const double& transit_time)
    {
        desired_sv_ = sv_in;
        T_transmit_ = transmit_time;
        T_transit_ = transit_time;

        calcSvPVStates(sv_state_);

        return sv_state_;
    }

    void Common::setCurrentEphem(vec_1_27& ephem,const int& sv)
    {
        ephem = sv_ephem.block<1,27>(sv-1,0);
    }

    void Common::calcSvPVStates(vec_7_1& sv_state)
    {
        vec_1_27 current_ephem;
        setCurrentEphem(current_ephem,desired_sv_);
        
        // Making current ephemerides their own variables
        double T_GD      = current_ephem(23);
        double t_oc      = current_ephem(22);
        double a_f2      = current_ephem(26);
        double a_f1      = current_ephem(25);
        double a_f0      = current_ephem(24);
        double C_rc      = current_ephem(13);
        double C_rs      = current_ephem(14);
        double C_uc      = current_ephem(11);
        double C_us      = current_ephem(12);
        double C_ic      = current_ephem(15);
        double C_is      = current_ephem(16);
        double Delta_n   = current_ephem(7);
        double M_0       = current_ephem(8);
        double e         = current_ephem(9);
        double A         = current_ephem(6);
        double t_oe      = current_ephem(5);
        double Omega_0   = current_ephem(19);
        double i_0       = current_ephem(17);
        double omega     = current_ephem(10);
        double dot_Omega = current_ephem(20);
        double I_dot     = current_ephem(18);
                
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

    void Common::sendUnitVectors(vec_8_1& X_hat,Eigen::MatrixXd& SvPVT,Eigen::MatrixXd& H)
    {
        // sending position states to state estimate
        x_hat_ = X_hat;

        num_sv_ = SvPVT.rows();
        
        sv_pos_.resize(num_sv_,7);
        H_.resize(num_sv_,4); // resize H_ to be the size it needs to be

        sv_pos_ = SvPVT;
                
        calcUnitVectors();

        H = H_;
    }

    void Common::sendMeasEst(vec_8_1& X_hat,Eigen::MatrixXd& SvPVT,Eigen::MatrixXd& Yhat)
    {
        x_hat_ = X_hat;

        num_sv_ = SvPVT.rows();

        H_.resize(num_sv_,4);
        sv_pos_.resize(num_sv_,7);
        Yhat_.resize(num_sv_*2,1);

        sv_pos_ = SvPVT;

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
        x_comp = sv_pos_(sv_id,0) - x_hat_[0];
        y_comp = sv_pos_(sv_id,1) - x_hat_[1];
        z_comp = sv_pos_(sv_id,2) - x_hat_[2];

        psr_hat = sqrt(pow(x_comp,2) + pow(y_comp,2) + pow(z_comp,2)) + x_hat_[3];
    }

    void Common::calcPsrRate(double sv_id)
    {
        vec_3_1 u; // unit vector to the sv_id satellite
        vec_3_1 relative_velocity;

        x_comp_vel = sv_pos_(sv_id,3) - x_hat_[4];
        y_comp_vel = sv_pos_(sv_id,4) - x_hat_[5];
        z_comp_vel = sv_pos_(sv_id,5) - x_hat_[6];
        
        relative_velocity << x_comp_vel,y_comp_vel,z_comp_vel;

        u = H_.block<1,3>(sv_id,0);

        psr_rate_hat = (u.dot(relative_velocity) + x_hat_[7]);        
    }

}// end of namespace
