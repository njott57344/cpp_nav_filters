#include "sv_manager/sv_manager.h"

namespace cpp_nav_filt
{
    SvManager::SvManager()
    {
        ones_32_1.setOnes();
        ones_3_1.setOnes();

        sv_state_.setZero();

        nanEphemerisMap();
        ephem_vect.resize(32);
    }

    SvManager::~SvManager()
    {

    }    

    vec_7_1 SvManager::sendSvStates(const int& sv_in,const double& transmit_time,const double& transit_time)
    {
        desired_sv_ = sv_in;
        T_transmit_ = transmit_time;
        T_transit_ = transit_time;

        calcSvPVStates(sv_state_);

        return sv_state_;
    }

    // ============= SV Ephemeris and SV PVT State Calc ============= //
    void SvManager::setCurrentEphem(const int& sv)
    {
        for(int i = 0;i<ephem_vect.size();i++)
        {
            if(ephem_vect[i].sv == sv)
            {
                current_ephem_ = ephem_vect[i].ephem_map;
            }
        }
    }

    void SvManager::calcSvPVStates(vec_7_1& sv_state)
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

    double SvManager::checkT(double time)
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

    void SvManager::nanEphemerisMap()
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
}// end of namespace
