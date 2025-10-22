#include "tightly_coupled_ins.h"

namespace cpp_nav_filt
{

    tightlyCoupledNavigator::tightlyCoupledNavigator()
    {
        P_.setZero();

        P_.diagonal() << 9.0, 9.0, 9.0, 0.01, 0.01, 0.01, 0.0025, 0.0025, 0.0025, 0.04, 0.04, 0.04,
                         0.0025, 0.0025, 0.0025, 1,0.1; 

        clk_b_ = 0;
        clk_d_ = 0;

        bg_.setZero();
        ba_.setZero();
    }

    tightlyCoupledNavigator::~tightlyCoupledNavigator()
    {

    }

    bool tightlyCoupledNavigator::Dynamics(const vec_3_1 & wb, const vec_3_1& fb, const vec_3_1& dt)
    {
        vec_3_1 w_ib_b = wb - bg_;
        vec_3_1 f_ib_b = fb - ba_;

        // Functions of Latitude
        double sL = std::sin(phi_);
        double cL = std::cos(phi_);
        double tL = sL / cL;
        double sL2 = sL * sL;
        double cL2 = cL * cL;

        // Radii of curvature
        double X1ME2 = 1.0 - e2;
        double tmp = 1.0 - e2 * sL2;
        double sqtmp = std::sqrt(tmp);
        double Re = Ro / sqtmp;
        double Rn = Ro * X1ME2 / (tmp * tmp / sqtmp);
        double R_es_e = Re * std::sqrt(cL2 + sL2 * X1ME2 * X1ME2);

        // Adjusted heights
        double hn = h_ + Rn;
        double he = h_ + Re;
        double hn2 = hn * hn;
        double he2 = he * he;

        // Somigliana gravity model
        double g0 = 9.7803253359 * ((1.0 + 0.001931853 * sL2) / std::sqrt(1.0 - e2 * sL2));
        Eigen::Vector3d g{-8.08e-9 * h_ * std::sin(2.0 * phi_), 0.0,
                        g0 * (1.0 -
                                (2.0 / Ro * h_) * (1.0 + flattening * (1.0 - 2.0 * sL2) +
                                                (std::pow(w_e * Ro, 2.0) * Rp / mu_g)) +
                                (3.0 * std::pow(h_ / Ro, 2.0)))};

        // Coriolis effects (Groves 5.41 & 5.44)
        Eigen::Vector3d w_en_n = {ve_ / he, -vn_ / hn, -ve_ * tL / he};
        Eigen::Vector3d w_ie_n = {w_e * cL, 0.0, w_e * sL};        
    }

    NavState tightlyCoupledNavigator::returnNavState()
    {

    }

    bool tightlyCoupledNavigator::tightlyCoupledCorrection(const Eigen::MatrixXd & psr, const Eigen::MatrixXd & psr_rate, const vec_3_1 & la, const Eigen::MatrixXd & SvPVT)
    {

    }

    bool tightlyCoupledNavigator::setInitialPosition(const vec_3_1 & r_nb_n)
    {
        try
        {
            phi_  = r_nb_n(0);
            lamb_ = r_nb_n(1);
            h_    = r_nb_n(2);

            return true;
        }
        catch(std::exception &e)
        {
            std::cerr << e.what() << std::endl;
            return false;
        }
    }

    bool tightlyCoupledNavigator::setInitialVelocity(const vec_3_1 & v_nb_n)
    {
        try
        {
            vn_ = v_nb_n(0);
            ve_ = v_nb_n(1);
            vd_ = v_nb_n(2);
            
            return true;
        }
        catch(std::exception &e)
        {
            std::cerr << e.what() << std::endl;
            return false;
        }
    }

    bool tightlyCoupledNavigator::setInitialAttitude(const vec_3_1 & a_nb)
    {
        try
        {
            r_ = a_nb(0);
            p_ = a_nb(1);
            y_ = a_nb(2);
            
            q_nb = eul2q(a_nb);
            C_nb_ = eul2Rotm(a_nb);
            return true;
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << std::endl;
            return false;
        }
        

    }
} // end of namespace