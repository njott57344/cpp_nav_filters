#include "tightly_coupled_ins.h"

namespace cpp_nav_filt
{

    tightlyCoupledNavigator::tightlyCoupledNavigator()
    {
        P_.setZero();
        A_.setZero();
        Q_.setZero();

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

    bool tightlyCoupledNavigator::Dynamics(const vec_3_1 & wb, const vec_3_1& fb, const double & dt)
    {
        try
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
            double Re = transverseRadiusOfCurvature(phi_);
            double Rn = meridianRadiusOfCurvature(phi_);
            double R_es_e = geocentricRadius(phi_);

            // Adjusted heights
            double hn = h_ + Rn;
            double he = h_ + Re;
            double hn2 = hn * hn;
            double he2 = he * he;

            // TODO: Make this part of cpp_nav_filt_lib
            // Somigliana gravity model
            double g0 = 9.7803253359 * ((1.0 + 0.001931853 * sL2) / std::sqrt(1.0 - e2 * sL2));
            vec_3_1 g{-8.08e-9 * h_ * std::sin(2.0 * phi_), 0.0,
                            g0 * (1.0 -
                                    (2.0 / Ro * h_) * (1.0 + flattening * (1.0 - 2.0 * sL2) +
                                                    (std::pow(w_e * Ro, 2.0) * Rp / mu_g)) +
                                    (3.0 * std::pow(h_ / Ro, 2.0)))};


            // Coriolis effects (Groves 5.41 & 5.44)
            vec_3_1 w_en_n = transportRate(ve_,vn_,h_,phi_);
            vec_3_1 w_ie_n = navFrameRotationRate(phi_);;  

            //* ===== MECHANIZATION =====

            // Groves E.38
            vec_3_1 a_ib_b = w_ib_b * dt;
            double mag_a_ib_b = a_ib_b.norm();
            vec_4_1 q_new_old;
            if (mag_a_ib_b > 1.0e-8) {
            q_new_old(0) = std::cos(0.5 * mag_a_ib_b);
            q_new_old.tail(3) << std::sin(0.5 * mag_a_ib_b) * a_ib_b / mag_a_ib_b;
            } else {
            q_new_old << 1.0, 0.0, 0.0, 0.0;
            }

            vec_4_1 b;
            b(0) = 0.0;
            b.tail<3>() = 0.5 * (w_en_n + w_ie_n) * dt;
            vec_4_1 q_new = qMult(q_bn_, q_new_old) - qMult(b, q_bn_);
            q_new = qNormalize(q_new);  // <- possibly unnecessary

            // Groves 5.48
            vec_4_1 q_avg = 0.5 * (q_bn_ + q_new);
            mat_3_3 C_avg = q2DCM(q_avg);
            vec_3_1 f_ib_n = C_avg * f_ib_b;
            // std::cout << "f_ib_n=[" << f_ib_n(0) << "," << f_ib_n(1) << "," << f_ib_n(2) << "]\n";

            // Groves 5.54
            vec_3_1 coriolis = w_en_n + 2.0 * w_ie_n;
            vec_3_1 v_old = {vn_, ve_, vd_};
            vec_3_1 v_new = v_old + (f_ib_n + g - makeSkewSymmetic(coriolis) * v_old) * dt;

            double h_new = h_ - 0.5 * (vd_ + v_new(2)) * dt;
            double phi_new = phi_ + 0.5 * (vn_ / hn + v_new(0) / (Rn + h_new)) * dt;
            double Re_new = transverseRadiusOfCurvature(phi_new);
            double lamb_new = lamb_ + 0.5 * (ve_ / (he * cL) + v_new(1) / ((Re_new + h_new) * std::cos(phi_new))) * dt;

            double new_clock_bias = clk_b_ + clk_d_ * dt;

            //* ===== SAVE PROPAGETED VALUES =====
            q_bn_ = q_new;
            vec_3_1 rpy = q2eul(q_bn_);
            r_ = rpy(0);
            p_ = rpy(1);
            y_ = rpy(2);
            vn_ = v_new(0);
            ve_ = v_new(1);
            vd_ = v_new(2);
            phi_ = phi_new;
            lamb_ = lamb_new;
            h_ = h_new;
            clk_b_ = new_clock_bias;

            // * ========= Propagate State Covariance (INS States) ======
            // Groves 14.64
            Eigen::Vector3d w_ib_n = C_avg * w_ib_b;
            Eigen::Matrix3d F11 = -makeSkewSymmetic(w_ib_n);

            // Groves 14.65
            Eigen::Matrix3d F12, F13, F22, F23, F33;
            
            F12 << 0.0, -1.0/he, 0.0, 
                1.0/he,     0.0, 0.0, 
                0.0,   tL/he, 0.0;

            // Groves I.14
            F13 << w_e*sL/hn, 0.0,    -ve_/he2,
                0.0, 0.0,     vn_/hn2,
                w_e*cL/hn + ve_/(hn*he*cL2), 0.0,  ve_*tL/he2;

            // Groves 14.67
            Eigen::Matrix3d F21 = -makeSkewSymmetic(f_ib_n);

            // Groves 14.68
            F22 << vd_/hn, -2.0*(ve_*tL/he + w_e*sL), vn_/hn,
                ve_*tL/he + 2.0*w_e*sL, (vn_*tL + vd_)/he, ve_/he + 2.0*w_e*cL,
                -2.0*vn_/hn,    -2.0*(ve_/he - w_e*cL), 0.0;

            // Groves I.15
            F23 << -ve_*ve_/(hn*he*cL2) - 2.0*ve_*w_e*cL/hn, 0.0,                      -std::pow(ve_/he,2.0)*tL + vn_*vd_/hn2, 
                    vn_*ve_/(hn*he*cL2) + 2.0*(w_e*vn_*cL - vd_*sL)/hn, 0.0,                                   (vn_*ve_*tL + ve_*vd_)/he2, 
                    2.0*ve_*w_e*sL/hn, 0.0, -std::pow(ve_/he,2.0) - std::pow(vn_/hn,2.0) + 2*g0/R_es_e;

            // Groves I.16
            F33 << 0.0, 0.0, vn_/hn, 
                ve_*tL/hn, 0.0, ve_/he, 
                0.0, 0.0,   0.0;

            // Groves I.17
            A_.block<3,3>(0,0)   = I3 + F33 * dt;
            A_.block<3,3>(0,3)   = I3 * dt;
            A_.block<3,3>(3,0)   = F23 * dt;
            A_.block<3,3>(3,3)   = I3 + F22 * dt;
            A_.block<3,3>(3,6)   = F21 * dt;
            A_.block<3,3>(3,9)   = C_avg * dt;
            A_.block<3,3>(6,0)   = F13 * dt;
            A_.block<3,3>(6,3)   = F12 * dt;
            A_.block<3,3>(6,6)   = I3 + F11 * dt;
            A_.block<3,3>(6,12)  = C_avg * dt;
            A_.block<3,3>(9,9)   = I3;
            A_.block<3,3>(12,12) = I3;

            // Groves 14.81 (but divided in half for Q/2)
            double T2 = dt * dt;
            double T3 = dt * T2;
            double T4 = dt * T3;
            double T5 = dt * T4;
            double T6 = dt * T5;
            double T7 = dt * T6;
            Eigen::Matrix3d F21_F21t = F21 * F21.transpose();
            Eigen::Matrix3d F21_C = F21 * C_avg;
            Eigen::Matrix3d Q11 = 0.5 * ((Srg_ * dt + (Sbgd_ * T3) / 3.0) * I3);
            Eigen::Matrix3d Q15 = 0.5 * ((Sbgd_ * T2) / 2.0 * C_avg);
            Eigen::Matrix3d Q21 = 0.5 * (((Srg_ * T2) / 2.0 + (Sbgd_ * T4) / 4.0) * F21);
            Eigen::Matrix3d Q22 = 0.5 * ((Sra_ * dt + (Sbad_ * T3) / 3.0) * I3 + ((Srg_ * T3) / 3.0 + (Sbgd_ * T5) / 5.0) * F21_F21t);
            Eigen::Matrix3d Q24 = 0.5 * ((Sbad_ * T2) / 2.0 * C_avg);
            Eigen::Matrix3d Q31 = 0.5 * (((Srg_ * T3) / 3.0 + (Sbgd_ * T5) / 5.0) * F21);
            Eigen::Matrix3d Q32 = 0.5 * (((Sra_ * T2) / 2.0 + (Sbad_ * T4) / 4.0) * I3 + ((Srg_ * T4) / 4.0 + (Sbgd_ * T6) / 6.0) * F21_F21t);
            Eigen::Matrix3d Q33 = 0.5 * (((Sra_ * T3) / 3.0 + (Sbad_ * T5) / 5.0) * I3 + ((Srg_ * T5) / 5.0 + (Sbgd_ * T7) / 7.0) * F21_F21t);
            Eigen::Matrix3d Q34 = 0.5 * ((Sbad_ * T3) / 3.0 * C_avg);
            Eigen::Matrix3d Q35 = 0.5 * ((Sbgd_ * T4) / 4.0 * F21_C);
            double q25 = 0.5 * ((Sbgd_ * T3) / 3.0);

            // Groves 14.80
            Q_.block<3,3>(0,0)   = Q33;
            Q_.block<3,3>(0,3)   = Q32.transpose();
            Q_.block<3,3>(0,6)   = Q31.transpose();
            Q_.block<3,3>(0,9)   = Q34;
            Q_.block<3,3>(0,12)  = Q35;
            Q_.block<3,3>(3,0)   = Q32.transpose();
            Q_.block<3,3>(3,3)   = Q22;
            Q_.block<3,3>(3,6)   = Q21;
            Q_.block<3,3>(3,9)   = Q24;
            Q_.block<3,3>(3,12)  = q25 * F21 * C_avg;
            Q_.block<3,3>(6,0)   = Q31.transpose();
            Q_.block<3,3>(6,3)   = Q21.transpose();
            Q_.block<3,3>(6,6)   = Q11;
            Q_.block<3,3>(6,12)  = Q15;
            Q_.block<3,3>(9,0)   = Q34.transpose();
            Q_.block<3,3>(9,3)   = Q24;
            Q_.block<3,3>(9,9)   = 0.5 * Sbad_ * dt * I3;
            Q_.block<3,3>(12,0)  = Q32.transpose();
            Q_.block<3,3>(12,3)  = q25 * F21.transpose() * C_avg;
            Q_.block<3,3>(12,6)  = Q15;
            Q_.block<3,3>(12,12) = 0.5 * Sbgd_ * dt * I3;

            //* ====== Propagate State Covariance (Clock States) ====== //
            mat_2_2 A_clock = mat_2_2::Identity();
            mat_2_2 Q_clock = mat_2_2::Identity();

            A_clock(0,1) = dt;

            double q11,q12,q22;
            double h11,h12,h22;

            q11 = Scd_ * (std::pow(dt,3) / 3);
            q12 = Scd_ * std::pow(dt,2) * 0.5;
            q22 = Scd_ * dt;

            h11 = 0.5 * q11 - 0.5 * dt * q12;
            h12 = 0.5 * q12 - 0.25 * dt * q22;
            h22 = 0.5 * q22;

            Q_clock(0,0) = h11;
            Q_clock(0,1) = h12;
            Q_clock(1,0) = Q_clock(0,1);
            Q_clock(2,2) = h22;

            Q_.block<2,2>(15,15) = Q_clock;
            A_.block<2,2>(15,15) = A_clock;

            P_ = A_ * (P_ + Q_) * A_.transpose() + Q_;

            return true;      
        }
        catch(std::exception &e)
        {
            std::cerr << e.what() << std::endl;
            return false;
        }
    }

    bool tightlyCoupledNavigator::returnNavState(NavState & out)
    {
        try
        {
            out.lat = phi_;
            out.lat_stdev = std::sqrt(P_(0,0));
            
            out.lon = lamb_;
            out.lon_stdev = std::sqrt(P_(1,1));

            out.alt = h_;
            out.alt_stdev = std::sqrt(P_(2,2));

            out.vn = vn_;
            out.vn_stdev = std::sqrt(P_(3,3));

            out.ve = ve_;
            out.ve_stdev = std::sqrt(P_(4,4));

            out.vd = vd_;
            out.vd_stdev = std::sqrt(P_(5,5));

            out.r = r_;
            out.r_stdev = std::sqrt(P_(6,6));

            out.p = p_;
            out.p_stdev = std::sqrt(P_(7,7));

            out.y = y_;
            out.y_stdev = std::sqrt(P_(8,8));

            out.clk_b = clk_b_;
            out.clk_b_stdev = std::sqrt(P_(15,15));

            out.clk_d = clk_d_;
            out.clk_d_stdev = std::sqrt(P_(16,16));

            return true;
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << std::endl;

            return false;
        }
        
    }

    bool tightlyCoupledNavigator::tightlyCoupledCorrection(const Eigen::MatrixXd & psr, const Eigen::MatrixXd & psr_rate, const vec_3_1 & la,
                                                           const vec_3_1 & w_ib_b, const Eigen::MatrixXd & SvPVT, const Eigen::MatrixXd & R)
    {
        /* Order of Operations
        1) Remove biases from angular rates [done]
        2) Put lla estimate in ECEF frame [done]
        3) Put ned velocity in ECEF frame [done]
        4) Estimate Pseudoranges and Pseudorange Rates [done]
        5) Form innovation and jacobian
        6) Kalman Update
        */
        try
        {
            int num_sv = SvPVT.rows();

            vec_3_1 wb = w_ib_b - bg_;

            vec_3_1 r_nb_n{phi_ * cpp_nav_filt::R2D,lamb_ * cpp_nav_filt::R2D ,h_};
            vec_3_1 v_nb_n{vn_,ve_,vd_};
            vec_3_1 a_nb{r_,p_,y_};

            vec_3_1 r_nb_e = lla2ecefPos(r_nb_e); // position of antenna in ecef frame
            vec_3_1 v_nb_e = ned2ecefVel(v_nb_n,r_nb_n); // velocity of antenna in ecef frame
            mat_3_3 C_be = eul2EcefDCM(a_nb,r_nb_n); // DCM body to ecef frame

            vec_3_1 r_no_e = r_nb_e + C_be * la;
            vec_3_1 v_no_e = v_nb_e + C_be * (makeSkewSymmetic(wb) * la);

            mat_3_3 C_ne = ned2ecefDCM(r_nb_n);
            mat_3_3 C_en = C_ne.transpose();

            Eigen::MatrixXd Y_hat = calcMeasEst(SvPVT,r_no_e,v_no_e,clk_b_,clk_d_);
            Eigen::MatrixXd Y,dz;

            Y.resize(num_sv,1);
            dz.resize(num_sv,1);

            Y_hat.block(0,0,num_sv,1) = psr;
            Y_hat.block(num_sv,0,num_sv,1) = psr_rate;

            dz = Y - Y_hat;

            Eigen::MatrixXd U = calcUnitVectors(SvPVT,r_no_e,clk_b_);
            
            Eigen::MatrixXd H;
            H.resize(num_sv * 2,17);

            H.block(0,15,num_sv,1) = U.block(0,3,num_sv,1);
            H.block(num_sv,16,num_sv,1) = U.block(0,3,num_sv,1);

            for (int i = 0;i<num_sv;i++)
            {
                vec_1_3 ned_U = (C_en * (U.block<1,3>(i,0).transpose())).transpose();
                H.block<1,3>(i,0) = ned_U;
                H.block<1,3>(i+num_sv,3) = ned_U;
            }


            return kalmanUpdate(H,R,dz);
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << std::endl;
            return false;
        }
        
    }

    bool tightlyCoupledNavigator::kalmanUpdate(const Eigen::MatrixXd & H, const Eigen::MatrixXd & R, const Eigen::MatrixXd & dz)
    {
        // TODO: Normalized Innovation Filtering (make sure normal works first)
        try
        {
            Eigen::MatrixXd S = H + P_ + H.transpose() + R;

            Eigen::MatrixXd K = P_ * H.transpose() * S.inverse();

            Eigen::MatrixXd L = I17 - K * H;

            P_ = L * P_ * L.transpose() + K * R * K.transpose();

            dx_ = K * dz;

            return closedLoopCorrection();
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << std::endl;
            return false;
        }
        
    }

    bool tightlyCoupledNavigator::closedLoopCorrection()
    {
        try
        {
            double cL = std::cos(phi_);
            double sL = std::sin(phi_);

            double Re = transverseRadiusOfCurvature(phi_);
            double Rn = meridianRadiusOfCurvature(phi_);

            vec_3_1 Tpr{Rn + h_,cL * (Re + h_),-1.0};

            vec_4_1 att_err{1.0, -0.5 * dx_(6), -0.5 * dx_(7), -0.5 * dx_(8)};

            // correct atttitude
            q_bn_ = qMult(att_err,q_bn_);
            q_bn_ = qNormalize(q_bn_);
            mat_3_3 C_bn = q2DCM(q_bn_);
            r_ = std::atan2(C_bn(2, 1), C_bn(2, 2));
            p_ = -std::asin(C_bn(2, 0));
            y_ = std::atan2(C_bn(1, 0), C_bn(0, 0));

            // correct velocity 
            vn_ -= dx_(3);
            ve_ -= dx_(4);
            vd_ -= dx_(5);
            
            // correct position
            phi_ -= dx_(0) / Tpr(0);
            lamb_ -= dx_(1) / Tpr(1);
            h_ -= dx_(2) / Tpr(2);

            // update bias estimates
            ba_(0) += dx_(9);
            ba_(1) += dx_(10);
            ba_(2) += dx_(11);

            bg_(0) += dx_(12);
            bg_(1) += dx_(13);
            bg_(2) += dx_(14);

            clk_b_ += dx_(15);
            clk_d_ += dx_(16);

            dx_.setZero();

            return true;
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << std::endl;
            return false;
        }
        
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
            
            q_bn_ = eul2q(a_nb);
            return true;
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << std::endl;
            return false;
        }
    }

    bool tightlyCoupledNavigator::setImuClockErrors(const double & ba, const double & na, const double & ka, 
                                                    const double & bg, const double & ng, const double & kg, const double & nd)
    {
        try
        {
            Sra_ = 2 * na * na;
            Sbad_ = 2 * ba * ba;
            Srg_ = 2 * ng * ng;
            Sbgd_ = 2 * bg * bg;
            Scd_ = 2 * nd * nd;
            return true;
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << std::endl;
            return false;
        }
    }
} // end of namespace