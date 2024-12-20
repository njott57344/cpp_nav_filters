#include "gps_least_squares.h"

namespace cpp_nav_filt
{

    GpsLeastSquares::GpsLeastSquares(GpsLeastSquaresSettings& settings_in)
    {
        ls_settings_ = settings_in;
        ones_8_1.setOnes();
    }

    GpsLeastSquares::~GpsLeastSquares()
    {

    }

    void GpsLeastSquares::sendStateEstimate(Eigen::MatrixXd& Y,Eigen::MatrixXd& SvPVT,vec_8_1& X)
    {
        num_measurements_ = Y.rows();
        num_svs_ = SvPVT.rows();

        // Resizing
        Y_.resize(num_measurements_,1);
        Yhat_.resize(num_measurements_,1);
        deltaY_.resize(num_measurements_,1);
        SvPVT_.resize(num_svs_,7);        
        G_.resize(num_svs_,4);
        H_.resize(num_measurements_,8);

        // Setting Internal Variables to the right things
        Y_ = Y;
        SvPVT_ = SvPVT;
        
        // applying clock corrections to pseudoranges
        Y_.block(0,0,num_svs_,1) = Y_.block(0,0,num_svs_,1)  + c*SvPVT_.col(6);

        Y_.block(num_svs_,0,num_svs_,1) = -Y_.block(num_svs_,0,num_svs_,1)*(c/f_l1);
        
        x_ = X;

        pos_ = x_.block(0,0,3,1);
        clk_ = x_[3];
        vel_ = x_.block(4,0,3,1);
        clk_drift_ = x_[7];
                
        delta_x_.setOnes();
        H_.setZero();

        int rows_G_,cols_G_;
        if(num_measurements_>=8)
        {
            while(ctr_<100 && delta_x_.norm()>0.0000001)
            {
                G_ = cpp_nav_filt::calcUnitVectors(SvPVT_,pos_,clk_);

                Yhat_ = cpp_nav_filt::calcMeasEst(SvPVT_,pos_,vel_,clk_,clk_drift_);
                deltaY_ = Y_ - Yhat_;

                H_.block(0,0,num_svs_,4) = G_;
                H_.block(num_svs_,4,num_svs_,4) = G_;

                // std::cout<<deltaY_<<std::endl<<std::endl;
    
                // std::cout<<Yhat_<<std::endl;

                if(!ls_settings_.weighted_least_squares)
                {
                    delta_x_ = ((H_.transpose()*H_).inverse())*H_.transpose()*deltaY_;
                }
                else
                {

                }

                x_ = x_ + delta_x_;
                
                pos_ = x_.block(0,0,3,1);
                clk_ = x_[3];
                vel_ = x_.block(4,0,3,1);
                clk_drift_ = x_[7];
                
                ctr_++;
            }
        
            ctr_ = 0;
            delta_x_.setOnes();
            X = x_; // setting output state to estimated state         
        }
        else
        {
            std::cout<<"WARNING: NOT ENOUGH SVS TO COMPUTE SOLUTION SKIPPING EPOCH!!"<<std::endl;
            X = NAN*ones_8_1;
        }
    }

    void GpsLeastSquares::sendDOPEstimate(vec_8_1& x_hat,Eigen::MatrixXd& SvPVT,mat_4_4& DOP)
    {
        SvPVT_ = SvPVT;
        x_ = x_hat;
        
        pos_ = x_.block(0,0,3,1);
        clk_ = x_[3];
        vel_ = x_.block(4,0,3,1);
        clk_drift_ = x_[7];
        
        num_svs_ = SvPVT_.rows();

        H_.resize(num_svs_,4);

        H_ = cpp_nav_filt::calcUnitVectors(SvPVT_,pos_,clk_);

        DOP = (H_.transpose()*H_).inverse();
    }

}
