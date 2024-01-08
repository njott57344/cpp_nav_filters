#include "gps_least_squares/gps_least_squares.h"

namespace cpp_nav_filt
{

    GpsLeastSquares::GpsLeastSquares()
    {
    
    }

    GpsLeastSquares::~GpsLeastSquares()
    {

    }

    void GpsLeastSquares::sendStateEstimate(Eigen::MatrixXd& Y,Eigen::MatrixXd& SvPVT,Common& common,vec_8_1& X)
    {
        int num_measurements = Y.rows();
        int num_svs = SvPVT.rows();

        // Resizing
        Y_.resize(num_measurements,1);
        Yhat_.resize(num_measurements,1);
        deltaY_.resize(num_measurements,1);
        SvPVT_.resize(num_svs,7);        
        G_.resize(num_svs,4);
        H_.resize(num_measurements,8);

        // Setting Internal Variables to the right things
        Y_ = Y;
        SvPVT_ = SvPVT;
        
        // applying clock corrections to pseudoranges
        Y_.block(0,0,num_svs,1) = Y_.block(0,0,num_svs,1) + common.c*SvPVT_.col(6);

        Y_.block(num_svs,0,num_svs,0) = Y_.block(num_svs,0,num_svs,0)*(common.c/common.f_l1);

        x_ = X;
        
        delta_x_.setOnes();
        H_.setZero();

        int rows_G_,cols_G_;
        
        while(ctr_<10 && delta_x_.norm()>0.0001)
        {
            common.sendUnitVectors(x_,SvPVT_,G_);
            common.sendMeasEst(x_,SvPVT_,Yhat_);

            deltaY_ = Y_ - Yhat_;

            H_.block(0,0,num_svs,4) = G_;
            H_.block(num_svs,4,num_svs,4) = G_;
            
            // std::cout<<Yhat_<<std::endl;

            delta_x_ = ((H_.transpose()*H_).inverse())*H_.transpose()*deltaY_;

            x_ = x_ + delta_x_;

            ctr_++;
        }
    
        ctr_ = 0;
        delta_x_.setOnes();
        X = x_; // setting output state to estimated state
        std::cout<<X<<std::endl;
    }

    void GpsLeastSquares::calcStateEstimate()
    {
        
    }
}
