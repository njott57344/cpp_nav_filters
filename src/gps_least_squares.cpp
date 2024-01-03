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

        x_.setZero();
        delta_x_.setOnes();
        H_.setZero();

        while(ctr_<100 && delta_x_.norm()>0.0001)
        {
            common.sendUnitVectors(x_,SvPVT_,G_);
            common.sendMeasEst(x_,SvPVT_,Yhat_);

            deltaY_ = Y_ - Yhat_;

        }

        X = x_; // setting output state to estimated state
    }

    void GpsLeastSquares::calcStateEstimate()
    {
        
    }
}
