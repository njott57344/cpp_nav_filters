#ifndef GPS_LEAST_SQUARES_H
#define GPS_LEAST_SQUARES_H

#include "common/common.h"

namespace cpp_nav_filt
{
    class GpsLeastSquares
    {
        public:
            
            GpsLeastSquares();
            ~GpsLeastSquares();

            void sendStateEstimate(Eigen::MatrixXd& Y,Eigen::MatrixXd& SvPVT,Common& common,vec_8_1& X);

        private:

            // internal variables
            vec_8_1 x_; // state estimate
            Eigen::MatrixXd H_; // Measurement Model
            Eigen::MatrixXd SvPVT_;
            Eigen::MatrixXd Y_;
            Eigen::MatrixXd Yhat_;
            Eigen::MatrixXd deltaY_;
            Eigen::MatrixXd G_; // Unit Vectors
            vec_8_1 delta_x_; // update to state estimate (error state)

            vec_8_1 ones_8_1;

            int ctr_ = 0; // counter for solver while loop

            // internal functions
            void calcStateEstimate();

    }; // end of class

} // end of namespace

#endif