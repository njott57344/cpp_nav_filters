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
            void sendDOPEstimate(vec_8_1& x_hat,Eigen::MatrixXd& SvPVT,Common& common,mat_4_4& DOP);

        private:

            // internal variables
            vec_3_1 pos_;
            vec_3_1 vel_;
            double clk_;
            double clk_drift_;
            vec_8_1 x_;

            Eigen::MatrixXd H_; // Measurement Model
            Eigen::MatrixXd SvPVT_;
            Eigen::MatrixXd Y_;
            Eigen::MatrixXd Yhat_;
            Eigen::MatrixXd deltaY_;
            Eigen::MatrixXd G_; // Unit Vectors
            vec_8_1 delta_x_; // update to state estimate (error state)

            vec_8_1 ones_8_1;

            double num_svs_;
            double num_measurements_;
            
            int ctr_ = 0; // counter for solver while loop

    }; // end of class

} // end of namespace

#endif