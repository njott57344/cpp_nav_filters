#ifndef GPS_LEAST_SQUARES_H
#define GPS_LEAST_SQUARES_H

#include "cpp_nav_filt_lib/cpp_nav_filt_lib.h"

namespace cpp_nav_filt
{
    typedef struct
    {
        bool weighted_least_squares;
    }GpsLeastSquaresSettings;

    class GpsLeastSquares
    {
        public:
            
            GpsLeastSquares(GpsLeastSquaresSettings& settings_in);
            ~GpsLeastSquares();

            void sendStateEstimate(Eigen::MatrixXd& Y,Eigen::MatrixXd& SvPVT,vec_8_1& X);
            void sendDOPEstimate(vec_8_1& x_hat,Eigen::MatrixXd& SvPVT,mat_4_4& DOP);

        private:

            GpsLeastSquaresSettings ls_settings_;

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
            Eigen::MatrixXd R_; // weighting matrix

            vec_8_1 delta_x_; // update to state estimate (error state)

            vec_8_1 ones_8_1;

            double num_svs_;
            double num_measurements_;
            
            int ctr_ = 0; // counter for solver while loop

    }; // end of class

} // end of namespace

#endif