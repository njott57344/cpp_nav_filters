#ifndef COMMON_H
#define COMMON_H

// Eigen
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"

// CPP Headers
#include <iostream>
#include <cmath>

typedef Eigen::Matrix<double,3,1> vec_3_1;
typedef Eigen::Matrix<double,32,3> mat_32_3;
typedef Eigen::Matrix<double,32,27> mat_32_27;
typedef Eigen::Matrix<double,1,27> vec_1_27;

namespace cpp_nav_filt
{
    class Common
    {
        public:

            Common();
            ~Common();

            void passSvEphem(vec_1_27& ephem_in,const int& sv_in);

        private:

            mat_32_27 sv_ephem;

        protected:


    }; // end of class

} // end of namespace

#endif