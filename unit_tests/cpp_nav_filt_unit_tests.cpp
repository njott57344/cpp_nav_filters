#include "gtest/gtest.h"
#include "cpp_nav_filt_lib.h"
#include "lc_ins.h"

/*
    I am confident that functions like calcPsr,calcUnitVectors,calcPsrRate,calcMeasEst already work from
    gps least squares demo. I am going to focus on testing stuff relating to the Loosely Coupled INS for now
*/

TEST(eul2Rotm,RowsColsMagOne)
{
    double row_mag,col_mag;

    vec_3_1 euler_angles;
    mat_3_3 C;

    euler_angles << 10,10,10; // roll,pitch,yaw
    euler_angles = cpp_nav_filt::D2R*euler_angles;

    C = cpp_nav_filt::eul2Rotm(euler_angles);

    for(int j = 0;j<3;j++)
    {
        col_mag = sqrt(pow(C(0,j),2) + pow(C(1,j),2) + pow(C(2,j),2));

        row_mag = sqrt(pow(C(j,0),2) + pow(C(j,1),2) + pow(C(j,2),2));

        std::cout<<col_mag<<" "<<row_mag<<std::endl;

        EXPECT_LT(std::abs(1-col_mag),0.00001) << "Column not Mag 1";
        EXPECT_LT(std::abs(1-row_mag),0.00001) << "Row not Mag 1";
    }


}

TEST(eul2Rotm,C_CtranposeEqualIdenity)
{
    vec_3_1 euler_angles;
    mat_3_3 C;
    
    euler_angles<<10,10,10; // roll,pitch,yaw
    euler_angles = euler_angles*cpp_nav_filt::D2R;
    
    C = cpp_nav_filt::eul2Rotm(euler_angles);

    
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc,argv);
    return RUN_ALL_TESTS();
}