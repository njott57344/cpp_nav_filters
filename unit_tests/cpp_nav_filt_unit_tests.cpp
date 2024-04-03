#include "gtest/gtest.h"
#include "cpp_nav_filt_lib.h"

namespace {
    int GetMeaningOfLife() {return 42;}
}

TEST(TestTopic,TrivialEquality)
{
    EXPECT_EQ(GetMeaningOfLife(),42);
}

TEST(TestTopic,MoreEqualityTests)
{
    ASSERT_EQ(GetMeaningOfLife(),42) << "Oh no, a mistake!";
    EXPECT_FLOAT_EQ(23.23F,23.23F);
}