#include <gtest/gtest.h>
#include <stdexcept>
#include "structures.h"


TEST(ParamTest, ConstructorParsesNoOrderSuffix) {
    Param<double> p({}, {}, "X", ParamRole::Input, {});

    EXPECT_TRUE(p.activeOrders.empty());
    EXPECT_EQ(p.highestOrder, 0);
}

TEST(ParamTest, ConstructorParsesSingleOrderSuffix) {
    Param<double> p({}, {}, "X_1", ParamRole::Input, {});

    std::vector<int> expectedOrders{1};
    ASSERT_EQ(p.activeOrders, expectedOrders);
    ASSERT_DOUBLE_EQ(p.highestOrder, 1);
}

TEST(ParamTest, ConstructorParsesMultipleOrderSuffixes) {
    Param<double> p({}, {}, "X_1_2", ParamRole::Input, {});

    std::vector<int> expectedOrders{1, 2};
    ASSERT_EQ(p.activeOrders, expectedOrders);
    ASSERT_DOUBLE_EQ(p.highestOrder, 2);
}

TEST(ParamTest, ConstructorTracksHighestOrderRegardlessOfSuffixOrder) {
    Param<double> p({}, {}, "X_3_1_2", ParamRole::Input, {});

    std::vector<int> expectedOrders{3, 1, 2};
    ASSERT_EQ(p.activeOrders, expectedOrders);
    ASSERT_DOUBLE_EQ(p.highestOrder, 3);
}

TEST(ParamTest, IsSeedReflectsRole) {
    Param<double> input({}, {}, "X", ParamRole::Input, {});
    Param<double> output({}, {}, "Y", ParamRole::Output, {});
    EXPECT_TRUE(input.isSeed());
    EXPECT_FALSE(output.isSeed());
}

TEST(ParamTest, IsActiveChecksOrderMembership) {
    Param<double> p({}, {}, "X_1_2", ParamRole::Input, {});

    EXPECT_TRUE(p.isActive(1)); 
    EXPECT_TRUE(p.isActive(2));
    EXPECT_FALSE(p.isActive(3));
}

TEST(ParamTest, GetShapeAtReturnsConfiguredExtent) {
    Param<double> p({{1, 2}, {2, 3}}, {}, "X_1_2", ParamRole::Input, {});

    EXPECT_EQ(p.getShapeAt(1), 2); 
    EXPECT_EQ(p.getShapeAt(2), 3);
}

TEST(ParamTest, GetShapeAtBehaviorForMissingOrder) {
    Param<double> p({{1, 2}}, {}, "X_1_2", ParamRole::Input, {});

    EXPECT_THROW(p.getShapeAt(99), std::out_of_range);
}
