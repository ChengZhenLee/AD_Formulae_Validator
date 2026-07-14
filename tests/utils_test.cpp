#include <gtest/gtest.h>
#include "utils.hpp"
#include "configManager.h"


class UtilsTest : public ::testing::Test {
protected:
    void SetUp() override {
        auto& cm = ConfigManager::getInstance();
        cm.set("T", "double");
        cm.set("V", "3");
        cm.set("U", "2");
        cm.set("x", "2");
        cm.set("y", "2");
        cm.set("f", "f");
    }
};

TEST_F(UtilsTest, GenerateParametersBaseCase) {
    auto params = generateParameters<double>("");

    // Assert
    // TODO: params.size() == 2; find "X" and "Y" by name; check role,
    // tensor.shape, and indexNames match SetUp's x=2/y=2.
    ASSERT_TRUE(params.size() == 2);

    Param<double> x;
    std::deque<size_t> expectedShapeX{2};
    std::deque<std::string> indexNamesX{"i"};
    Param<double> y;
    std::deque<size_t> expectedShapeY{2};
    std::deque<std::string> indexNamesY{"j"};

    for (const auto& p : params) {
        if (p.name == "X") x = p;
        if (p.name == "Y") y = p;
    }

    ASSERT_TRUE(x.role == ParamRole::Input);
    ASSERT_TRUE(x.tensor.shape == expectedShapeX);
    ASSERT_TRUE(x.indexNames == indexNamesX);

    ASSERT_TRUE(y.role == ParamRole::Output);
    ASSERT_TRUE(y.tensor.shape == expectedShapeY);
    ASSERT_TRUE(y.indexNames == indexNamesY);
}

TEST_F(UtilsTest, GenerateParametersSingleTangentOrder) {
    auto params = generateParameters<double>("t");

    ASSERT_TRUE(params.size() == 4);

    Param<double> x_1;
    std::deque<size_t> expectedShapeX1{2, 3};
    std::deque<std::string> indexNamesX1{"i", "v_1"};

    Param<double> y_1;
    std::deque<size_t> expectedShapeY1{2, 3};
    std::deque<std::string> indexNamesY1{"j", "v_1"};

    for (const auto& p : params) {
        if (p.name == "X_1") x_1 = p;
        if (p.name == "Y_1") y_1 = p;
    }

    ASSERT_TRUE(x_1.role == ParamRole::Input);
    ASSERT_TRUE(x_1.tensor.shape == expectedShapeX1);
    ASSERT_TRUE(x_1.indexNames == indexNamesX1);

    ASSERT_TRUE(y_1.role == ParamRole::Output);
    ASSERT_TRUE(y_1.tensor.shape == expectedShapeY1);
    ASSERT_TRUE(y_1.indexNames == indexNamesY1);
}

TEST_F(UtilsTest, GenerateParametersSingleAdjointOrder) {
    auto params = generateParameters<double>("a");

    ASSERT_TRUE(params.size() == 4);

    Param<double> x_1;
    std::deque<size_t> expectedShapeX1{2, 2};
    std::deque<std::string> indexNamesX1{"u_1", "i"};

    Param<double> y_1;
    std::deque<size_t> expectedShapeY1{2, 2};
    std::deque<std::string> indexNamesY1{"u_1", "j"};

    for (const auto& p : params) {
        if (p.name == "X_1") x_1 = p;
        if (p.name == "Y_1") y_1 = p;
    }

    ASSERT_TRUE(x_1.role == ParamRole::Output);
    ASSERT_TRUE(x_1.tensor.shape == expectedShapeX1);
    ASSERT_TRUE(x_1.indexNames == indexNamesX1);

    ASSERT_TRUE(y_1.role == ParamRole::Input);
    ASSERT_TRUE(y_1.tensor.shape == expectedShapeY1);
    ASSERT_TRUE(y_1.indexNames == indexNamesY1);
}

TEST_F(UtilsTest, GenerateParametersTwoOrdersCompound) {
    auto params = generateParameters<double>("at");

    ASSERT_TRUE(params.size() == 8);

    Param<double> x_1_2;
    std::deque<size_t> expectedShapeX12{2, 2, 3};
    std::deque<std::string> indexNamesX12{"u_1", "i", "v_2"};

    Param<double> y_1_2;
    std::deque<size_t> expectedShapeY12{2, 2, 3};
    std::deque<std::string> indexNamesY12{"u_1", "j", "v_2"};

    for (const auto& p : params) {
        if (p.name == "X_1_2") x_1_2 = p;
        if (p.name == "Y_1_2") y_1_2 = p;
    }

    ASSERT_TRUE(x_1_2.role == ParamRole::Output);
    ASSERT_TRUE(x_1_2.tensor.shape == expectedShapeX12);
    ASSERT_TRUE(x_1_2.indexNames == indexNamesX12);

    ASSERT_TRUE(y_1_2.role == ParamRole::Input);
    ASSERT_TRUE(y_1_2.tensor.shape == expectedShapeY12);
    ASSERT_TRUE(y_1_2.indexNames == indexNamesY12);
}

TEST(FindParamByNameTest, FindsExistingParam) {
    std::deque<Param<double>> params;
    Param<double> param = Param<double>({}, {}, "name", ParamRole::Input, {});
    params.push_back(param);

    ASSERT_EQ(findParamByName("name", params).name, "name");
}

TEST(FindParamByNameTest, ThrowsOnMissingParam) {
    std::deque<Param<double>> empty;

    EXPECT_THROW(findParamByName<double>("nonexistent", empty), std::runtime_error);
}

TEST(IsADNestedTypeTest, TrueForADTypes) {
    EXPECT_TRUE((isADNestedType<ad::tangent_t<double,3>>));
    EXPECT_TRUE((isADNestedType<ad::adjoint_t<double,2>>));
}

TEST(IsADNestedTypeTest, FalseForPlainNumericTypes) {
    EXPECT_FALSE(isADNestedType<double>);
}

TEST(SeedElementExtractElementTest, RoundTripsThroughOneLayer) {
    Param<double> p({}, {1}, "X", ParamRole::Input, {"i"});
    p.tensor.data = {42.0};
    ad::tangent_t<double,3> elem;

    seedElement(elem, p, 0);
    EXPECT_DOUBLE_EQ(elem.v, 42.0);

    elem.v = 99.0;
    extractElement(elem, p, 0);
    EXPECT_DOUBLE_EQ(p.tensor.data[0], 99.0);
}

TEST(SeedADExtractADTest, RoundTripsASingleTangentOrder) {
    Param<double> p({{1, 3}}, {2, 3}, "X_1", ParamRole::Input, {"i", "v_1"});
    p.tensor.data = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0};

    ad::tangent_t<double,3> x;
    std::deque<size_t> leftCoords, rightCoords;
    size_t primalIndex = 0;

    seedAD(x, p, 1, "t", leftCoords, rightCoords, primalIndex);

    EXPECT_DOUBLE_EQ(x.tangent(0), 10.0);
    EXPECT_DOUBLE_EQ(x.tangent(1), 20.0);
    EXPECT_DOUBLE_EQ(x.tangent(2), 30.0);

    x.tangent(0) = 100.0;
    x.tangent(1) = 200.0;
    x.tangent(2) = 300.0;

    extractAD(x, p, 1, "t", leftCoords, rightCoords, primalIndex);

    EXPECT_DOUBLE_EQ(p.tensor.data[0], 100.0);
    EXPECT_DOUBLE_EQ(p.tensor.data[1], 200.0);
    EXPECT_DOUBLE_EQ(p.tensor.data[2], 300.0);
}
