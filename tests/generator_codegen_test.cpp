#include <gtest/gtest.h>
#include "codegen/generator.h"
#include "configManager.h"


class GeneratorCodegenTest : public ::testing::Test {
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

TEST_F(GeneratorCodegenTest, GenerateNestedADTypeSingleTangent) {
    auto result = generateNestedADType("t");

    ASSERT_EQ(result, "T_t<double,3>");
}

TEST_F(GeneratorCodegenTest, GenerateNestedADTypeSingleAdjoint) {
    auto result = generateNestedADType("a");

    // Assert
    ASSERT_EQ(result, "A_t<double,2>");
}

TEST_F(GeneratorCodegenTest, GenerateNestedADTypeTangentOverAdjoint) {
    auto result = generateNestedADType("at");

    ASSERT_EQ(result, "T_t<A_t<double,2>,3>");
}

TEST_F(GeneratorCodegenTest, GenerateXNestedADTypeWrapsInXT) {
    auto result = generateXNestedADType("t");

    ASSERT_EQ(result, "X_t<T_t<double,3>>");
}

TEST_F(GeneratorCodegenTest, GenerateYNestedADTypeWrapsInYT) {
    auto result = generateYNestedADType("t");

    ASSERT_EQ(result, "Y_t<T_t<double,3>>");
}

TEST_F(GeneratorCodegenTest, GetCurrentLayerADTypeUsesPrefixOnly) {
    std::string result = getCurrentLayerADType(1, "at");
    
    ASSERT_EQ(result, "A_t<double,2>");
}

TEST_F(GeneratorCodegenTest, GetCurrentLayerFunctionNamePrimalAtOrderZero) {
    std::string result = getCurrentLayerFunctionName(0);
    
    ASSERT_EQ(result, "f");
}

TEST_F(GeneratorCodegenTest, GetCurrentLayerFunctionNameADLayers) {
    std::string result = getCurrentLayerFunctionName(2);

    ASSERT_EQ(result, "AD_F_2");
}

TEST_F(GeneratorCodegenTest, GenerateRegisterInputStringSkipsTangentOrders) {
    std::string result = generateRegisterInputString("tt");

    ASSERT_TRUE(result.find(".register_input()") == std::string::npos);
}

TEST_F(GeneratorCodegenTest, GenerateRegisterInputStringUnwrapsPerOrder) {
    // sequence "at": the 'a' is at position 0, with one layer ('t') above
    // it, so reaching it from the outermost x needs exactly one .value()
    // unwrap before .register_input().
    std::string result = generateRegisterInputString("at");

    ASSERT_TRUE(result.find("x[i].value().register_input();") != std::string::npos);
}

TEST_F(GeneratorCodegenTest, GenerateResetTapeStringOneCallPerAdjointOrder) {
    std::string resultTangent = generateResetTapeString("tt");
    std::string resultAdjoint = generateResetTapeString("at");

    ASSERT_TRUE(resultTangent.find("tape::reset();") == std::string::npos);
    ASSERT_TRUE(resultAdjoint.find("A_t<double,2>::tape::reset();") != std::string::npos);
}

TEST_F(GeneratorCodegenTest, GenerateTangentEmitsPrimalCallAtOrderOne) {
    auto xType = generateXNestedADType("t");
    auto yType = generateYNestedADType("t");

    auto body = generateTangent(1, "t", xType, yType);

    ASSERT_TRUE(body.find("f") != std::string::npos);
    ASSERT_TRUE(body.find("extractPrimal(") != std::string::npos);
}

TEST_F(GeneratorCodegenTest, GenerateTangentCallsPreviousLayerAtHigherOrder) {
    auto xType = generateXNestedADType("tt");
    auto yType = generateYNestedADType("tt");

    auto body = generateTangent(2, "tt", xType, yType);

    ASSERT_TRUE(body.find("AD_F_1") != std::string::npos);
    ASSERT_TRUE(body.find("extractPrimal(") == std::string::npos);
}

TEST_F(GeneratorCodegenTest, GenerateAdjointEmitsTapeInitAndInterpret) {
    auto xType = generateXNestedADType("a");
    auto yType = generateYNestedADType("a");

    auto body = generateAdjoint(1, "a", xType, yType);

    ASSERT_TRUE(body.find("::tape::init_adjoints();") != std::string::npos);
    ASSERT_TRUE(body.find("::tape::interpret();") != std::string::npos);
}