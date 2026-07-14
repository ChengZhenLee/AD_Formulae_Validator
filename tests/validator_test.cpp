#include <gtest/gtest.h>
#include <filesystem>
#include <limits>
#include "generator/symbolic/validator.h"
#include "generator/readWrite.h"

TEST(CompareTensorsTest, ValidWithinTolerance) {
    Tensor<double> left({2, 2});
    left.data = {1, 2, 3, 4};
    Tensor<double> right({2, 2});
    right.data = {1, 2, 3, 4};

    auto status = compareTensors(
        left, right, 
        std::sqrt(std::numeric_limits<double>::epsilon())
    );

    EXPECT_EQ(status, "valid");
}

TEST(CompareTensorsTest, InvalidOutsideTolerance) {
    Tensor<double> left({2, 2});
    left.data = {1, 2, 3, 4};
    Tensor<double> right({2, 2});
    right.data = {1, 2, 3, 5};

    auto status = compareTensors(
        left, right, 
        std::sqrt(std::numeric_limits<double>::epsilon())
    );

    EXPECT_EQ(status, "Incorrect data values");
}

TEST(CompareTensorsTest, InvalidOnShapeSizeMismatch) {
    Tensor<double> left({2, 2});
    left.data = {1, 2, 3, 4};
    Tensor<double> right({1});
    right.data = {1, 2};

    auto status = compareTensors(
        left, right, 
        std::sqrt(std::numeric_limits<double>::epsilon())
    );
    
    EXPECT_EQ(status, "Incorrect tensor shape");
}

TEST(CompareTensorsTest, InvalidOnShapeMismatch) {
    Tensor<double> left({2, 2});
    left.data = {1, 2, 3, 4};
    Tensor<double> right({2, 3});
    right.data = {1, 2, 3, 4, 5, 6};

    auto status = compareTensors(
        left, right, 
        std::sqrt(std::numeric_limits<double>::epsilon())
    );
    
    EXPECT_EQ(status, "Incorrect tensor shape");
}

TEST(CompareTensorsTest, InvalidOnDataSizeMismatch) {
    Tensor<double> left({2, 2});
    left.data = {1, 2, 3, 4};
    Tensor<double> right({2, 2});
    right.data = {1, 2, 3, 4, 5};

    auto status = compareTensors(
        left, right, 
        std::sqrt(std::numeric_limits<double>::epsilon())
    );
    
    EXPECT_EQ(status, "Incorrect tensor data size");
}


class ValidatorTest : public ::testing::Test {
protected:
    static Param<double> makeParam(
        const std::string& name,
        std::deque<double> data,
        std::deque<size_t> shape,
        std::deque<std::string> indexNames)
    {
        Param<double> p;
        p.name = name;
        p.indexNames = indexNames;
        p.tensor = Tensor<double>(shape);
        p.tensor.data = data;
        return p;
    }

    void SetUp() override {
        inputPath = (std::filesystem::temp_directory_path() / "validator_test_input.bin").string();
        outputPath = (std::filesystem::temp_directory_path() / "validator_test_output.txt").string();
    }

    void TearDown() override {
        std::filesystem::remove(inputPath);
        std::filesystem::remove(outputPath);
    }

    std::string inputPath;
    std::string outputPath;
};

TEST_F(ValidatorTest, ValidateValidParameters) {
    std::vector<Param<double>> adResults{
        makeParam("first", {1, 2, 3, 4}, {2, 2}, {"i", "j"})
    };
    writeParameters<double>(adResults, inputPath);
    
    std::vector<Param<double>> formulaResults{
        makeParam("first", {1, 2, 3, 4}, {2, 2}, {"i", "j"})
    };

    bool allValid = validateParameters(formulaResults, inputPath, outputPath);

    EXPECT_TRUE(allValid);
}

TEST_F(ValidatorTest, ValidateMismatchedParameters) {
    std::vector<Param<double>> adResults{
        makeParam("first", {1, 2, 3, 4}, {2, 2}, {"i", "j"})
    };
    writeParameters<double>(adResults, inputPath);
    
    std::vector<Param<double>> formulaResults{
        makeParam("first", {1, 2, 3, 5}, {2, 2}, {"i", "j"})
    };

    bool allValid = validateParameters(formulaResults, inputPath, outputPath);

    EXPECT_FALSE(allValid);
}
