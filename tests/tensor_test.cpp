#include <gtest/gtest.h>
#include <stdexcept>
#include "structures.h"


class TensorTest : public ::testing::Test {
protected:
    void SetUp() override {
        left = Tensor<double>({2, 2});
        right = Tensor<double>({2, 2});
    }

    Tensor<double> left;
    Tensor<double> right;
};

TEST_F(TensorTest, Addition) {
    left.data = {1, 2, 3, 4};
    right.data = {1, 2, 3, 4};
    auto result = left + right;
    std::deque<size_t> expectedShape{2, 2};

    EXPECT_DOUBLE_EQ(result.data[0], 2);
    EXPECT_DOUBLE_EQ(result.data[1], 4);
    EXPECT_DOUBLE_EQ(result.data[2], 6);
    EXPECT_DOUBLE_EQ(result.data[3], 8);
    EXPECT_EQ(result.shape, expectedShape);
}

TEST_F(TensorTest, AdditionThrowsOnShapeMismatch) {
    Tensor<double> mismatched({3, 3});

    EXPECT_THROW(left + mismatched, std::runtime_error);
}

TEST_F(TensorTest, Contraction) {
    left.data = {1, 2, 3, 4};
    right.data = {1, 2, 3, 4};
    std::deque<size_t> expectedShape{2, 2};

    auto result = Tensor<double>::productGeneralContraction(
        left, right, {1}, {0});

    EXPECT_DOUBLE_EQ(result.data[0], 7);
    EXPECT_DOUBLE_EQ(result.data[1], 10);
    EXPECT_DOUBLE_EQ(result.data[2], 15);
    EXPECT_DOUBLE_EQ(result.data[3], 22);
    EXPECT_EQ(result.shape, expectedShape);
}

TEST_F(TensorTest, ContractionAxisOutOfRange) {
    EXPECT_THROW(Tensor<double>::productGeneralContraction(
        left, right, {2}, {1}),
        std::out_of_range
    );
}

TEST_F(TensorTest, ContractionMismatchedDimensions) {
    Tensor<double> mismatched({3, 3});

    EXPECT_THROW(Tensor<double>::productGeneralContraction(
        left, mismatched, {1}, {0}),
        std::invalid_argument
    );
}

TEST_F(TensorTest, PermuteAxis) {
    Tensor<double> tensor({2, 3});
    tensor.data = {1, 2, 3, 4, 5, 6};

    // Transpose a 2x3 matrix
    auto result = Tensor<double>::permuteAxes(tensor, {1, 0});

    std::deque<size_t> expectedShape{3, 2};

    EXPECT_DOUBLE_EQ(result.data[0], 1);
    EXPECT_DOUBLE_EQ(result.data[1], 4);
    EXPECT_DOUBLE_EQ(result.data[2], 2);
    EXPECT_DOUBLE_EQ(result.data[3], 5);
    EXPECT_DOUBLE_EQ(result.data[4], 3);
    EXPECT_DOUBLE_EQ(result.data[5], 6);
    EXPECT_EQ(result.shape, expectedShape);
}

TEST_F(TensorTest, PermuteAxisThrowsOnRankMismatch) {
    EXPECT_THROW(Tensor<double>::permuteAxes(left, {0, 1, 2}),
    std::invalid_argument);
}