#include <gtest/gtest.h>
#include "symbolic/formulaDriver.hpp"


class FormulaDriverTest : public ::testing::Test {
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
};

TEST_F(FormulaDriverTest, ContractByMetadataSharedIndex) {
    auto left = makeParam("left", {1, 2, 3, 4}, {2, 2}, {"i", "j"});
    auto right = makeParam("right", {1, 2, 3, 4}, {2, 2}, {"j", "v"});

    std::deque<std::string> outIndexNames;
    auto result = contractByMetadata(left, right, outIndexNames);

    std::deque<std::string> expectedIndexNames({"i", "v"});
    ASSERT_DOUBLE_EQ(result.data[0], 7);
    ASSERT_DOUBLE_EQ(result.data[1], 10);
    ASSERT_DOUBLE_EQ(result.data[2], 15);
    ASSERT_DOUBLE_EQ(result.data[3], 22);
    ASSERT_EQ(outIndexNames, expectedIndexNames);
}

TEST_F(FormulaDriverTest, ContractByMetadataNoSharedIndex) {
    auto left = makeParam("left", {1, 2, 3, 4}, {2, 2}, {"i", "j"});
    auto right = makeParam("right", {1, 2, 3, 4}, {2, 2}, {"u", "v"});

    std::deque<std::string> outIndexNames;
    ASSERT_THROW(contractByMetadata(left, right, outIndexNames), std::invalid_argument);
}

TEST_F(FormulaDriverTest, AlignToTargetNoOpWhenAlreadyAligned) {
    Tensor<double> tensor({2, 2});
    tensor.data = {1, 2, 3, 4};
    std::deque<std::string> source{"i", "j"};
    std::deque<std::string> target{"i", "j"};

    auto result = alignToTarget(tensor, source, target);
    
    ASSERT_DOUBLE_EQ(tensor.data[0], 1);
    ASSERT_DOUBLE_EQ(tensor.data[1], 2);
    ASSERT_DOUBLE_EQ(tensor.data[2], 3);
    ASSERT_DOUBLE_EQ(tensor.data[3], 4);
}

TEST_F(FormulaDriverTest, AlignToTargetPermutes) {
        Tensor<double> tensor({2, 2});
    tensor.data = {1, 2, 3, 4};
    std::deque<std::string> source{"i", "j"};
    std::deque<std::string> target{"j", "i"};

    auto result = alignToTarget(tensor, source, target);

    ASSERT_DOUBLE_EQ(result.data[0], 1);
    ASSERT_DOUBLE_EQ(result.data[1], 3);
    ASSERT_DOUBLE_EQ(result.data[2], 2);
    ASSERT_DOUBLE_EQ(result.data[3], 4);
}

TEST_F(FormulaDriverTest, AlignToTargetThrowsOnMissingIndex) {
    Tensor<double> tensor({2, 2});
    tensor.data = {1, 2, 3, 4};
    std::deque<std::string> source{"i", "j"};
    std::deque<std::string> target{"v", "j"};

    ASSERT_THROW(alignToTarget(tensor, source, target), std::runtime_error);
}

TEST_F(FormulaDriverTest, EvaluateMonomialSequenceChainsContractions) {
    std::deque<Param<double>> chain = {
        makeParam("first", {1, 2, 3, 4}, {2, 2}, {"i", "j"}),
        makeParam("second", {5, 6, 7, 8}, {2, 2}, {"k", "v"}),
        makeParam("third", {9, 10, 11, 12}, {2, 2}, {"j", "k"}),
    };

    std::deque<std::string> targetIndexNames = {"i", "v"};
    auto result = evaluateMonomialSequence(chain, targetIndexNames);

    ASSERT_DOUBLE_EQ(result.data[0], 393);
    ASSERT_DOUBLE_EQ(result.data[1], 458);
    ASSERT_DOUBLE_EQ(result.data[2], 901);
    ASSERT_DOUBLE_EQ(result.data[3], 1050);
}

TEST_F(FormulaDriverTest, EvaluateMonomialSequenceThrowsOnEmptyInput) {
    std::deque<Param<double>> empty;

    EXPECT_THROW(evaluateMonomialSequence(empty, {}), std::invalid_argument);
}
