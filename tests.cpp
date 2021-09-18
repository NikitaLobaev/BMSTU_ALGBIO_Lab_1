#include "gtest/gtest.h"
#include "Matrix.h"
#include "NeedlemanWunsch.cpp"

using Lobaev::Math::Matrix;

const std::map<char, size_t> default_matrix_map{
        {'A', 0},
        {'C', 1},
        {'G', 2},
        {'T', 3},
};

const Lobaev::Math::Matrix<int> default_matrix({
    {5, -4, -4, -4},
    {-4, 5, -4, -4},
    {-4, -4, 5, -4},
    {-4, -4, -4, 5},
});

TEST(basic, test1) {
    const std::vector<char> seq1{'A'};
    const std::vector<char> seq2{'A'};
    const int gap = -1;

    auto result_pair = needleman_wunsch(default_matrix_map, default_matrix, gap, seq1, seq2);

    auto expected_pair = std::make_pair(std::vector<char>{'A'}, 5);

    ASSERT_EQ(expected_pair, result_pair);
}

TEST(basic, test2) {
    const std::vector<char> seq1{'A', 'T'};
    const std::vector<char> seq2{'A'};
    const int gap = -1;

    auto result_pair = needleman_wunsch(default_matrix_map, default_matrix, gap, seq1, seq2);

    auto expected_pair = std::make_pair(std::vector<char>{'A'}, 5);

    ASSERT_EQ(expected_pair, result_pair);
}

TEST(basic, test3) {
    const std::vector<char> seq1{'A', 'C', 'T'};
    const std::vector<char> seq2{'A', 'G', 'T'};
    const int gap = -1;

    auto result_pair = needleman_wunsch(default_matrix_map, default_matrix, gap, seq1, seq2);

    auto expected_pair = std::make_pair(std::vector<char>{'A', 'T'}, 8);

    ASSERT_EQ(expected_pair, result_pair);
}

TEST(basic, test4) {
    const std::vector<char> seq1{'A', 'A', 'C', 'T'};
    const std::vector<char> seq2{'C', 'G', 'A', 'T'};
    const int gap = -1;

    auto result_pair = needleman_wunsch(default_matrix_map, default_matrix, gap, seq1, seq2);

    auto expected_pair = std::make_pair(std::vector<char>{'A', 'T'}, 8);

    ASSERT_EQ(expected_pair, result_pair);
}
