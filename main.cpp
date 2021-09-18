#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include "IO.cpp"
#include "Matrix.cpp"

typedef std::pair<std::string, std::string> DNA; //<dna_id, dna_string>

void read_dnas(std::istream&, std::vector<DNA>&);

template <class T, class L>
std::pair<std::vector<T>, L> needleman_wunsch(const std::map<T, size_t>&, const Lobaev::Math::Matrix<L>&, L,
                                              const std::vector<T>&, const std::vector<T>&);

const std::map<char, size_t> default_matrix_map{
        {'A', 0},
        {'R', 1},
        {'N', 2},
        {'D', 3},
        {'C', 4},
        {'Q', 5},
        {'E', 6},
        {'G', 7},
        {'H', 8},
        {'I', 9},
        {'L', 10},
        {'K', 11},
        {'M', 12},
        {'F', 13},
        {'P', 14},
        {'S', 15},
        {'T', 16},
        {'W', 17},
        {'Y', 18},
        {'V', 19},
        {'B', 20},
        {'Z', 21},
        {'X', 22},
        {'*', 23},
};

const std::string usage = "Usage: lab1 (-m | --matrix)=<input matrix filename> [(-i | --input)=<input filename>] [(-o | --output)=<output filename>] [(-g | --gap)=<gap>]";

int main(int argc, char **argv) {
    std::string input_filename, output_filename, matrix_filename;
    long gap = -2;
    for (size_t i = 1; i + 1 < argc; i += 2) {
        std::string arg(argv[i]);

        if (arg == "-i" || arg == "--input") {
            input_filename = argv[i + 1];
        } else if (arg == "-o" || arg == "--output") {
            output_filename = argv[i + 1];
        } else if (arg == "-g" || arg == "--gap") {
            gap = std::stol(argv[i + 1]);
        } else if (arg == "-m" || arg == "--matrix") {
            matrix_filename = argv[i + 1];
        } else {
            std::cerr << usage << std::endl;
            return 1;
        }
    }

    if (matrix_filename.empty()) {
        std::cerr << usage << std::endl;
        return 1;
    }

    if (!input_filename.empty()) {
        auto file = std::freopen(input_filename.c_str(), "r", stdin);
        if (!file) {
            std::cerr << "Input file doesn't exist." << std::endl << usage << std::endl;
            return 1;
        }
    }

    if (!output_filename.empty()) {
        auto file = std::freopen(output_filename.c_str(), "w", stdout);
        if (!file) {
            std::cerr << "Output file doesn't exist." << std::endl << usage << std::endl;
            return 1;
        }
    }

    std::ifstream matrix_ifstream(matrix_filename);
    if (!matrix_ifstream) {
        std::cerr << "Input matrix file doesn't exist." << std::endl << usage << std::endl;
        return 1;
    }

    auto matrix = Lobaev::IO::read_matrix<long>(matrix_ifstream);

    matrix_ifstream.close();

    if (matrix.rows_count() != default_matrix_map.size() || matrix.columns_count() != default_matrix_map.size()) {
        std::cerr << "Input matrix is invalid. Size should be (" <<
        default_matrix_map.size() << 'x' << default_matrix_map.size() << ')' << std::endl << usage << std::endl;
        return 1;
    }

    std::vector<DNA> dnas;
    read_dnas(std::cin, dnas);

    if (dnas.size() != 2) {
        std::cerr << "Only 2 dna's in input file are allowed." << std::endl << usage << std::endl;
        return 1;
    }

    const std::pair<std::vector<char>, long> result_pair = needleman_wunsch(default_matrix_map, matrix, gap,
                                                                            std::vector<char>(dnas[0].second.begin(), dnas[0].second.end()),
                                                                                    std::vector<char>(dnas[1].second.begin(), dnas[1].second.end()));
    const std::string result_string = std::string(result_pair.first.begin(), result_pair.first.end());
    const long result_score = result_pair.second;

    std::cout << result_string << std::endl;
    std::cout << "Score: " << result_score << std::endl;

    return 0;
}

template <class T, class L>
std::pair<std::vector<T>, L> needleman_wunsch(const std::map<T, size_t> &matrix_map,
                                              const Lobaev::Math::Matrix<L> &matrix, L gap, const std::vector<T> &seq1,
                                              const std::vector<T> &seq2) {
    std::vector<std::vector<L>> dynamic(seq1.size() + 1, std::vector<L>(seq2.size() + 1));

    for (size_t i = 1; i <= seq1.size(); i++) {
        for (size_t j = 1; j <= seq2.size(); j++) {
            const size_t map_index_i = matrix_map.at(seq1[i - 1]);
            const size_t map_index_j = matrix_map.at(seq1[i - 1]);

            dynamic[i][j] = dynamic[i - 1][j - 1] + matrix(map_index_i, map_index_j);
            dynamic[i][j] = std::max(dynamic[i][j], dynamic[i - 1][j] + gap);
            dynamic[i][j] = std::max(dynamic[i][j], dynamic[i][j - 1] + gap);
        }
    }

    size_t best_index_i = 0;
    for (size_t i = 1; i <= seq1.size(); i++) {
        if (dynamic[best_index_i][seq2.size()] < dynamic[i][seq2.size()]) {
            best_index_i = i;
        }
    }

    size_t best_index_j = 0;
    for (size_t j = 1; j <= seq2.size(); j++) {
        if (dynamic[seq1.size()][best_index_j] < dynamic[seq1.size()][j]) {
            best_index_j = j;
        }
    }

    size_t i = seq1.size(), j = seq2.size();
    if (dynamic[best_index_i][seq2.size()] > dynamic[seq1.size()][best_index_j]) {
        i = best_index_i;
    } else {
        j = best_index_j;
    }

    std::pair<std::vector<T>, L> result;
    result.second = dynamic[i][j];

    while (i > 0 && j > 0) {
        if (dynamic[i][j] == dynamic[i - 1][j] + gap) {
            i--;
        } else if (dynamic[i][j] == dynamic[i][j - 1] + gap) {
            j--;
        } else {
            i--;
            j--;
            result.first.emplace_back(seq1[i]);
        }
    }
    return result;
}

void read_dnas(std::istream &in, std::vector<DNA> &dnas) {
    std::string cur_line, cur_dna_buf;
    while (std::getline(in, cur_line)) {
        if (cur_line.empty()) {
            continue;
        }

        if (cur_line[0] == '>') {
            if (!cur_dna_buf.empty()) {
                dnas.back().second = cur_dna_buf;
                cur_dna_buf.clear();
            }

            const size_t id_start_index = cur_line.find('|');
            const size_t id_length = cur_line.substr(id_start_index + 1).find('|');
            const std::string id = cur_line.substr(id_start_index + 1, id_length);
            dnas.emplace_back(id, "");
        } else {
            cur_dna_buf += cur_line;
        }
    }

    if (!cur_dna_buf.empty()) {
        dnas.back().second = cur_dna_buf;
    }
}
