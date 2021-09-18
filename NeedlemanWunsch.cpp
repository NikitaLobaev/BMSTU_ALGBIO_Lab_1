#include "Matrix.cpp"

template <class T, class L> //T - тип элементов последовательности (и матриц; в частном случае char), L - числовой тип метрики
//функция возвращает пару <наилучшая последовательность, метрика>
std::pair<std::vector<T>, L> needleman_wunsch(const std::map<T, size_t> &matrix_map, //отображение элемента (в частном случае символа) в его индекс в матрице
                                              const Lobaev::Math::Matrix<L> &matrix, //матрица схожести
                                              L gap, //стоимость замены элементов последовательностей на пробел (стоимость вставки пробела)
                                              const std::vector<T> &seq1, //первая входная последовательность (в частном случае строка ДНК)
                                              const std::vector<T> &seq2  //вторая входная последовательность (в частном случае строка ДНК)
                                              ) {
    if (matrix.rows_count() != matrix_map.size() || matrix.columns_count() != matrix_map.size()) {
        throw "Weights matrix is invalid, it's size should be " +
        std::to_string(matrix_map.size()) + "*" + std::to_string(matrix_map.size());
    }

    std::vector<std::vector<L>> dynamic(seq1.size() + 1, std::vector<L>(seq2.size() + 1)); //матрица динамики

    for (size_t i = 1; i <= seq1.size(); i++) {
        for (size_t j = 1; j <= seq2.size(); j++) {
            const size_t map_index_i = matrix_map.at(seq1[i - 1]); //определение индексов в матрице схожести
            const size_t map_index_j = matrix_map.at(seq2[j - 1]); // -/-

            dynamic[i][j] = dynamic[i - 1][j - 1] + matrix(map_index_i, map_index_j); //выполняется замена символов (если они не равны)
            dynamic[i][j] = std::max(dynamic[i][j], dynamic[i - 1][j] + gap); //вставка пробела в первую последовательность
            dynamic[i][j] = std::max(dynamic[i][j], dynamic[i][j - 1] + gap); //вставка пробела во вторую последовательность
        }
    }

    size_t best_index_i = 0;
    for (size_t i = 1; i <= seq1.size(); i++) { //определение наилучшего результата в крайнем правом столбце матрицы динамики
        if (dynamic[best_index_i][seq2.size()] < dynamic[i][seq2.size()]) {
            best_index_i = i;
        }
    }

    size_t best_index_j = 0;
    for (size_t j = 1; j <= seq2.size(); j++) { //определение наилучшего результата в крайней нижней строке матрицы динамики
        if (dynamic[seq1.size()][best_index_j] < dynamic[seq1.size()][j]) {
            best_index_j = j;
        }
    }

    size_t i = seq1.size(), j = seq2.size();
    if (dynamic[best_index_i][seq2.size()] > dynamic[seq1.size()][best_index_j]) { //определение итогового наилучшего результата
        i = best_index_i;
    } else {
        j = best_index_j;
    }

    std::pair<std::vector<T>, L> result;
    result.second = dynamic[i][j];

    while (i > 0 && j > 0) { //вычисление ответа путём перемещения влево вверх по матрице динамики
        if (dynamic[i][j] == dynamic[i - 1][j] + gap) { //получили этот ответ путём вставки пробела в первую последовательность
            i--;
        } else if (dynamic[i][j] == dynamic[i][j - 1] + gap) { // -/- во вторую последовательность
            j--;
        } else { //выполняли замену
            i--;
            j--;
            result.first.emplace_back(seq1[i]); //p.s. можно брать и из второй последовательности, поскольку по умолчанию считается, что входная матрица схожести симметричная
        }
    }

    std::reverse(result.first.begin(), result.first.end()); //поскольку ответ вычислялся с конца

    return result;
}
