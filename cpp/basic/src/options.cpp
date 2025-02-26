#include "basic/include/options.hpp"
#include "basic/include/log.hpp"

namespace phys_sim_complementary_dynamics {

const integer Options::GetMatrixOptionRows(const std::string& key) const {
    return static_cast<integer>(matrix_option_.at(key).rows());
}

const integer Options::GetMatrixOptionCols(const std::string& key) const {
    return static_cast<integer>(matrix_option_.at(key).cols());
}

const integer Options::GetVectorOptionSize(const std::string& key) const {
    return static_cast<integer>(vector_option_.at(key).size());
}

void Options::Clear() {
    integer_option_.Clear();
    real_option_.Clear();
    string_option_.Clear();
    bool_option_.Clear();
    matrix_option_.Clear();
    vector_option_.Clear();
}

}