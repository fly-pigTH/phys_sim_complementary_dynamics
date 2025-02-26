#ifndef PHYS_SIM_COMPLEMENTARY_DYNAMICS_BASIC_OPTIONS
#define PHYS_SIM_COMPLEMENTARY_DYNAMICS_BASIC_OPTIONS

#include "basic/include/option.hpp"

namespace phys_sim_complementary_dynamics {

class Options {
public:
    Options() {}
    ~Options() {}

    const Option<integer>& integer_option() const { return integer_option_; }
    Option<integer>& integer_option() { return integer_option_; }
    const Option<real>& real_option() const { return real_option_; }
    Option<real>& real_option() { return real_option_; }
    const Option<std::string>& string_option() const { return string_option_; }
    Option<std::string>& string_option() { return string_option_; }
    const Option<bool>& bool_option() const { return bool_option_; }
    Option<bool>& bool_option() { return bool_option_; }
    const Option<MatrixXr>& matrix_option() const { return matrix_option_; }
    Option<MatrixXr>& matrix_option() { return matrix_option_; }
    const Option<VectorXr>& vector_option() const { return vector_option_; }
    Option<VectorXr>& vector_option() { return vector_option_; }
    const integer GetMatrixOptionRows(const std::string& key) const;
    const integer GetMatrixOptionCols(const std::string& key) const;
    const integer GetVectorOptionSize(const std::string& key) const;

    void Clear();

private:
    Option<integer> integer_option_;
    Option<real> real_option_;
    Option<std::string> string_option_;
    Option<bool> bool_option_;
    Option<MatrixXr> matrix_option_;
    Option<VectorXr> vector_option_;
};

}

#endif