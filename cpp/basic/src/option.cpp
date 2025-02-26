#include "basic/include/option.hpp"
#include "basic/include/log.hpp"

namespace phys_sim_complementary_dynamics {

template<typename DataType>
DataType& Option<DataType>::operator[](const std::string& key) {
    return option_[key];
}

template<typename DataType>
const DataType& Option<DataType>::at(const std::string& key) const {
    CheckCondition(HasKey(key), "basic::Option::at",
        "Missing key " + key + ".");
    return option_.at(key);
}

template<typename DataType>
const bool Option<DataType>::HasKey(const std::string& key) const {
    return option_.find(key) != option_.end();
}

template<typename DataType>
std::ostream& operator<<(std::ostream& out, const Option<DataType>& option) {
    bool leading_endl = false;
    for (const auto& pair : option()) {
        if (leading_endl) out << std::endl;
        out << pair.first << ": " << pair.second;
        leading_endl = true;
    }
    return out;
}

template class Option<integer>;
template class Option<real>;
template class Option<std::string>;
template class Option<bool>;
template class Option<MatrixXr>;
template class Option<VectorXr>;

template
std::ostream& operator<<<integer>(std::ostream& out,
    const Option<integer>& option);

template
std::ostream& operator<<<real>(std::ostream& out, const Option<real>& option);

template
std::ostream& operator<<<std::string>(std::ostream& out,
    const Option<std::string>& option);

template
std::ostream& operator<<<bool>(std::ostream& out, const Option<bool>& option);

template
std::ostream& operator<<<MatrixXr>(std::ostream& out,
    const Option<MatrixXr>& option);

template
std::ostream& operator<<<VectorXr>(std::ostream& out,
    const Option<VectorXr>& option);

}