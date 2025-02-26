#ifndef PHYS_SIM_COMPLEMENTARY_DYNAMICS_BASIC_OPTION
#define PHYS_SIM_COMPLEMENTARY_DYNAMICS_BASIC_OPTION

#include "basic/include/config.hpp"

namespace phys_sim_complementary_dynamics {

template<typename DataType>
class Option {
public:
    Option() : option_() {}
    ~Option() {}

    DataType& operator[](const std::string& key);
    const DataType& at(const std::string& key) const;
    const std::map<std::string, DataType>& operator()() const {
        return option_;
    }
    void Clear() { return option_.clear(); }

private:
    const bool HasKey(const std::string& key) const;
    std::map<std::string, DataType> option_;
};

template<typename DataType>
std::ostream& operator<<(std::ostream& out, const Option<DataType>& option);

}

#endif