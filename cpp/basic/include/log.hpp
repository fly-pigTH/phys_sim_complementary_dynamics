#ifndef PHYS_SIM_COMPLEMENTARY_DYNAMICS_BASIC_LOG
#define PHYS_SIM_COMPLEMENTARY_DYNAMICS_BASIC_LOG

#include "basic/include/config.hpp"

namespace phys_sim_complementary_dynamics {

void PrintInfo(const std::string& location, const std::string& message);
void PrintWarning(const std::string& location, const std::string& message);
void PrintError(const std::string& location, const std::string& message);
void PrintSuccess(const std::string& location, const std::string& message);

void CheckCondition(const bool condition, const std::string& location,
    const std::string& message);

const bool StartWith(const std::string& full, const std::string& prefix);
const bool EndWith(const std::string& full, const std::string& suffix);

// Timing.
void Tic();
void Toc(const std::string& location, const std::string& message);
const real Toc();

// Load and save sparse matrices.
template<typename DataType>
void Save(std::ofstream& f, const DataType& val);
template<typename DataType>
const DataType Load(std::ifstream& f);

void SaveVectorXr(const std::string& file_name, const VectorXr& vec);
const VectorXr LoadVectorXr(const std::string& file_name);

template<typename DataType>
const Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>
LoadMatrix(const std::string& file_name) {
    std::ifstream file(file_name);

    integer row = 0, col = 0;
    file >> row >> col;

    Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> matrix(row, col);
    for (integer i = 0; i < row; ++i)
        for (integer j = 0; j < col; ++j)
            file >> matrix(i, j);

    return matrix;
}

const std::pair<Matrix3Xr, Matrix3Xi> LoadObjMesh(const std::string& file_name);
// Parse the .mesh file from Complementary Dynamics.
// https://libigl.github.io/file-formats/.
struct MeditMesh {
public:
    integer vertex_num;
    Matrix3Xr vertices;

    integer triangle_num;
    Matrix3Xi triangles;

    integer tetrahedron_num;
    Matrix4Xi tetrahedrons;
};

const MeditMesh LoadMeditMesh(const std::string& file_name);

}

#endif