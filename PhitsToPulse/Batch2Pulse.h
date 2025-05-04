#pragma once
#include <vector>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "struct.h"

struct InputParameters;

// パルスに関する構造体
struct PulseParameters {
    const double k_b = 1.381 * 1.0e-23;
    const double ptfn_Flink = 0.5;
    const double e_const = 1.602e-19;

    const double C_abs;
    const double G_abs_abs;
    const double I;
    const double t_el;
    const double L_I;
    const double t_I;

    PulseParameters(const InputParameters& InputPara);
};

// 等差数列を作成する関数
std::vector<double> linspace(double start, double stop, int num);

// 入射したブロックを見つける関数（だと思う）
int AssignBlock(const InputParameters& InputPara, const std::vector<double>& Block, const double& x_deposit, const double& y_deposit, const double& z_deposit);

Eigen::MatrixXd MakeMatrix_A(const PulseParameters& PulsePara, const InputParameters& InputPara);

Eigen::MatrixXd MakeMatrix_X(const PulseParameters& PulsePara, const InputParameters& InputPara, const std::vector<int>& pixel, const std::vector<double>& Positions);

// 比較関数：ベクトルの最初の要素の絶対値に基づいて比較
bool compareByMagnitude(int i, int j, const Eigen::VectorXcd& eigenvalues);

void SortEigen(Eigen::VectorXcd& EigenValues, Eigen::MatrixXcd& EigenVectors);

void checkImaginaryPart(const Eigen::VectorXcd& vec);
void checkImaginaryPart(const Eigen::MatrixXcd& mat);