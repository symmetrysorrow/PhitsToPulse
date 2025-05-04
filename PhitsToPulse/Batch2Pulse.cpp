#pragma once

#include "Batch2Pulse.h"
#include<filesystem>
#include <iostream>
#include <numeric>

#include"Dump2Batch.h"

PulseParameters::PulseParameters(const InputParameters& InputPara)
	: C_abs(InputPara.C_abs / InputPara.n_abs),
	G_abs_abs(InputPara.G_abs_abs* (InputPara.n_abs - 1)),
	I(std::sqrt((InputPara.G_tes_bath* InputPara.T_c* (1 - std::pow(InputPara.T_bath / InputPara.T_c, InputPara.n))) / (InputPara.n * InputPara.R))),
	t_el(InputPara.L / (InputPara.R_l + InputPara.R * (1 + InputPara.beta))),
	L_I((InputPara.alpha* std::pow(I, 2)* InputPara.R) / (InputPara.G_tes_bath * InputPara.T_c)),
	t_I(InputPara.C_tes / ((1 - L_I) * InputPara.G_tes_bath)) {
}

// 等差数列を作成する関数
std::vector<double> linspace(double start, double stop, int num) {
	std::vector<double> values;
	if (num <= 0) return values;

	double step = (stop - start) / (num - 1);

	for (int i = 0; i < num; ++i) {
		values.push_back(start + i * step);
	}

	return values;
}

// なぞ関数1
int AssignBlock(const InputParameters& InputPara, const std::vector<double>& Block, const double& x_deposit, const double& y_deposit, const double& z_deposit) {
	if ((-InputPara.length / 20 <= x_deposit || x_deposit <= InputPara.length / 20) && (-InputPara.depth / 20 <= y_deposit || y_deposit <= InputPara.depth / 20) && (-InputPara.height / 10 <= z_deposit || z_deposit <= 0.0)) {
		for (int i = 0; i < Block.size(); ++i) {
			if (Block[i] >= x_deposit) {
				return i;
			}
		}
	}
	return 0;
}

Eigen::MatrixXd MakeMatrix_A(const PulseParameters& PulsePara, const InputParameters& InputPara)
{
	const int n_abs = InputPara.n_abs;
	const int n_abs_1 = n_abs + 1;
	const int n_abs_2 = n_abs + 2;
	const int n_abs_3 = n_abs + 3;
	const int n_abs_4 = n_abs + 4;

	Eigen::MatrixXd Matrix_A(n_abs_4, n_abs_4);
	Matrix_A = Eigen::MatrixXd::Zero(n_abs_4, n_abs_4);

	for (int i = 0; i < n_abs_4; i++) {
		for (int j = 0; j < n_abs_4; j++) {
			if (j == i - 1) {
				Matrix_A(i, j) = -PulsePara.G_abs_abs / PulsePara.C_abs;
			}
			if (j == i) {
				Matrix_A(i, j) = 2 * PulsePara.G_abs_abs / PulsePara.C_abs;
			}
			if (j == i + 1) {
				Matrix_A(i, j) = -PulsePara.G_abs_abs / PulsePara.C_abs;
			}
		}
	}

	Matrix_A(0, 0) = 1 / PulsePara.t_el;
	Matrix_A(0, 1) = PulsePara.L_I * InputPara.G_tes_bath / (PulsePara.I * InputPara.L);

	Matrix_A(1, 0) = -PulsePara.I * InputPara.R * (2 + InputPara.beta) / InputPara.C_tes;
	Matrix_A(1, 1) = 1 / PulsePara.t_I + (InputPara.G_abs_tes / InputPara.C_tes);
	Matrix_A(1, 2) = -InputPara.G_abs_tes / InputPara.C_tes;

	Matrix_A(2, 1) = -InputPara.G_abs_tes / PulsePara.C_abs;
	Matrix_A(2, 2) = InputPara.G_abs_tes / PulsePara.C_abs + PulsePara.G_abs_abs / PulsePara.C_abs;
	Matrix_A(2, 3) = -PulsePara.G_abs_abs / PulsePara.C_abs;

	Matrix_A(n_abs_1, n_abs) = -PulsePara.G_abs_abs / PulsePara.C_abs;
	Matrix_A(n_abs_1, n_abs_1) = InputPara.G_abs_tes / PulsePara.C_abs + PulsePara.G_abs_abs / PulsePara.C_abs;
	Matrix_A(n_abs_1, n_abs_2) = -InputPara.G_abs_tes / PulsePara.C_abs;

	Matrix_A(n_abs_2, n_abs_1) = -InputPara.G_abs_tes / InputPara.C_tes;
	Matrix_A(n_abs_2, n_abs_2) = 1 / PulsePara.t_I + InputPara.G_abs_tes / InputPara.C_tes;
	Matrix_A(n_abs_2, n_abs_3) = -PulsePara.I * InputPara.R * (2 + InputPara.beta) / InputPara.C_tes;

	Matrix_A(n_abs_3, n_abs_2) = PulsePara.L_I * InputPara.G_tes_bath / (PulsePara.I * InputPara.L);
	Matrix_A(n_abs_3, n_abs_3) = 1 / PulsePara.t_el;

	Matrix_A *= -1;

	return Matrix_A;
}

Eigen::MatrixXd MakeMatrix_X(const PulseParameters& PulsePara, const InputParameters& InputPara, const std::vector<int>& pixel, const std::vector<double>& Positions)
{
	const int n_abs = InputPara.n_abs;
	const int n_abs_1 = n_abs + 1;
	const int n_abs_2 = n_abs + 2;
	const int n_abs_3 = n_abs + 3;
	const int n_abs_4 = n_abs + 4;

	Eigen::MatrixXd Matrix_X(n_abs_2, n_abs_4);
	Matrix_X = Eigen::MatrixXd::Zero(n_abs_2, n_abs_4);

	const double Same = InputPara.E * 1e3 * PulsePara.e_const / InputPara.C_tes;

	for (const int& pix : pixel){
		if (pix <= n_abs_2){
			Matrix_X(pix, pix + 1) = Positions[pix - 1] * PulsePara.e_const / PulsePara.C_abs;
		}
	}

	Matrix_X(0, 1) = Same;
	Matrix_X(n_abs_1, n_abs_2) = Same;

	return Matrix_X;
}

bool compareByMagnitude(int i, int j, const Eigen::VectorXcd& eigenvalues)
{
	return std::abs(eigenvalues(i)) < std::abs(eigenvalues(j));
}

void SortEigen(Eigen::VectorXcd& EigenValues, Eigen::MatrixXcd& EigenVectors)
{
	std::vector<int> indices(EigenValues.size());
	std::iota(indices.begin(), indices.end(), 0);

	std::sort(indices.begin(), indices.end(),
		[&EigenValues](int i, int j) { return compareByMagnitude(i, j, EigenValues); });

	// ソートされた固有値と固有ベクトルを作成
	Eigen::VectorXcd sorted_eigenvalues(EigenValues.size());
	Eigen::MatrixXcd sorted_eigenvectors(EigenVectors.rows(), EigenVectors.cols());

	for (int i = 0; i < indices.size(); ++i) {
		sorted_eigenvalues(i) = EigenValues(indices[i]);
		sorted_eigenvectors.col(i) = EigenVectors.col(indices[i]);
	}

	EigenValues = sorted_eigenvalues;
	EigenVectors = sorted_eigenvectors;

	// 各列（固有ベクトル）を正規化
	for (int i = 0; i < EigenVectors.cols(); ++i) {
		EigenVectors.col(i).normalize();
	}

	// 固有ベクトルを調整し、最初の要素が正になるようにする
	for (int i = 0; i < EigenVectors.cols(); ++i) {
		if (EigenVectors(0, i).real() < 0) {
			EigenVectors.col(i) = -EigenVectors.col(i); // ベクトル全体を反転
		}
	}

	return;
}

void checkImaginaryPart(const Eigen::VectorXcd& vec)
{
	bool hasImaginaryPart = (vec.imag().array() != 0).any();

	if (hasImaginaryPart) {
		std::cout << "Warning: Vector contains non-zero imaginary parts!" << std::endl;
	}
	else {
		std::cout << "No non-zero imaginary parts detected in the vector." << std::endl;
	}
}

void checkImaginaryPart(const Eigen::MatrixXcd& mat)
{
	bool hasImaginaryPart = (mat.imag().array() != 0).any();

	if (hasImaginaryPart) {
		std::cout << "Warning: Matrix contains non-zero imaginary parts!" << std::endl;
	}
	else {
		std::cout << "No non-zero imaginary parts detected in the matrix." << std::endl;
	}
}