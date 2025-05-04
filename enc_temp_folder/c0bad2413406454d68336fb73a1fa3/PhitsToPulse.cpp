#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <nlohmann/json.hpp>
#include <iomanip>
#include <thread>
#include <ppl.h> 
#include<algorithm>
#include <mutex>
#include "Batch2Pulse.h"
#include"Dump2Batch.h"
#include <concurrent_vector.h>
#include "SpinProgress.hpp"

void to_json(nlohmann::json& j, const EventInfo& e) {
	j = nlohmann::json{
		{"ityp", e.ityp},
		{"x", e.x}, {"y", e.y}, {"z", e.z}, {"E", e.E},
		{"x_deposit", e.x_deposit}, {"y_deposit", e.y_deposit},
		{"z_deposit", e.z_deposit}, {"E_deposit", e.E_deposit}
	};
}

int main()
{
	std::string DataPath;
	std::cout << "Data Path:";
	std::cin >> DataPath;
	std::string InputPath = DataPath + "/input.json";

	InputParameters InputPara = ReadInputJson(InputPath);

	std::vector<std::string> DumpPathes;

	for (const auto& posi : InputPara.positions) {
		DumpPathes.push_back(DataPath + "/" + std::to_string(static_cast<int>(InputPara.E)) + "keV_" + std::to_string(posi));
	}

	PulseParameters PulsePara(InputPara);

	int TotalCounter = 1;

	Eigen::MatrixXd Matrix_A = MakeMatrix_A(PulsePara, InputPara);

	Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(Matrix_A);
	// 固有値
	Eigen::VectorXcd EigenValues = eigensolver.eigenvalues();
	// 固有ベクトル
	Eigen::MatrixXcd EigenVectors = eigensolver.eigenvectors();

	// 各固有ベクトルを正規化する
	for (int i = 0; i < EigenVectors.cols(); ++i) {
		std::complex<double> norm = EigenVectors.col(i).norm();  // ノルム計算
		if (std::abs(norm) > 0) {  // ノルムがゼロでないか確認
			EigenVectors.col(i) /= norm;  // 固有ベクトルをそのノルムで割って正規化
		}
	}

	const int n_abs = InputPara.n_abs;
	const int n_abs_1 = n_abs + 1;
	const int n_abs_2 = n_abs + 2;
	const int n_abs_3 = n_abs + 3;
	const int n_abs_4 = n_abs + 4;

	const std::vector<double> Block = linspace(-InputPara.length / 20, InputPara.length / 20, n_abs + 1);

	std::vector<double> time = linspace(0, InputPara.samples / InputPara.rate, static_cast<int>(InputPara.samples));

	for (const auto& posi : InputPara.positions) {
		Eigen::VectorXd VecInit(n_abs_4);
		VecInit.setZero();
		VecInit[posi + 1] = InputPara.E * 1e3 * PulsePara.e_const / PulsePara.C_abs;
		std::cout << "Energy:" << InputPara.E * 1e3 << ",e_const:" << PulsePara.e_const << ",C_abs:" << PulsePara.C_abs << "\n";
		Eigen::MatrixXd EiVec = EigenVectors.real();
		Eigen::VectorXd Consts = EiVec.colPivHouseholderQr().solve(VecInit);

		Eigen::VectorXd Pulse0(static_cast<int>(InputPara.samples));
		Pulse0.setZero();

		std::vector<double> values(Consts.data(), Consts.data() + Consts.size());

		for (int i = 0; i < Consts.size(); i++) {
			for (int j = 0; j < time.size(); j++) {
				Pulse0[j] += Consts[i] * EiVec(0, i) * std::exp(EigenValues[i].real() * time[j]);
			}
		}

		Eigen::VectorXd Pulse1(static_cast<int>(InputPara.samples));
		Pulse1.setZero();

		for (int i = 0; i < Consts.size(); i++) {
			for (int j = 0; j < time.size(); j++) {
				Pulse1[j] += Consts[i] * EiVec(n_abs_3, i) * std::exp(EigenValues[i].real() * time[j]);
			}
		}

		std::string PulsePath = DataPath + "/" + std::to_string(static_cast<int>(InputPara.E)) + "keV_" + std::to_string(posi) + "/Pulse";

		std::string PulseFile_0 = PulsePath + "/Ch0/CH0_0.dat";
		std::string PulseFile_1 = PulsePath + "/Ch1/CH1_0.dat";

		std::ofstream PulseoutFile_0(PulseFile_0);
		if (!PulseoutFile_0) {
			std::cerr << "Failed to open file:" << PulseFile_0 << std::endl;
			return -1;
		}
		for (int i = 0; i < Pulse0.size(); ++i) {
			PulseoutFile_0 << Pulse0[i] << std::endl;
		}
		PulseoutFile_0.close();

		std::ofstream PulseoutFile_1(PulseFile_1);
		if (!PulseoutFile_1) {
			std::cerr << "Failed to open file:" << PulseFile_1 << std::endl;
			return -1;
		}
		for (int i = 0; i < Pulse0.size(); ++i) {
			PulseoutFile_1 << Pulse0[i] << std::endl;
		}
		PulseoutFile_1.close();
	}

	for (const auto& DumpPath : DumpPathes) {
		std::cout << TotalCounter << "/" << DumpPathes.size() << "\n";

		std::map<int, std::map<int, EventInfo>> batch;
		std::string jsonPath = DumpPath + "/batch.json";
		std::string dumpPath = DumpPath + "/dumpall.dat";

		
		std::filesystem::create_directories(DumpPath + "/Pulse/Ch0");
		std::filesystem::create_directories(DumpPath + "/Pulse/CH1");
		std::filesystem::create_directories(DumpPath + "/Pulse_ms/Ch0");
		std::filesystem::create_directories(DumpPath + "/Pulse_ms/CH1");

		int Counter = 0;

		concurrency::concurrent_vector<std::tuple<int, double, double>> PulseInfo_Ch0;
		concurrency::concurrent_vector<std::tuple<int, double, double>> PulseInfo_Ch1;

		int Input = -1;
		while (true) {
			std::cout << "Continue to analyze phits result? [0]yes [1]no\n";
			std::cin >> Input;

			if (std::cin.fail()) {
				std::cin.clear(); // エラー状態クリア
				std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // 入力バッファクリア
				std::cout << "Invalid (non-numeric) input. Please enter 0 or 1.\n";
				continue;
			}

			if (Input == 0 || Input == 1) {
				break; // 正常値なのでループを抜ける
			}
			else {
				std::cout << "Invalid Input. Please enter 0 or 1.\n";
			}
		}

		if (Input == 0) {
			// JSONファイルが存在するか確認
			if (std::filesystem::exists(jsonPath)) {
				SpinProgress spinner;
				spinner.set_message("Processing batch file...");
				// JSONから読み込み
				std::ifstream file(jsonPath);
				nlohmann::json json_batch;
				file >> json_batch;

				for (auto& [key1, inner] : json_batch.items()) {
					int int_key1 = std::stoi(key1);
					for (auto& [key2, event_json] : inner.items()) {
						int int_key2 = std::stoi(key2);
						EventInfo event = event_json.get<EventInfo>();
						batch[int_key1][int_key2] = event;
					}
				}
				spinner.complete("Loaded from batch");

			}
			else {
				// dumpファイルから読み込み
				SpinProgress spinner;
				spinner.set_message("Processing dumpall file...");
				int ReadReturn = ReadDump(dumpPath, batch, InputPara.E / 1000);
				spinner.complete("Loaded from dump");

				if (ReadReturn == -1) {
					std::cout << "Error in dump file\n";
					return -1;
				}

				// JSONとして保存
				nlohmann::json json_batch;
				for (const auto& [key1, map_inner] : batch) {
					for (const auto& [key2, event] : map_inner) {
						json_batch[std::to_string(key1)][std::to_string(key2)] = event;
					}
				}
				std::ofstream file(jsonPath);
				file << json_batch.dump(4);
				file.close();
			}
		}
		else if (Input == 1) {
			return 2;
		}

		size_t total_items = batch.size(); // 全体の要素数を取得
		std::atomic<size_t> completed_items(0); // 完了したアイテム数
		std::mutex output_mutex; // 出力用ミューテックス

		concurrency::parallel_for_each(batch.begin(), batch.end(), [&](const std::pair<const int, std::map<int, EventInfo>>& outer_pair) {

			std::vector<double> BlockDeposit(n_abs, 0.0);

			for (const auto& inner_pair : outer_pair.second) {
				for (int i = 0; i < inner_pair.second.x_deposit.size(); i++) {
					int EnergyedBlocks = AssignBlock(InputPara, Block, inner_pair.second.x_deposit[i], inner_pair.second.y_deposit[i], inner_pair.second.z_deposit[i]);
					BlockDeposit[EnergyedBlocks - 1] += inner_pair.second.E_deposit[i]*1e6;
				}
			}
			
			std::vector<int> EnergyedBlocks;
			int count = 1;

			for (const double& pos : BlockDeposit) {
				if (pos > 0) {
					EnergyedBlocks.push_back(count);
				}
				count++;
			}

			Eigen::VectorXd VecInit(n_abs_4);
			VecInit.setZero();

			for (int i=0; i < BlockDeposit.size(); i++) {
				VecInit[i + 1] = BlockDeposit[i] * PulsePara.e_const / PulsePara.C_abs;
			}
			
			Eigen::MatrixXd EiVec = EigenVectors.real();

			Eigen::VectorXd Consts = EiVec.colPivHouseholderQr().solve(VecInit);

			Eigen::VectorXd Pulse0(static_cast<int>(InputPara.samples));
			Pulse0.setZero();
			
			for (int i = 0; i < Consts.size();i++) {
				for (int j = 0; j < time.size(); j++) {
					Pulse0[j] += Consts[i] * EiVec(0,i) *std::exp(EigenValues[i].real() * time[j]);
				}
			}
			
			Eigen::VectorXd Pulse1(static_cast<int>(InputPara.samples));
			Pulse1.setZero();

			for (int i = 0; i < Consts.size(); i++) {
				for (int j = 0; j < time.size(); j++) {
					Pulse1[j] += Consts[i] * EiVec(n_abs_3, i) * std::exp(EigenValues[i].real() * time[j]);
				}
			}

			std::string PulseFile_0 = DumpPath + "/Pulse_ms/Ch0/CH0_" + std::to_string(outer_pair.first) + ".dat";
			std::string PulseFile_1 = DumpPath + "/Pulse_ms/Ch1/CH1_" + std::to_string(outer_pair.first) + ".dat";

			std::ofstream PulseoutFile_0(PulseFile_0);
			if (!PulseoutFile_0) {
				std::cerr << "Failed to open file:" << PulseFile_0 << std::endl;
				return -1;
			}
			for (int i = 0; i < Pulse0.size(); ++i) {
				PulseoutFile_0 << Pulse0[i] << std::endl;
			}
			PulseoutFile_0.close();

			std::ofstream PulseoutFile_1(PulseFile_1);
			if (!PulseoutFile_1) {
				std::cerr << "Failed to open file:" << PulseFile_1 << std::endl;
				return -1;
			}
			for (int i = 0; i < Pulse0.size(); ++i) {
				PulseoutFile_1 << Pulse0[i] << std::endl;
			}
			PulseoutFile_1.close();

			size_t completed = completed_items.fetch_add(1);

			// 進捗を更新（1%ごと）
			size_t progress = static_cast<size_t>((completed + 1) * 100 / total_items); // +1は現在のアイテムを含むため
			if (progress > 0 && progress % 1 == 0) { // 1%ごとに更新
				std::lock_guard<std::mutex> lock(output_mutex); // スレッドセーフな出力
				std::cout << "\rConverting to Pulse: " << progress << "%" << std::flush; // プログレスを表示
			}
			});

		{
			std::lock_guard<std::mutex> lock(output_mutex);
			std::cout << "\nFinished\n";
		}
		TotalCounter++;
	}
	std::string Input;
	std::cout << "Key input to close...\n";
	std::cin >> Input;
	return 0;
}