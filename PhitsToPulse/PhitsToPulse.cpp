#pragma once

#include <algorithm>
#include <atomic>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <string>
#include <tuple>
#include <vector>

#include <nlohmann/json.hpp>
#include <ppl.h>

#include "Batch2Pulse.h"
#include "Dump2Batch.h"
#include "SpinProgress.hpp"

namespace fs = std::filesystem;

// ------------------------------------------------------------
// JSON
// ------------------------------------------------------------
void to_json(nlohmann::json& j, const EventInfo& e) {
    j = nlohmann::json{
        {"ityp", e.ityp},
        {"x", e.x}, {"y", e.y}, {"z", e.z}, {"E", e.E},
        {"x_deposit", e.x_deposit}, {"y_deposit", e.y_deposit},
        {"z_deposit", e.z_deposit}, {"E_deposit", e.E_deposit}
    };
}

// ------------------------------------------------------------
// Utility
// ------------------------------------------------------------
static bool WriteVectorToFile(const std::string& filename, const Eigen::Ref<const Eigen::VectorXd>& vec) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << filename << std::endl;
        return false;
    }

    file << std::setprecision(16);
    for (Eigen::Index i = 0; i < vec.size(); ++i) {
        file << vec[i] << '\n';
    }
    return true;
}

static void exportTempsCSV(
    const std::string& filename,
    const std::vector<double>& Block,
    const std::vector<double>& time,
    const Eigen::MatrixXd& Temps // [samples x n_abs]
) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }

    file << std::setprecision(16);

    file << "time";
    for (Eigen::Index i = 0; i < Temps.cols(); ++i) {
        file << "," << Block[static_cast<size_t>(i)];
    }
    file << "\n";

    for (Eigen::Index j = 0; j < Temps.rows(); ++j) {
        file << time[static_cast<size_t>(j)];
        for (Eigen::Index i = 0; i < Temps.cols(); ++i) {
            file << "," << Temps(j, i);
        }
        file << "\n";
    }

    std::cout << "Exported to " << filename << std::endl;
}

static Eigen::MatrixXd BuildExpTable(
    const Eigen::VectorXcd& eigenValues,
    const std::vector<double>& time
) {
    const int modes = static_cast<int>(eigenValues.size());
    const int samples = static_cast<int>(time.size());

    Eigen::MatrixXd expTable(modes, samples);

    for (int k = 0; k < modes; ++k) {
        const double lambda = eigenValues[k].real();
        for (int j = 0; j < samples; ++j) {
            expTable(k, j) = std::exp(lambda * time[static_cast<size_t>(j)]);
        }
    }

    return expTable;
}

static Eigen::VectorXd BuildPulseFromConsts(
    const Eigen::VectorXd& consts,
    const Eigen::VectorXd& channelEigenRow,
    const Eigen::MatrixXd& expTable
) {
    // coeff[k] = consts[k] * channelEigenRow[k]
    const Eigen::RowVectorXd coeff =
        (consts.array() * channelEigenRow.array()).matrix().transpose();

    // [1 x modes] * [modes x samples] -> [1 x samples]
    return (coeff * expTable).transpose();
}

static Eigen::MatrixXd BuildAllBlockPulses(
    const Eigen::MatrixXd& constsAll,      // [n_abs_4 x n_abs]
    const Eigen::VectorXd& channelEigenRow, // [n_abs_4]
    const Eigen::MatrixXd& expTable         // [n_abs_4 x samples]
) {
    // 各列 i について pulse_i = ((constsAll.col(i) .* channelEigenRow)^T * expTable)^T
    // まとめると:
    // coeffMat = channelEigenRow.asDiagonal() * constsAll   [modes x n_abs]
    // pulses^T = coeffMat^T * expTable                      [n_abs x samples]
    // pulses   = (coeffMat^T * expTable)^T                 [samples x n_abs]

    const Eigen::MatrixXd coeffMat = channelEigenRow.asDiagonal() * constsAll;
    return (coeffMat.transpose() * expTable).transpose();
}

static bool LoadOrCreateBatch(
    const std::string& jsonPath,
    const std::string& dumpPath,
    const InputParameters& inputPara,
    std::map<int, std::map<int, EventInfo>>& batch
) {
    if (fs::exists(jsonPath)) {
        SpinProgress spinner;
        spinner.set_message("Processing batch file...");

        std::ifstream file(jsonPath);
        if (!file.is_open()) {
            std::cerr << "Failed to open " << jsonPath << std::endl;
            return false;
        }

        nlohmann::json json_batch;
        file >> json_batch;

        for (auto& [key1, inner] : json_batch.items()) {
            const int int_key1 = std::stoi(key1);
            for (auto& [key2, event_json] : inner.items()) {
                const int int_key2 = std::stoi(key2);
                batch[int_key1][int_key2] = event_json.get<EventInfo>();
            }
        }

        spinner.complete("Loaded from batch");
        return true;
    }

    SpinProgress spinner;
    spinner.set_message("Processing dumpall file...");

    std::vector<int> FullEnergyList;
    const int readReturn = ReadDump(
        dumpPath,
        batch,
        inputPara.E / 1000,
        inputPara.SaveAll,
        FullEnergyList
    );

    if (readReturn == -1) {
        std::cout << "Error in dump file\n";
        return false;
    }

    if (inputPara.SaveAll) {
        const std::string indexListFile = fs::path(dumpPath).parent_path().string() + "/FullEnergyList.dat";
        std::ofstream indexListS(indexListFile);
        if (!indexListS.is_open()) {
            std::cerr << "Failed to open " << indexListFile << std::endl;
            return false;
        }
        for (const int& index : FullEnergyList) {
            indexListS << index << "\n";
        }
    }

    spinner.complete("Loaded from dump");

    nlohmann::json json_batch;
    for (const auto& [key1, map_inner] : batch) {
        for (const auto& [key2, event] : map_inner) {
            json_batch[std::to_string(key1)][std::to_string(key2)] = event;
        }
    }

    std::ofstream file(jsonPath);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << jsonPath << " for writing" << std::endl;
        return false;
    }
    file << json_batch.dump(4);

    return true;
}

// ------------------------------------------------------------
// main
// ------------------------------------------------------------
int main() {
    std::string DataPath;
    std::cout << "Data Path:";
    std::cin >> DataPath;

    const std::string InputPath = DataPath + "/input.json";
    InputParameters InputPara = ReadInputJson(InputPath);
    PulseParameters PulsePara(InputPara);

    std::vector<std::string> DumpPathes;
    DumpPathes.reserve(InputPara.positions.size());
    for (const auto& posi : InputPara.positions) {
        DumpPathes.push_back(
            DataPath + "/" +
            std::to_string(static_cast<int>(InputPara.E)) + "keV_" +
            std::to_string(posi)
        );
    }

    const int n_abs = InputPara.n_abs;
    const int n_abs_1 = n_abs + 1;
    const int n_abs_2 = n_abs + 2;
    const int n_abs_3 = n_abs + 3;
    const int n_abs_4 = n_abs + 4;
    const int samples = static_cast<int>(InputPara.samples);

    const std::vector<double> Block =
        linspace(-InputPara.length / 20, InputPara.length / 20, n_abs + 1);

    const std::vector<double> time =
        linspace(0, InputPara.samples / InputPara.rate, samples);

    // --------------------------------------------------------
    // 固有値解析: 1回だけ
    // --------------------------------------------------------
    const Eigen::MatrixXd Matrix_A = MakeMatrix_A(PulsePara, InputPara);
    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(Matrix_A);

    Eigen::VectorXcd EigenValues = eigensolver.eigenvalues();
    Eigen::MatrixXcd EigenVectors = eigensolver.eigenvectors();

    for (int i = 0; i < EigenVectors.cols(); ++i) {
        const double norm = EigenVectors.col(i).norm();
        if (norm > 0.0) {
            EigenVectors.col(i) /= norm;
        }
    }

    // 実部行列は毎回作らない
    const Eigen::MatrixXd EiVec = EigenVectors.real();

    // 指数テーブルを前計算
    const Eigen::MatrixXd ExpTable = BuildExpTable(EigenValues, time);

    // 線形方程式ソルバも1回だけ
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(EiVec);

    // --------------------------------------------------------
    // positions ごとのパルス/温度出力
    // --------------------------------------------------------
    for (const auto& posi : InputPara.positions) {
        Eigen::VectorXd VecInit = Eigen::VectorXd::Zero(n_abs_4);
        VecInit[posi + 1] = InputPara.E * 1e3 * PulsePara.e_const / PulsePara.C_abs;

        const Eigen::VectorXd Consts = solver.solve(VecInit);

        const Eigen::VectorXd Pulse0 =
            BuildPulseFromConsts(Consts, EiVec.row(0).transpose(), ExpTable);

        const Eigen::VectorXd Pulse1 =
            BuildPulseFromConsts(Consts, EiVec.row(n_abs_3).transpose(), ExpTable);

        const std::string PulsePath =
            DataPath + "/" +
            std::to_string(static_cast<int>(InputPara.E)) + "keV_" +
            std::to_string(posi) + "/Pulse";

        fs::create_directories(PulsePath + "/Ch0");
        fs::create_directories(PulsePath + "/Ch1");

        const std::string PulseFile_0 = PulsePath + "/Ch0/CH0_0.dat";
        const std::string PulseFile_1 = PulsePath + "/Ch1/CH1_0.dat";

        if (!WriteVectorToFile(PulseFile_0, Pulse0)) {
            return -1;
        }
        if (!WriteVectorToFile(PulseFile_1, Pulse1)) {
            return -1;
        }

        // Temps を [samples x n_abs] 行列でまとめて構築
        Eigen::MatrixXd Temps(samples, n_abs);
        for (int blockIdx = 0; blockIdx < n_abs; ++blockIdx) {
            const int rowIdx = blockIdx + 2;
            Temps.col(blockIdx) =
                BuildPulseFromConsts(Consts, EiVec.row(rowIdx).transpose(), ExpTable);
        }

        const std::string TempFile =
            DataPath + "/" +
            std::to_string(static_cast<int>(InputPara.E)) + "keV_" +
            std::to_string(posi) + "/Temps.csv";

        exportTempsCSV(TempFile, Block, time, Temps);
    }

    // --------------------------------------------------------
    // dump ごとの処理
    // --------------------------------------------------------
    int TotalCounter = 1;

    for (const auto& DumpPath : DumpPathes) {
        std::cout << TotalCounter << "/" << DumpPathes.size() << "\n";

        std::map<int, std::map<int, EventInfo>> batch;
        const std::string jsonPath = DumpPath + "/batch.json";
        const std::string dumpPath = DumpPath + "/dumpall.dat";

        fs::create_directories(DumpPath + "/Pulse_ms/Ch0");
        fs::create_directories(DumpPath + "/Pulse_ms/Ch1");

        if (!LoadOrCreateBatch(jsonPath, dumpPath, InputPara, batch)) {
            return -1;
        }

        // ----------------------------------------------------
        // ブロックごとの基底パルスをまとめて計算
        // ----------------------------------------------------
        // RHS(:, i) = block i への単位入力
        Eigen::MatrixXd RHS = Eigen::MatrixXd::Zero(n_abs_4, n_abs);
        for (int i = 0; i < n_abs; ++i) {
            RHS(i + 2, i) = PulsePara.e_const / PulsePara.C_abs;
        }

        // まとめて解く
        const Eigen::MatrixXd ConstsAll = solver.solve(RHS);

        // 各ブロックの Ch0 / Ch1 パルス基底を [samples x n_abs] で持つ
        const Eigen::MatrixXd Pulse_Blocks_ch0 =
            BuildAllBlockPulses(ConstsAll, EiVec.row(0).transpose(), ExpTable);

        const Eigen::MatrixXd Pulse_Blocks_ch1 =
            BuildAllBlockPulses(ConstsAll, EiVec.row(n_abs_3).transpose(), ExpTable);

        const size_t total_items = batch.size();
        std::atomic<size_t> completed_items{ 0 };
        std::mutex output_mutex;

        concurrency::parallel_for_each(
            batch.begin(),
            batch.end(),
            [&](const std::pair<const int, std::map<int, EventInfo>>& outer_pair) {
                Eigen::VectorXd BlockDeposit = Eigen::VectorXd::Zero(n_abs);

                for (const auto& inner_pair : outer_pair.second) {
                    const auto& event = inner_pair.second;
                    const int depositCount = static_cast<int>(event.x_deposit.size());

                    for (int i = 0; i < depositCount; ++i) {
                        const int EnergyedBlocks = AssignBlock(
                            InputPara,
                            Block,
                            event.x_deposit[static_cast<size_t>(i)],
                            event.y_deposit[static_cast<size_t>(i)],
                            event.z_deposit[static_cast<size_t>(i)]
                        );

                        const int blockIdx = EnergyedBlocks - 1;
                        if (0 <= blockIdx && blockIdx < n_abs) {
                            BlockDeposit[blockIdx] +=
                                event.E_deposit[static_cast<size_t>(i)] * 1e6;
                        }
                    }
                }

                // 元コードの対称加算を保持
                Eigen::VectorXd DepositSym = BlockDeposit;
                for (int i = 0; i < n_abs; ++i) {
                    DepositSym[i] += BlockDeposit[n_abs - 1 - i];
                }

                // [samples x n_abs] * [n_abs] = [samples]
                const Eigen::VectorXd Pulse0 = Pulse_Blocks_ch0 * DepositSym;
                const Eigen::VectorXd Pulse1 = Pulse_Blocks_ch1 * DepositSym;

                const std::string PulseFile_0 =
                    DumpPath + "/Pulse_ms/Ch0/CH0_" + std::to_string(outer_pair.first) + ".dat";
                const std::string PulseFile_1 =
                    DumpPath + "/Pulse_ms/Ch1/CH1_" + std::to_string(outer_pair.first) + ".dat";

                if (!WriteVectorToFile(PulseFile_0, Pulse0)) {
                    std::lock_guard<std::mutex> lock(output_mutex);
                    std::cerr << "Failed to write " << PulseFile_0 << std::endl;
                    return;
                }

                if (!WriteVectorToFile(PulseFile_1, Pulse1)) {
                    std::lock_guard<std::mutex> lock(output_mutex);
                    std::cerr << "Failed to write " << PulseFile_1 << std::endl;
                    return;
                }

                const size_t completed = completed_items.fetch_add(1) + 1;
                const size_t progress =
                    static_cast<size_t>(completed * 100 / (total_items == 0 ? 1 : total_items));

                {
                    std::lock_guard<std::mutex> lock(output_mutex);
                    std::cout << "\rConverting to Pulse: " << progress << "%" << std::flush;
                }
            }
        );

        {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "\nFinished\n";
        }

        ++TotalCounter;
    }

    std::string Input;
    std::cout << "Key input to close...\n";
    std::cin >> Input;
    return 0;
}