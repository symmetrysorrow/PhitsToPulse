#pragma once

#include <filesystem>
#include <map>
#include <vector>
#include "struct.h"
#include <string>

// intと粒子の種類の対応マップ
extern std::map<int, std::string> itype;

// 数字を受け取り対応する粒子を返す関数
std::string GetItype(const double& ityp);

// 空白で文章を分割する関数
std::vector<double> split_line(const std::string& line);

// Input.jsonを読み出す関数
InputParameters ReadInputJson(const std::string& InputPath);

// dumpall.datをbatchにする関数
int ReadDump(const std::string& DumpPath, std::map<int, std::map<int, EventInfo>>& batch, const double& InputEnergy);

// batchをoutput.jsonに書き出す関数
void WriteOutput(const std::map<int, std::map<int, EventInfo>>& batch, const std::string& output_file);