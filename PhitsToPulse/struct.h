#pragma once
#include <nlohmann/json.hpp>
#include<vector>

struct EventInfo {
    int ityp;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> E;
    std::vector<double> x_deposit;
    std::vector<double> y_deposit;
    std::vector<double> z_deposit;
    std::vector<double> E_deposit;
};

inline void from_json(const nlohmann::json& j, EventInfo& e) {
    j.at("ityp").get_to(e.ityp);
    j.at("x").get_to(e.x);
    j.at("y").get_to(e.y);
    j.at("z").get_to(e.z);
    j.at("E").get_to(e.E);
    j.at("x_deposit").get_to(e.x_deposit);
    j.at("y_deposit").get_to(e.y_deposit);
    j.at("z_deposit").get_to(e.z_deposit);
    j.at("E_deposit").get_to(e.E_deposit);
}

struct InputParameters
{
    double C_abs;
    double C_tes;
    double G_abs_abs;
    double G_abs_tes;
    double G_tes_bath;
    double R;
    double R_l;
    double T_c;
    double T_bath;
    double alpha;
    double beta;
    double L;
    double n;
    double E;
    double length;
    double height;
    double depth;
    int n_abs;
    double rate;
    double samples;
    std::vector<int> positions;
    int cutoff;
    int history;
    bool SaveAll;
};