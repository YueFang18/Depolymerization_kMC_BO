// =================================================================
// 性能优化版本：一体化聚合-解聚KMC模拟主程序
// 新增功能：解聚过程中分子量分布的详细追踪和输出
// 主要优化：减少缺陷分析频率、算法优化、内存优化
// 预期性能提升：50-80%
// =================================================================

#include <cstdio>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <utility>
#include <ctime>
#include <map>
#include <string>
#include <sstream>
#include <unordered_map>  // 性能优化：使用hash map
#ifdef _WIN32
#include <direct.h>
    #define mkdir _mkdir
#else
#include <sys/stat.h>
#include <unistd.h>
#endif

// Include all headers
#include "my_rand.h"
#include "DataStructureAndConstant.h"
#include "efficient_explicit_sequence_record_number.h"
#include "depolymerization_explicit_sequence_record_number.h"
#include "ode.h"
#include "MolecularWeightDistribution.h"
#include "MolecularWeightOutput.h"
#include "GlassEffectPropagationRate.h"
#include "PolymerFeathure.h"
#include "InitiatorEfficiency.h"
#include "LengthDependentTerminationRate.h"
#include "ReadPolymer.h"

using namespace std;

// Global mode flags
bool g_is_bo_mode = false;
bool g_is_sweep_mode = false;

// Configuration structure
struct SimulationConfig {
    double poly_target_conversion = 0.50;
    double poly_temp = 343.0;
    double poly_initiator_ratio = 0.003;
    string poly_output_suffix = "poly_";

    double depoly_target_conversion = 0.80;
    double depoly_temp = 363.0;
    string depoly_output_suffix = "depoly_";

    double scaling_factor = 1.0;
    string base_output_dir;
    string input_file;

    bool enable_defect_analysis = false;
    bool enable_mmd_tracking = false;  // 新增：分子量分布追踪开关
    int sweep_index = 0;

    // 新增：解聚过程中需要记录分子量分布的转化率点
    vector<double> mmd_tracking_conversions = {0.0, 0.35, 0.70};
};

// 新增：分子量分布数据结构
struct MolecularWeightDistribution {
    double conversion_rate;
    double time_point;
    vector<int> chain_lengths;
    vector<int> chain_counts;
    double Mn = 0.0;  // 数均分子量
    double Mw = 0.0;  // 重均分子量
    double PDI = 0.0; // 分散度

    void calculateMolecularWeights() {
        if (chain_lengths.empty() || chain_counts.empty()) return;

        double sum_ni = 0.0;      // 总链数
        double sum_ni_Mi = 0.0;   // Σ(ni * Mi)
        double sum_ni_Mi2 = 0.0;  // Σ(ni * Mi^2)

        const double monomer_mw = 100.0;  // MMA单体分子量

        for (size_t i = 0; i < chain_lengths.size(); i++) {
            int ni = chain_counts[i];
            double Mi = chain_lengths[i] * monomer_mw;

            sum_ni += ni;
            sum_ni_Mi += ni * Mi;
            sum_ni_Mi2 += ni * Mi * Mi;
        }

        if (sum_ni > 0) {
            Mn = sum_ni_Mi / sum_ni;
            Mw = sum_ni_Mi2 / sum_ni_Mi;
            PDI = (Mn > 0) ? (Mw / Mn) : 0.0;
        }
    }
};

// 性能优化版本的缺陷分析结构
struct DefectAnalysis {
    double unsaturated_fraction = 0.0;
    double saturated_fraction = 0.0;
    double hh_bond_fraction = 0.0;

    double avg_chain_length = 0.0;
    double chain_length_std = 0.0;
    double max_chain_length = 0.0;

    // 性能优化：使用unordered_map
    unordered_map<int, int> defects_by_length;

    double final_conversion_rate = 0.0;
    double time_to_target = 0.0;
    double depoly_time_to_target = 0.0;
    double depoly_conversion_rate = 0.0;
    double recyclability_index = 0.0;

    void clear() {
        defects_by_length.clear();
        defects_by_length.reserve(100);
    }
};

// Data transition structure
struct TransitionData {
    double final_conversion;
    double final_time;
    double final_volume;
    double M0_total;
    double final_temperature;

    vector<int> all_chain1;
    vector<int> all_chain2;
    vector<int> all_num_M1;
    vector<int> all_chain_feature;
    vector<int> all_chain1_feature;
    vector<int> all_chain2_feature;

    DefectAnalysis defect_stats;

    // 新增：聚合结束时的分子量分布
    MolecularWeightDistribution initial_mmd;

    void clear() {
        all_chain1.clear();
        all_chain2.clear();
        all_num_M1.clear();
        all_chain_feature.clear();
        all_chain1_feature.clear();
        all_chain2_feature.clear();
        defect_stats.clear();
    }
};

// 新增：解聚过程分子量分布追踪结构
struct DepolymerizationMMDTracker {
    vector<MolecularWeightDistribution> mmd_snapshots;
    vector<double> target_conversions;
    size_t current_target_index = 0;

    void initialize(const vector<double>& conversions) {
        target_conversions = conversions;
        mmd_snapshots.clear();
        mmd_snapshots.reserve(conversions.size());
        current_target_index = 0;
    }

    bool shouldRecord(double current_conversion) {
        if (current_target_index >= target_conversions.size()) return false;
        return current_conversion >= target_conversions[current_target_index];
    }

    void recordMMD(double conversion, double time) {
        if (current_target_index >= target_conversions.size()) return;

        MolecularWeightDistribution mmd;
        mmd.conversion_rate = conversion;
        mmd.time_point = time;

        // 统计当前所有聚合物的链长分布
        unordered_map<int, int> length_count;

        // 统计未饱和聚合物
        for (int length : unsaturated_Polymers.num_M1) {
            length_count[length]++;
        }

        // 统计饱和聚合物
        for (int length : saturated_Polymers.num_M1) {
            length_count[length]++;
        }

        // 统计HH键聚合物
        for (int length : HH_Polymers.num_M1) {
            length_count[length]++;
        }

        // 转换为向量格式
        for (const auto& pair : length_count) {
            mmd.chain_lengths.push_back(pair.first);
            mmd.chain_counts.push_back(pair.second);
        }

        // 计算分子量
        mmd.calculateMolecularWeights();

        mmd_snapshots.push_back(mmd);
        current_target_index++;
    }
};

// Utility functions
bool createDirectory(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
#ifdef _WIN32
        return (_mkdir(path.c_str()) == 0);
#else
        return (mkdir(path.c_str(), 0755) == 0);
#endif
    } else if (S_ISDIR(info.st_mode)) {
        return true;
    } else {
        return false;
    }
}

string ensureTrailingSlash(const string& path) {
    if (path.empty()) return "/";
    if (path.back() != '/' && path.back() != '\\') {
        return path + "/";
    }
    return path;
}

template <typename T>
void shrinkVector(std::vector<T>& vec, double scaling_factor) {
    size_t newSize = static_cast<size_t>(std::ceil(vec.size() * scaling_factor));
    if (newSize < vec.size()) {
        vec.erase(vec.begin() + newSize, vec.end());
    }
}

void clearGlobalState() {
    poly = poly_chain();
    unsaturated_Polymers = poly_chain();
    saturated_Polymers = poly_chain();
    HH_Polymers = poly_chain();
    Pr.clear();
    for (int i = 0; i < 12; i++) {
        species[i] = 0;
    }
    for (int i = 0; i < 2; i++) {
        num_monomer[i] = 0;
        M_initial[i] = 0;
    }
    M1M1.clear();
    M1M2.clear();
    M2M1.clear();
    M2M2.clear();
}

// 新增：计算聚合物链长分布的函数
MolecularWeightDistribution calculateMMDFromPolymer(const poly_chain& polymer_data, double conversion, double time) {
    MolecularWeightDistribution mmd;
    mmd.conversion_rate = conversion;
    mmd.time_point = time;

    if (polymer_data.num_M1.empty()) {
        return mmd;
    }

    // 统计链长分布
    unordered_map<int, int> length_count;
    for (int length : polymer_data.num_M1) {
        length_count[length]++;
    }

    // 转换为向量格式
    for (const auto& pair : length_count) {
        mmd.chain_lengths.push_back(pair.first);
        mmd.chain_counts.push_back(pair.second);
    }

    // 计算分子量
    mmd.calculateMolecularWeights();

    return mmd;
}

// 新增：输出分子量分布到文件
void exportMMDToFile(const string& filename, const MolecularWeightDistribution& mmd, const string& description) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Warning: Could not create MMD file: " << filename << endl;
        return;
    }

    file << "# Molecular Weight Distribution - " << description << "\n";
    file << "# Conversion: " << mmd.conversion_rate << "\n";
    file << "# Time: " << mmd.time_point << " s\n";
    file << "# Mn: " << mmd.Mn << " g/mol\n";
    file << "# Mw: " << mmd.Mw << " g/mol\n";
    file << "# PDI: " << mmd.PDI << "\n";
    file << "# Format: chain_length count weight_fraction log10_molecular_weight\n";
    file << "# ========================================\n";

    const double monomer_mw = 100.0;
    double total_weight = 0.0;

    // 先计算总重量
    for (size_t i = 0; i < mmd.chain_lengths.size(); i++) {
        total_weight += mmd.chain_counts[i] * mmd.chain_lengths[i] * monomer_mw;
    }

    // 输出数据
    for (size_t i = 0; i < mmd.chain_lengths.size(); i++) {
        int length = mmd.chain_lengths[i];
        int count = mmd.chain_counts[i];
        double molecular_weight = length * monomer_mw;
        double weight_fraction = (total_weight > 0) ? (count * molecular_weight) / total_weight : 0.0;
        double log_mw = (molecular_weight > 0) ? log10(molecular_weight) : 0.0;

        file << length << " " << count << " " << weight_fraction << " " << log_mw << "\n";
    }

    file.close();
}

// 新增：输出所有解聚过程的分子量分布
void exportDepolymerizationMMD(const string& output_dir, const DepolymerizationMMDTracker& tracker,
                               const MolecularWeightDistribution& initial_mmd) {
    string mmd_dir = ensureTrailingSlash(output_dir) + "molecular_weight_distributions/";
    createDirectory(mmd_dir);

    // 输出初始分子量分布（聚合后）
    exportMMDToFile(mmd_dir + "mmd_initial_after_polymerization.dat", initial_mmd, "After Polymerization");

    // 输出解聚过程中的分子量分布
    for (size_t i = 0; i < tracker.mmd_snapshots.size(); i++) {
        const auto& mmd = tracker.mmd_snapshots[i];
        stringstream ss;
        ss << "mmd_depolymerization_" << fixed << setprecision(0) << (mmd.conversion_rate * 100) << "percent.dat";
        string filename = mmd_dir + ss.str();

        stringstream desc;
        desc << "Depolymerization at " << (mmd.conversion_rate * 100) << "% conversion";
        exportMMDToFile(filename, mmd, desc.str());
    }

    // 输出汇总信息
    ofstream summary(mmd_dir + "mmd_summary.dat");
    if (summary.is_open()) {
        summary << "# Molecular Weight Distribution Summary\n";
        summary << "# Format: stage conversion time Mn Mw PDI\n";
        summary << "# =====================================\n";

        summary << "Initial " << initial_mmd.conversion_rate << " " << initial_mmd.time_point << " "
                << initial_mmd.Mn << " " << initial_mmd.Mw << " " << initial_mmd.PDI << "\n";

        for (const auto& mmd : tracker.mmd_snapshots) {
            summary << "Depolymerization " << mmd.conversion_rate << " " << mmd.time_point << " "
                    << mmd.Mn << " " << mmd.Mw << " " << mmd.PDI << "\n";
        }

        summary.close();
    }
}

// 性能优化版本的缺陷分析函数
DefectAnalysis analyzePolymerDefectsOptimized(const poly_chain& polymer_data) {
    DefectAnalysis analysis;

    if (polymer_data.num_M1.empty()) {
        return analysis;
    }

    const size_t total_chains = polymer_data.num_M1.size();
    int unsaturated_count = 0, saturated_count = 0, hh_bond_count = 0;

    // 性能优化：预分配内存
    analysis.defects_by_length.reserve(total_chains / 10);

    double total_length = 0.0;
    double max_length = 0.0;
    double sum_of_squares = 0.0;  // 避免二次遍历

    // 关键优化：单次遍历完成所有计算
    for (size_t i = 0; i < total_chains; i++) {
        const int feature = polymer_data.chain_feature[i];
        const int length = polymer_data.num_M1[i];

        total_length += length;
        max_length = max(max_length, (double)length);

        // 性能优化：避免pow()函数
        const double length_d = static_cast<double>(length);
        sum_of_squares += length_d * length_d;

        // 分类缺陷
        if (feature == 0 || feature == 00000) {
            unsaturated_count++;
        } else if (feature == 1 || feature == 11111) {
            saturated_count++;
        } else if (feature == 2) {
            hh_bond_count++;
        }

        // 使用unordered_map提升性能
        ++analysis.defects_by_length[length];
    }

    // 计算统计量
    const double total_chains_d = static_cast<double>(total_chains);
    analysis.unsaturated_fraction = unsaturated_count / total_chains_d;
    analysis.saturated_fraction = saturated_count / total_chains_d;
    analysis.hh_bond_fraction = hh_bond_count / total_chains_d;

    analysis.avg_chain_length = total_length / total_chains_d;
    analysis.max_chain_length = max_length;

    // 性能优化：单次计算标准差，避免sqrt(pow())
    const double mean = analysis.avg_chain_length;
    const double variance = (sum_of_squares / total_chains_d) - (mean * mean);
    analysis.chain_length_std = (variance > 0) ? sqrt(variance) : 0.0;

    return analysis;
}

// 快速缺陷分析（仅用于实时跟踪）
void quickDefectAnalysis(const poly_chain& polymer_data, double& unsaturated_frac,
                         double& saturated_frac, double& hh_frac, double& avg_length) {
    if (polymer_data.num_M1.empty()) {
        unsaturated_frac = saturated_frac = hh_frac = avg_length = 0.0;
        return;
    }

    int unsaturated = 0, saturated = 0, hh_bond = 0;
    double total_length = 0.0;
    const size_t total_chains = polymer_data.num_M1.size();

    for (size_t i = 0; i < total_chains; i++) {
        const int feature = polymer_data.chain_feature[i];
        total_length += polymer_data.num_M1[i];

        if (feature == 0 || feature == 00000) unsaturated++;
        else if (feature == 1 || feature == 11111) saturated++;
        else if (feature == 2) hh_bond++;
    }

    const double total_chains_d = static_cast<double>(total_chains);
    unsaturated_frac = unsaturated / total_chains_d;
    saturated_frac = saturated / total_chains_d;
    hh_frac = hh_bond / total_chains_d;
    avg_length = total_length / total_chains_d;
}

// Export detailed defect analysis
void exportDefectAnalysis(const string& output_dir, const DefectAnalysis& analysis,
                          double temp, double ratio) {
    string defect_file = ensureTrailingSlash(output_dir) + "defect_analysis.out";
    ofstream defect_out(defect_file);

    if (!defect_out.is_open()) {
        cerr << "Warning: Could not create defect analysis file" << endl;
        return;
    }

    defect_out << "# Defect Analysis Report\n";
    defect_out << "# Polymerization Conditions: T=" << temp << "K, Ratio=" << ratio << "\n";
    defect_out << "# ================================================\n\n";

    defect_out << "DEFECT_FRACTIONS\n";
    defect_out << "unsaturated_fraction " << analysis.unsaturated_fraction << "\n";
    defect_out << "saturated_fraction " << analysis.saturated_fraction << "\n";
    defect_out << "hh_bond_fraction " << analysis.hh_bond_fraction << "\n\n";

    defect_out << "CHAIN_STATISTICS\n";
    defect_out << "avg_chain_length " << analysis.avg_chain_length << "\n";
    defect_out << "chain_length_std " << analysis.chain_length_std << "\n";
    defect_out << "max_chain_length " << analysis.max_chain_length << "\n\n";

    defect_out << "KINETIC_PROPERTIES\n";
    defect_out << "time_to_target " << analysis.time_to_target << "\n";
    defect_out << "depoly_time_to_target " << analysis.depoly_time_to_target << "\n";
    defect_out << "recyclability_index " << analysis.recyclability_index << "\n\n";

    defect_out << "DEFECTS_BY_LENGTH\n";
    defect_out << "chain_length count\n";
    for (const auto& pair : analysis.defects_by_length) {
        defect_out << pair.first << " " << pair.second << "\n";
    }

    defect_out.close();
}

// 性能优化版本的聚合阶段
TransitionData runPolymerizationPhase(const SimulationConfig& config) {
    if (!g_is_bo_mode && !g_is_sweep_mode) {
        cout << "\n=== Starting Polymerization Phase ===" << endl;
        cout << "  Temperature: " << config.poly_temp << " K, Initiator Ratio: " << config.poly_initiator_ratio << endl;
    }

    clock_t start = clock();
    clearGlobalState();

    // 性能优化：智能分析模式控制
    const bool detailed_analysis = config.enable_defect_analysis && (!g_is_bo_mode && !g_is_sweep_mode);
    const bool quick_analysis = config.enable_defect_analysis && g_is_sweep_mode;
    const bool mmd_tracking = config.enable_mmd_tracking && (!g_is_bo_mode && !g_is_sweep_mode);

    // 变量声明（保持与原程序一致）
    double mol_den[2], mol_wt[2], mol_frac[2];
    double temp, init_conc, solvent_conc, time_limit, M0;
    double total_rate = 0.0, reaction = 0.0, time = 0.0, out_time = 0.0;
    double d_time = 0.0, Aver_d_time = 0.0, Aver_d_time_1 = 0.0;
    double rate[7] = {0.0}, A[7] = {0.0}, E[7] = {0.0}, c[7] = {0.0};
    double conversion[2] = {0.0};
    double volume = 0.0, scale_fac = 1.0;
    double radical_ode = 0.0, Aver_add_spe = 0.0;
    double all_species_counter = 0, species_counter = 0;
    double all_species = 0.0;
    int num_chain = -1, reaction_index;
    double rand_11, rand_22, tau = 0.0, reaction_counter = 0.0;
    int fl_aver = 0, flag1 = 0;

    stiff_system stiff_sys;
    stiff_system_jacobi system_jacobi;

    // 读取输入参数
    ifstream input0(config.input_file);
    if (!input0.is_open()) {
        throw runtime_error("无法打开输入文件: " + config.input_file);
    }

    double temp_from_file;
    input0 >> scale_fac >> temp_from_file;
    temp = config.poly_temp;

    for (int i = 0; i < 2; i++) input0 >> mol_den[i];
    for (int i = 0; i < 2; i++) input0 >> mol_wt[i];
    for (int i = 0; i < 2; i++) input0 >> mol_frac[i];
    input0 >> init_conc >> solvent_conc >> time_limit >> M0;

    for (int i = 0; i < 7; i++) {
        input0 >> A[i] >> E[i];
        c[i] = A[i] * exp(-E[i] / (R * temp));
    }
    input0.close();

    double ratio_AIBN_to_MMA = config.poly_initiator_ratio;
    species[2] = round(M0 / (1 + ratio_AIBN_to_MMA) + 0.5);
    species[0] = M0 - species[2];
    volume = (species[2] / Na) * (mol_wt[0] / mol_den[0]);
    M_initial[0] = species[2];

    // 设置输出文件（性能优化：减少不必要的文件I/O）
    ofstream conv, molwt, poly_feature, defect_tracker;
    const bool create_standard_output = !g_is_bo_mode && !g_is_sweep_mode;

    if (create_standard_output || detailed_analysis) {
        string base_dir = ensureTrailingSlash(config.base_output_dir);
        createDirectory(base_dir);

        string poly_dir = base_dir + "poly_";
        createDirectory(poly_dir);
        poly_dir = ensureTrailingSlash(poly_dir);

        if (create_standard_output) {
            conv.open(poly_dir + "conversion.out");
            molwt.open(poly_dir + "molecularweight.out");
            poly_feature.open(poly_dir + "polymer_feature.out");

            conv << "time\tconversion\toverall_conversion" << endl;
            molwt << "conversion\tMn\tMw" << endl;
            poly_feature << "conversion num_saturated num_unsaturated num_HHbond" << endl;
        }

        if (detailed_analysis) {
            defect_tracker.open(poly_dir + "defect_evolution.out");
            defect_tracker << "conversion time unsaturated_frac saturated_frac hh_frac avg_length" << endl;
        }
    }

    // 初始化反应速率常数
    double k_pchem = c[2];
    double initiator_eff = 0.58;
    double kd_chem = c[0];

    c[0] = kd_chem * initiator_eff;
    c[3] = c[4] * (1 - DisproportionationCoefficient) / DisproportionationCoefficient;

    double kt_re = c[3];
    double kt_dis = c[4];
    double kt_i = c[0];
    double kt_p = c[2];

    vector_type r(6);
    r[0] = species[0] / (volume * Na);
    r[2] = species[2] / (volume * Na);
    vector_type Rad(6);
    Rad[0] = species[0] / (volume * Na);
    Rad[2] = species[2] / (volume * Na);

    r[1] = 0.0; r[3] = 0.0; r[4] = 0.0; r[5] = 0.0;
    Rad[1] = 0.0; Rad[3] = 0.0; Rad[4] = 0.0; Rad[5] = 0.0;

    // 初始化反应速率
    rate[0] = c[0] * species[0];
    rate[1] = c[1] / (Na * volume) * species[1] * species[2];
    rate[2] = c[2] / scale_fac / (Na * volume) * species[2] * species[3];
    rate[3] = c[3] / scale_fac / scale_fac / (Na * volume) * species[3] * (species[3] - 1) / 2;
    rate[4] = c[4] / scale_fac / scale_fac / (Na * volume) * species[3] * (species[3] - 1) / 2;
    rate[5] = c[5] / scale_fac / (Na * volume) * species[2] * species[3];
    rate[6] = c[6] / (Na * volume) * species[2] * species[5];

    // 关键性能优化：大幅减少缺陷分析频率
    double next_conversion = 0;
    double conversionsteps = 1e-3;
    // 从每5%改为每20%，减少4倍分析次数
    double defect_track_interval = detailed_analysis ? 0.20 : 0.50;
    double next_defect_track = defect_track_interval;

    // 聚合KMC主循环
    while (conversion[0] <= config.poly_target_conversion) {
        if (conversion[0] > 0.2) {
            conversionsteps = 1e-4;
        }

        // KMC步骤（保持原始逻辑）
        rand_11 = 1.0 * (my_rand(RAND_MAX - 1) + 1.0) / RAND_MAX;
        rand_22 = 1.0 * (my_rand(RAND_MAX - 1) + 1.0) / RAND_MAX;

        total_rate = 0.0;
        for (int i = 0; i < 7; i++) {
            total_rate += rate[i];
        }
        if (total_rate == 0.0) break;

        tau = (1 / total_rate) * log(1 / rand_11);
        reaction = rand_22 * total_rate;

        reaction_counter = rate[0];
        reaction_index = 0;
        while (reaction_counter < reaction) {
            reaction_index++;
            reaction_counter += rate[reaction_index];
        }

        time += tau;
        out_time += tau;
        d_time += tau;

        all_species_counter += species[3] * tau;
        species_counter++;
        all_species = species[3];
        Aver_add_spe += all_species * tau;
        Aver_d_time += tau;
        Aver_d_time_1 += tau;

        // 计算缩放因子（保持原始逻辑）
        if (scale_fac == 1 && d_time > 0.001 && flag1 == 0) {
            ode(stiff_sys, system_jacobi, r, Rad, d_time, kt_re, kt_dis, kt_i, kt_p);
            radical_ode = r[3];
            scale_fac = 2 / (Na * volume * radical_ode);
            flag1 = 1;
            d_time = 0;
        }

        // 动态缩放因子更新（保持原始逻辑）
        if (((conversion[0] <= 0.2 && Aver_d_time_1 > 1.0) ||
             (conversion[0] > 0.2 && Aver_d_time_1 > 0.1)) &&
            Aver_add_spe / Aver_d_time > 2 && fl_aver == 0) {
            fl_aver = 1;
            double aver_spe3 = Aver_add_spe / Aver_d_time;

            ode(stiff_sys, system_jacobi, r, Rad, Aver_d_time_1, kt_re, kt_dis, kt_i, kt_p);
            radical_ode = r[3];
            scale_fac = aver_spe3 / (Na * volume * radical_ode);

            Aver_add_spe = 0;
            Aver_d_time = 0;
            Aver_d_time_1 = 0;
            fl_aver = 0;
        }

        // 定期输出（仅在非BO模式下）
        if (create_standard_output && out_time >= 10) {
            double num_unsaturated, num_saturated, num_HHbond;
            countPolymerFeatures(poly, num_unsaturated, num_saturated, num_HHbond);
            poly_feature << conversion[0] << " " << num_saturated << " "
                         << num_unsaturated << " " << num_HHbond << endl;

            MolecularWeightOutput(molwt, Pr, poly.num_M1, conversion[0]);

            double over_num = num_monomer[0];
            double over_den = M_initial[0];
            double over_x = over_num / over_den;
            conv << time << "\t" << conversion[0] << "\t" << over_x << endl;

            out_time = 0.0;
        }

        // 关键性能优化：大幅减少缺陷跟踪频率
        if (detailed_analysis && conversion[0] >= next_defect_track) {
            DefectAnalysis current_defects = analyzePolymerDefectsOptimized(poly);
            if (defect_tracker.is_open()) {
                defect_tracker << conversion[0] << " " << time << " "
                               << current_defects.unsaturated_fraction << " "
                               << current_defects.saturated_fraction << " "
                               << current_defects.hh_bond_fraction << " "
                               << current_defects.avg_chain_length << endl;
            }
            next_defect_track += defect_track_interval;
        }

        // 调用反应处理函数（保持原始逻辑）
        efficient_explicit_sequence_record_number(reaction_index, num_monomer, species,
                                                  num_chain, M1M1, M1M2, M2M1, M2M2,
                                                  poly, conversion[0]);

        // 更新转化率和反应常数（保持原始逻辑）
        conversion[0] = num_monomer[0] / M_initial[0];
        if (conversion[0] > next_conversion) {
            initiator_eff = CalculateInitiatorEfficiency(temp, conversion[0]);
            c[0] = kd_chem * initiator_eff;
            c[4] = DisproportionationCoefficient * PopulationWeightedTerminationRate(Pr, conversion[0]);
            c[3] = c[4] * (1 - DisproportionationCoefficient) / DisproportionationCoefficient;
            if (conversion[0] > 0.80) {
                c[2] = CalculatePropagationRate(conversion[0], temp, k_pchem);
            }
            next_conversion += conversionsteps;
            kt_re = c[3];
            kt_dis = c[4];
            kt_i = c[0];
            kt_p = c[2];
        }

        // 更新反应速率（保持原始逻辑）
        rate[0] = c[0] * species[0];
        rate[1] = c[1] / (Na * volume) * species[1] * species[2];
        rate[2] = c[2] / scale_fac / (Na * volume) * species[2] * species[3];
        rate[3] = c[3] / scale_fac / scale_fac / (Na * volume) * species[3] * (species[3] - 1) / 2;
        rate[4] = c[4] / scale_fac / scale_fac / (Na * volume) * species[3] * (species[3] - 1) / 2;
        rate[5] = c[5] / scale_fac / (Na * volume) * species[2] * species[3];
        rate[6] = c[6] / (Na * volume) * species[2] * species[5];
    }

    // 关闭文件
    if (create_standard_output) {
        conv.close();
        molwt.close();
        poly_feature.close();
    }
    if (defect_tracker.is_open()) {
        defect_tracker.close();
    }

    if (!g_is_bo_mode && !g_is_sweep_mode) {
        clock_t end = clock();
        cout << "聚合阶段完成!" << endl;
        cout << "用时: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
        cout << "最终转化率: " << conversion[0] << endl;
        cout << "总聚合物链数: " << poly.num_M1.size() << endl;
    }

    // 准备传递数据
    TransitionData transition_data;
    transition_data.final_conversion = conversion[0];
    transition_data.final_time = time;
    transition_data.final_volume = volume;
    transition_data.M0_total = M0;
    transition_data.final_temperature = temp;

    // 完整复制聚合物数据
    transition_data.all_chain1 = poly.chain1;
    transition_data.all_chain2 = poly.chain2;
    transition_data.all_num_M1 = poly.num_M1;
    transition_data.all_chain_feature = poly.chain_feature;
    transition_data.all_chain1_feature = poly.chain1_feature;
    transition_data.all_chain2_feature = poly.chain2_feature;

    // 新增：保存聚合结束时的分子量分布
    if (mmd_tracking) {
        transition_data.initial_mmd = calculateMMDFromPolymer(poly, conversion[0], time);
    }

    // 性能优化：仅在需要时进行缺陷分析
    if (detailed_analysis || quick_analysis) {
        transition_data.defect_stats = analyzePolymerDefectsOptimized(poly);
        transition_data.defect_stats.time_to_target = time;
    } else if (g_is_sweep_mode) {
        // 快速分析版本
        quickDefectAnalysis(poly,
                            transition_data.defect_stats.unsaturated_fraction,
                            transition_data.defect_stats.saturated_fraction,
                            transition_data.defect_stats.hh_bond_fraction,
                            transition_data.defect_stats.avg_chain_length);
        transition_data.defect_stats.time_to_target = time;
    }

    return transition_data;
}

// 解聚阶段实现（新增分子量分布追踪）
double runDepolymerizationPhase(const TransitionData& transition_data, const SimulationConfig& config) {
    if (!g_is_bo_mode && !g_is_sweep_mode) cout << "\n=== 开始解聚阶段 ===" << endl;

    clock_t start = clock();
    clearGlobalState();

    // 新增：初始化分子量分布追踪器
    DepolymerizationMMDTracker mmd_tracker;
    if (config.enable_mmd_tracking && (!g_is_bo_mode && !g_is_sweep_mode)) {
        mmd_tracker.initialize(config.mmd_tracking_conversions);
        cout << "分子量分布追踪已启用，将在以下转化率记录: ";
        for (double conv : config.mmd_tracking_conversions) {
            cout << conv * 100 << "% ";
        }
        cout << endl;
    }

    // 缩放参数
    double scaling_times = config.scaling_factor;
    double M0 = transition_data.M0_total * scaling_times;

    // 重新初始化聚合物数据结构
    vector<int> temp_chain1 = transition_data.all_chain1;
    vector<int> temp_chain2 = transition_data.all_chain2;
    vector<int> temp_num_M1 = transition_data.all_num_M1;
    vector<int> temp_chain_feature = transition_data.all_chain_feature;
    vector<int> temp_chain1_feature = transition_data.all_chain1_feature;
    vector<int> temp_chain2_feature = transition_data.all_chain2_feature;

    // 应用缩放
    if (scaling_times != 1.0) {
        shrinkVector(temp_chain1, scaling_times);
        shrinkVector(temp_chain2, scaling_times);
        shrinkVector(temp_num_M1, scaling_times);
        shrinkVector(temp_chain_feature, scaling_times);
        shrinkVector(temp_chain1_feature, scaling_times);
        shrinkVector(temp_chain2_feature, scaling_times);
    }

    // 按照特征分类聚合物
    for (size_t i = 0; i < temp_num_M1.size(); i++) {
        int feature = temp_chain_feature[i];

        if (feature == 0 || feature == 00000) {
            unsaturated_Polymers.chain1.push_back(temp_chain1[i]);
            unsaturated_Polymers.chain2.push_back(temp_chain2[i]);
            unsaturated_Polymers.num_M1.push_back(temp_num_M1[i]);
            unsaturated_Polymers.chain_feature.push_back(temp_chain_feature[i]);
            unsaturated_Polymers.chain1_feature.push_back(temp_chain1_feature[i]);
            unsaturated_Polymers.chain2_feature.push_back(temp_chain2_feature[i]);
        } else if (feature == 1 || feature == 11111) {
            saturated_Polymers.chain1.push_back(temp_chain1[i]);
            saturated_Polymers.chain2.push_back(temp_chain2[i]);
            saturated_Polymers.num_M1.push_back(temp_num_M1[i]);
            saturated_Polymers.chain_feature.push_back(temp_chain_feature[i]);
            saturated_Polymers.chain1_feature.push_back(temp_chain1_feature[i]);
            saturated_Polymers.chain2_feature.push_back(temp_chain2_feature[i]);
        } else if (feature == 2) {
            HH_Polymers.chain1.push_back(temp_chain1[i]);
            HH_Polymers.chain2.push_back(temp_chain2[i]);
            HH_Polymers.num_M1.push_back(temp_num_M1[i]);
            HH_Polymers.chain_feature.push_back(temp_chain_feature[i]);
            HH_Polymers.chain1_feature.push_back(temp_chain1_feature[i]);
            HH_Polymers.chain2_feature.push_back(temp_chain2_feature[i]);
        }
    }

    // 正确初始化species数组
    species[0] = static_cast<double>(unsaturated_Polymers.num_M1.size());
    species[1] = static_cast<double>(saturated_Polymers.num_M1.size());
    species[2] = static_cast<double>(HH_Polymers.num_M1.size());
    species[3] = 0;
    species[4] = 0;
    species[5] = 0;

    // 计算新的限制转化率
    int new_M0 = 0;
    for (int length : temp_num_M1) {
        new_M0 += length;
    }
    if (new_M0 == 0) return 1e12;

    double new_limiting_conversion = static_cast<double>(new_M0) / M0;
    M_initial[0] = M0 * new_limiting_conversion;
    if (M_initial[0] == 0) return 1e12;

    // 解聚变量声明
    double time = 0.0;
    double conversion[2] = {0.0};
    double total_rate = 0.0, reaction = 0.0, tau = 0.0;
    double rate[7] = {0.0}, c[7] = {0.0};
    int reaction_index, num_chain = -1;
    double rand_11, rand_22, reaction_counter = 0.0;

    // 设置解聚反应常数
    c[0] = 8.6e-4;
    c[1] = 2e-7;
    c[2] = 1.5e-2;
    c[3] = 5e3;
    c[4] = 1e4;
    c[5] = 1e4;

    double density = 940;
    double monomer_weight = 100;
    double vol = 0;

    // 设置输出文件
    ofstream deconv, depoly_feature;
    if (!g_is_bo_mode && !g_is_sweep_mode) {
        string base_dir = ensureTrailingSlash(config.base_output_dir);
        createDirectory(base_dir);

        string depoly_dir = base_dir + "depoly_";
        createDirectory(depoly_dir);
        depoly_dir = ensureTrailingSlash(depoly_dir);

        deconv.open(depoly_dir + "conversion.out");
        depoly_feature.open(depoly_dir + "polymer_feature.out");

        deconv << "time\tconversion" << endl;
        depoly_feature << "monomer_yield num_unsaturated num_saturated num_HHbond time" << endl;
    }

    // 新增：记录初始分子量分布（0%转化率）
    if (config.enable_mmd_tracking && (!g_is_bo_mode && !g_is_sweep_mode)) {
        if (mmd_tracker.shouldRecord(0.0)) {
            mmd_tracker.recordMMD(0.0, 0.0);
        }
    }

    // 解聚KMC主循环
    double next_conversion = 0;
    double conversionsteps = 0.01;
    double scale_fac = 1.0;
    double tot_dyad[4] = {0}, tot_triad[8] = {0};

    while (conversion[0] < config.depoly_target_conversion) {
        vol = (M_initial[0] * (1 - conversion[0]) * monomer_weight) / (Na * density);
        if (vol <= 0) break;

        rand_11 = (double) rand() / RAND_MAX;
        rand_22 = (double) rand() / RAND_MAX;

        int total_bonds_saturated = 0;
        for (size_t i = 0; i < saturated_Polymers.num_M1.size(); ++i) {
            int chain_length = saturated_Polymers.num_M1[i];
            if (chain_length > 1) {
                total_bonds_saturated += (chain_length - 1);
            }
        }

        rate[0] = c[0] * species[0];
        rate[1] = c[1] * total_bonds_saturated;
        rate[2] = c[2] * species[2];
        rate[3] = c[3] * (species[3] - species[5]) / scale_fac;
        rate[4] = c[4] * 2 * species[3] * (species[3] - 1) / (Na * vol) / 2 / scale_fac / scale_fac;
        rate[5] = c[5] * 2 * species[3] * (species[3] - 1) / (Na * vol) / 2 / scale_fac / scale_fac;

        total_rate = 0.0;
        for (int i = 0; i < 6; i++) {
            total_rate += rate[i];
        }
        if (total_rate <= 1e-20) break;

        tau = (1 / total_rate) * log(1 / rand_11);
        reaction = rand_22 * total_rate;

        reaction_counter = rate[0];
        reaction_index = 0;
        while (reaction_counter < reaction) {
            reaction_index++;
            reaction_counter += rate[reaction_index];
        }

        time += tau;

        depolymerization_explicit_sequence_record_number(reaction_index, num_monomer, num_chain,
                                                         M1M1, M1M2, M2M1, M2M2, poly,
                                                         tot_dyad, tot_triad, conversion[0]);

        conversion[0] = num_monomer[0] / M_initial[0];

        // 新增：检查是否需要记录分子量分布
        if (config.enable_mmd_tracking && (!g_is_bo_mode && !g_is_sweep_mode)) {
            if (mmd_tracker.shouldRecord(conversion[0])) {
                mmd_tracker.recordMMD(conversion[0], time);
                cout << "已记录 " << (conversion[0] * 100) << "% 转化率时的分子量分布" << endl;
            }
        }

        if (!g_is_bo_mode && !g_is_sweep_mode && conversion[0] > next_conversion) {
            deconv << time << "\t" << conversion[0] << endl;
            depoly_feature << conversion[0] << " " << unsaturated_Polymers.num_M1.size() << " "
                           << saturated_Polymers.num_M1.size() << " " << HH_Polymers.num_M1.size()
                           << " " << time << endl;
            next_conversion += conversionsteps;
        }
    }

    // 新增：输出分子量分布数据
    if (config.enable_mmd_tracking && (!g_is_bo_mode && !g_is_sweep_mode)) {
        cout << "正在输出分子量分布数据..." << endl;
        exportDepolymerizationMMD(config.base_output_dir, mmd_tracker, transition_data.initial_mmd);
        cout << "分子量分布数据已保存到 molecular_weight_distributions/ 文件夹" << endl;
    }

    if (!g_is_bo_mode && !g_is_sweep_mode) {
        deconv.close();
        depoly_feature.close();

        clock_t end = clock();
        cout << "解聚阶段完成!" << endl;
        cout << "用时: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
        cout << "最终单体产率: " << conversion[0] << endl;
    }

    return time;
}

int main(int argc, char* argv[]) {
    // 默认配置
    SimulationConfig config;
    config.input_file = "src/input0.txt";
    config.base_output_dir = "/Users/yueyue/CLionProjects/untitled3/output/";
    config.poly_target_conversion = 0.50;
    config.depoly_target_conversion = 0.80;
    config.scaling_factor = 1.0;
    config.poly_temp = 343.0;
    config.poly_initiator_ratio = 0.003;
    config.enable_mmd_tracking = false;  // 默认关闭分子量分布追踪

    // 命令行参数解析
    if (argc > 1) {
        string arg1 = argv[1];
        if (arg1 == "--bo" && argc == 4) {
            g_is_bo_mode = true;
            config.poly_temp = stod(argv[2]);
            config.poly_initiator_ratio = stod(argv[3]);
            config.enable_defect_analysis = true;
        } else if (arg1 == "--analyze" && argc == 5) {
            g_is_bo_mode = false;
            config.poly_temp = stod(argv[2]);
            config.poly_initiator_ratio = stod(argv[3]);
            config.base_output_dir = argv[4];
            config.enable_defect_analysis = true;
            config.enable_mmd_tracking = true;  // 分析模式下启用分子量分布追踪
        } else if (arg1 == "--sweep" && argc == 6) {
            g_is_sweep_mode = true;
            config.poly_temp = stod(argv[2]);
            config.poly_initiator_ratio = stod(argv[3]);
            config.base_output_dir = argv[4];
            config.sweep_index = stoi(argv[5]);
            config.enable_defect_analysis = true;
        } else if (arg1 == "--mmd" && argc == 5) {
            // 新增：专门用于分子量分布分析的模式
            g_is_bo_mode = false;
            config.poly_temp = stod(argv[2]);
            config.poly_initiator_ratio = stod(argv[3]);
            config.base_output_dir = argv[4];
            config.enable_defect_analysis = true;
            config.enable_mmd_tracking = true;
            cout << "分子量分布分析模式已启用" << endl;
        } else {
            cerr << "错误：无效的参数！" << endl;
            cerr << "用法 1 (单次默认运行): " << argv[0] << endl;
            cerr << "用法 2 (优化器调用): " << argv[0] << " --bo <temp> <ratio>" << endl;
            cerr << "用法 3 (详细分析): " << argv[0] << " --analyze <temp> <ratio> <output_directory>" << endl;
            cerr << "用法 4 (参数扫描): " << argv[0] << " --sweep <temp> <ratio> <output_directory> <index>" << endl;
            cerr << "用法 5 (分子量分布分析): " << argv[0] << " --mmd <temp> <ratio> <output_directory>" << endl;
            return 1;
        }
    } else {
        // 默认运行模式下也启用分子量分布追踪
        config.enable_mmd_tracking = true;
    }

    try {
        if (!g_is_bo_mode && !g_is_sweep_mode) {
            if (argc > 1 && (string(argv[1]) == "--analyze" || string(argv[1]) == "--mmd")) {
                cout << "=== 开始详细分析（含分子量分布追踪） ===" << endl;
                cout << "  分析参数: Temp=" << config.poly_temp << " K, Ratio=" << config.poly_initiator_ratio << endl;
                cout << "  输出目录: " << config.base_output_dir << endl;
                cout << "  分子量分布追踪: " << (config.enable_mmd_tracking ? "启用" : "禁用") << endl;
            } else {
                cout << "=== 运行单次模拟（默认参数） ===" << endl;
                cout << "  分子量分布追踪: " << (config.enable_mmd_tracking ? "启用" : "禁用") << endl;
            }
            createDirectory(config.base_output_dir);
        } else if (g_is_sweep_mode) {
            stringstream ss;
            ss << config.base_output_dir << "/sweep_" << config.sweep_index << "_T" << config.poly_temp << "_R" << config.poly_initiator_ratio;
            config.base_output_dir = ss.str();
            createDirectory(config.base_output_dir);
        }

        // 第一阶段：聚合
        TransitionData transition_data = runPolymerizationPhase(config);

        // 第二阶段：解聚
        double final_time = runDepolymerizationPhase(transition_data, config);

        // 更新缺陷统计
        if (config.enable_defect_analysis) {
            transition_data.defect_stats.depoly_time_to_target = final_time;
            transition_data.defect_stats.recyclability_index =
                    1.0 / (1.0 + transition_data.defect_stats.time_to_target + final_time);
        }

        if (g_is_bo_mode) {
            // 贝叶斯优化模式：输出总时间（聚合+解聚）
            double total_time = transition_data.defect_stats.time_to_target + final_time;
            cout << fixed << setprecision(10) << total_time << endl;
        } else if (g_is_sweep_mode) {
            // 参数扫描模式：输出结构化数据
            cout << config.poly_temp << " " << config.poly_initiator_ratio << " "
                 << transition_data.defect_stats.unsaturated_fraction << " "
                 << transition_data.defect_stats.saturated_fraction << " "
                 << transition_data.defect_stats.hh_bond_fraction << " "
                 << transition_data.defect_stats.avg_chain_length << " "
                 << transition_data.defect_stats.time_to_target << " "
                 << final_time << " "
                 << transition_data.defect_stats.recyclability_index << endl;

            // 导出详细缺陷分析
            if (config.enable_defect_analysis) {
                exportDefectAnalysis(config.base_output_dir, transition_data.defect_stats,
                                     config.poly_temp, config.poly_initiator_ratio);
            }
        } else {
            cout << "\n=== 模拟完成 ===" << endl;
            cout << "达到 " << config.depoly_target_conversion * 100 << "% 解聚转化率所需时间为: " << final_time << " s" << endl;
            cout << "总循环时间为: " << (transition_data.defect_stats.time_to_target + final_time) << " s" << endl;
            cout << "所有输出文件已保存到: " << config.base_output_dir << endl;

            if (config.enable_defect_analysis) {
                cout << "\n=== 缺陷分析总结 ===" << endl;
                cout << "不饱和链: " << transition_data.defect_stats.unsaturated_fraction * 100 << "%" << endl;
                cout << "饱和链: " << transition_data.defect_stats.saturated_fraction * 100 << "%" << endl;
                cout << "HH键链: " << transition_data.defect_stats.hh_bond_fraction * 100 << "%" << endl;
                cout << "平均链长: " << transition_data.defect_stats.avg_chain_length << endl;
                cout << "可回收性指数: " << transition_data.defect_stats.recyclability_index << endl;

                exportDefectAnalysis(config.base_output_dir, transition_data.defect_stats,
                                     config.poly_temp, config.poly_initiator_ratio);
            }

            if (config.enable_mmd_tracking) {
                cout << "\n=== 分子量分布分析 ===" << endl;
                cout << "初始分子量分布 (聚合后):" << endl;
                cout << "  Mn: " << transition_data.initial_mmd.Mn << " g/mol" << endl;
                cout << "  Mw: " << transition_data.initial_mmd.Mw << " g/mol" << endl;
                cout << "  PDI: " << transition_data.initial_mmd.PDI << endl;
                cout << "分子量分布文件已保存到: " << config.base_output_dir << "/molecular_weight_distributions/" << endl;
            }
        }

    } catch (const exception& e) {
        cerr << "错误: " << e.what() << endl;
        return 1;
    }

    return 0;
}