import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import pearsonr, spearmanr
from scipy.interpolate import griddata
import warnings
warnings.filterwarnings('ignore')

# Try importing scikit-optimize libraries
try:
    from skopt import gp_minimize
    from skopt.space import Categorical
    from skopt.utils import use_named_args
    from skopt.plots import plot_convergence, plot_objective, plot_evaluations
except ImportError:
    print("Warning: scikit-optimize not found. Bayesian optimization will be disabled.")
    print("Install with: pip install scikit-optimize")

# =================================================================
# 1. CONFIGURATION SECTION
# =================================================================

# Path configuration
KMC_EXECUTABLE_PATH = "/Users/yueyue/CLionProjects/untitled3/cmake-build-debug/integrated_kmc"

# =================================================================
# 2. ANALYSIS MODES CONFIGURATION
# =================================================================

# Mode 1: Single detailed analysis
DO_SINGLE_ANALYSIS = False
ANALYSIS_TEMP = 343.0
ANALYSIS_RATIO = 0.003
ANALYSIS_OUTPUT_DIR = os.path.join("output", f"analysis_results_T{ANALYSIS_TEMP}_R{ANALYSIS_RATIO}")

# Mode 2: Parameter sweep for defect investigation
DO_PARAMETER_SWEEP = False
SWEEP_OUTPUT_DIR = os.path.join("output", "parameter_sweep_results")

# Temperature range for parameter sweep (K)
TEMP_RANGE = np.arange(323.0, 373.0 + 0.1, 10.0)  # Every 10K from 323K to 373K
# Initiator ratio range
RATIO_RANGE = np.arange(0.001, 0.005 + 1e-6, 0.001)  # Every 0.001 from 0.001 to 0.005

# Mode 3: Bayesian optimization with chain length constraint
DO_OPTIMIZATION = True

# Optimization parameter space (can be different from sweep)
#temp_values = list(np.arange(323.0, 373.0 + 1e-6, 2.0))
#ratio_values = list(np.arange(0.001, 0.005 + 1e-9, 0.0001))
temp_values = list(np.arange(323.0, 373.0 + 1e-6, 1.0))
ratio_values = list(np.arange(0.001, 0.005 + 1e-9, 0.0002))

if 'gp_minimize' in globals():
    dim_temp = Categorical(categories=temp_values, name='poly_temp')
    dim_ratio = Categorical(categories=ratio_values, name='poly_initiator_ratio')
    dimensions = [dim_temp, dim_ratio]

# Optimization settings
N_CALLS = 120
N_INITIAL_POINTS = 20

# =================================================================
# 3. CHAIN LENGTH CONSTRAINT CONFIGURATION
# =================================================================

# Chain length constraint parameters
MIN_CHAIN_LENGTH = 1000.0  # 最小平均链长要求
CHAIN_LENGTH_PENALTY_FACTOR = 1000.0  # 违反约束时的惩罚因子
USE_SOFT_CONSTRAINT = True  # True: 软约束(惩罚), False: 硬约束(拒绝)

# 用于记录优化历史的全局变量
optimization_history = {
    'temperatures': [],
    'ratios': [],
    'total_times': [],
    'chain_lengths': [],
    'constraint_violations': [],
    'objective_values': []
}

# =================================================================
# 4. FUNCTION DEFINITIONS
# =================================================================

def run_single_analysis(temp, ratio, output_dir):
    """Run single detailed analysis with comprehensive output."""
    print("="*60)
    print("🔬 RUNNING SINGLE DETAILED ANALYSIS")
    print(f"  📊 Temperature: {temp:.1f} K")
    print(f"  📊 Initiator Ratio: {ratio:.4f}")
    print(f"  📁 Output Directory: '{output_dir}/'")
    print("="*60)

    os.makedirs(output_dir, exist_ok=True)

    if not os.path.exists(KMC_EXECUTABLE_PATH):
        print(f"❌ Error: Executable not found at '{KMC_EXECUTABLE_PATH}'")
        return

    command = [
        KMC_EXECUTABLE_PATH,
        "--analyze",
        str(temp),
        str(ratio),
        output_dir
    ]

    try:
        print("\n📝 C++ Simulation Log:")
        print("-" * 40)
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)

        for line in iter(process.stdout.readline, ''):
            print(line, end='')
        process.stdout.close()
        return_code = process.wait()

        if return_code == 0:
            print("\n✅ Analysis completed successfully!")
            print(f"📁 Check detailed output files in '{output_dir}' folder.")
        else:
            print(f"\n❌ Analysis failed with exit code {return_code}")

    except Exception as e:
        print(f"❌ Error occurred during analysis: {e}")

    print("\n")

def run_parameter_sweep():
    """Run comprehensive parameter sweep to investigate defect formation."""
    print("="*60)
    print("🔍 PARAMETER SWEEP: INVESTIGATING DEFECT FORMATION")
    print(f"  🌡️  Temperature range: {TEMP_RANGE[0]:.1f} - {TEMP_RANGE[-1]:.1f} K ({len(TEMP_RANGE)} points)")
    print(f"  ⚗️  Ratio range: {RATIO_RANGE[0]:.4f} - {RATIO_RANGE[-1]:.4f} ({len(RATIO_RANGE)} points)")
    print(f"  📊 Total simulations: {len(TEMP_RANGE) * len(RATIO_RANGE)}")
    print("="*60)

    os.makedirs(SWEEP_OUTPUT_DIR, exist_ok=True)

    if not os.path.exists(KMC_EXECUTABLE_PATH):
        print(f"❌ Error: Executable not found at '{KMC_EXECUTABLE_PATH}'")
        return None

    results = []
    total_sims = len(TEMP_RANGE) * len(RATIO_RANGE)
    current_sim = 0

    for temp in TEMP_RANGE:
        for ratio in RATIO_RANGE:
            current_sim += 1
            print(f"\n⏳ Running simulation {current_sim}/{total_sims}: T={temp:.1f}K, R={ratio:.4f}")

            command = [
                KMC_EXECUTABLE_PATH,
                "--sweep",
                str(temp),
                str(ratio),
                SWEEP_OUTPUT_DIR,
                str(current_sim)
            ]

            try:
                result = subprocess.run(command, capture_output=True, text=True, check=True, timeout=600)

                # Parse the structured output from C++ program
                # Format: temp ratio unsaturated_frac saturated_frac hh_frac avg_length poly_time depoly_time recyclability
                output_line = result.stdout.strip()
                if output_line:
                    data = output_line.split()
                    if len(data) >= 9:
                        results.append({
                            'temperature': float(data[0]),
                            'initiator_ratio': float(data[1]),
                            'unsaturated_fraction': float(data[2]),
                            'saturated_fraction': float(data[3]),
                            'hh_bond_fraction': float(data[4]),
                            'avg_chain_length': float(data[5]),
                            'polymerization_time': float(data[6]),
                            'depolymerization_time': float(data[7]),
                            'recyclability_index': float(data[8])
                        })
                        print(f"✅ Completed: Recyclability = {float(data[8]):.6f}, Chain Length = {float(data[5]):.1f}")
                    else:
                        print(f"⚠️  Warning: Unexpected output format")
                else:
                    print(f"⚠️  Warning: No output received")

            except subprocess.TimeoutExpired:
                print(f"⏰ Timeout for T={temp:.1f}K, R={ratio:.4f}")
                results.append({
                    'temperature': temp,
                    'initiator_ratio': ratio,
                    'unsaturated_fraction': np.nan,
                    'saturated_fraction': np.nan,
                    'hh_bond_fraction': np.nan,
                    'avg_chain_length': np.nan,
                    'polymerization_time': np.nan,
                    'depolymerization_time': np.nan,
                    'recyclability_index': np.nan
                })
            except Exception as e:
                print(f"❌ Error for T={temp:.1f}K, R={ratio:.4f}: {e}")
                results.append({
                    'temperature': temp,
                    'initiator_ratio': ratio,
                    'unsaturated_fraction': np.nan,
                    'saturated_fraction': np.nan,
                    'hh_bond_fraction': np.nan,
                    'avg_chain_length': np.nan,
                    'polymerization_time': np.nan,
                    'depolymerization_time': np.nan,
                    'recyclability_index': np.nan
                })

    # Convert to DataFrame and save
    df = pd.DataFrame(results)
    csv_file = os.path.join(SWEEP_OUTPUT_DIR, "parameter_sweep_results.csv")
    df.to_csv(csv_file, index=False)

    print(f"\n✅ Parameter sweep completed!")
    print(f"📊 Results saved to: {csv_file}")
    print(f"📈 Valid results: {df.dropna().shape[0]}/{len(df)}")

    return df

def analyze_defect_correlations(df):
    """Analyze correlations between polymerization conditions and defect formation."""
    print("\n" + "="*60)
    print("📊 CORRELATION ANALYSIS: CONDITIONS vs DEFECTS")
    print("="*60)

    if df is None or df.empty:
        print("❌ No data available for correlation analysis")
        return

    # Remove rows with NaN values for correlation analysis
    df_clean = df.dropna()

    if df_clean.empty:
        print("❌ No valid data for correlation analysis")
        return

    # Define correlation pairs
    condition_vars = ['temperature', 'initiator_ratio']
    defect_vars = ['unsaturated_fraction', 'saturated_fraction', 'hh_bond_fraction', 'avg_chain_length']
    performance_vars = ['polymerization_time', 'depolymerization_time', 'recyclability_index']

    print("\n🔗 CORRELATIONS: Polymerization Conditions → Defect Formation")
    print("-" * 60)

    for condition in condition_vars:
        for defect in defect_vars:
            if condition in df_clean.columns and defect in df_clean.columns:
                pearson_r, pearson_p = pearsonr(df_clean[condition], df_clean[defect])
                spearman_r, spearman_p = spearmanr(df_clean[condition], df_clean[defect])

                print(f"{condition.replace('_', ' ').title()} → {defect.replace('_', ' ').title()}:")
                print(f"  📈 Pearson:  r = {pearson_r:+.3f} (p = {pearson_p:.3f})")
                print(f"  📊 Spearman: ρ = {spearman_r:+.3f} (p = {spearman_p:.3f})")

    print("\n🔗 CORRELATIONS: Defects → Recyclability Performance")
    print("-" * 60)

    for defect in defect_vars:
        for performance in performance_vars:
            if defect in df_clean.columns and performance in df_clean.columns:
                pearson_r, pearson_p = pearsonr(df_clean[defect], df_clean[performance])
                spearman_r, spearman_p = spearmanr(df_clean[defect], df_clean[performance])

                print(f"{defect.replace('_', ' ').title()} → {performance.replace('_', ' ').title()}:")
                print(f"  📈 Pearson:  r = {pearson_r:+.3f} (p = {pearson_p:.3f})")
                print(f"  📊 Spearman: ρ = {spearman_r:+.3f} (p = {spearman_p:.3f})")

def create_comprehensive_defect_analysis_plots(df):
    """
    Create a comprehensive 3x4 visualization grid of the parameter sweep results.
    Final version with custom legends for MMD plot.
    """
    if df is None or df.empty:
        print("❌ No data available for plotting")
        return

    print("\n" + "="*60)
    print("📊 CREATING RESTRUCTURED 3x4 ANALYSIS PLOT (FINAL LAYOUT)")
    print("="*60)

    # ... (数据预处理和网格创建部分保持不变) ...
    plt.style.use('default')
    df_clean = df.dropna().copy()
    if df_clean.empty: return
    df_clean['total_time'] = df_clean['polymerization_time'] + df_clean['depolymerization_time']
    df_clean.rename(columns={'recyclability_index': 'process_efficiency'}, inplace=True)
    new_metric_name = 'Process Efficiency Index'
    fig = plt.figure(figsize=(26, 18))
    gs = fig.add_gridspec(3, 4, wspace=0.4, hspace=0.45)

    # --- 绘制第一行 (热力图 - 保持不变) ---
    heatmap_vars = {
        'unsaturated_fraction': 'Unsaturated Chains Fraction',
        'saturated_fraction': 'Saturated Chains Fraction',
        'hh_bond_fraction': 'HH-Bond Fraction',
        'avg_chain_length': 'Average Chain Length'
    }
    for i, (var, label) in enumerate(heatmap_vars.items()):
        ax = fig.add_subplot(gs[0, i])
        pivot_data = df_clean.pivot_table(values=var, index='temperature', columns='initiator_ratio')
        sns.heatmap(pivot_data, annot=True, fmt='.3g', cmap='viridis', ax=ax, cbar_kws={'label': label})
        ax.set_xlabel('Initiator Ratio')
        ax.set_ylabel('Temperature (K)')

    # --- 绘制第二行 ---

    # 图 2-1 & 2-2: 箱线图 (保持不变)
    boxplot_palette = ["#e41a1c", "#377eb8", "#4daf4a"]
    for i, col in enumerate(['temperature', 'initiator_ratio']):
        ax_box = fig.add_subplot(gs[1, i])
        range_col = f'{col}_range'
        bins = min(5, len(df_clean[col].unique()))
        df_clean[range_col] = pd.cut(df_clean[col], bins=bins)
        melted_data = df_clean.melt(id_vars=[range_col], value_vars=['unsaturated_fraction', 'saturated_fraction', 'hh_bond_fraction'], var_name='Defect Type', value_name='Defect Fraction')
        sns.boxplot(data=melted_data, x=range_col, y='Defect Fraction', hue='Defect Type', ax=ax_box, palette=boxplot_palette, showfliers=False)
        if col == 'temperature':
            labels = [f"({cat.left:.0f}, {cat.right:.0f}]" for cat in df_clean[range_col].cat.categories]
            ax_box.set_xlabel('Temperature Range (K)')
        else:
            labels = [f"({cat.left:.4f}, {cat.right:.4f}]" for cat in df_clean[range_col].cat.categories]
            ax_box.set_xlabel('Initiator Ratio Range')
        ax_box.set_xticklabels(labels, rotation=45, ha="right")
        ax_box.legend(title='Defect Type')


    # *** --- START OF MODIFICATION --- ***
    # 图 2-3: Log Molar Mass Distribution by Temperature (新绘图逻辑)
    ax_mmd = fig.add_subplot(gs[1, 2])
    df_clean['log_avg_M'] = np.log10(df_clean['avg_chain_length'] * 100.0)

    # 1. 定义温度分组
    num_temp_groups = min(5, len(df_clean['temperature'].unique()))
    df_clean['temp_group_mmd'] = pd.cut(df_clean['temperature'], bins=num_temp_groups)

    # 2. 获取分组和颜色
    temp_groups = df_clean.groupby('temp_group_mmd')
    colors = plt.cm.get_cmap('viridis', temp_groups.ngroups)

    # 3. 遍历每个分组，手动绘图并指定label
    for i, (name, group) in enumerate(temp_groups):
        if group.empty:
            continue
        # 创建可读的标签
        label = f"{name.left:.0f}-{name.right:.0f}K"
        sns.kdeplot(data=group, x='log_avg_M', ax=ax_mmd, color=colors(i),
                    fill=True, common_norm=False, alpha=0.5, linewidth=1.5, label=label)

    ax_mmd.set_xlabel(r'$\log(M_{avg})$')
    ax_mmd.set_ylabel('Density')
    ax_mmd.legend(title='Temperature Range') # 现在图例将显示我们自定义的标签
    # *** --- END OF MODIFICATION --- ***


    # 图 2-4: Correlation Matrix
    ax_corr = fig.add_subplot(gs[1, 3])
    corr_cols = ['temperature', 'initiator_ratio', 'unsaturated_fraction', 'saturated_fraction', 'hh_bond_fraction', 'avg_chain_length', 'process_efficiency']
    corr_matrix = df_clean[corr_cols].corr()
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    sns.heatmap(corr_matrix, mask=mask, annot=True, fmt='.2f', cmap='coolwarm', center=0, ax=ax_corr)
    ax_corr.set_xticklabels(ax_corr.get_xticklabels(), rotation=45, ha="right")
    ax_corr.set_yticklabels(ax_corr.get_yticklabels(), rotation=30, style='italic', ha='right')

    # ... (第三行和结尾部分的代码保持不变) ...
    # --- 5. 绘制第三行 ---
    temp_grid, ratio_grid = np.meshgrid(np.linspace(df_clean['temperature'].min(), df_clean['temperature'].max(), 100),
                                        np.linspace(df_clean['initiator_ratio'].min(), df_clean['initiator_ratio'].max(), 100))
    ax_poly_depoly = fig.add_subplot(gs[2, 0])
    scatter_pd = ax_poly_depoly.scatter(df_clean['polymerization_time'], df_clean['depolymerization_time'], c=df_clean['process_efficiency'], cmap='viridis', s=50, edgecolors='k', linewidth=0.5)
    ax_poly_depoly.set_xlabel('Polymerization Time (s)')
    ax_poly_depoly.set_ylabel('Depolymerization Time (s)')
    fig.colorbar(scatter_pd, ax=ax_poly_depoly, label=new_metric_name)
    ax_eff_time = fig.add_subplot(gs[2, 1])
    scatter_eff = ax_eff_time.scatter(df_clean['total_time'], df_clean['process_efficiency'], c=df_clean['temperature'], cmap='coolwarm', s=50, edgecolors='k', linewidth=0.5)
    ax_eff_time.set_xlabel('Total Cycle Time (s)')
    ax_eff_time.set_ylabel(new_metric_name)
    fig.colorbar(scatter_eff, ax=ax_eff_time, label='Temperature (K)')
    ax_eff_landscape = fig.add_subplot(gs[2, 2])
    efficiency_grid = griddata((df_clean['temperature'], df_clean['initiator_ratio']), df_clean['process_efficiency'], (temp_grid, ratio_grid), method='cubic')
    if not np.all(np.isnan(efficiency_grid)):
        contour_eff = ax_eff_landscape.contourf(temp_grid, ratio_grid, efficiency_grid, levels=20, cmap='RdYlBu')
        fig.colorbar(contour_eff, ax=ax_eff_landscape, label=new_metric_name)
    ax_eff_landscape.set_xlabel('Temperature (K)')
    ax_eff_landscape.set_ylabel('Initiator Ratio')
    opt_idx = df_clean['process_efficiency'].idxmax()
    ax_eff_landscape.plot(df_clean.loc[opt_idx, 'temperature'], df_clean.loc[opt_idx, 'initiator_ratio'], 'r*', markersize=15, markeredgecolor='white', label='Max Efficiency')
    ax_eff_landscape.legend()
    ax_time_landscape = fig.add_subplot(gs[2, 3])
    time_grid = griddata((df_clean['temperature'], df_clean['initiator_ratio']), df_clean['total_time'], (temp_grid, ratio_grid), method='cubic')
    if not np.all(np.isnan(time_grid)):
        contour_time = ax_time_landscape.contourf(temp_grid, ratio_grid, time_grid, levels=20, cmap='viridis_r')
        fig.colorbar(contour_time, ax=ax_time_landscape, label='Total Time (s)')
    ax_time_landscape.set_xlabel('Temperature (K)')
    ax_time_landscape.set_ylabel('Initiator Ratio')
    min_time_idx = df_clean['total_time'].idxmin()
    ax_time_landscape.plot(df_clean.loc[min_time_idx, 'temperature'], df_clean.loc[min_time_idx, 'initiator_ratio'], 'r*', markersize=15, markeredgecolor='white', label='Min Time')
    ax_time_landscape.legend()

    # --- 6. 添加总标题并保存 ---
    fig.suptitle('Comprehensive Analysis of Polymerization-Depolymerization Cycle', fontsize=28, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    output_filename = 'structured_analysis_dashboard_final_v3.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"✅ Final restructured analysis plot saved: '{output_filename}'")

def get_chain_length_from_simulation(poly_temp, poly_initiator_ratio):
    """获取仿真结果中的平均链长信息"""
    if not os.path.exists(KMC_EXECUTABLE_PATH):
        print(f"❌ Error: Executable not found at '{KMC_EXECUTABLE_PATH}'")
        return None, None

    command = [KMC_EXECUTABLE_PATH, "--sweep", str(poly_temp), str(poly_initiator_ratio), "/tmp", "0"]

    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True, timeout=600)
        output_line = result.stdout.strip()
        if output_line:
            data = output_line.split()
            if len(data) >= 9:
                total_time = float(data[6]) + float(data[7])  # poly_time + depoly_time
                avg_chain_length = float(data[5])
                return total_time, avg_chain_length
    except Exception as e:
        print(f"❌ [Chain Length Check] Error occurred: {e}")

    return None, None

@use_named_args(dimensions=dimensions if 'dimensions' in globals() else [])
def objective_with_chain_length_constraint(poly_temp, poly_initiator_ratio):
    """带链长约束的目标函数 - 现在同时优化总时间并确保链长满足要求"""
    global optimization_history

    print(f"🎯 [Optimizer] Evaluating: Temp = {poly_temp:.1f} K, Ratio = {poly_initiator_ratio:.4f}")

    if not os.path.exists(KMC_EXECUTABLE_PATH):
        print(f"❌ Error: Executable not found at '{KMC_EXECUTABLE_PATH}'")
        return 1e10

    command = [KMC_EXECUTABLE_PATH, "--sweep", str(poly_temp), str(poly_initiator_ratio), "/tmp", "0"]

    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True, timeout=600)
        output_line = result.stdout.strip()

        if output_line:
            data = output_line.split()
            if len(data) >= 9:
                total_time = float(data[6]) + float(data[7])  # poly_time + depoly_time
                avg_chain_length = float(data[5])

                # 记录历史数据
                optimization_history['temperatures'].append(poly_temp)
                optimization_history['ratios'].append(poly_initiator_ratio)
                optimization_history['total_times'].append(total_time)
                optimization_history['chain_lengths'].append(avg_chain_length)

                # 检查链长约束
                constraint_violation = max(0, MIN_CHAIN_LENGTH - avg_chain_length)
                optimization_history['constraint_violations'].append(constraint_violation)

                if USE_SOFT_CONSTRAINT:
                    # 软约束：违反约束时添加惩罚
                    if constraint_violation > 0:
                        penalty = CHAIN_LENGTH_PENALTY_FACTOR * (constraint_violation / MIN_CHAIN_LENGTH) ** 2
                        objective_value = total_time + penalty
                        print(f"⚠️  [Optimizer] Chain length {avg_chain_length:.1f} < {MIN_CHAIN_LENGTH:.1f}, penalty = {penalty:.2f}")
                        print(f"✅ [Optimizer] Total Time: {total_time:.4f} s, Penalized Objective: {objective_value:.4f}")
                    else:
                        objective_value = total_time
                        print(f"✅ [Optimizer] Total Time: {total_time:.4f} s, Chain Length: {avg_chain_length:.1f} (OK)")
                else:
                    # 硬约束：违反约束时直接返回大值
                    if constraint_violation > 0:
                        objective_value = 1e8
                        print(f"❌ [Optimizer] Chain length {avg_chain_length:.1f} < {MIN_CHAIN_LENGTH:.1f}, REJECTED")
                    else:
                        objective_value = total_time
                        print(f"✅ [Optimizer] Total Time: {total_time:.4f} s, Chain Length: {avg_chain_length:.1f} (OK)")

                optimization_history['objective_values'].append(objective_value)
                print()
                return objective_value

    except Exception as e:
        print(f"❌ [Optimizer] Error occurred: {e}")

    # 失败情况
    optimization_history['temperatures'].append(poly_temp)
    optimization_history['ratios'].append(poly_initiator_ratio)
    optimization_history['total_times'].append(1e10)
    optimization_history['chain_lengths'].append(0)
    optimization_history['constraint_violations'].append(MIN_CHAIN_LENGTH)
    optimization_history['objective_values'].append(1e10)
    return 1e10

def run_bayesian_optimization():
    """运行带链长约束的贝叶斯优化"""
    if 'gp_minimize' not in globals():
        print("❌ Scikit-optimize not available. Skipping Bayesian optimization.")
        return None

    print("="*80)
    print("🎯 CONSTRAINED BAYESIAN OPTIMIZATION: FINDING OPTIMAL CONDITIONS")
    print(f"Objective: Minimize total cycle time (polymerization + depolymerization)")
    print(f"Constraint: Average chain length ≥ {MIN_CHAIN_LENGTH:.1f}")
    print(f"Constraint Type: {'Soft (penalty)' if USE_SOFT_CONSTRAINT else 'Hard (rejection)'}")
    if USE_SOFT_CONSTRAINT:
        print(f"Penalty Factor: {CHAIN_LENGTH_PENALTY_FACTOR:.0f}")
    print(f"Parameter Space:")
    print(f"  🌡️  Temperature (K): {temp_values}")
    print(f"  ⚗️  Initiator Ratio: {[f'{r:.4f}' for r in ratio_values]}")
    print(f"Running {N_CALLS} evaluations...")
    print("="*80)

    # 清空历史记录
    global optimization_history
    optimization_history = {
        'temperatures': [],
        'ratios': [],
        'total_times': [],
        'chain_lengths': [],
        'constraint_violations': [],
        'objective_values': []
    }

    # 运行优化
    result_gp = gp_minimize(
        func=objective_with_chain_length_constraint,
        dimensions=dimensions,
        n_calls=N_CALLS,
        n_initial_points=N_INITIAL_POINTS,
        #acq_func='EI',
        acq_func='LCB',      # 使用 LCB 采集函数
        kappa=5.0,
        random_state=123
    )

    # 提取结果
    best_temp = result_gp.x[0]
    best_ratio = result_gp.x[1]
    best_objective = result_gp.fun

    # 获取最优点的实际性能指标
    best_total_time, best_chain_length = get_chain_length_from_simulation(best_temp, best_ratio)

    print(f"\n🏆 CONSTRAINED OPTIMIZATION RESULTS")
    print(f"  🌡️  Best Temperature: {best_temp:.1f} K")
    print(f"  ⚗️  Best Ratio: {best_ratio:.4f}")
    if best_total_time is not None and best_chain_length is not None:
        print(f"  ⏱️  Total Time: {best_total_time:.2f} s")
        print(f"  🔗 Chain Length: {best_chain_length:.1f}")
        print(f"  ✅ Constraint: {'SATISFIED' if best_chain_length >= MIN_CHAIN_LENGTH else 'VIOLATED'}")
    print(f"  🎯 Best Objective: {best_objective:.2f}")
    print("-" * 60)

    # 分析约束满足情况
    history_df = pd.DataFrame(optimization_history)
    satisfied_count = sum(1 for cl in history_df['chain_lengths'] if cl >= MIN_CHAIN_LENGTH)
    print(f"\n📊 CONSTRAINT SATISFACTION ANALYSIS:")
    print(f"  ✅ Satisfied evaluations: {satisfied_count}/{len(history_df)} ({satisfied_count/len(history_df)*100:.1f}%)")
    print(f"  📈 Average chain length: {history_df['chain_lengths'].mean():.1f}")
    print(f"  📊 Chain length range: {history_df['chain_lengths'].min():.1f} - {history_df['chain_lengths'].max():.1f}")

    # 创建优化图表
    create_constrained_optimization_plots(result_gp, best_temp, best_ratio, best_objective, history_df)

    return result_gp

def create_constrained_optimization_plots(result_gp, best_temp, best_ratio, best_objective, history_df):
    """
    Create a single, publication-quality dashboard visualizing the constrained optimization process.
    Titles on subplots are removed.
    """
    print("\n📊 Generating Enhanced Optimization Dashboard...")

    # --- 1. 设置高级绘图风格 ---
    plt.style.use('seaborn-v0_8-ticks')
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 14,
        'axes.labelsize': 16,
        'axes.titlesize': 18,
        'axes.labelweight': 'bold',
        'axes.titleweight': 'bold',
        'legend.fontsize': 12,
        'xtick.major.size': 7,
        'ytick.major.size': 7,
        'xtick.minor.size': 4,
        'ytick.minor.size': 4,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'axes.linewidth': 1.5,
    })

    # --- 2. 创建图形和布局 ---
    fig = plt.figure(figsize=(20, 10))
    gs = fig.add_gridspec(2, 2, width_ratios=(1.8, 1), height_ratios=(1, 1), hspace=0.35, wspace=0.25)

    ax_main = fig.add_subplot(gs[:, 0])
    ax_conv = fig.add_subplot(gs[0, 1])
    ax_constraint = fig.add_subplot(gs[1, 1])

    # --- 3. 主图：参数空间探索 ---
    satisfied_mask = history_df['chain_lengths'] >= MIN_CHAIN_LENGTH

    temp_grid, ratio_grid = np.meshgrid(
        np.linspace(min(temp_values), max(temp_values), 100),
        np.linspace(min(ratio_values), max(ratio_values), 100)
    )
    time_grid = griddata((history_df['temperatures'], history_df['ratios']),
                         history_df['total_times'],
                         (temp_grid, ratio_grid), method='cubic')

    contour = ax_main.contourf(temp_grid, ratio_grid, time_grid, levels=20, cmap='viridis_r', alpha=0.7)
    cbar = fig.colorbar(contour, ax=ax_main, orientation='vertical', pad=0.02)
    cbar.set_label('Interpolated Total Cycle Time (s)', rotation=270, labelpad=20)

    ax_main.scatter(
        history_df.loc[satisfied_mask, 'temperatures'],
        history_df.loc[satisfied_mask, 'ratios'],
        c='lime', s=100,
        marker='o', edgecolors='k', linewidth=1.5,
        alpha=1.0, label='Constraint Satisfied'
    )

    ax_main.scatter(
        history_df.loc[~satisfied_mask, 'temperatures'],
        history_df.loc[~satisfied_mask, 'ratios'],
        c='red', s=120,
        marker='x', linewidth=3,
        alpha=1.0, label='Constraint Violated'
    )

    ax_main.plot(best_temp, best_ratio, marker='*', color='gold', markersize=25,
                 markeredgecolor='black', markeredgewidth=1.5,
                 label=f'Optimum ({best_temp:.1f}K, {best_ratio:.4f})')

    ax_main.set_xlabel("Temperature (K)")
    ax_main.set_ylabel("Initiator Ratio")
    # ax_main.set_title("Parameter Space Exploration") # <-- 已移除
    ax_main.legend(loc='best')
    ax_main.minorticks_on()

    # --- 4. 右上：收敛图 ---
    from skopt.plots import plot_convergence
    plot_convergence(result_gp, ax=ax_conv)
    # ax_conv.set_title("Convergence of Objective Function") # <-- 标题已被移除
    ax_conv.set_xlabel("Number of Calls")
    ax_conv.set_ylabel("Best Objective Value Found")
    ax_conv.grid(True, which='major', linestyle='--', alpha=0.6)
    ax_conv.minorticks_on()

    # --- 5. 右下：约束满足情况 ---
    iterations = np.arange(1, len(history_df) + 1)
    chain_lengths = history_df['chain_lengths'].values

    ax_constraint.plot(iterations, chain_lengths, marker='', ls='-', color='lightgray', lw=2)

    ax_constraint.scatter(iterations[satisfied_mask], chain_lengths[satisfied_mask],
                          color='green', s=80, zorder=5, label='Satisfied', edgecolors='k', linewidth=1)
    ax_constraint.scatter(iterations[~satisfied_mask], chain_lengths[~satisfied_mask],
                          color='red', s=80, zorder=5, label='Violated', edgecolors='k', linewidth=1)

    ax_constraint.axhline(y=MIN_CHAIN_LENGTH, color='blue', linestyle='--', linewidth=2.5,
                          label=f'Constraint (≥ {MIN_CHAIN_LENGTH:.0f})')

    ax_constraint.set_xlabel("Number of Calls")
    ax_constraint.set_ylabel("Average Chain Length")
    # ax_constraint.set_title("Chain Length Constraint Fulfillment") # <-- 已移除
    ax_constraint.legend()
    ax_constraint.set_ylim(bottom=0)
    ax_constraint.set_xlim(0, len(iterations) + 1)
    ax_constraint.minorticks_on()

    # --- 6. 添加总标题并保存 ---
    fig.suptitle(f"Bayesian Optimization Dashboard (Target: Minimize Time, Constraint: Chain Length ≥ {MIN_CHAIN_LENGTH:.0f})",
                 fontsize=24, y=0.98)

    output_filename = 'new_1000__bayesian_optimization_dashboard_notitle.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"✅ Optimization dashboard (no subtitles) saved: '{output_filename}'")

# =================================================================
# 5. MAIN EXECUTION
# =================================================================

if __name__ == "__main__":
    print("🚀 PMMA LIFECYCLE SIMULATION FRAMEWORK WITH CHAIN LENGTH CONSTRAINTS")
    print("Investigating polymerization conditions → structural defects → recyclability")
    print(f"🔗 Chain Length Constraint: Average length ≥ {MIN_CHAIN_LENGTH}")
    print("="*80)

    # Step 1: Single detailed analysis
    if DO_SINGLE_ANALYSIS:
        run_single_analysis(ANALYSIS_TEMP, ANALYSIS_RATIO, ANALYSIS_OUTPUT_DIR)

    # Step 2: Parameter sweep for defect investigation
    sweep_df = None
    if DO_PARAMETER_SWEEP:
        sweep_df = run_parameter_sweep()

        if sweep_df is not None:
            # Analyze correlations
            analyze_defect_correlations(sweep_df)

            # Create comprehensive plots
            create_comprehensive_defect_analysis_plots(sweep_df)

    # Step 3: Constrained Bayesian optimization
    if DO_OPTIMIZATION:
        optimization_result = run_bayesian_optimization()

    # Final summary
    print("\n" + "="*80)
    print("🎉 CONSTRAINED PMMA LIFECYCLE INVESTIGATION COMPLETE!")
    print("="*80)

    summary_parts = []
    if DO_SINGLE_ANALYSIS:
        summary_parts.append(f"✅ Single analysis completed (T={ANALYSIS_TEMP}K, R={ANALYSIS_RATIO})")
    if DO_PARAMETER_SWEEP and sweep_df is not None:
        summary_parts.append(f"✅ Parameter sweep completed ({len(sweep_df)} simulations)")
    if DO_OPTIMIZATION:
        summary_parts.append(f"✅ Constrained Bayesian optimization completed (Chain Length ≥ {MIN_CHAIN_LENGTH})")

    for part in summary_parts:
        print(part)

    print("\n📊 Generated Analysis Files:")
    if DO_SINGLE_ANALYSIS:
        print(f"  📁 Single analysis: {ANALYSIS_OUTPUT_DIR}/")
    if DO_PARAMETER_SWEEP:
        print(f"  📁 Parameter sweep: {SWEEP_OUTPUT_DIR}/")
        print(f"  📊 Comprehensive plot: comprehensive_defect_analysis.png")
    if DO_OPTIMIZATION:
        print(f"  📊 Constrained optimization plots: constrained_optimization_*.png")

    print("\n🔬 This enhanced framework provides comprehensive insights into:")
    print("  1️⃣  How temperature and initiator ratio affect defect formation")
    print("  2️⃣  Correlations between structural defects and recyclability")
    print("  3️⃣  Optimal synthesis conditions with chain length constraints")
    print("  4️⃣  Complete PMMA lifecycle modeling from synthesis to recycling")
    print("  5️⃣  Trade-offs between process efficiency, product quality, and chain length")
    print("  6️⃣  Constraint satisfaction analysis in optimization")

    if not any([DO_SINGLE_ANALYSIS, DO_PARAMETER_SWEEP, DO_OPTIMIZATION]):
        print("⚠️  All analysis modes are disabled. Enable at least one mode in the configuration section.")