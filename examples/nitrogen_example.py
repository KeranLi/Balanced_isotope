"""
氮同位素体系使用示例

基于 Kang et al. (2023) 和 Ma et al. (2025) 的氮循环箱式模型

本示例展示如何：
1. 使用正向模型计算沉积物氮同位素
2. 使用反向模型反演硝酸盐占比
3. 进行蒙特卡洛不确定性分析
4. 模拟不同地质时期的氮循环

可选：设置 SAVE_PLOTS = True 保存图表
"""

import numpy as np

# 可选：启用图表保存
SAVE_PLOTS = False
PLOT_DIR = "examples/output"

# 导入氮同位素体系
from systems.n import NIsotopeSystem, get_scenario_info
from toolkit.math.statistics import Bootstrap, LOWESS, ChangepointDetector

# 导入绘图库
import matplotlib.pyplot as plt
import os


def example_1_forward_model():
    """
    示例 1: 正向模型 - 从硝酸盐占比计算沉积物氮同位素
    """
    print("=" * 70)
    print("示例 1: 正向模型计算")
    print("=" * 70)
    
    # 创建模型实例
    n_system = NIsotopeSystem(scenario='modern')
    
    # 不同硝酸盐占比下的沉积物同位素
    f_values = np.linspace(0, 1, 11)
    
    print(f"{'f_assimilator':<15} {'δ¹⁵N_sed (‰)':<15}")
    print("-" * 30)
    
    results = []
    for f in f_values:
        delta15N = n_system.forward_model(f_assimilator=f)
        results.append((f, delta15N))
        print(f"{f:<15.2f} {delta15N:<15.2f}")
    
    print("\n解释:")
    print("- f = 0: 完全由固氮生物贡献，δ¹⁵N 接近 0‰")
    print("- f ≈ 0.48: δ¹⁵N 达到最大值")
    print("- f = 1: 完全由硝酸盐同化生物贡献，δ¹⁵N 再次降低")
    
    # 可选：绘图
    if SAVE_PLOTS:
        os.makedirs(PLOT_DIR, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 5))
        f_range = np.linspace(0, 1, 100)
        delta_range = [n_system.forward_model(f) for f in f_range]
        ax.plot(f_range, delta_range, 'b-', linewidth=2)
        ax.scatter([r[0] for r in results], [r[1] for r in results], 
                  c='red', s=50, zorder=5)
        ax.set_xlabel('f_assimilator')
        ax.set_ylabel('δ¹⁵N_sed (‰)')
        ax.set_title('Forward Model: f_assimilator vs δ¹⁵N_sed')
        ax.grid(True, alpha=0.3)
        plt.savefig(f"{PLOT_DIR}/example_1_forward_model.png", dpi=200)
        plt.close()
        print(f"\n图表已保存: {PLOT_DIR}/example_1_forward_model.png")
    
    return results


def example_2_inverse_model():
    """
    示例 2: 反向模型 - 从沉积物氮同位素反演硝酸盐占比
    """
    print("\n" + "=" * 70)
    print("示例 2: 反向模型反演")
    print("=" * 70)
    
    n_system = NIsotopeSystem(scenario='early_triassic')
    
    # 模拟早三叠世不同阶段的观测数据
    observed_delta15N = [0.5, 1.5, 3.0, 4.5, 2.0]
    stages = ['极端缺氧', '缺氧', '弱氧化', '氧化', '再缺氧']
    
    print(f"{'阶段':<12} {'观测δ¹⁵N':<12} {'f_assimilator':<15} {'解释'}")
    print("-" * 70)
    
    results = []
    for stage, delta15N in zip(stages, observed_delta15N):
        # 反演计算 (限制在合理范围 0-0.48 对应缺氧环境)
        result = n_system.inverse_model(delta15N_sed=delta15N, 
                                       f_range=(0.0, 0.48))
        f = result['f_assimilator']
        
        # 解释
        if f < 0.1:
            interp = "硝酸盐极度匮乏"
        elif f < 0.2:
            interp = "硝酸盐受限"
        else:
            interp = "硝酸盐相对充足"
        
        print(f"{stage:<12} {delta15N:<12.2f} {f:<15.3f} {interp}")
        results.append((delta15N, f, stage))
    
    return results


def example_3_monte_carlo():
    """
    示例 3: 蒙特卡洛不确定性分析
    """
    print("\n" + "=" * 70)
    print("示例 3: 蒙特卡洛不确定性分析")
    print("=" * 70)
    
    n_system = NIsotopeSystem(scenario='neoproterozoic')
    
    # 新元古代早期假设的硝酸盐占比
    f_assimilator = 0.11  # Kang et al. 2023 估计值
    
    print(f"假设硝酸盐占比 f_assimilator = {f_assimilator}")
    print("\n进行 10,000 次蒙特卡洛模拟...")
    
    # 蒙特卡洛模拟
    mc_result = n_system.monte_carlo_simulation(
        f_assimilator=f_assimilator,
        n_samples=10000,
        epsilon_fix_range=(-2.0, 1.0),
        epsilon_wcd_range=(-30.0, -22.0)
    )
    
    print(f"\n结果:")
    print(f"  δ¹⁵N_sed 均值: {mc_result['delta15N_sed_mean']:.2f}‰")
    print(f"  δ¹⁵N_sed 中位数: {mc_result['delta15N_sed_median']:.2f}‰")
    print(f"  δ¹⁵N_sed 标准差: {mc_result['delta15N_sed_std']:.2f}‰")
    print(f"  68% 置信区间: [{mc_result['delta15N_sed_ci68'][0]:.2f}, "
          f"{mc_result['delta15N_sed_ci68'][1]:.2f}]‰")
    print(f"  95% 置信区间: [{mc_result['delta15N_sed_ci95'][0]:.2f}, "
          f"{mc_result['delta15N_sed_ci95'][1]:.2f}]‰")
    
    # 可选：绘图
    if SAVE_PLOTS:
        os.makedirs(PLOT_DIR, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 5))
        samples = mc_result['samples']
        ax.hist(samples, bins=50, density=True, alpha=0.7, edgecolor='black')
        ax.axvline(mc_result['delta15N_sed_mean'], color='red', 
                  linestyle='--', linewidth=2, label=f"Mean: {mc_result['delta15N_sed_mean']:.2f}‰")
        ax.axvline(mc_result['delta15N_sed_median'], color='green',
                  linestyle='-.', linewidth=2, label=f"Median: {mc_result['delta15N_sed_median']:.2f}‰")
        ax.axvspan(mc_result['delta15N_sed_ci68'][0], mc_result['delta15N_sed_ci68'][1],
                  alpha=0.2, color='yellow', label='68% CI')
        ax.set_xlabel('δ¹⁵N_sed (‰)')
        ax.set_ylabel('Probability Density')
        ax.set_title(f'Monte Carlo Analysis (f_assimilator={f_assimilator})')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.savefig(f"{PLOT_DIR}/example_3_monte_carlo.png", dpi=200)
        plt.close()
        print(f"\n图表已保存: {PLOT_DIR}/example_3_monte_carlo.png")
    
    return mc_result


def example_4_scenario_comparison():
    """
    示例 4: 不同情景对比
    """
    print("\n" + "=" * 70)
    print("示例 4: 不同地质情景对比")
    print("=" * 70)
    
    scenarios = [
        ('modern', '现代氧化海洋'),
        ('early_triassic', '早三叠世'),
        ('neoproterozoic', '新元古代')
    ]
    
    print(f"{'情景':<20} {'f=0.1':<12} {'f=0.3':<12} {'f=0.5':<12}")
    print("-" * 60)
    
    results = []
    for scenario_key, scenario_name in scenarios:
        n_system = NIsotopeSystem(scenario=scenario_key)
        
        delta_at_01 = n_system.forward_model(f_assimilator=0.1)
        delta_at_03 = n_system.forward_model(f_assimilator=0.3)
        delta_at_05 = n_system.forward_model(f_assimilator=0.5)
        
        print(f"{scenario_name:<20} {delta_at_01:<12.2f} "
              f"{delta_at_03:<12.2f} {delta_at_05:<12.2f}")
        results.append((scenario_name, delta_at_01, delta_at_03, delta_at_05))
    
    return results


def example_5_curve_calculation():
    """
    示例 5: 计算 f_assimilator 与 δ¹⁵N_sed 的关系曲线
    """
    print("\n" + "=" * 70)
    print("示例 5: 关系曲线计算 (用于绘图)")
    print("=" * 70)
    
    n_system = NIsotopeSystem(scenario='modern')
    
    # 计算曲线
    curve = n_system.calculate_f_assimilator_curve(
        f_range=(0.0, 1.0), n_points=50, n_monte_carlo=500
    )
    
    print(f"计算了 {len(curve['f_assimilator'])} 个点")
    print(f"f_assimilator 范围: [{curve['f_assimilator'][0]:.2f}, "
          f"{curve['f_assimilator'][-1]:.2f}]")
    print(f"δ¹⁵N_sed 范围: [{curve['delta15N_sed_mean'].min():.2f}, "
          f"{curve['delta15N_sed_mean'].max():.2f}]‰")
    
    # 找到最大值点
    max_idx = np.argmax(curve['delta15N_sed_mean'])
    print(f"\nδ¹⁵N_sed 最大值 {curve['delta15N_sed_mean'][max_idx]:.2f}‰ "
          f"出现在 f = {curve['f_assimilator'][max_idx]:.2f}")
    
    # 可选：绘图
    if SAVE_PLOTS:
        os.makedirs(PLOT_DIR, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(curve['f_assimilator'], curve['delta15N_sed_mean'], 
               'b-', linewidth=2, label='Mean')
        ax.fill_between(curve['f_assimilator'],
                       curve['delta15N_sed_ci95_lower'],
                       curve['delta15N_sed_ci95_upper'],
                       alpha=0.2, color='blue', label='95% CI')
        ax.fill_between(curve['f_assimilator'],
                       curve['delta15N_sed_ci68_lower'],
                       curve['delta15N_sed_ci68_upper'],
                       alpha=0.3, color='blue', label='68% CI')
        ax.plot(curve['f_assimilator'][max_idx], 
               curve['delta15N_sed_mean'][max_idx],
               'r*', markersize=15, label=f'Peak (f={curve["f_assimilator"][max_idx]:.2f})')
        ax.set_xlabel('f_assimilator')
        ax.set_ylabel('δ¹⁵N_sed (‰)')
        ax.set_title('Nitrogen Isotope Forward Model Curve')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.savefig(f"{PLOT_DIR}/example_5_curve.png", dpi=200)
        plt.close()
        print(f"\n图表已保存: {PLOT_DIR}/example_5_curve.png")
    
    return curve


def example_6_temporal_analysis():
    """
    示例 6: 时间序列分析 (模拟早三叠世数据)
    """
    print("\n" + "=" * 70)
    print("示例 6: 时间序列分析 (Bootstrap + LOWESS)")
    print("=" * 70)
    
    # 模拟早三叠世氮同位素时间序列数据
    np.random.seed(42)
    
    # 年龄 (Ma)
    ages = np.array([251.9, 251.5, 251.2, 250.8, 250.5, 250.2, 249.8, 
                     249.5, 249.2, 248.8, 248.5, 248.2, 247.8, 247.5])
    
    # 模拟 δ¹⁵N 数据 (加入噪声)
    # 模拟从低值(~1‰)升高到高值(~4‰)再降低的趋势
    true_trend = 1.0 + 3.0 * np.exp(-((252-ages)**2) / 2.0)
    noise = np.random.normal(0, 0.3, len(ages))
    delta15N_data = true_trend + noise
    
    print(f"模拟了 {len(ages)} 个数据点")
    print(f"年龄范围: {ages.min():.1f} - {ages.max():.1f} Ma")
    print(f"δ¹⁵N 范围: {delta15N_data.min():.2f} - {delta15N_data.max():.2f}‰")
    
    # Bootstrap 均值估计
    print("\nBootstrap 分析...")
    bootstrap_result = Bootstrap.confidence_interval(
        delta15N_data, statistic_func=np.mean, n_bootstrap=5000
    )
    print(f"  均值: {bootstrap_result['statistic']:.2f}‰")
    print(f"  95% 置信区间: [{bootstrap_result['ci_lower']:.2f}, "
          f"{bootstrap_result['ci_upper']:.2f}]‰")
    
    # LOWESS 平滑
    print("\nLOWESS 平滑...")
    lowess_result = LOWESS.fit(ages, delta15N_data, frac=0.4)
    print(f"  残差标准差: {np.std(lowess_result['residuals']):.2f}‰")
    
    # 反演硝酸盐占比
    print("\n反演硝酸盐占比...")
    n_system = NIsotopeSystem(scenario='early_triassic')
    
    f_values = []
    for delta in delta15N_data:
        result = n_system.inverse_model(delta, f_range=(0.0, 0.48))
        f_values.append(result['f_assimilator'])
    
    f_values = np.array(f_values)
    print(f"  f_assimilator 范围: {f_values.min():.3f} - {f_values.max():.3f}")
    print(f"  平均 f_assimilator: {np.mean(f_values):.3f}")
    
    # 可选：绘图
    if SAVE_PLOTS:
        os.makedirs(PLOT_DIR, exist_ok=True)
        fig, axes = plt.subplots(2, 1, figsize=(10, 8))
        
        # 上图: δ¹⁵N 时间序列
        ax1 = axes[0]
        ax1.scatter(ages, delta15N_data, c='blue', s=50, label='Observed')
        ax1.plot(ages, lowess_result['y_smoothed'], 'r-', linewidth=2, label='LOWESS')
        ax1.fill_between(ages, delta15N_data - 0.5, delta15N_data + 0.5,
                        alpha=0.2, color='gray', label='Uncertainty')
        ax1.set_ylabel('δ¹⁵N_sed (‰)')
        ax1.set_title('Simulated Early Triassic Nitrogen Isotope Record')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.invert_xaxis()
        
        # 下图: 反演的 f_assimilator
        ax2 = axes[1]
        ax2.scatter(ages, f_values, c='green', s=50)
        ax2.plot(ages, f_values, 'g-', linewidth=1.5)
        ax2.axhspan(0, 0.1, alpha=0.2, color='red', label='Stage I (Anoxic)')
        ax2.axhspan(0.15, 0.25, alpha=0.2, color='green', label='Stage II (Oxic)')
        ax2.set_xlabel('Age (Ma)')
        ax2.set_ylabel('f_assimilator')
        ax2.set_title('Inferred Nitrate Availability')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.invert_xaxis()
        
        plt.tight_layout()
        plt.savefig(f"{PLOT_DIR}/example_6_temporal.png", dpi=200)
        plt.close()
        print(f"\n图表已保存: {PLOT_DIR}/example_6_temporal.png")
    
    return ages, delta15N_data, f_values


def example_7_changepoint_detection():
    """
    示例 7: 变点检测 (模拟新元古代 800Ma 变化)
    """
    print("\n" + "=" * 70)
    print("示例 7: 变点检测")
    print("=" * 70)
    
    # 模拟新元古代 δ¹⁵N 数据 (1000-700 Ma)
    # 在 ~800 Ma 处有一个阶跃上升
    np.random.seed(123)
    
    ages = np.arange(1000, 700, -10)  # 每 10 Ma 一个点
    
    # 生成带阶跃的数据
    delta15N = np.where(ages > 800, 
                       np.random.normal(1.5, 0.8, len(ages)),  # 800 Ma 前
                       np.random.normal(4.0, 1.0, len(ages)))  # 800 Ma 后
    
    print(f"时间序列: {ages.max()} - {ages.min()} Ma, {len(ages)} 个点")
    print(f"800 Ma 前平均 δ¹⁵N: {np.mean(delta15N[ages > 800]):.2f}‰")
    print(f"800 Ma 后平均 δ¹⁵N: {np.mean(delta15N[ages <= 800]):.2f}‰")
    
    # 变点检测
    print("\n变点检测 (AMOC)...")
    cp_result = ChangepointDetector.detect_mean_variance(delta15N)
    
    detected_age = ages[cp_result['changepoint']]
    print(f"  检测到变点位置: {cp_result['changepoint']} (年龄: {detected_age:.0f} Ma)")
    print(f"  变点前均值: {cp_result['before_mean']:.2f}‰")
    print(f"  变点后均值: {cp_result['after_mean']:.2f}‰")
    print(f"  似然比统计量: {cp_result['likelihood']:.2f}")
    
    # 可选：绘图
    if SAVE_PLOTS:
        os.makedirs(PLOT_DIR, exist_ok=True)
        fig, ax = plt.subplots(figsize=(10, 5))
        
        colors = ['blue' if a > 800 else 'red' for a in ages]
        ax.scatter(ages, delta15N, c=colors, s=50)
        ax.axvline(800, color='green', linestyle='--', linewidth=2, 
                  label='Expected Change (800 Ma)')
        ax.axvline(detected_age, color='orange', linestyle='-.', linewidth=2,
                  label=f'Detected Change ({detected_age:.0f} Ma)')
        ax.axhline(cp_result['before_mean'], xmin=0, xmax=0.5, 
                  color='blue', linestyle=':', alpha=0.7)
        ax.axhline(cp_result['after_mean'], xmin=0.5, xmax=1,
                  color='red', linestyle=':', alpha=0.7)
        
        ax.set_xlabel('Age (Ma)')
        ax.set_ylabel('δ¹⁵N_sed (‰)')
        ax.set_title('Changepoint Detection: Neoproterozoic Nitrogen Isotope Shift')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.invert_xaxis()
        
        plt.tight_layout()
        plt.savefig(f"{PLOT_DIR}/example_7_changepoint.png", dpi=200)
        plt.close()
        print(f"\n图表已保存: {PLOT_DIR}/example_7_changepoint.png")
    
    return ages, delta15N, cp_result


def main():
    """运行所有示例"""
    print("\n" + "=" * 70)
    print("氮同位素体系使用示例")
    print("基于 Kang et al. (2023) 和 Ma et al. (2025)")
    print("=" * 70)
    
    if SAVE_PLOTS:
        print(f"\n图表保存功能已启用，输出目录: {PLOT_DIR}/")
        print("如需禁用图表保存，请设置 SAVE_PLOTS = False")
    else:
        print("\n图表保存功能已禁用。")
        print("如需启用图表保存，请设置 SAVE_PLOTS = True")
    
    try:
        example_1_forward_model()
        example_2_inverse_model()
        example_3_monte_carlo()
        example_4_scenario_comparison()
        example_5_curve_calculation()
        example_6_temporal_analysis()
        example_7_changepoint_detection()
        
        print("\n" + "=" * 70)
        print("所有示例运行完成!")
        if SAVE_PLOTS:
            print(f"图表已保存到: {PLOT_DIR}/")
        print("=" * 70)
        
    except Exception as e:
        print(f"\n错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
