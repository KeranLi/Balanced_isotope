#!/usr/bin/env python3
"""
铀同位素模型 - 不确定度分析示例

演示如何使用蒙特卡洛、Bootstrap和敏感性分析评估模型结果的不确定度
"""

import numpy as np
import sys
from pathlib import Path

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from systems.u import UIsotopeSystem, UncertaintyAnalyzer, ParameterUncertainty


def example_1_monte_carlo():
    """
    示例1: 蒙特卡洛不确定度分析
    
    评估在考虑所有参数不确定度的情况下，f_anox的不确定度范围
    """
    print("=" * 70)
    print("示例1: 蒙特卡洛不确定度分析")
    print("=" * 70)
    
    # 创建铀同位素体系
    u_system = UIsotopeSystem(scenario='modern')
    analyzer = UncertaintyAnalyzer(u_system)
    
    # 观测数据
    delta238_carb = -0.65  # ‰
    
    print(f"\n观测数据: δ²³⁸U_carb = {delta238_carb:.2f}‰")
    print("\n参数不确定度设置:")
    print("  - Δ_sw-anox: 0.77 ± 0.04 ‰ (正态分布)")
    print("  - Δ_diag: 0.40 (范围 0.30-0.50) ‰ (均匀分布)")
    print("  - δ_river: -0.29 ± 0.16 ‰ (正态分布)")
    print("  - 测量误差: ±0.05 ‰ (正态分布)")
    
    # 运行蒙特卡洛模拟
    print("\n运行蒙特卡洛模拟 (n=50,000)...")
    result = analyzer.monte_carlo_steady_state(
        delta238_carb=delta238_carb,
        n_samples=50000,
        confidence_level=0.95,
        random_seed=42  # 为了结果可复现
    )
    
    print("\n" + "-" * 50)
    print("结果:")
    print("-" * 50)
    print(f"  f_anox 均值:      {result['f_anox_mean']:.1%}")
    print(f"  f_anox 中位数:    {result['f_anox_median']:.1%}")
    print(f"  标准差:           {result['f_anox_std']:.1%}")
    print(f"  68% 置信区间:     [{result['f_anox_ci'][0]:.1%}, {result['f_anox_ci'][1]:.1%}]")
    print(f"  95% 置信区间:     [{result['f_anox_ci'][0]:.1%}, {result['f_anox_ci'][1]:.1%}]")
    
    # 收敛性诊断
    conv = result['convergence']
    print(f"\n  收敛性诊断:")
    print(f"    R̂ (Gelman-Rubin统计量): {conv['r_hat']:.3f}")
    print(f"    有效样本数: {conv['n_effective']}")
    print(f"    收敛状态: {'✓ 已收敛' if conv['converged'] else '⚠ 未收敛'}")
    
    # 解释
    print("\n解释:")
    print(f"  - 最可能的 f_anox 值约为 {result['f_anox_median']:.0%}")
    print(f"  - 有95%的置信度认为 f_anox 在 [{result['f_anox_ci'][0]:.0%}, {result['f_anox_ci'][1]:.0%}] 范围内")
    print(f"  - 相对不确定度: {result['f_anox_std']/result['f_anox_mean']*100:.1f}%")


def example_2_bootstrap():
    """
    示例2: Bootstrap不确定度分析
    
    适用于有多个重复测量样本的情况
    """
    print("\n" + "=" * 70)
    print("示例2: Bootstrap不确定度分析")
    print("=" * 70)
    
    u_system = UIsotopeSystem(scenario='modern')
    analyzer = UncertaintyAnalyzer(u_system)
    
    # 模拟地层中的多个测量值（同一层位多个样品）
    # 例如：某个缺氧事件层位的多个碳酸盐样品
    delta_measurements = np.array([
        -0.62, -0.58, -0.65, -0.60, -0.63, -0.59, -0.61, -0.64
    ])
    
    print(f"\n测量数据 (n={len(delta_measurements)}):")
    print(f"  δ²³⁸U = {delta_measurements}")
    print(f"  均值 = {np.mean(delta_measurements):.2f}‰")
    print(f"  标准差 = {np.std(delta_measurements, ddof=1):.2f}‰")
    
    # 运行Bootstrap分析
    print("\n运行Bootstrap重采样 (n=10,000)...")
    result = analyzer.bootstrap_analysis(
        delta238_measurements=delta_measurements,
        n_bootstrap=10000,
        random_seed=42
    )
    
    print("\n" + "-" * 50)
    print("结果:")
    print("-" * 50)
    print(f"  原始估计值:       {result['original_estimate']:.1%}")
    print(f"  Bootstrap均值:    {result['f_anox_mean']:.1%}")
    print(f"  Bootstrap标准差:  {result['f_anox_std']:.1%}")
    print(f"  95% 置信区间:     [{result['f_anox_ci'][0]:.1%}, {result['f_anox_ci'][1]:.1%}]")
    print(f"  偏差:             {result['bias']:.1%}")


def example_3_sensitivity():
    """
    示例3: 敏感性分析
    
    识别对f_anox计算影响最大的参数
    """
    print("\n" + "=" * 70)
    print("示例3: 敏感性分析 (Tornado分析)")
    print("=" * 70)
    
    u_system = UIsotopeSystem(scenario='modern')
    analyzer = UncertaintyAnalyzer(u_system)
    
    delta238_carb = -0.50  # ‰
    
    print(f"\n输入: δ²³⁸U_carb = {delta238_carb:.2f}‰")
    
    # 运行敏感性分析
    result = analyzer.sensitivity_analysis(
        delta238_carb=delta238_carb,
        n_points=50
    )
    
    print(f"\n基准 f_anox: {result['baseline_f_anox']:.1%}")
    print("\n" + "-" * 50)
    print("参数敏感性排名 (按影响大小):")
    print("-" * 50)
    
    for i, item in enumerate(result['tornado_data'], 1):
        param = item['parameter']
        min_eff = item['min_effect']
        max_eff = item['max_effect']
        range_val = item['range']
        
        print(f"{i}. {param:15}")
        print(f"   影响范围: [{min_eff:+.1%}, {max_eff:+.1%}]")
        print(f"   归一化敏感性: {range_val:.1%}")
    
    print("\n结论:")
    print(f"  最关键参数: {result['most_important']}")
    print("  建议: 为提高f_anox估计精度，应优先改进该参数的测定")


def example_4_custom_uncertainty():
    """
    示例4: 自定义参数不确定度
    
    根据具体研究情况调整参数不确定度
    """
    print("\n" + "=" * 70)
    print("示例4: 自定义参数不确定度")
    print("=" * 70)
    
    u_system = UIsotopeSystem(scenario='modern')
    analyzer = UncertaintyAnalyzer(u_system)
    
    # 自定义参数不确定度
    # 假设我们对自己的测量更自信，但成岩校正不太确定
    analyzer.config.delta_sw_anox = ParameterUncertainty(
        mean=0.77, std=0.04, distribution='normal'
    )
    analyzer.config.delta_diag = ParameterUncertainty(
        mean=0.40, lower=0.20, upper=0.60, distribution='uniform'
    )
    analyzer.config.delta_river = ParameterUncertainty(
        mean=-0.29, std=0.10, distribution='normal'  # 更精确的河流值
    )
    analyzer.config.delta_measurement = ParameterUncertainty(
        mean=0.0, std=0.02, distribution='normal'  # 高精度测量
    )
    
    delta238_carb = -0.55
    
    print(f"\n自定义不确定度设置:")
    print(f"  - δ_river: -0.29 ± 0.10 ‰ (更精确)")
    print(f"  - Δ_diag: 0.40 (范围 0.20-0.60) ‰ (更保守)")
    print(f"  - 测量误差: ±0.02 ‰ (高精度)")
    print(f"\n输入: δ²³⁸U_carb = {delta238_carb:.2f}‰")
    
    result = analyzer.monte_carlo_steady_state(
        delta238_carb=delta238_carb,
        n_samples=20000,
        random_seed=42
    )
    
    print("\n" + "-" * 50)
    print("结果:")
    print("-" * 50)
    print(f"  f_anox = {result['f_anox_mean']:.1%} ± {result['f_anox_std']:.1%}")
    print(f"  95% CI: [{result['f_anox_ci'][0]:.1%}, {result['f_anox_ci'][1]:.1%}]")


def example_5_measurement_uncertainty_only():
    """
    示例5: 仅考虑测量不确定度
    
    当其他参数确定时，仅评估测量误差的影响
    """
    print("\n" + "=" * 70)
    print("示例5: 仅测量不确定度评估")
    print("=" * 70)
    
    u_system = UIsotopeSystem(scenario='modern')
    analyzer = UncertaintyAnalyzer(u_system)
    
    delta238_carb = -0.50
    measurement_std = 0.10  # 较大的测量误差
    
    print(f"\n输入: δ²³⁸U_carb = {delta238_carb:.2f} ± {measurement_std:.2f}‰")
    print("假设: 所有模型参数完全确定")
    
    result = analyzer.measurement_uncertainty_only(
        delta238_carb=delta238_carb,
        measurement_std=measurement_std,
        n_samples=10000
    )
    
    print("\n" + "-" * 50)
    print("结果:")
    print("-" * 50)
    print(f"  f_anox = {result['f_anox_mean']:.1%} ± {result['f_anox_std']:.1%}")
    print(f"  95% CI: [{result['f_anox_ci'][0]:.1%}, {result['f_anox_ci'][1]:.1%}]")
    
    print("\n与全参数不确定度对比:")
    full_result = analyzer.monte_carlo_steady_state(
        delta238_carb=delta238_carb,
        n_samples=10000,
        random_seed=42
    )
    print(f"  仅测量不确定度: ±{result['f_anox_std']:.1%}")
    print(f"  全参数不确定度: ±{full_result['f_anox_std']:.1%}")
    print(f"  测量贡献: {(result['f_anox_std']/full_result['f_anox_std'])**2*100:.1f}%")


def example_6_likelihood_analysis():
    """
    示例6: 似然分析
    
    计算不同f_anox值的似然度
    """
    print("\n" + "=" * 70)
    print("示例6: 似然分析")
    print("=" * 70)
    
    from systems.u.uncertainty import analyze_likelihood
    
    u_system = UIsotopeSystem(scenario='modern')
    
    observed_delta = -0.60  # 观测值
    measurement_std = 0.05
    
    # f_anox 搜索网格
    f_anox_grid = np.linspace(0, 1, 101)
    
    print(f"\n观测: δ²³⁸U_carb = {observed_delta:.2f} ± {measurement_std:.2f}‰")
    print("计算不同 f_anox 值的似然度...")
    
    result = analyze_likelihood(
        u_system=u_system,
        observed_delta238=observed_delta,
        f_anox_grid=f_anox_grid,
        measurement_std=measurement_std
    )
    
    print("\n" + "-" * 50)
    print("结果:")
    print("-" * 50)
    print(f"  最大似然估计: f_anox = {result['f_anox_mle']:.1%}")
    print(f"  95% 可信区间: [{result['credible_interval'][0]:.1%}, "
          f"{result['credible_interval'][1]:.1%}]")
    
    # 显示似然分布的一些特征
    likelihood = result['likelihood']
    print(f"\n  似然分布特征:")
    print(f"    最大值: {np.max(likelihood):.4f}")
    print(f"    半高宽范围: {np.sum(likelihood > 0.5*np.max(likelihood))/len(likelihood)*100:.1f}% of grid")


def main():
    """运行所有示例"""
    print("\n" + "=" * 70)
    print("铀同位素模型 - 不确定度分析完整示例")
    print("=" * 70)
    
    example_1_monte_carlo()
    example_2_bootstrap()
    example_3_sensitivity()
    example_4_custom_uncertainty()
    example_5_measurement_uncertainty_only()
    example_6_likelihood_analysis()
    
    print("\n" + "=" * 70)
    print("所有示例运行完成!")
    print("=" * 70)
    print("\nCLI使用示例:")
    print("  python cli.py u --delta-carb -0.65 --steady-state --uncertainty mc")
    print("  python cli.py u --delta-carb -0.65 --steady-state --sensitivity-analysis")
    print()


if __name__ == '__main__':
    main()
