#!/usr/bin/env python3
"""
Sr同位素模拟示例
基于Wang et al. (2021) Earth-Science Reviews 218: 103679

演示内容：
1. 基本质量平衡计算
2. 蒙特卡洛随机模拟
3. 不同情景分析
4. 参数敏感性分析
"""

import sys
from pathlib import Path

# 添加项目根目录到Python路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import numpy as np


def demo_basic_calculation():
    """演示基本质量平衡计算"""
    print("\n" + "="*70)
    print("示例1: 基本质量平衡计算")
    print("="*70)
    
    from systems.sr import SrIsotopeSystem, SrFluxConfig
    
    # 创建Sr同位素体系
    sr = SrIsotopeSystem(scenario='modern')
    
    # 示例1a: 使用现代参数
    print("\n[1a] 现代参数下的海水⁸⁷Sr/⁸⁶Sr:")
    modern_ratio = sr.calculate_seawater_ratio()
    print(f"  计算值: {modern_ratio:.5f}")
    print(f"  观测值: 0.70917")
    
    # 示例1b: 改变河流输入
    print("\n[1b] 不同河流同位素组成的影响:")
    river_ratios = [0.708, 0.711, 0.714, 0.720]
    for r_riv in river_ratios:
        ratio = sr.calculate_seawater_ratio(R_river=r_riv)
        print(f"  R_river = {r_riv:.3f}: R_seawater = {ratio:.5f}")
    
    # 示例1c: 改变热液输入
    print("\n[1c] 高温热液通量变化的影响:")
    f_highT_values = [2e9, 8e9, 20e9, 35e9]
    for f in f_highT_values:
        ratio = sr.calculate_seawater_ratio(F_hydrothermal_highT=f)
        change = (ratio - modern_ratio) * 10000  # 转换为ppm
        print(f"  F_highT = {f/1e9:.0f}×10⁹ mol/yr: "
              f"R_seawater = {ratio:.5f} ({change:+.1f} ppm)")


def demo_monte_carlo():
    """演示蒙特卡洛随机模拟"""
    print("\n" + "="*70)
    print("示例2: 蒙特卡洛随机模拟")
    print("="*70)
    print("基于 Wang et al. (2021) 的随机海洋箱模型")
    
    from systems.sr import SrIsotopeSystem
    
    # 创建二叠纪情景的Sr体系
    sr = SrIsotopeSystem(scenario='permian')
    
    # 运行蒙特卡洛模拟（使用较少数量的运行以加快演示）
    print("\n运行蒙特卡洛模拟 (n_runs=2000)...")
    result = sr.monte_carlo_simulation(
        time_span=(299, 252),  # 二叠纪
        n_runs=2000,
        n_time_points=30,
        verbose=False
    )
    
    if result.success:
        print(f"\n模拟结果:")
        print(f"  成功率: {result.data['success_rate']:.2%}")
        print(f"  成功运行次数: {result.data['n_successful']}")
        
        print(f"\n关键参数统计:")
        for param_name in ['F_river', 'R_river', 'F_hydrothermal_highT']:
            if param_name in result.statistics:
                stats = result.statistics[param_name]
                print(f"  {param_name}:")
                print(f"    平均值: {stats['mean']:.3e}")
                print(f"    中位数: {stats['median']:.3e}")
                print(f"    范围: [{stats['min']:.3e}, {stats['max']:.3e}]")
        
        print(f"\n海水Sr同位素范围:")
        mean_ratio = result.data['mean_ratio']
        print(f"  平均值范围: [{mean_ratio.min():.5f}, {mean_ratio.max():.5f}]")
        print(f"  2.5%分位数: {result.data['percentile_2.5'][-1]:.5f}")
        print(f"  97.5%分位数: {result.data['percentile_97.5'][-1]:.5f}")


def demo_scenarios():
    """演示不同情景分析"""
    print("\n" + "="*70)
    print("示例3: 情景分析")
    print("="*70)
    print("测试Wang et al. (2021)论文中的不同约束情景")
    
    from systems.sr import SrIsotopeSystem
    
    sr = SrIsotopeSystem(scenario='permian')
    
    # 运行几个关键情景（使用较少的运行次数）
    scenarios = ['scenario1', 'scenario2', 'scenario7']
    
    for scenario_name in scenarios:
        print(f"\n{'-'*50}")
        result = sr.scenario_analysis(
            scenario_name=scenario_name,
            n_runs=1000,
            filter_by_data=False,
            verbose=False
        )
        
        if result.success:
            print(f"\n情景 {scenario_name} 结果:")
            print(f"  成功运行次数: {result.data['n_successful']}")
            
            # 显示该情景下的关键参数
            if 'F_river' in result.statistics:
                f_riv_stats = result.statistics['F_river']
                print(f"  河流通量: {f_riv_stats['mean']:.2e} ± "
                      f"{f_riv_stats['std']:.2e} mol/yr")
            
            if 'R_river' in result.statistics:
                r_riv_stats = result.statistics['R_river']
                print(f"  河流同位素: {r_riv_stats['mean']:.5f} ± "
                      f"{r_riv_stats['std']:.5f}")


def demo_sensitivity():
    """演示参数敏感性分析"""
    print("\n" + "="*70)
    print("示例4: 参数敏感性分析")
    print("="*70)
    
    from systems.sr import SrIsotopeSystem
    
    sr = SrIsotopeSystem(scenario='modern')
    
    # 分析河流同位素组成的敏感性
    print("\n[4a] 河流同位素组成 (R_river) 的敏感性:")
    result = sr.sensitivity_analysis(
        param_name='R_river',
        param_range=(0.705, 0.725),
        n_steps=20
    )
    
    if result.success:
        r_vals = result.data['param_values']
        sr_vals = result.data['sr_ratios']
        sens = result.data['sensitivity']
        
        print(f"  R_river范围: [{r_vals.min():.3f}, {r_vals.max():.3f}]")
        print(f"  R_seawater范围: [{sr_vals.min():.5f}, {sr_vals.max():.5f}]")
        print(f"  平均敏感性: {np.mean(np.abs(sens)):.4f} (ΔR_sw/ΔR_riv)")
        
        # 显示几个关键点
        print(f"\n  关键点:")
        for i in [0, len(r_vals)//2, -1]:
            print(f"    R_river={r_vals[i]:.3f}: R_sw={sr_vals[i]:.5f}")
    
    # 分析高温热液通量的敏感性
    print("\n[4b] 高温热液通量 (F_highT) 的敏感性:")
    result = sr.sensitivity_analysis(
        param_name='F_hydrothermal_highT',
        param_range=(2e9, 35e9),
        n_steps=20
    )
    
    if result.success:
        f_vals = result.data['param_values']
        sr_vals = result.data['sr_ratios']
        
        print(f"  F_highT范围: [{f_vals.min()/1e9:.0f}, {f_vals.max()/1e9:.0f}] ×10⁹ mol/yr")
        print(f"  R_seawater范围: [{sr_vals.min():.5f}, {sr_vals.max():.5f}]")
        
        # 计算变化幅度
        sr_change = (sr_vals.max() - sr_vals.min()) * 10000  # ppm
        print(f"  Sr同位素变化幅度: {sr_change:.1f} ppm")


def demo_permian_evolution():
    """演示二叠纪Sr同位素演化"""
    print("\n" + "="*70)
    print("示例5: 二叠纪Sr同位素演化趋势")
    print("="*70)
    
    from systems.sr import SrIsotopeSystem, SrFluxConfig
    
    sr = SrIsotopeSystem(scenario='permian')
    
    # 模拟二叠纪Sr同位素演化
    # 基于Wang et al. (2021)的观察
    print("\n二叠纪Sr同位素演化阶段:")
    
    ages = [299, 273, 265, 259, 252]  # Ma
    stages = [
        "Asselian (299 Ma) - 早二叠世",
        "Kungurian (273 Ma) - 早二叠世末",
        "Wordian (265 Ma) - 中二叠世",
        "Capitanian (259 Ma) - 中二叠世末",
        "Changhsingian (252 Ma) - 晚二叠世"
    ]
    
    # 典型观测值
    observed_ratios = [0.70827, 0.7070, 0.70683, 0.7069, 0.70717]
    
    print(f"\n  {'时期':<40} {'观测⁸⁷Sr/⁸⁶Sr':<15} {'趋势':<10}")
    print("  " + "-"*65)
    
    for i, (age, stage, ratio) in enumerate(zip(ages, stages, observed_ratios)):
        if i > 0:
            change = ratio - observed_ratios[i-1]
            trend = "↓ 下降" if change < 0 else "↑ 上升"
        else:
            trend = "→ 开始"
        print(f"  {stage:<40} {ratio:<15.5f} {trend:<10}")
    
    print("\n  趋势解释 (基于Wang et al., 2021):")
    print("  1. 早-中二叠世持续下降: 可能反映热液活动增强")
    print("  2. Capitanian达到最低值: 最低⁸⁷Sr/⁸⁶Sr ~0.70683")
    print("  3. 晚二叠世上升: 可能与大陆风化增强有关")


def demo_flux_contributions():
    """演示不同端元的贡献"""
    print("\n" + "="*70)
    print("示例6: 各端元对海水Sr的贡献")
    print("="*70)
    
    from systems.sr import SrIsotopeSystem
    
    sr = SrIsotopeSystem(scenario='modern')
    
    # 现代各端元参数
    fluxes = {
        '河流': (47.6e9, 0.71107),
        '高温热液': (8.04e9, 0.7037),
        '低温热液': (10e9, 0.7084),
        '成岩作用': (3.4e9, 0.7084),
    }
    
    print("\n现代各端元参数:")
    print(f"{'端元':<15} {'通量 (10⁹ mol/yr)':<20} {'⁸⁷Sr/⁸⁶Sr':<15} {'通量占比':<10}")
    print("-"*70)
    
    total_flux = sum(f for f, _ in fluxes.values())
    
    for name, (flux, ratio) in fluxes.items():
        percentage = flux / total_flux * 100
        print(f"{name:<15} {flux/1e9:<20.2f} {ratio:<15.5f} {percentage:<10.1f}%")
    
    print(f"{'总计':<15} {total_flux/1e9:<20.2f}")
    
    # 计算混合结果
    print("\n质量平衡计算:")
    weighted_sum = sum(f * r for f, r in fluxes.values())
    seawater_ratio = weighted_sum / total_flux
    print(f"  加权平均: {seawater_ratio:.5f}")
    print(f"  现代海水观测: 0.70917")
    
    # 各端元对Sr同位素的相对贡献
    print("\n各端元对海水Sr同位素的贡献 (以通量×比值计算):")
    contributions = {}
    for name, (flux, ratio) in fluxes.items():
        contrib = flux * ratio / weighted_sum * 100
        contributions[name] = contrib
        print(f"  {name}: {contrib:.1f}%")


def main():
    """主函数"""
    print("\n" + "#"*70)
    print("#" + " "*68 + "#")
    print("#" + " Sr同位素海洋箱模型演示".center(68) + "#")
    print("#" + " 基于Wang et al. (2021) Earth-Science Reviews".center(68) + "#")
    print("#" + " "*68 + "#")
    print("#"*70)
    
    # 运行所有演示
    demo_basic_calculation()
    demo_monte_carlo()
    demo_scenarios()
    demo_sensitivity()
    demo_permian_evolution()
    demo_flux_contributions()
    
    print("\n" + "#"*70)
    print("演示完成!")
    print("#"*70 + "\n")


if __name__ == "__main__":
    main()
