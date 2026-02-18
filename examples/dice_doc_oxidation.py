#!/usr/bin/env python3
"""
DICE事件DOC氧化与氧化剂消耗分析
基于 Li et al. 2020, Precambrian Research

本示例展示如何使用新架构的C同位素体系来复现原Carbon_Modeling/Oxi_Est.py的功能
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import matplotlib.pyplot as plt

from systems.c import CIsotopeSystem, OXIDANT_SCENARIOS


def reproduce_li_2020_figure7():
    """
    复现 Li et al. 2020 论文中的图7
    展示不同氧化剂情景下，DOC氧化与碳同位素负漂的关系
    """
    print("="*80)
    print("DICE事件DOC氧化分析")
    print("基于: Li et al. 2020, Precambrian Research")
    print("="*80)
    
    # 创建C同位素体系（DICE情景）
    c_system = CIsotopeSystem(scenario='dice')
    
    # 运行完整模型
    print("\n[1] 计算不同DOC通量下的碳同位素偏移...")
    result = c_system.doc_excursion_model(
        F_odoc_range=(0, 10e18),
        n_points=300
    )
    
    if not result.success:
        print(f"模型运行失败: {result.message}")
        return
    
    F_odoc = result.get('F_odoc')
    delta_delta13C = result.get('delta_delta13C')
    oxidant_scenarios = result.get('oxidant_scenarios')
    
    # 关键验证点
    print("\n[2] 模型验证（与论文对比）:")
    print("-"*60)
    
    # 找到2‰和4‰偏移对应的DOC通量
    idx_2 = np.argmin(np.abs(delta_delta13C + 2.0))
    idx_4 = np.argmin(np.abs(delta_delta13C + 4.0))
    
    F_odoc_at_2 = F_odoc[idx_2] / 1e18
    F_odoc_at_4 = F_odoc[idx_4] / 1e18
    
    print(f"产生 ~2‰ 偏移需要 F_odoc ≈ {F_odoc_at_2:.2f}×10¹⁸ mol/Ma (论文: 2.1×10¹⁸)")
    print(f"产生 ~4‰ 偏移需要 F_odoc ≈ {F_odoc_at_4:.2f}×10¹⁸ mol/Ma (论文: 4.6×10¹⁸)")
    
    # 氧化剂消耗计算
    print("\n[3] 氧化剂消耗计算（论文图7对应范围 -4‰到-2‰）:")
    print("-"*60)
    
    for scenario_key, scenario_data in oxidant_scenarios.items():
        consumption = scenario_data['consumption']
        
        # 找到-4‰到-2‰范围内的消耗
        mask = (delta_delta13C >= -4.0) & (delta_delta13C <= -2.0)
        if np.any(mask):
            F_min = np.min(consumption[mask]) / 1e18
            F_max = np.max(consumption[mask]) / 1e18
            print(f"  {scenario_data['name']:30}: {F_min:.2f} - {F_max:.2f} ×10¹⁸ mol/Ma")
    
    # 绘图
    print("\n[4] 生成图表...")
    plot_results(F_odoc, delta_delta13C, oxidant_scenarios)
    
    print("\n" + "="*80)
    print("分析完成！图表已保存: DICE_Figure7_Reproduction.png")
    print("="*80)


def plot_results(F_odoc, delta_delta13C, oxidant_scenarios):
    """
    绘制结果图表（复现论文图7）
    """
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['axes.unicode_minus'] = False
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 绘制各情景曲线
    for scenario_key, scenario_data in oxidant_scenarios.items():
        ax.plot(
            delta_delta13C,
            scenario_data['consumption'] / 1e18,
            color=scenario_data['color'],
            linestyle=scenario_data['linestyle'],
            linewidth=2,
            label=scenario_data['name']
        )
    
    # 标注DICE偏移范围
    mask = (delta_delta13C >= -4.0) & (delta_delta13C <= -2.0)
    if np.any(mask):
        ax.fill_between(
            delta_delta13C[mask],
            0, 8,
            color="yellow",
            alpha=0.3,
            label="Global DICE Excursion (2‰~4‰)"
        )
    
    # 添加标注
    ax.text(-3.0, 7.5, "2‰~4‰ δ¹³C Excursion", 
            fontsize=11, ha="center", fontweight="bold")
    ax.annotate('', xy=(-2.2, 7.2), xytext=(-3.8, 7.2),
                arrowprops=dict(arrowstyle='->', color='black'))
    ax.annotate('', xy=(-3.8, 7.2), xytext=(-2.2, 7.2),
                arrowprops=dict(arrowstyle='->', color='black'))
    
    # 坐标轴设置
    ax.set_xlabel(r"$\delta^{13}$C Excursion (‰)", fontsize=12)
    ax.set_ylabel(r"Oxidant Consumption Flux ($\times10^{18}$ mol/Ma)", fontsize=12)
    ax.set_ylim(0, 8)
    ax.set_xlim(-6, 0)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)
    
    plt.tight_layout()
    plt.savefig("DICE_Figure7_Reproduction.png", dpi=300, bbox_inches="tight")
    plt.show()


def detailed_scenario_analysis():
    """
    详细分析特定氧化剂情景
    """
    print("\n" + "="*80)
    print("详细情景分析")
    print("="*80)
    
    c_system = CIsotopeSystem(scenario='dice')
    
    scenarios_to_analyze = [
        'modern_o2_high_sulfate',
        'low_o2_high_sulfate',
        'complete_oxidation'
    ]
    
    target_excursions = [-2, -3, -4]
    
    for scenario_key in scenarios_to_analyze:
        print(f"\n情景: {OXIDANT_SCENARIOS[scenario_key]['name']}")
        print(f"描述: {OXIDANT_SCENARIOS[scenario_key]['description']}")
        print("-"*60)
        
        for target in target_excursions:
            result = c_system.find_doc_for_excursion(target)
            F_odoc = result['F_odoc']
            
            ox = c_system.calculate_oxidant_consumption(F_odoc, scenario_key)
            
            print(f"  目标偏移: {target:+.0f}‰")
            print(f"    所需 F_odoc: {F_odoc/1e18:.2f}×10¹⁸ mol/Ma")
            print(f"    O₂消耗: {ox.o2_consumption/1e18:.2f}×10¹⁸ mol/Ma")
            print(f"    硫酸盐消耗: {ox.sulfate_consumption/1e18:.2f}×10¹⁸ mol/Ma")
            print(f"    总氧化剂: {ox.total_consumption/1e18:.2f}×10¹⁸ mol/Ma")


def compare_architectures():
    """
    对比新旧架构的实现方式
    """
    print("\n" + "="*80)
    print("新旧架构对比")
    print("="*80)
    
    print("""
旧架构 (Carbon_Modeling/Oxi_Est.py):
  - 所有代码在一个文件中
  - 参数硬编码
  - 难以复用和扩展
  - 无法方便地切换情景

新架构 (本示例):
  - 使用 systems.c.CIsotopeSystem 类
  - 参数通过 parameters.py 配置
  - 易于扩展到其他同位素体系
  - 支持多种氧化剂情景的灵活切换
  - 可与其他体系（如Mg、S）组合分析
""")
    
    # 展示新架构的灵活性
    print("\n新架构优势演示:")
    print("-"*60)
    
    # 1. 轻松切换情景
    print("\n1. 轻松切换情景:")
    for scenario in ['dice', 'modern']:
        c = CIsotopeSystem(scenario=scenario)
        result = c.solve_steady_state(F_odoc=3e18)
        print(f"   {scenario:10}: δ¹³C_carb = {result.get('delta13C_carb'):+.2f}‰")
    
    # 2. 访问基础公式
    print("\n2. 直接使用底层公式:")
    from toolkit.isotope.formulas import MassBalance
    delta_mix = MassBalance.two_component_mixing(-4, -30, 0.7)
    print(f"   混合计算: -4‰和-30‰按7:3混合 = {delta_mix:.2f}‰")
    
    # 3. 访问物理常数
    print("\n3. 访问物理常数:")
    from toolkit.physics.constants import ReservoirConstants
    print(f"   现代海水DIC储量: {ReservoirConstants.SEAWATER_DIC_MOL:.2e} mol")


if __name__ == '__main__':
    # 主分析：复现论文图7
    reproduce_li_2020_figure7()
    
    # 详细情景分析
    detailed_scenario_analysis()
    
    # 架构对比
    compare_architectures()
