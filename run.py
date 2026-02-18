#!/usr/bin/env python3
"""
同位素质量平衡模型 - 统一入口
整合所有功能，展示新架构的使用
"""

import sys
from pathlib import Path

# 添加项目根目录到Python路径
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))


def run_mg_analysis(data_file=None):
    """
    Mg同位素分析示例
    """
    print("\n" + "="*80)
    print("Mg同位素风化分析")
    print("="*80)
    
    from systems.mg import MgIsotopeSystem
    
    # 创建Mg同位素体系
    mg = MgIsotopeSystem()
    
    # 示例1：风化比例计算
    print("\n[1] 风化端元比例计算")
    print("-"*60)
    
    test_samples = [
        (-1.5, "富集碳酸盐风化"),
        (-2.5, "混合风化"),
        (-0.5, "富集硅酸盐风化")
    ]
    
    for delta_sample, description in test_samples:
        ratios = mg.calculate_weathering_ratio(
            delta_sample=delta_sample,
            delta_seawater=-0.83
        )
        print(f"  δ²⁶Mg = {delta_sample:+.2f}‰ ({description}):")
        print(f"    碳酸盐风化: {ratios['f_carbonate']:.1%}")
        print(f"    硅酸盐风化: {ratios['f_silicate']:.1%}")
    
    # 示例2：海水演化模拟
    print("\n[2] 海水Mg同位素演化模拟")
    print("-"*60)
    
    for scenario in ['modern', 'high_weathering', 'low_weathering']:
        result = mg.seawater_evolution(
            time_span=(0, 100),
            initial_delta=-0.5,
            flux_scenario=scenario
        )
        
        if result.success:
            values = result.values
            delta_final = values[-1, 1]  # δ²⁶Mg
            print(f"  {scenario:20}: 100Ma后 δ²⁶Mg = {delta_final:+.3f}‰")
    
    # 示例3：如果有数据文件，进行完整分析
    if data_file and Path(data_file).exists():
        print(f"\n[3] 数据文件分析: {data_file}")
        print("-"*60)
        print("  (数据加载功能待实现)")
    
    print()


def run_c_analysis():
    """
    C同位素分析示例
    """
    print("\n" + "="*80)
    print("C同位素DOC氧化分析")
    print("="*80)
    
    from systems.c import CIsotopeSystem
    
    # 创建C同位素体系（DICE情景）
    c = CIsotopeSystem(scenario='dice')
    
    # 示例1：特定DOC通量的稳态
    print("\n[1] 特定DOC通量下的稳态同位素")
    print("-"*60)
    
    F_odoc_values = [1e18, 2.1e18, 4.6e18, 10e18]
    for F_odoc in F_odoc_values:
        result = c.solve_steady_state(F_odoc=F_odoc)
        if result.success:
            delta_carb = result.get('delta13C_carb')
            print(f"  F_doc = {F_odoc/1e18:5.2f}×10¹⁸: δ¹³C_carb = {delta_carb:+.2f}‰")
    
    # 示例2：氧化剂消耗
    print("\n[2] 氧化剂消耗计算 (F_doc = 4×10¹⁸ mol/Ma)")
    print("-"*60)
    
    F_odoc = 4e18
    scenarios = ['modern_o2_high_sulfate', 'low_o2_high_sulfate', 'complete_oxidation']
    
    for scenario_key in scenarios:
        ox = c.calculate_oxidant_consumption(F_odoc, scenario_key)
        print(f"  {scenario_key:30}:")
        print(f"    O₂消耗: {ox.o2_consumption/1e18:.2f}×10¹⁸ mol/Ma")
        print(f"    硫酸盐消耗: {ox.sulfate_consumption/1e18:.2f}×10¹⁸ mol/Ma")
        print(f"    总计: {ox.total_consumption/1e18:.2f}×10¹⁸ mol/Ma")
    
    # 示例3：反演计算
    print("\n[3] 从观测偏移反演DOC通量")
    print("-"*60)
    
    target_excursions = [-2, -4, -6]
    for target in target_excursions:
        result = c.find_doc_for_excursion(target)
        print(f"  目标偏移 {target:+.0f}‰: 需要 F_doc = {result['F_odoc']/1e18:.2f}×10¹⁸ mol/Ma")
    
    print()


def run_core_tools_demo():
    """
    核心工具演示
    """
    print("\n" + "="*80)
    print("核心工具演示")
    print("="*80)
    
    import numpy as np
    
    # 1. ODE求解
    print("\n[1] ODE求解器")
    print("-"*60)
    
    from toolkit.math.numerical import ODESolver
    
    # 指数衰减模型
    def decay_model(t, y, k):
        return -k * y
    
    result = ODESolver.solve(
        decay_model,
        y0=100.0,
        t_span=(0, 10),
        args=(0.5,)
    )
    
    if result.success:
        y_arr = np.array(result.y).flatten()
        print(f"  指数衰减模型: y' = -0.5*y")
        print(f"  初始: y(0) = {y_arr[0]:.2f}")
        print(f"  终值: y(10) = {y_arr[-1]:.2f}")
    
    # 2. 同位素公式
    print("\n[2] 同位素公式")
    print("-"*60)
    
    from toolkit.isotope.formulas import (
        DeltaCalculator, MassBalance, RayleighFractionation
    )
    
    # Delta转换
    ratio = DeltaCalculator.delta_to_ratio(-5, 0.0112372)
    print(f"  δ = -5‰ 对应比值: {ratio:.6f}")
    
    # 二元混合
    delta_mix = MassBalance.two_component_mixing(-4, -30, 0.7)
    print(f"  -4‰和-30‰按7:3混合: {delta_mix:.2f}‰")
    
    # 瑞利分馏
    delta_res = RayleighFractionation.residual_fraction(0.5, 1.025)
    print(f"  瑞利分馏(50%剩余, α=1.025): {delta_res:.2f}‰")
    
    # 3. 物理常数
    print("\n[3] 物理常数")
    print("-"*60)
    
    from toolkit.physics.constants import (
        ReservoirConstants, get_element_info, FractionationTheory
    )
    
    print(f"  海水DIC储量: {ReservoirConstants.SEAWATER_DIC_MOL:.2e} mol")
    print(f"  海水Mg储量: {ReservoirConstants.SEAWATER_MG_MOL:.2e} mol")
    
    mg_info = get_element_info('mg')
    print(f"  Mg同位素: {', '.join(mg_info.get('isotopes', []))}")
    
    # 温度依赖分馏
    epsilon = FractionationTheory.equilibrium_fractionation(298, 0, -0.8)
    print(f"  298K平衡分馏(简化模型): {epsilon:.3f}‰")
    
    print()


def run_comparison():
    """
    多体系对比分析
    """
    print("\n" + "="*80)
    print("多体系对比分析")
    print("="*80)
    
    from systems.mg import MgIsotopeSystem
    from systems.c import CIsotopeSystem
    
    print("\n不同同位素体系的特征对比:")
    print("-"*60)
    
    # Mg体系
    mg = MgIsotopeSystem()
    mg_info = mg.get_info()
    print(f"\nMg同位素体系:")
    print(f"  参考标准: {mg_info['parameters']['reference_standard']}")
    print(f"  端元数量: {len(mg_info['parameters']['end_members'])}")
    print(f"  主要应用: 风化来源示踪")
    
    # C体系
    c = CIsotopeSystem()
    c_info = c.get_info()
    print(f"\nC同位素体系:")
    print(f"  参考标准: {c_info['parameters']['reference_standard']}")
    print(f"  储库规模: {c_info['parameters']['reservoir_mass']:.2e} mol")
    print(f"  主要应用: 碳循环与氧化还原")
    
    print("\n潜在的多体系联合分析:")
    print("-"*60)
    print("  • Mg + C: 风化强度与碳循环耦合")
    print("  • C + S: 碳-硫循环与氧化还原")
    print("  • Sr + Os: 风化通量与气候演化")
    
    print()


def main():
    """
    主函数
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description='同位素质量平衡模型 - 统一入口',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python run.py all              # 运行所有演示
  python run.py mg               # 仅运行Mg同位素分析
  python run.py c                # 仅运行C同位素分析
  python run.py tools            # 仅运行核心工具演示
  python run.py mg --data FILE   # 使用数据文件运行Mg分析
        """
    )
    
    parser.add_argument(
        'command',
        choices=['all', 'mg', 'c', 'tools', 'compare'],
        help='要运行的分析类型'
    )
    parser.add_argument(
        '--data',
        type=str,
        help='数据文件路径（用于Mg分析）'
    )
    
    args = parser.parse_args()
    
    print("\n" + "="*80)
    print("同位素质量平衡模型框架")
    print("Isotope Mass Balance Modeling Framework")
    print("="*80)
    
    if args.command == 'all':
        run_core_tools_demo()
        run_mg_analysis(args.data)
        run_c_analysis()
        run_comparison()
    elif args.command == 'mg':
        run_mg_analysis(args.data)
    elif args.command == 'c':
        run_c_analysis()
    elif args.command == 'tools':
        run_core_tools_demo()
    elif args.command == 'compare':
        run_comparison()
    
    print("="*80)
    print("分析完成！")
    print("="*80 + "\n")


if __name__ == '__main__':
    main()
