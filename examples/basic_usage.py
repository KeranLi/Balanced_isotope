"""
基础使用示例
展示如何使用新架构的各层功能
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np


def example_1_core_math():
    """示例1: 使用核心数学工具"""
    print("=== Example 1: Core Math Tools ===\n")
    
    from toolkit.math.numerical import ODESolver, Interpolator
    
    # ODE求解示例：指数衰减
    def exponential_decay(t, y, k):
        return -k * y
    
    result = ODESolver.solve(
        exponential_decay,
        y0=1.0,
        t_span=(0, 10),
        args=(0.5,),  # k = 0.5
        n_points=100
    )
    
    print(f"ODE solved: {len(result.t)} points")
    y_arr = np.array(result.y)
    if y_arr.ndim > 1:
        y_arr = y_arr.flatten()
    print(f"Initial: y(0) = {y_arr[0]:.4f}")
    print(f"Final: y(10) = {y_arr[-1]:.4f}")
    print()
    
    # 插值示例
    x = np.array([0, 1, 2, 3, 4])
    y = np.array([0, 1, 4, 9, 16])  # y = x^2
    x_new = np.linspace(0, 4, 20)
    
    y_interp = Interpolator.interpolate(x, y, x_new, method='cubic')
    print(f"Interpolation: {len(x)} points -> {len(x_new)} points")
    print()


def example_2_core_isotope():
    """示例2: 使用同位素公式"""
    print("=== Example 2: Isotope Formulas ===\n")
    
    from toolkit.isotope.formulas import (
        DeltaCalculator, MassBalance, RayleighFractionation
    )
    
    # Delta转换
    standard_ratio = 0.0112372  # VPDB
    sample_ratio = 0.0111500
    delta = DeltaCalculator.ratio_to_delta(sample_ratio, standard_ratio)
    print(f"Sample ratio: {sample_ratio}")
    print(f"Delta value: {delta:.2f}‰")
    
    # 混合计算
    delta_a, delta_b = -4.0, -30.0  # 碳酸盐 vs 有机碳
    f_a = 0.7
    delta_mix = MassBalance.two_component_mixing(delta_a, delta_b, f_a)
    print(f"\nTwo-component mixing:")
    print(f"  End A: {delta_a}‰ ({f_a*100:.0f}%)")
    print(f"  End B: {delta_b}‰ ({(1-f_a)*100:.0f}%)")
    print(f"  Mixture: {delta_mix:.2f}‰")
    
    # 瑞利分馏
    f_residual = 0.5  # 50%剩余
    alpha = 1.025  # 分馏系数
    delta_residual = RayleighFractionation.residual_fraction(f_residual, alpha)
    print(f"\nRayleigh fractionation:")
    print(f"  Residual fraction: {f_residual}")
    print(f"  Residual delta: {delta_residual:.2f}‰")
    print()


def example_3_mg_system():
    """示例3: Mg同位素体系"""
    print("=== Example 3: Mg Isotope System ===\n")
    
    from systems.mg import MgIsotopeSystem
    
    # 创建体系实例
    mg_system = MgIsotopeSystem()
    
    # 获取参数信息
    info = mg_system.get_info()
    print(f"Element: {info['name']} ({info['element']})")
    print(f"Reference: {info['parameters']['reference_standard']}")
    
    # 风化比例计算
    delta_sample = -2.5  # 样品值
    delta_sw = -0.83     # 海水值
    
    ratios = mg_system.calculate_weathering_ratio(delta_sample, delta_sw)
    print(f"\nWeathering ratio calculation:")
    print(f"  Sample δ²⁶Mg: {delta_sample}‰")
    print(f"  Seawater δ²⁶Mg: {delta_sw}‰")
    print(f"  Carbonate fraction: {ratios['f_carbonate']:.1%}")
    print(f"  Silicate fraction: {ratios['f_silicate']:.1%}")
    
    # 海水演化模拟
    print(f"\nSeawater evolution (100 Ma):")
    result = mg_system.seawater_evolution(
        time_span=(0, 100),
        initial_delta=-0.5,
        flux_scenario='modern'
    )
    
    if result.success:
        values = np.atleast_2d(result.values)
        delta_mg = values[:, 1]  # δ²⁶Mg
        print(f"  Initial δ²⁶Mg_sw: {delta_mg[0]:.3f}‰")
        print(f"  Final δ²⁶Mg_sw: {delta_mg[-1]:.3f}‰")
    print()


def example_4_c_system():
    """示例4: C同位素体系"""
    print("=== Example 4: C Isotope System ===\n")
    
    from systems.c import CIsotopeSystem
    
    # 创建DICE情景的C体系
    c_system = CIsotopeSystem(scenario='dice')
    
    print(f"Scenario: DICE (Doushantuo Carbon Isotope Excursion)")
    print(f"Reference: Li et al. 2020, Precambrian Research")
    
    # 计算特定DOC通量下的稳态
    F_odoc = 4.0e18  # mol/Ma
    result = c_system.solve_steady_state(F_odoc=F_odoc)
    
    if result.success:
        print(f"\nSteady-state with F_odoc = {F_odoc:.2e} mol/Ma:")
        print(f"  δ¹³C_carb: {result.get('delta13C_carb'):.2f}‰")
        print(f"  δ¹³C_org: {result.get('delta13C_org'):.2f}‰")
    
    # 查找产生4‰负漂的DOC通量
    print(f"\nFinding DOC flux for -4‰ excursion...")
    result = c_system.find_doc_for_excursion(target_excursion=-4.0)
    print(f"  Required F_odoc: {result['F_odoc']:.2e} mol/Ma")
    
    # 氧化剂消耗
    from systems.c.parameters import OXIDANT_SCENARIOS
    print(f"\nOxidant consumption (F_odoc = {F_odoc:.2e} mol/Ma):")
    for scenario_key in ['modern_o2_high_sulfate', 'low_o2_high_sulfate']:
        ox = c_system.calculate_oxidant_consumption(F_odoc, scenario_key)
        name = OXIDANT_SCENARIOS[scenario_key]['name']
        print(f"  {name}:")
        print(f"    Total: {ox.total_consumption/1e18:.2f}×10¹⁸ mol/Ma")
    print()


def example_5_custom_model():
    """示例5: 使用底层工具构建自定义模型"""
    print("=== Example 5: Custom Model ===\n")
    
    from toolkit.math.numerical import ODESolver
    from toolkit.physics.constants import FractionationTheory
    
    # 自定义：温度依赖的分馏演化模型
    
    def isotope_evolution_with_cooling(t, delta, params):
        """
        模拟随温度降低的同位素分馏演化
        
        dδ/dt = k * (δ_eq(T(t)) - δ)
        """
        T_initial, cooling_rate, k, A, B = params
        
        # 温度随时间降低
        T = T_initial - cooling_rate * t
        
        # 平衡分馏（温度依赖）
        epsilon_eq = FractionationTheory.equilibrium_fractionation(T, A, B)
        
        # 向平衡态弛豫
        return k * (epsilon_eq - delta)
    
    # 参数：初始温度300K，冷却速率10K/Ma，A=0，B=-0.8
    params = (300, 10, 0.1, 0, -0.8)
    
    result = ODESolver.solve(
        isotope_evolution_with_cooling,
        y0=0.0,  # 初始与标准相同
        t_span=(0, 20),  # 20 Ma
        args=(params,),
        n_points=200
    )
    
    if result.success:
        y_arr = np.array(result.y)
        if y_arr.ndim > 1:
            y_arr = y_arr.flatten()
        print("Custom cooling model:")
        print(f"  Simulation: 20 Ma cooling from 300K")
        print(f"  Initial δ: {y_arr[0]:.3f}‰")
        print(f"  Final δ: {y_arr[-1]:.3f}‰")
        print(f"  Total shift: {y_arr[-1] - y_arr[0]:.3f}‰")
    print()


def run_all_examples():
    """运行所有示例"""
    print("\n" + "="*60)
    print("Isotope Mass Balance Framework - Usage Examples")
    print("="*60 + "\n")
    
    example_1_core_math()
    example_2_core_isotope()
    example_3_mg_system()
    example_4_c_system()
    example_5_custom_model()
    
    print("="*60)
    print("All examples completed!")
    print("="*60)


if __name__ == '__main__':
    run_all_examples()
