#!/usr/bin/env python3
"""
同位素质量平衡模型 - 统一入口
基于Kasemann等(2014)论文的Mg同位素风化模型
"""

import sys
from pathlib import Path

# 添加项目根目录到Python路径
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

import numpy as np


def run_mg_weathering_analysis():
    """
    Mg同位素风化通量分析 - 基于Kasemann等(2014)论文
    """
    print("\n" + "="*80)
    print("Mg同位素风化通量模型 (Kasemann et al., 2014)")
    print("="*80)
    
    from systems.mg import MgIsotopeSystem
    
    # 创建Mg同位素体系
    mg = MgIsotopeSystem(scenario='modern')
    
    # 示例1：端元混合计算
    print("\n[1] 风化端元混合计算")
    print("-"*60)
    print(f"  硅酸盐端元: δ²⁶Mg = {mg._delta_sil:+.2f}‰")
    print(f"  碳酸盐端元: δ²⁶Mg = {mg._delta_carb:+.2f}‰")
    
    test_ratios = [0.0, 0.25, 0.5, 0.75, 1.0]
    print("\n  不同风化比例下的河流δ²⁶Mg:")
    for f_sil in test_ratios:
        delta_riv = mg.weathering_model.river_input(
            type('C', (), {
                'f_silicate': f_sil,
                'F_riv_multiplier': 1.0,
                'delta_silicate': mg._delta_sil,
                'delta_carbonate': mg._delta_carb,
                'Delta_carb': -2.7
            })()
        )[1]
        print(f"    f_silicate={f_sil:.2f}: δ²⁶Mg_riv = {delta_riv:+.2f}‰")
    
    # 示例2：稳态海水δ²⁶Mg
    print("\n[2] 稳态海水δ²⁶Mg计算")
    print("-"*60)
    from systems.mg.carbonate import WeatheringFluxConfig
    
    configs = [
        (0.2, "碳酸盐主导风化"),
        (0.5, "混合风化"),
        (0.8, "硅酸盐主导风化")
    ]
    
    for f_sil, desc in configs:
        config = WeatheringFluxConfig(f_silicate=f_sil)
        delta_sw = mg.weathering_model.steady_state_seawater(config)
        print(f"  {desc} (f_sil={f_sil:.1f}): δ²⁶Mg_sw = {delta_sw:+.2f}‰")
    
    # 示例3：风化转变演化模拟
    print("\n[3] 风化转变演化模拟 (0-5 Myr)")
    print("-"*60)
    print("  情景：从碳酸盐主导(f_sil=0.2) → 硅酸盐主导(f_sil=0.8)")
    
    transition = {
        't_start': 0.0,
        't_end': 2.0,
        'f_initial': 0.2,
        'f_final': 0.8,
        'mode': 'linear'
    }
    
    result = mg.simulate_weathering_transition(
        time_span=(0, 5.0),
        transition=transition,
        initial_delta=-2.0,
        n_points=500
    )
    
    if result.success:
        times_ma = result.time / 1e6
        delta_sw = result.data['delta_sw']
        f_sil = result.data['f_silicate']
        
        print(f"\n  时间点对比:")
        for t_target in [0, 1, 2, 3, 5]:
            idx = np.argmin(np.abs(times_ma - t_target))
            print(f"    t={t_target} Myr: f_sil={f_sil[idx]:.2f}, δ²⁶Mg_sw={delta_sw[idx]:+.2f}‰")
    
    # 示例4：Cryogenian冰期后情景
    print("\n[4] Cryogenian冰期后情景模拟 (Kasemann et al., 2014)")
    print("-"*60)
    print("  论文参数：")
    print("    - 阶段1 (0-0.5 Myr): 混合风化，9×现代通量")
    print("    - 阶段2 (0.5-1.5 Myr): 硅酸盐主导，6×现代通量")
    
    result_cryo = mg.simulate_cryogenian_scenario(duration_ma=3.0, n_points=300)
    
    if result_cryo.success:
        times_ma = result_cryo.time / 1e6
        delta_sw = result_cryo.data['delta_sw']
        f_sil = result_cryo.data['f_silicate']
        flux_mult = result_cryo.data['flux_multiplier']
        
        print(f"\n  演化结果:")
        key_times = [0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
        for t_target in key_times:
            if t_target <= times_ma[-1]:
                idx = np.argmin(np.abs(times_ma - t_target))
                print(f"    t={t_target:4.1f} Myr: "
                      f"flux={flux_mult[idx]:.1f}×, "
                      f"f_sil={f_sil[idx]:.2f}, "
                      f"δ²⁶Mg_sw={delta_sw[idx]:+.2f}‰")
    
    # 示例5：反演计算
    print("\n[5] 风化通量反演示例")
    print("-"*60)
    print("  假设观测数据（模拟碳酸盐δ²⁶Mg）:")
    
    # 模拟观测数据
    age_obs = np.array([0.5, 1.0, 1.5, 2.0, 2.5])
    # 模拟从-2.5‰到-1.0‰的碳酸盐同位素变化
    delta_carb_obs = np.array([-2.5, -2.2, -1.8, -1.4, -1.0])
    
    result_inv = mg.inverse_weathering_flux(
        age_data=age_obs,
        delta_carb_data=delta_carb_obs,
        assume_steady_state=True
    )
    
    if result_inv.success:
        f_sil_inv = result_inv.get('f_silicate')
        delta_riv_inv = result_inv.get('delta_river_inferred')
        
        print(f"\n  反演结果:")
        print(f"  {'Age (Myr)':<12} {'δ²⁶Mg_carb':<12} {'δ²⁶Mg_riv':<12} {'f_silicate':<12}")
        print("  " + "-"*48)
        for i in range(len(age_obs)):
            print(f"  {age_obs[i]:<12.1f} {delta_carb_obs[i]:<12.2f} "
                  f"{delta_riv_inv[i]:<12.2f} {f_sil_inv[i]:<12.2f}")
    
    print("\n" + "="*80)
    print("Mg同位素分析完成")
    print("="*80)


def run_c_analysis():
    """C同位素分析示例"""
    print("\n" + "="*80)
    print("C同位素DOC氧化分析")
    print("="*80)
    
    from systems.c import CIsotopeSystem
    
    c = CIsotopeSystem(scenario='dice')
    
    print("\n[1] 特定DOC通量下的稳态同位素")
    print("-"*60)
    
    F_odoc_values = [1e18, 2.1e18, 4.6e18, 10e18]
    for F_odoc in F_odoc_values:
        result = c.solve_steady_state(F_odoc=F_odoc)
        if result.success:
            delta_carb = result.get('delta13C_carb')
            print(f"  F_doc = {F_odoc/1e18:5.2f}×10¹⁸: δ¹³C_carb = {delta_carb:+.2f}‰")


def run_sr_analysis():
    """Sr同位素分析示例 - 基于Wang et al. (2021)"""
    print("\n" + "="*80)
    print("Sr同位素海洋箱模型 (Wang et al., 2021)")
    print("="*80)
    
    from systems.sr import SrIsotopeSystem, SrFluxConfig
    
    # 创建Sr同位素体系
    sr = SrIsotopeSystem(scenario='modern')
    
    # 示例1: 基本质量平衡
    print("\n[1] 现代海水⁸⁷Sr/⁸⁶Sr计算")
    print("-"*60)
    modern_ratio = sr.calculate_seawater_ratio()
    print(f"  计算值: {modern_ratio:.5f}")
    print(f"  观测值: 0.70917")
    print(f"  偏差: {(modern_ratio - 0.70917)*10000:+.2f} ppm")
    
    # 示例2: 不同端元混合
    print("\n[2] 四端元混合示例")
    print("-"*60)
    configs = [
        (47.6e9, 0.71107, "现代河流"),
        (100e9, 0.71107, "增强风化通量"),
        (47.6e9, 0.705, "低放射成因河流（玄武岩主导）"),
        (47.6e9, 0.720, "高放射成因河流（古老地壳主导）"),
    ]
    
    for f_riv, r_riv, desc in configs:
        ratio = sr.calculate_seawater_ratio(F_river=f_riv, R_river=r_riv)
        print(f"  {desc}: {ratio:.5f}")
    
    # 示例3: 热液活动影响
    print("\n[3] 热液活动对海水Sr的影响")
    print("-"*60)
    print("  高温热液通量变化:")
    f_highT_values = [2e9, 8e9, 20e9, 35e9]
    for f in f_highT_values:
        ratio = sr.calculate_seawater_ratio(F_hydrothermal_highT=f)
        print(f"    F_highT = {f/1e9:.0f}×10⁹ mol/yr: R_sw = {ratio:.5f}")
    
    # 示例4: 二叠纪演化趋势（简化演示）
    print("\n[4] 二叠纪Sr同位素演化趋势")
    print("-"*60)
    sr_permian = SrIsotopeSystem(scenario='permian')
    
    ages_ma = [299, 273, 265, 259, 252]
    stages = ["Asselian", "Kungurian", "Wordian", "Capitanian", "Changhsingian"]
    
    print(f"  {'时期':<15} {'年龄(Ma)':<10} {'⁸⁷Sr/⁸⁶Sr':<12} {'趋势':<10}")
    print("  " + "-"*50)
    
    # 典型观测值（基于Wang et al., 2021）
    observed = [0.70827, 0.7070, 0.70683, 0.7069, 0.70717]
    
    for i, (age, stage, ratio) in enumerate(zip(ages_ma, stages, observed)):
        if i > 0:
            change = "↓" if ratio < observed[i-1] else "↑"
        else:
            change = "→"
        print(f"  {stage:<15} {age:<10} {ratio:<12.5f} {change:<10}")
    
    print("\n  趋势解释:")
    print("    - 早-中二叠世: Sr同位素持续下降，反映热液活动增强")
    print("    - Capitanian: 达到最低值(~0.70683)")
    print("    - 晚二叠世: 上升，可能与大陆风化增强有关")
    
    print("\n" + "="*80)
    print("Sr同位素分析完成")
    print("="*80)


def main():
    """主函数"""
    print("\n" + "#"*80)
    print("#" + " "*78 + "#")
    print("#" + "  同位素质量平衡模型 v2.0".center(78) + "#")
    print("#" + "  基于Kasemann等(2014)Mg同位素风化模型".center(78) + "#")
    print("#" + "  以及Wang等(2021)Sr同位素海洋箱模型".center(78) + "#")
    print("#" + " "*78 + "#")
    print("#"*80)
    
    # 运行Mg同位素分析
    run_mg_weathering_analysis()
    
    # 运行C同位素分析
    run_c_analysis()
    
    # 运行Sr同位素分析
    run_sr_analysis()
    
    print("\n" + "#"*80)
    print("所有分析完成!")
    print("#"*80 + "\n")


if __name__ == "__main__":
    main()
