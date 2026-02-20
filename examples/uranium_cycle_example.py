#!/usr/bin/env python3
"""
铀同位素海洋循环模型 - 使用示例

本示例演示:
1. 稳态模型: 从碳酸盐δ²³⁸U计算缺氧汇比例
2. 非稳态模型: 模拟缺氧事件的时间演化
3. 批量计算和反演
"""

import numpy as np
import sys
from pathlib import Path

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from systems.u import UIsotopeSystem, get_scenario_info


def example_1_steady_state():
    """
    示例1: 稳态模型计算
    
    从测量的碳酸盐δ²³⁸U值计算缺氧汇比例 f_anox
    """
    print("=" * 60)
    print("示例1: 稳态模型 - 从δ²³⁸U计算 f_anox")
    print("=" * 60)
    
    # 创建铀同位素体系实例
    u_system = UIsotopeSystem(scenario='modern')
    
    # 示例数据: Frasnian-Famennian 过渡期的碳酸盐δ²³⁸U
    delta_measurements = [
        (-0.32, "BS-1 (background)"),
        (-0.42, "BS-12 (early anoxia)"),
        (-0.54, "BS-13 (peak anoxia)"),
        (-0.18, "BS-3 (recovery)"),
    ]
    
    print("\n计算结果 (使用成岩校正 Δ_diag = 0.4‰):")
    print("-" * 60)
    print(f"{'Sample':<20} {'δ²³⁸U_meas':<12} {'f_anox':<10} {'Anoxic Area':<15}")
    print("-" * 60)
    
    for delta_meas, sample_name in delta_measurements:
        result = u_system.calculate_f_anox_steady_state(
            delta238_carb=delta_meas,
            apply_diagenetic_correction=True,
            delta_diag=0.4
        )
        
        # 估算缺氧面积
        anoxic_area = u_system.estimate_anoxic_area(result['f_anox'])
        
        print(f"{sample_name:<20} {delta_meas:+.2f}‰       {result['f_anox']:.0%}        ~{anoxic_area:.1f}%")
    
    print("-" * 60)
    print("\n解释:")
    print("  - f_anox = 0.2 (20%): 现代氧化海洋条件")
    print("  - f_anox = 0.6 (60%): 中等缺氧条件")
    print("  - f_anox = 0.8 (80%): 严重缺氧条件")


def example_2_transient():
    """
    示例2: 非稳态模型 - 模拟缺氧事件
    
    模拟一次持续1 Myr的缺氧事件，观察海水δ²³⁸U的响应
    """
    print("\n" + "=" * 60)
    print("示例2: 非稳态模型 - 缺氧事件模拟")
    print("=" * 60)
    
    u_system = UIsotopeSystem(scenario='modern', model_type='transient')
    
    # 模拟参数
    event_duration = 1.0  # 事件持续 1 Myr
    peak_f_anox = 0.8     # 峰值缺氧汇比例 80%
    background_f_anox = 0.2  # 背景值 20%
    
    print(f"\n模拟参数:")
    print(f"  事件持续时间: {event_duration} Myr")
    print(f"  背景 f_anox: {background_f_anox} (20%)")
    print(f"  峰值 f_anox: {peak_f_anox} (80%)")
    
    # 运行模拟
    result = u_system.simulate_anoxic_event(
        event_duration=event_duration,
        peak_f_anox=peak_f_anox,
        background_f_anox=background_f_anox,
        n_points=500
    )
    
    if result.success:
        data = result.data
        time = data['time_myr']
        delta_sw = data['delta_seawater']
        f_anox = data['f_anox']
        
        # 关键时间点
        bg_idx = 0  # 背景
        min_idx = np.argmin(delta_sw)  # 最小值
        end_idx = -1  # 结束
        
        print(f"\n模拟结果:")
        print(f"  背景 δ²³⁸U_sw:     {delta_sw[bg_idx]:+.2f}‰")
        print(f"  最小 δ²³⁸U_sw:     {delta_sw[min_idx]:+.2f}‰ (t={time[min_idx]:.2f} Myr)")
        print(f"  结束 δ²³⁸U_sw:     {delta_sw[end_idx]:+.2f}‰")
        print(f"  偏移幅度:          {delta_sw[min_idx] - delta_sw[bg_idx]:.2f}‰")
        
        # 检查系统是否达到平衡
        equil = u_system.solve_equilibration_time(
            target_f_anox=peak_f_anox,
            initial_f_anox=background_f_anox
        )
        
        print(f"\n系统响应特征:")
        print(f"  铀停留时间 (τ):    {equil['residence_time']:.2f} Myr")
        print(f"  平衡特征时间 (3τ): {equil['equilibration_time']:.2f} Myr")
        print(f"  理论平衡态 δ²³⁸U:  {equil['final_delta']:+.2f}‰")
        
        if event_duration < equil['equilibration_time']:
            print(f"\n  注: 事件持续时间 < 平衡时间")
            print(f"     系统未达到新的稳态，呈现瞬态特征")
        else:
            print(f"\n  注: 事件持续时间 > 平衡时间")
            print(f"     系统有足够时间接近新的稳态")


def example_3_comparison():
    """
    示例3: 稳态 vs 非稳态对比
    
    比较相同 f_anox 变化下，稳态预测与非稳态实际演化的差异
    """
    print("\n" + "=" * 60)
    print("示例3: 稳态 vs 非稳态对比")
    print("=" * 60)
    
    u_system = UIsotopeSystem(scenario='modern')
    
    # 场景: f_anox 从 20% 快速变化到 80%
    f_anox_values = [0.2, 0.4, 0.6, 0.8]
    
    print("\n稳态预测:")
    print("-" * 40)
    print(f"{'f_anox':<10} {'δ²³⁸U_sw (稳态)':<15}")
    print("-" * 40)
    
    steady_results = {}
    for f_anox in f_anox_values:
        result = u_system.calculate_seawater_delta_steady_state(f_anox=f_anox)
        steady_results[f_anox] = result['delta238_seawater']
        print(f"{f_anox:<10.1f} {result['delta238_seawater']:+.2f}‰")
    
    # 非稳态: 观察突变后的响应
    print("\n非稳态响应 (突变到 f_anox=0.8):")
    print("-" * 40)
    
    # 使用阶跃函数模拟突变
    def step_function(t):
        return 0.8 if t >= 0 else 0.2
    
    result = u_system.solve_transient(
        initial_delta238=steady_results[0.2],
        time_span=(0, 2.0),  # 2 Myr
        f_anox_function=step_function,
        n_points=100
    )
    
    if result.success:
        data = result.data
        time = data['time_myr']
        delta_sw = data['delta_seawater']
        
        # 显示关键时间点
        print(f"{'Time (Myr)':<12} {'δ²³⁸U_sw':<12} {'vs 稳态':<12}")
        print("-" * 40)
        
        for t_check in [0, 0.1, 0.25, 0.5, 1.0, 2.0]:
            idx = np.argmin(np.abs(time - t_check))
            deviation = delta_sw[idx] - steady_results[0.8]
            print(f"{t_check:<12.2f} {delta_sw[idx]:+.2f}‰      {deviation:+.2f}‰")


def example_4_inverse_model():
    """
    示例4: 反演模型与不确定性分析
    
    从观测数据推断 f_anox，并评估不确定性
    """
    print("\n" + "=" * 60)
    print("示例4: 反演模型 - 不确定性分析")
    print("=" * 60)
    
    u_system = UIsotopeSystem(scenario='modern')
    
    # 模拟观测数据
    observed_delta238 = -0.65  # ‰
    uncertainty = 0.05  # 1σ = 0.05‰
    
    print(f"\n观测数据:")
    print(f"  δ²³⁸U_carb = {observed_delta238:.2f} ± {uncertainty:.2f}‰")
    
    # 反演计算
    print(f"\n蒙特卡洛反演 (n=10000):")
    print("  考虑分馏系数的不确定性...")
    
    result = u_system.inverse_model_steady_state(
        observed_delta238_carb=observed_delta238,
        uncertainty=uncertainty,
        n_monte_carlo=10000
    )
    
    print(f"\n结果:")
    print(f"  f_anox 均值:    {result['f_anox_mean']:.1%}")
    print(f"  f_anox 中位数:  {result['f_anox_median']:.1%}")
    print(f"  标准差:         {result['f_anox_std']:.1%}")
    print(f"  95% 置信区间:   [{result['f_anox_ci95'][0]:.1%}, {result['f_anox_ci95'][1]:.1%}]")


def example_5_batch_processing():
    """
    示例5: 批量处理地层数据
    
    处理一列地层样品的δ²³⁸U数据
    """
    print("\n" + "=" * 60)
    print("示例5: 批量处理 - 地层剖面分析")
    print("=" * 60)
    
    u_system = UIsotopeSystem(scenario='frasnian_famennian')
    
    # 模拟地层数据 (高程, δ²³⁸U)
    stratigraphic_data = [
        (0.05, -0.32),   # 背景
        (0.85, -0.37),
        (1.65, -0.41),
        (3.15, -0.18),
        (5.35, -0.22),
        (9.55, -0.21),
        (11.55, -0.12),
        (17.55, -0.15),
        (19.55, -0.42),  # 缺氧开始
        (21.05, -0.54),  # 缺氧峰值
        (23.05, -0.38),
        (25.05, -0.29),
        (27.05, -0.25),
        (29.05, -0.18),
    ]
    
    elevations = np.array([d[0] for d in stratigraphic_data])
    delta_values = np.array([d[1] for d in stratigraphic_data])
    
    print(f"\n处理 {len(stratigraphic_data)} 个样品...")
    
    # 批量计算
    results = u_system.batch_steady_state_calculation(
        delta238_carb_array=delta_values,
        apply_diagenetic_correction=True,
        delta_diag=0.4
    )
    
    print("\n结果:")
    print("-" * 70)
    print(f"{'Elev (m)':<10} {'δ²³⁸U':<10} {'f_anox':<10} {'Interpretation':<30}")
    print("-" * 70)
    
    for i, (elev, delta) in enumerate(stratigraphic_data):
        f_anox = results['f_anox'][i]
        
        # 解释
        if f_anox < 0.3:
            interp = "Oxic conditions"
        elif f_anox < 0.5:
            interp = "Moderate anoxia"
        elif f_anox < 0.7:
            interp = "Strong anoxia"
        else:
            interp = "Severe anoxia/euxinia"
        
        print(f"{elev:<10.2f} {delta:<+10.2f} {f_anox:<10.0%} {interp:<30}")
    
    print("-" * 70)
    
    # 识别缺氧事件
    max_anoxia_idx = np.argmax(results['f_anox'])
    print(f"\n缺氧峰值:")
    print(f"  位置: {elevations[max_anoxia_idx]:.2f} m")
    print(f"  f_anox: {results['f_anox'][max_anoxia_idx]:.0%}")


def main():
    """运行所有示例"""
    print("\n" + "=" * 60)
    print("铀同位素海洋循环模型 - 完整示例")
    print("=" * 60)
    
    # 运行所有示例
    example_1_steady_state()
    example_2_transient()
    example_3_comparison()
    example_4_inverse_model()
    example_5_batch_processing()
    
    print("\n" + "=" * 60)
    print("示例运行完成!")
    print("=" * 60)
    print("\n更多信息:")
    print("  - CLI使用: python cli.py u --help")
    print("  - 查看源码: systems/u/model.py")
    print()


if __name__ == '__main__':
    main()
