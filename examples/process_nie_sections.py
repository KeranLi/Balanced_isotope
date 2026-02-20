#!/usr/bin/env python3
"""
Nie剖面Mg同位素数据处理示例
展示如何分析Nie_Section_A和Nie_Section_B两个剖面

参考标准: D3MS (Dead Sea Magnesium Standard)
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np
from systems.mg import MgIsotopeSystem, calculate_river_delta26, solve_f_silicate


def analyze_section(file_path, section_name, delta_silicate, delta_carbonate):
    """
    分析单个剖面数据
    
    Parameters
    ----------
    file_path : str
        Excel文件路径
    section_name : str
        剖面名称
    delta_silicate : float
        硅酸盐端元δ²⁶Mg
    delta_carbonate : float
        碳酸盐端元δ²⁶Mg
    
    Returns
    -------
    dict
        分析结果
    """
    print(f"\n{'='*60}")
    print(f"Analyzing: {section_name}")
    print(f"{'='*60}")
    
    # 读取数据
    df = pd.read_excel(file_path)
    df = df.dropna(axis=1, how='all')  # 移除空列
    
    print(f"Samples: {len(df)}")
    print(f"Reference: D3MS")
    print(f"Silicate end-member: {delta_silicate:+.2f}‰")
    print(f"Carbonate end-member: {delta_carbonate:+.2f}‰")
    
    # 提取数据
    delta_26 = df['delta_26_Mg_iso'].values
    delta_26_err = df['delta_26_Mg_iso_2sd'].values / 2  # 1σ
    
    print(f"\nδ²⁶Mg range: {delta_26.min():+.2f} to {delta_26.max():+.2f}‰")
    print(f"Mean: {delta_26.mean():+.2f} ± {delta_26.std():.2f}‰")
    
    # 计算风化比例
    f_silicate = []
    f_carbonate = []
    delta_river = []
    
    for d in delta_26:
        f_sil = solve_f_silicate(d, delta_silicate, delta_carbonate)
        f_silicate.append(f_sil)
        f_carbonate.append(1 - f_sil)
        delta_river.append(d)  # 假设数据代表河流/源区组成
    
    f_silicate = np.array(f_silicate)
    f_carbonate = np.array(f_carbonate)
    
    # 统计
    print(f"\nWeathering Ratios:")
    print(f"  Mean f_silicate: {f_silicate.mean():.1%} (range: {f_silicate.min():.0%} - {f_silicate.max():.0%})")
    print(f"  Mean f_carbonate: {f_carbonate.mean():.1%}")
    
    # 趋势分析
    x = np.arange(len(f_silicate))
    slope, intercept = np.polyfit(x, f_silicate, 1)
    print(f"\nTrend:")
    print(f"  Slope: {slope:+.4f} per sample")
    
    return {
        'name': section_name,
        'n_samples': len(df),
        'delta_26_mean': delta_26.mean(),
        'delta_26_std': delta_26.std(),
        'delta_26_range': (delta_26.min(), delta_26.max()),
        'f_silicate_mean': f_silicate.mean(),
        'f_silicate_std': f_silicate.std(),
        'f_silicate_range': (f_silicate.min(), f_silicate.max()),
        'trend_slope': slope,
        'data': pd.DataFrame({
            'delta_26Mg': delta_26,
            'delta_26Mg_err': delta_26_err,
            'f_silicate': f_silicate,
            'f_carbonate': f_carbonate
        })
    }


def main():
    """主函数"""
    print("\n" + "#"*70)
    print("#" + " Nie Section Mg Isotope Analysis ".center(68) + "#")
    print("#" + " Reference: D3MS ".center(68) + "#")
    print("#"*70)
    
    # 定义文件路径
    data_dir = Path(__file__).parent.parent / "data"
    results_dir = Path(__file__).parent.parent / "results"
    results_dir.mkdir(exist_ok=True)
    
    # 根据数据范围设置端元值
    # Section A: 范围 +0.44 到 +0.81‰
    # Section B: 范围 -0.38 到 +0.15‰
    
    config_a = {
        'file': data_dir / "Nie_Section_A.xlsx",
        'name': "Nie Section A",
        'delta_silicate': 0.70,
        'delta_carbonate': -1.20
    }
    
    config_b = {
        'file': data_dir / "Nie_Section_B.xlsx",
        'name': "Nie Section B",
        'delta_silicate': 0.10,
        'delta_carbonate': -1.50
    }
    
    # 分析两个剖面
    result_a = analyze_section(
        file_path=config_a['file'],
        section_name=config_a['name'],
        delta_silicate=config_a['delta_silicate'],
        delta_carbonate=config_a['delta_carbonate']
    )
    result_b = analyze_section(
        file_path=config_b['file'],
        section_name=config_b['name'],
        delta_silicate=config_b['delta_silicate'],
        delta_carbonate=config_b['delta_carbonate']
    )
    
    # 对比分析
    print(f"\n{'='*60}")
    print("Comparison: Section A vs Section B")
    print(f"{'='*60}")
    print(f"{'Metric':<30} {'Section A':<15} {'Section B':<15}")
    print("-"*60)
    print(f"{'Samples':<30} {result_a['n_samples']:<15} {result_b['n_samples']:<15}")
    print(f"{'Mean δ²⁶Mg':<30} {result_a['delta_26_mean']:+.2f}‰{'':<10} {result_b['delta_26_mean']:+.2f}‰")
    print(f"{'δ²⁶Mg range':<30} {result_a['delta_26_range'][0]:+.2f} to {result_a['delta_26_range'][1]:+.2f}{'':<3} {result_b['delta_26_range'][0]:+.2f} to {result_b['delta_26_range'][1]:+.2f}")
    print(f"{'Mean f_silicate':<30} {result_a['f_silicate_mean']:.1%}{'':<10} {result_b['f_silicate_mean']:.1%}")
    print(f"{'Trend slope':<30} {result_a['trend_slope']:+.4f}{'':<10} {result_b['trend_slope']:+.4f}")
    
    # 保存结果
    print(f"\n{'='*60}")
    print("Saving Results")
    print(f"{'='*60}")
    
    result_a['data'].to_csv(results_dir / "Nie_Section_A_results.csv", index=False)
    print(f"  Saved: Nie_Section_A_results.csv")
    
    result_b['data'].to_csv(results_dir / "Nie_Section_B_results.csv", index=False)
    print(f"  Saved: Nie_Section_B_results.csv")
    
    # 合并对比表
    comparison = pd.DataFrame({
        'Section': ['A', 'B'],
        'N_Samples': [result_a['n_samples'], result_b['n_samples']],
        'Delta26_Mean': [result_a['delta_26_mean'], result_b['delta_26_mean']],
        'Delta26_Std': [result_a['delta_26_std'], result_b['delta_26_std']],
        'F_Silicate_Mean': [result_a['f_silicate_mean'], result_b['f_silicate_mean']],
        'F_Silicate_Std': [result_a['f_silicate_std'], result_b['f_silicate_std']],
        'Trend_Slope': [result_a['trend_slope'], result_b['trend_slope']]
    })
    comparison.to_csv(results_dir / "Nie_Sections_comparison.csv", index=False)
    print(f"  Saved: Nie_Sections_comparison.csv")
    
    print("\n" + "#"*70)
    print("#" + " Analysis Complete ".center(68) + "#")
    print("#"*70 + "\n")


if __name__ == "__main__":
    main()
