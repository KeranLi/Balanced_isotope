#!/usr/bin/env python3
"""
基于Nie Section A和B剖面的Mg同位素风化分析

本脚本实现：
1. 数据加载和可视化
2. 风化端元比例计算（碳酸盐 vs 硅酸盐）
3. 风化通量反演
4. 剖面对比分析
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, Tuple

from systems.mg import MgIsotopeSystem


def load_section_data(file_path: str) -> pd.DataFrame:
    """加载剖面数据"""
    df = pd.read_excel(file_path)
    
    # 添加样品编号（假设等间距采样）
    df['sample_id'] = range(1, len(df) + 1)
    
    return df


def calculate_weathering_ratios(df: pd.DataFrame, 
                                 delta_seawater: float = -0.83) -> pd.DataFrame:
    """
    计算每个样品的风化端元比例
    
    Parameters
    ----------
    df : DataFrame
        包含Mg同位素数据
    delta_seawater : float
        同期海水Mg同位素值（现代海水约-0.83‰）
        
    Returns
    -------
    DataFrame
        添加风化比例列
    """
    mg_system = MgIsotopeSystem()
    
    # 计算每个样品的风化比例
    f_carb_list = []
    f_sil_list = []
    
    for _, row in df.iterrows():
        delta_sample = row['delta_26_Mg_iso']
        
        # 使用二元混合模型计算风化比例
        ratios = mg_system.calculate_weathering_ratio(
            delta_sample=delta_sample,
            delta_seawater=delta_seawater
        )
        
        f_carb_list.append(ratios['f_carbonate'])
        f_sil_list.append(ratios['f_silicate'])
    
    df['f_carbonate'] = f_carb_list
    df['f_silicate'] = f_sil_list
    
    return df


def inverse_weathering_flux(df: pd.DataFrame,
                           delta_seawater: float = -0.83,
                           total_flux_range: Tuple[float, float] = (1e18, 10e18)) -> pd.DataFrame:
    """
    反演风化通量
    
    基于观测的Mg同位素组成，反演碳酸盐和硅酸盐风化通量
    
    Parameters
    ----------
    df : DataFrame
        包含Mg同位素数据
    delta_seawater : float
        海水Mg同位素值
    total_flux_range : tuple
        总风化通量范围 (mol/Ma)
        
    Returns
    -------
    DataFrame
        添加通量估算列
    """
    mg_system = MgIsotopeSystem()
    params = mg_system.params
    
    # 端元值
    delta_carb = params.end_members['carbonate']['delta26']
    delta_sil = params.end_members['silicate']['delta26']
    
    # 假设总风化通量（需要根据地质背景调整）
    # 这里使用简化的估算：基于偏离海水的程度
    F_carb_list = []
    F_sil_list = []
    F_total_list = []
    
    for _, row in df.iterrows():
        delta_sample = row['delta_26_Mg_iso']
        f_carb = row['f_carbonate']
        f_sil = row['f_silicate']
        
        # 简化的通量估算模型
        # 假设与现代河流输入成比例
        # 可以根据具体地质背景调整此模型
        base_flux = 5e18  # 现代全球河流Mg通量约5e18 mol/Ma
        
        # 根据同位素偏移调整（偏移越大，风化通量可能越高）
        offset = abs(delta_sample - delta_seawater)
        flux_factor = 1 + offset  # 简化线性关系
        
        F_total = base_flux * flux_factor
        F_carb = F_total * f_carb
        F_sil = F_total * f_sil
        
        F_total_list.append(F_total)
        F_carb_list.append(F_carb)
        F_sil_list.append(F_sil)
    
    df['F_total'] = F_total_list
    df['F_carbonate'] = F_carb_list
    df['F_silicate'] = F_sil_list
    
    return df


def plot_section_comparison(df_a: pd.DataFrame, df_b: pd.DataFrame):
    """
    绘制两个剖面的对比图
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # 1. Mg同位素原始数据
    ax = axes[0, 0]
    ax.errorbar(df_a['sample_id'], df_a['delta_26_Mg_iso'], 
                yerr=df_a['delta_26_Mg_iso_2sd']/2, 
                fmt='o', label='Section A', capsize=3, color='blue')
    ax.errorbar(df_b['sample_id'], df_b['delta_26_Mg_iso'], 
                yerr=df_b['delta_26_Mg_iso_2sd']/2, 
                fmt='s', label='Section B', capsize=3, color='red')
    ax.axhline(y=-0.83, color='gray', linestyle='--', label='Modern Seawater')
    ax.axhline(y=-4.3, color='green', linestyle=':', alpha=0.5, label='Carbonate End-member')
    ax.axhline(y=-0.3, color='orange', linestyle=':', alpha=0.5, label='Silicate End-member')
    ax.set_xlabel('Sample ID')
    ax.set_ylabel('δ²⁶Mg (‰)')
    ax.set_title('Mg Isotope Compositions')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. 风化比例 - Section A
    ax = axes[0, 1]
    x = df_a['sample_id']
    ax.stackplot(x, df_a['f_carbonate'], df_a['f_silicate'], 
                 labels=['Carbonate Weathering', 'Silicate Weathering'],
                 colors=['lightblue', 'lightcoral'], alpha=0.7)
    ax.plot(x, df_a['f_carbonate'], 'o-', color='blue', label='f_carbonate')
    ax.plot(x, df_a['f_silicate'], 's-', color='red', label='f_silicate')
    ax.set_xlabel('Sample ID')
    ax.set_ylabel('Weathering Fraction')
    ax.set_title('Section A: Weathering End-member Proportions')
    ax.set_ylim(0, 1)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 3. 风化比例 - Section B
    ax = axes[0, 2]
    x = df_b['sample_id']
    ax.stackplot(x, df_b['f_carbonate'], df_b['f_silicate'], 
                 labels=['Carbonate Weathering', 'Silicate Weathering'],
                 colors=['lightblue', 'lightcoral'], alpha=0.7)
    ax.plot(x, df_b['f_carbonate'], 'o-', color='blue')
    ax.plot(x, df_b['f_silicate'], 's-', color='red')
    ax.set_xlabel('Sample ID')
    ax.set_ylabel('Weathering Fraction')
    ax.set_title('Section B: Weathering End-member Proportions')
    ax.set_ylim(0, 1)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 4. 通量估算 - Section A
    ax = axes[1, 0]
    x = df_a['sample_id']
    ax.fill_between(x, 0, df_a['F_carbonate']/1e18, 
                     alpha=0.5, color='blue', label='Carbonate Flux')
    ax.fill_between(x, df_a['F_carbonate']/1e18, 
                     (df_a['F_carbonate'] + df_a['F_silicate'])/1e18, 
                     alpha=0.5, color='red', label='Silicate Flux')
    ax.set_xlabel('Sample ID')
    ax.set_ylabel('Weathering Flux (×10¹⁸ mol/Ma)')
    ax.set_title('Section A: Estimated Weathering Fluxes')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 5. 通量估算 - Section B
    ax = axes[1, 1]
    x = df_b['sample_id']
    ax.fill_between(x, 0, df_b['F_carbonate']/1e18, 
                     alpha=0.5, color='blue', label='Carbonate Flux')
    ax.fill_between(x, df_b['F_carbonate']/1e18, 
                     (df_b['F_carbonate'] + df_b['F_silicate'])/1e18, 
                     alpha=0.5, color='red', label='Silicate Flux')
    ax.set_xlabel('Sample ID')
    ax.set_ylabel('Weathering Flux (×10¹⁸ mol/Ma)')
    ax.set_title('Section B: Estimated Weathering Fluxes')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 6. 剖面对比统计
    ax = axes[1, 2]
    stats = ['Mean δ²⁶Mg', 'Mean f_carb', 'Mean f_sil', 'Mean F_total']
    section_a_values = [
        df_a['delta_26_Mg_iso'].mean(),
        df_a['f_carbonate'].mean(),
        df_a['f_silicate'].mean(),
        df_a['F_total'].mean() / 1e18
    ]
    section_b_values = [
        df_b['delta_26_Mg_iso'].mean(),
        df_b['f_carbonate'].mean(),
        df_b['f_silicate'].mean(),
        df_b['F_total'].mean() / 1e18
    ]
    
    x = np.arange(len(stats))
    width = 0.35
    ax.bar(x - width/2, section_a_values, width, label='Section A', color='blue', alpha=0.7)
    ax.bar(x + width/2, section_b_values, width, label='Section B', color='red', alpha=0.7)
    ax.set_ylabel('Value')
    ax.set_title('Section Comparison Statistics')
    ax.set_xticks(x)
    ax.set_xticklabels(stats, rotation=15, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('Section_Mg_Isotope_Analysis.png', dpi=300, bbox_inches='tight')
    print("图表已保存: Section_Mg_Isotope_Analysis.png")
    plt.close()  # 关闭图形，不显示


def print_statistics(df_a: pd.DataFrame, df_b: pd.DataFrame):
    """打印统计分析结果"""
    print("\n" + "="*80)
    print("剖面统计分析")
    print("="*80)
    
    print("\n【Section A】")
    print(f"样品数量: {len(df_a)}")
    print(f"δ²⁶Mg 范围: {df_a['delta_26_Mg_iso'].min():.3f} ~ {df_a['delta_26_Mg_iso'].max():.3f}‰")
    print(f"δ²⁶Mg 均值: {df_a['delta_26_Mg_iso'].mean():.3f} ± {df_a['delta_26_Mg_iso'].std():.3f}‰")
    print(f"\n风化比例:")
    print(f"  碳酸盐风化: {df_a['f_carbonate'].mean():.1%} ± {df_a['f_carbonate'].std():.1%}")
    print(f"  硅酸盐风化: {df_a['f_silicate'].mean():.1%} ± {df_a['f_silicate'].std():.1%}")
    print(f"\n估算风化通量:")
    print(f"  总通量: {df_a['F_total'].mean()/1e18:.2f} ± {df_a['F_total'].std()/1e18:.2f} ×10¹⁸ mol/Ma")
    print(f"  碳酸盐通量: {df_a['F_carbonate'].mean()/1e18:.2f} ×10¹⁸ mol/Ma")
    print(f"  硅酸盐通量: {df_a['F_silicate'].mean()/1e18:.2f} ×10¹⁸ mol/Ma")
    
    print("\n【Section B】")
    print(f"样品数量: {len(df_b)}")
    print(f"δ²⁶Mg 范围: {df_b['delta_26_Mg_iso'].min():.3f} ~ {df_b['delta_26_Mg_iso'].max():.3f}‰")
    print(f"δ²⁶Mg 均值: {df_b['delta_26_Mg_iso'].mean():.3f} ± {df_b['delta_26_Mg_iso'].std():.3f}‰")
    print(f"\n风化比例:")
    print(f"  碳酸盐风化: {df_b['f_carbonate'].mean():.1%} ± {df_b['f_carbonate'].std():.1%}")
    print(f"  硅酸盐风化: {df_b['f_silicate'].mean():.1%} ± {df_b['f_silicate'].std():.1%}")
    print(f"\n估算风化通量:")
    print(f"  总通量: {df_b['F_total'].mean()/1e18:.2f} ± {df_b['F_total'].std()/1e18:.2f} ×10¹⁸ mol/Ma")
    print(f"  碳酸盐通量: {df_b['F_carbonate'].mean()/1e18:.2f} ×10¹⁸ mol/Ma")
    print(f"  硅酸盐通量: {df_b['F_silicate'].mean()/1e18:.2f} ×10¹⁸ mol/Ma")
    
    print("\n【剖面对比】")
    print(f"δ²⁶Mg 差异: {df_a['delta_26_Mg_iso'].mean() - df_b['delta_26_Mg_iso'].mean():.3f}‰")
    print(f"碳酸盐风化比例差异: {(df_a['f_carbonate'].mean() - df_b['f_carbonate'].mean())*100:.1f}%")
    print(f"总风化通量差异: {(df_a['F_total'].mean() - df_b['F_total'].mean())/1e18:.2f} ×10¹⁸ mol/Ma")


def analyze_temporal_trends(df_a: pd.DataFrame, df_b: pd.DataFrame):
    """分析时间趋势（假设样品顺序代表时间顺序）"""
    print("\n" + "="*80)
    print("时间趋势分析（基于样品顺序）")
    print("="*80)
    
    # 计算线性趋势
    from scipy import stats
    
    # Section A趋势
    slope_a, intercept_a, r_value_a, p_value_a, std_err_a = stats.linregress(
        df_a['sample_id'], df_a['delta_26_Mg_iso']
    )
    
    print("\n【Section A 趋势】")
    print(f"δ²⁶Mg 趋势: {'增加' if slope_a > 0 else '减少'} ({slope_a:.4f}‰/sample)")
    print(f"R² = {r_value_a**2:.3f}, p = {p_value_a:.4f}")
    
    if slope_a > 0:
        print("  → 向上变重，可能指示碳酸盐风化增强或硅酸盐风化减弱")
    else:
        print("  → 向上变轻，可能指示硅酸盐风化增强或碳酸盐风化减弱")
    
    # Section B趋势
    slope_b, intercept_b, r_value_b, p_value_b, std_err_b = stats.linregress(
        df_b['sample_id'], df_b['delta_26_Mg_iso']
    )
    
    print("\n【Section B 趋势】")
    print(f"δ²⁶Mg 趋势: {'增加' if slope_b > 0 else '减少'} ({slope_b:.4f}‰/sample)")
    print(f"R² = {r_value_b**2:.3f}, p = {p_value_b:.4f}")
    
    if slope_b > 0:
        print("  → 向上变重，可能指示碳酸盐风化增强或硅酸盐风化减弱")
    else:
        print("  → 向上变轻，可能指示硅酸盐风化增强或碳酸盐风化减弱")


def main():
    """主函数"""
    print("="*80)
    print("Nie Section A & B - Mg同位素风化分析")
    print("="*80)
    
    # 1. 加载数据
    print("\n[1] 加载剖面数据...")
    df_a = load_section_data('data/Nie_Section_A.xlsx')
    df_b = load_section_data('data/Nie_Section_B.xlsx')
    print(f"  Section A: {len(df_a)} 个样品")
    print(f"  Section B: {len(df_b)} 个样品")
    
    # 2. 计算风化比例
    print("\n[2] 计算风化端元比例...")
    # 注意：这里假设同期海水δ²⁶Mg = -0.83‰（现代值）
    # 实际研究中可能需要根据地质背景调整
    delta_seawater = -0.83
    df_a = calculate_weathering_ratios(df_a, delta_seawater)
    df_b = calculate_weathering_ratios(df_b, delta_seawater)
    print("  ✓ 风化比例计算完成")
    
    # 3. 反演风化通量
    print("\n[3] 反演风化通量...")
    df_a = inverse_weathering_flux(df_a, delta_seawater)
    df_b = inverse_weathering_flux(df_b, delta_seawater)
    print("  ✓ 通量估算完成")
    
    # 4. 打印统计结果
    print_statistics(df_a, df_b)
    
    # 5. 分析时间趋势
    analyze_temporal_trends(df_a, df_b)
    
    # 6. 绘制图表
    print("\n[4] 生成可视化图表...")
    plot_section_comparison(df_a, df_b)
    
    # 7. 保存结果
    print("\n[5] 保存分析结果...")
    df_a.to_csv('Section_A_weathering_results.csv', index=False)
    df_b.to_csv('Section_B_weathering_results.csv', index=False)
    print("  ✓ Section_A_weathering_results.csv")
    print("  ✓ Section_B_weathering_results.csv")
    
    print("\n" + "="*80)
    print("分析完成！")
    print("="*80)
    
    return df_a, df_b


if __name__ == '__main__':
    df_a, df_b = main()
