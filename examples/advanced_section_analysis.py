#!/usr/bin/env python3
"""
Nie Section A & B - 高级Mg同位素分析

考虑因素：
1. 白云石化作用（dolomitization）作为第三端元
2. 可变的古海水Mg同位素组成
3. 蒸发岩沉积环境
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from systems.mg import MgIsotopeSystem


def three_endmember_mixing(delta_sample: float, 
                          delta_carb: float = -4.3,
                          delta_sil: float = -0.3,
                          delta_dolo: float = -2.5) -> dict:
    """
    三端元混合模型（碳酸盐 + 硅酸盐 + 白云石）
    
    Parameters
    ----------
    delta_sample : float
        样品δ²⁶Mg
    delta_carb, delta_sil, delta_dolo : float
        各端元值
        
    Returns
    -------
    dict
        各端元比例
    """
    # 使用最小二乘法求解三端元混合
    def objective(f):
        f_carb, f_sil, f_dolo = f
        if f_carb < 0 or f_sil < 0 or f_dolo < 0:
            return 1e10
        if abs(f_carb + f_sil + f_dolo - 1) > 0.01:
            return 1e10
        
        delta_calc = f_carb * delta_carb + f_sil * delta_sil + f_dolo * delta_dolo
        return (delta_calc - delta_sample)**2
    
    # 初始猜测
    result = minimize(objective, [0.33, 0.33, 0.34], 
                     method='L-BFGS-B',
                     bounds=[(0, 1), (0, 1), (0, 1)])
    
    if result.success:
        f_carb, f_sil, f_dolo = result.x
        # 归一化
        total = f_carb + f_sil + f_dolo
        return {
            'f_carbonate': f_carb / total,
            'f_silicate': f_sil / total,
            'f_dolomite': f_dolo / total
        }
    else:
        return {'f_carbonate': 0, 'f_silicate': 1, 'f_dolomite': 0}


def analyze_with_variable_seawater(df: pd.DataFrame, 
                                   seawater_range=(-2.0, 0.0)) -> pd.DataFrame:
    """
    尝试不同的古海水Mg同位素值来拟合数据
    
    假设：样品是海水与沉积物端元的混合
    """
    results = []
    
    for delta_sw in np.linspace(seawater_range[0], seawater_range[1], 21):
        # 计算在此海水值下的风化比例
        mg = MgIsotopeSystem()
        
        f_carb_list = []
        for _, row in df.iterrows():
            delta_sample = row['delta_26_Mg_iso']
            
            # 简化的二元混合：样品 = 海水 + 沉积物输入
            # 假设沉积物输入介于碳酸盐和硅酸盐之间
            delta_sediment = -2.0  # 平均沉积物
            
            # 质量平衡：delta_sample = f_sw * delta_sw + (1-f_sw) * delta_sediment
            if abs(delta_sediment - delta_sw) > 0.001:
                f_sw = (delta_sample - delta_sediment) / (delta_sw - delta_sediment)
                f_sw = np.clip(f_sw, 0, 1)
            else:
                f_sw = 0.5
            
            f_carb_list.append(1 - f_sw)  # 沉积物比例
        
        # 计算残差
        residuals = np.array(f_carb_list) - 0.5  # 假设理想情况下沉积物占50%
        rmse = np.sqrt(np.mean(residuals**2))
        
        results.append({
            'delta_seawater': delta_sw,
            'mean_f_sediment': np.mean(f_carb_list),
            'rmse': rmse
        })
    
    results_df = pd.DataFrame(results)
    best_fit = results_df.loc[results_df['rmse'].idxmin()]
    
    return results_df, best_fit


def analyze_dolomitization_effect(df: pd.DataFrame) -> pd.DataFrame:
    """
    分析白云石化作用的影响
    
    白云石化分馏系数约为 -0.5 ~ -1.5‰
    """
    mg_system = MgIsotopeSystem()
    
    # 白云石化分馏参数
    epsilon_dolo = -1.0  # 白云石化分馏（‰）
    delta_fluid = -0.83  # 流体Mg同位素
    
    # 白云石Mg同位素 = 流体Mg + 分馏
    delta_dolomite = delta_fluid + epsilon_dolo
    
    print(f"\n白云石化分析:")
    print(f"  假设流体δ²⁶Mg = {delta_fluid}‰ (现代海水)")
    print(f"  白云石化分馏ε = {epsilon_dolo}‰")
    print(f"  理论白云石δ²⁶Mg = {delta_dolomite}‰")
    
    # 计算每个样品需要的白云石化程度
    df['dolomitization_degree'] = np.nan
    
    for i, row in df.iterrows():
        delta_sample = row['delta_26_Mg_iso']
        
        # 假设原始灰岩δ²⁶Mg ≈ -4‰
        delta_limestone = -4.0
        
        # 混合模型：样品 = f_dolo * delta_dolomite + (1-f_dolo) * delta_limestone
        if abs(delta_dolomite - delta_limestone) > 0.001:
            f_dolo = (delta_sample - delta_limestone) / (delta_dolomite - delta_limestone)
            f_dolo = np.clip(f_dolo, 0, 1)
            df.loc[i, 'dolomitization_degree'] = f_dolo
    
    return df


def create_comprehensive_plot(df_a: pd.DataFrame, df_b: pd.DataFrame, 
                              best_sw_a, best_sw_b):
    """创建综合分析图"""
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # 1. 原始数据对比
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(df_a['delta_26_Mg_iso'], bins=10, alpha=0.5, label='Section A', color='blue')
    ax1.hist(df_b['delta_26_Mg_iso'], bins=10, alpha=0.5, label='Section B', color='red')
    ax1.axvline(x=-0.83, color='black', linestyle='--', label='Modern Seawater')
    ax1.axvline(x=-4.3, color='green', linestyle=':', alpha=0.5, label='Carbonate')
    ax1.axvline(x=-0.3, color='orange', linestyle=':', alpha=0.5, label='Silicate')
    ax1.axvline(x=-2.5, color='purple', linestyle=':', alpha=0.5, label='Dolomite')
    ax1.set_xlabel('δ²⁶Mg (‰)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Mg Isotope Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. 三端元混合 - Section A
    ax2 = fig.add_subplot(gs[0, 1])
    # 计算三端元比例
    f_carb_a, f_sil_a, f_dolo_a = [], [], []
    for _, row in df_a.iterrows():
        ratios = three_endmember_mixing(row['delta_26_Mg_iso'])
        f_carb_a.append(ratios['f_carbonate'])
        f_sil_a.append(ratios['f_silicate'])
        f_dolo_a.append(ratios['f_dolomite'])
    
    x = range(len(df_a))
    ax2.stackplot(x, f_carb_a, f_sil_a, f_dolo_a,
                  labels=['Carbonate', 'Silicate', 'Dolomite'],
                  colors=['lightblue', 'lightcoral', 'plum'], alpha=0.7)
    ax2.set_xlabel('Sample ID')
    ax2.set_ylabel('Proportion')
    ax2.set_title('Section A: Three-endmember Mixing')
    ax2.set_ylim(0, 1)
    ax2.legend(loc='upper left')
    ax2.grid(True, alpha=0.3)
    
    # 3. 三端元混合 - Section B
    ax3 = fig.add_subplot(gs[0, 2])
    f_carb_b, f_sil_b, f_dolo_b = [], [], []
    for _, row in df_b.iterrows():
        ratios = three_endmember_mixing(row['delta_26_Mg_iso'])
        f_carb_b.append(ratios['f_carbonate'])
        f_sil_b.append(ratios['f_silicate'])
        f_dolo_b.append(ratios['f_dolomite'])
    
    x = range(len(df_b))
    ax3.stackplot(x, f_carb_b, f_sil_b, f_dolo_b,
                  labels=['Carbonate', 'Silicate', 'Dolomite'],
                  colors=['lightblue', 'lightcoral', 'plum'], alpha=0.7)
    ax3.set_xlabel('Sample ID')
    ax3.set_ylabel('Proportion')
    ax3.set_title('Section B: Three-endmember Mixing')
    ax3.set_ylim(0, 1)
    ax3.legend(loc='upper left')
    ax3.grid(True, alpha=0.3)
    
    # 4. 白云石化程度
    ax4 = fig.add_subplot(gs[1, 0])
    df_a_dolo = analyze_dolomitization_effect(df_a.copy())
    df_b_dolo = analyze_dolomitization_effect(df_b.copy())
    
    ax4.scatter(df_a['sample_id'], df_a_dolo['dolomitization_degree'], 
               c='blue', label='Section A', alpha=0.6)
    ax4.scatter(df_b['sample_id'], df_b_dolo['dolomitization_degree'], 
               c='red', label='Section B', alpha=0.6)
    ax4.set_xlabel('Sample ID')
    ax4.set_ylabel('Dolomitization Degree')
    ax4.set_title('Estimated Dolomitization Degree')
    ax4.set_ylim(0, 1)
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # 5. 古海水Mg同位素敏感性分析 - Section A
    ax5 = fig.add_subplot(gs[1, 1])
    results_a, _ = analyze_with_variable_seawater(df_a)
    ax5.plot(results_a['delta_seawater'], results_a['rmse'], 'b-', label='RMSE')
    best_sw_a_val = best_sw_a['delta_seawater']
    ax5.axvline(x=best_sw_a_val, color='red', linestyle='--', 
               label=f'Best fit: {best_sw_a_val:.2f}‰')
    ax5.set_xlabel('Assumed Seawater δ²⁶Mg (‰)')
    ax5.set_ylabel('RMSE')
    ax5.set_title('Section A: Seawater Mg Sensitivity')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # 6. 古海水Mg同位素敏感性分析 - Section B
    ax6 = fig.add_subplot(gs[1, 2])
    results_b, _ = analyze_with_variable_seawater(df_b)
    ax6.plot(results_b['delta_seawater'], results_b['rmse'], 'b-', label='RMSE')
    best_sw_b_val = best_sw_b['delta_seawater']
    ax6.axvline(x=best_sw_b_val, color='red', linestyle='--',
               label=f'Best fit: {best_sw_b_val:.2f}‰')
    ax6.set_xlabel('Assumed Seawater δ²⁶Mg (‰)')
    ax6.set_ylabel('RMSE')
    ax6.set_title('Section B: Seawater Mg Sensitivity')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    # 7. 散点图：Mg vs Sample
    ax7 = fig.add_subplot(gs[2, :])
    ax7.errorbar(df_a['sample_id'], df_a['delta_26_Mg_iso'],
                yerr=df_a['delta_26_Mg_iso_2sd']/2,
                fmt='o', label='Section A', capsize=3, color='blue', alpha=0.7)
    ax7.errorbar(df_b['sample_id'] + 35, df_b['delta_26_Mg_iso'],
                yerr=df_b['delta_26_Mg_iso_2sd']/2,
                fmt='s', label='Section B', capsize=3, color='red', alpha=0.7)
    
    # 添加端元参考线
    ax7.axhline(y=-4.3, color='green', linestyle='--', alpha=0.5, label='Carbonate end-member')
    ax7.axhline(y=-0.3, color='orange', linestyle='--', alpha=0.5, label='Silicate end-member')
    ax7.axhline(y=-2.5, color='purple', linestyle='--', alpha=0.5, label='Dolomite')
    ax7.axhline(y=-0.83, color='gray', linestyle=':', alpha=0.5, label='Modern seawater')
    
    ax7.set_xlabel('Sample ID')
    ax7.set_ylabel('δ²⁶Mg (‰)')
    ax7.set_title('Mg Isotope Compositions with End-members')
    ax7.legend(loc='upper right', ncol=2)
    ax7.grid(True, alpha=0.3)
    
    plt.savefig('Advanced_Mg_Analysis.png', dpi=300, bbox_inches='tight')
    print("\n高级分析图表已保存: Advanced_Mg_Analysis.png")
    plt.close()


def print_advanced_statistics(df_a: pd.DataFrame, df_b: pd.DataFrame,
                             best_sw_a, best_sw_b):
    """打印高级统计结果"""
    print("\n" + "="*80)
    print("高级分析结果")
    print("="*80)
    
    # 三端元混合结果
    print("\n【三端元混合分析】")
    print("假设端元:")
    print("  碳酸盐: δ²⁶Mg = -4.3‰")
    print("  硅酸盐: δ²⁶Mg = -0.3‰")
    print("  白云石: δ²⁶Mg = -2.5‰")
    
    for name, df in [('Section A', df_a), ('Section B', df_b)]:
        f_carb, f_sil, f_dolo = [], [], []
        for _, row in df.iterrows():
            ratios = three_endmember_mixing(row['delta_26_Mg_iso'])
            f_carb.append(ratios['f_carbonate'])
            f_sil.append(ratios['f_silicate'])
            f_dolo.append(ratios['f_dolomite'])
        
        print(f"\n{name}:")
        print(f"  碳酸盐比例: {np.mean(f_carb):.1%} ± {np.std(f_carb):.1%}")
        print(f"  硅酸盐比例: {np.mean(f_sil):.1%} ± {np.std(f_sil):.1%}")
        print(f"  白云石比例: {np.mean(f_dolo):.1%} ± {np.std(f_dolo):.1%}")
    
    # 白云石化分析
    print("\n【白云石化分析】")
    print("假设白云石δ²⁶Mg = -1.83‰ (海水-0.83‰ + 分馏-1.0‰)")
    
    df_a_dolo = analyze_dolomitization_effect(df_a.copy())
    df_b_dolo = analyze_dolomitization_effect(df_b.copy())
    
    print(f"\nSection A 估算白云石化程度: {df_a_dolo['dolomitization_degree'].mean():.1%}")
    print(f"Section B 估算白云石化程度: {df_b_dolo['dolomitization_degree'].mean():.1%}")
    
    # 古海水Mg同位素
    print("\n【古海水Mg同位素估算】")
    print(f"Section A 最佳拟合海水δ²⁶Mg: {best_sw_a['delta_seawater']:.2f}‰")
    print(f"Section B 最佳拟合海水δ²⁶Mg: {best_sw_b['delta_seawater']:.2f}‰")
    print("\n注：如果最佳拟合值与现代海水(-0.83‰)差异大，")
    print("    可能指示古海水Mg同位素组成与现今不同")


def main():
    """主函数"""
    print("="*80)
    print("Nie Section A & B - 高级Mg同位素分析")
    print("="*80)
    
    # 加载数据
    print("\n[1] 加载数据...")
    df_a = pd.read_excel('data/Nie_Section_A.xlsx')
    df_b = pd.read_excel('data/Nie_Section_B.xlsx')
    df_a['sample_id'] = range(1, len(df_a) + 1)
    df_b['sample_id'] = range(1, len(df_b) + 1)
    print(f"  Section A: {len(df_a)} 个样品")
    print(f"  Section B: {len(df_b)} 个样品")
    
    # 古海水Mg同位素敏感性分析
    print("\n[2] 分析古海水Mg同位素敏感性...")
    results_a, best_sw_a = analyze_with_variable_seawater(df_a)
    results_b, best_sw_b = analyze_with_variable_seawater(df_b)
    print(f"  Section A 最佳拟合: δ²⁶Mg_seawater = {best_sw_a['delta_seawater']:.2f}‰")
    print(f"  Section B 最佳拟合: δ²⁶Mg_seawater = {best_sw_b['delta_seawater']:.2f}‰")
    
    # 打印统计结果
    print_advanced_statistics(df_a, df_b, best_sw_a, best_sw_b)
    
    # 生成图表
    print("\n[3] 生成综合分析图表...")
    create_comprehensive_plot(df_a, df_b, best_sw_a, best_sw_b)
    
    print("\n" + "="*80)
    print("高级分析完成！")
    print("="*80)


if __name__ == '__main__':
    main()
