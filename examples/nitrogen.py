"""
氮同位素体系可视化示例

创建各种图表展示氮循环模型的结果:
1. f_assimilator vs δ¹⁵N_sed 关系曲线
2. 蒙特卡洛不确定性图
3. 地质时期对比图
4. 时间序列演化图
5. 箱式模型示意图
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.gridspec import GridSpec
import os

# 导入氮同位素体系
from systems.n import NIsotopeSystem, get_scenario_info, SCENARIO_PARAMETERS
from toolkit.math.statistics import Bootstrap, LOWESS, ChangepointDetector

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


def plot_1_forward_model_curve(save_path=None):
    """
    图1: 正向模型 - f_assimilator 与 δ¹⁵N_sed 关系曲线
    展示关键的非线性关系
    """
    print("\n正在生成图1: f_assimilator vs δ¹⁵N_sed 关系曲线...")
    
    n_system = NIsotopeSystem(scenario='modern')
    
    # 计算曲线
    curve = n_system.calculate_f_assimilator_curve(
        f_range=(0.0, 1.0), n_points=100, n_monte_carlo=2000
    )
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 绘制均值曲线
    ax.plot(curve['f_assimilator'], curve['delta15N_sed_mean'], 
            'b-', linewidth=2.5, label='Mean δ¹⁵N$_{sed}$')
    
    # 绘制置信区间
    ax.fill_between(curve['f_assimilator'], 
                    curve['delta15N_sed_ci95_lower'],
                    curve['delta15N_sed_ci95_upper'],
                    alpha=0.2, color='blue', label='95% CI')
    ax.fill_between(curve['f_assimilator'],
                    curve['delta15N_sed_ci68_lower'],
                    curve['delta15N_sed_ci68_upper'],
                    alpha=0.3, color='blue', label='68% CI')
    
    # 标记关键点
    max_idx = np.argmax(curve['delta15N_sed_mean'])
    f_max = curve['f_assimilator'][max_idx]
    delta_max = curve['delta15N_sed_mean'][max_idx]
    ax.plot(f_max, delta_max, 'r*', markersize=20, 
            label=f'Peak (f={f_max:.2f}, δ¹⁵N={delta_max:.1f}‰)')
    
    # 标记现代海洋值
    modern_f = 0.7
    modern_delta = n_system.forward_model(modern_f)
    ax.plot(modern_f, modern_delta, 'go', markersize=12,
            label=f'Modern Ocean (f={modern_f}, δ¹⁵N={modern_delta:.1f}‰)')
    
    # 标记早三叠世范围
    ax.axvspan(0.0, 0.1, alpha=0.2, color='red', label='Early Triassic Stage I')
    ax.axvspan(0.15, 0.25, alpha=0.2, color='green', label='Early Triassic Stage II')
    
    ax.set_xlabel('Nitrate Assimilator Fraction (f$_{assimilator}$)', fontsize=12)
    ax.set_ylabel('Sediment δ¹⁵N (‰)', fontsize=12)
    ax.set_title('Nitrogen Isotope Model: Forward Relationship\n' + 
                 r'Based on Kang et al. (2023) and Ma et al. (2025)', fontsize=13)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)
    ax.set_ylim(-2, 8)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"  已保存: {save_path}")
    
    plt.close()


def plot_2_monte_carlo_analysis(save_path=None):
    """
    图2: 蒙特卡洛不确定性分析
    展示不同 f_assimilator 值下的 δ¹⁵N 分布
    """
    print("\n正在生成图2: 蒙特卡洛不确定性分析...")
    
    n_system = NIsotopeSystem(scenario='neoproterozoic')
    
    # 测试不同的 f 值
    f_values = [0.05, 0.11, 0.20, 0.35, 0.50]
    f_labels = ['0.05\n(Anoxic)', '0.11\n(Neo-proterozoic\npre-800Ma)', 
                '0.20\n(Early Triassic\nStage II)', 
                '0.35\n(Neo-proterozoic\npost-800Ma)', 
                '0.50\n(Modern\nOxygenated)']
    colors = ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4', '#9467bd']
    
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.flatten()
    
    for idx, (f, label, color) in enumerate(zip(f_values, f_labels, colors)):
        ax = axes[idx]
        
        # 蒙特卡洛模拟
        mc_result = n_system.monte_carlo_simulation(f, n_samples=5000)
        samples = mc_result['samples']
        
        # 绘制直方图
        n_bins, bins, patches = ax.hist(samples, bins=50, density=True, 
                                        alpha=0.7, color=color, edgecolor='black')
        
        # 标记统计量
        mean_val = mc_result['delta15N_sed_mean']
        median_val = mc_result['delta15N_sed_median']
        ci68 = mc_result['delta15N_sed_ci68']
        
        ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, 
                  label=f'Mean: {mean_val:.2f}‰')
        ax.axvline(median_val, color='green', linestyle='-.', linewidth=2,
                  label=f'Median: {median_val:.2f}‰')
        ax.axvspan(ci68[0], ci68[1], alpha=0.2, color='yellow',
                  label=f'68% CI: [{ci68[0]:.1f}, {ci68[1]:.1f}]')
        
        ax.set_title(f'f = {label}', fontsize=10)
        ax.set_xlabel('δ¹⁵N$_{sed}$ (‰)', fontsize=9)
        ax.set_ylabel('Probability Density', fontsize=9)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
    
    # 最后一个子图显示汇总
    ax = axes[-1]
    ax.axis('off')
    
    summary_text = "Monte Carlo Analysis Summary\n" + "="*40 + "\n\n"
    for f in f_values:
        mc = n_system.monte_carlo_simulation(f, n_samples=1000)
        summary_text += f"f = {f:.2f}:\n"
        summary_text += f"  Mean: {mc['delta15N_sed_mean']:.2f}‰\n"
        summary_text += f"  Std:  {mc['delta15N_sed_std']:.2f}‰\n"
        summary_text += f"  95% CI: [{mc['delta15N_sed_ci95'][0]:.1f}, "
        summary_text += f"{mc['delta15N_sed_ci95'][1]:.1f}]‰\n\n"
    
    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes, fontsize=9,
           verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.suptitle('Monte Carlo Uncertainty Analysis\n' + 
                 r'($\varepsilon_{fix}$: -2‰ to +1‰, $\varepsilon_{wcd}$: -30‰ to -22‰)',
                 fontsize=14, y=1.02)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"  已保存: {save_path}")
    
    plt.close()


def plot_3_geological_comparison(save_path=None):
    """
    图3: 不同地质时期氮循环对比
    """
    print("\n正在生成图3: 地质时期氮循环对比...")
    
    scenarios = [
        ('anoxic_nitrate_depleted', 'Anoxic Ocean\n(Nitrate Depleted)', '#d62728'),
        ('early_triassic_stage_I', 'Early Triassic\nStage I', '#ff7f0e'),
        ('early_triassic_stage_II', 'Early Triassic\nStage II', '#2ca02c'),
        ('early_triassic_stage_III', 'Early Triassic\nStage III', '#ffbb78'),
        ('neoproterozoic_pre_800Ma', 'Neoproterozoic\npre-800 Ma', '#1f77b4'),
        ('neoproterozoic_post_800Ma', 'Neoproterozoic\npost-800 Ma', '#17becf'),
        ('modern_oxic', 'Modern Ocean\n(Oxic)', '#9467bd'),
    ]
    
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)
    
    # 子图1: f_assimilator 范围对比
    ax1 = fig.add_subplot(gs[0, 0])
    
    y_pos = np.arange(len(scenarios))
    f_ranges = [get_scenario_info(s[0])['f_assimilator_range'] for s in scenarios]
    f_means = [(r[0] + r[1])/2 for r in f_ranges]
    f_errors = [[(r[0]+r[1])/2 - r[0], r[1] - (r[0]+r[1])/2] for r in f_ranges]
    colors_list = [s[2] for s in scenarios]
    
    ax1.barh(y_pos, f_means, xerr=np.array(f_errors).T, 
            color=colors_list, alpha=0.7, edgecolor='black', capsize=5)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels([s[1] for s in scenarios], fontsize=9)
    ax1.set_xlabel('Nitrate Assimilator Fraction (f$_{assimilator}$)', fontsize=11)
    ax1.set_title('Nitrate Availability Comparison\nAcross Geological Periods', fontsize=12)
    ax1.axvline(x=0.7, color='red', linestyle='--', linewidth=2, label='Modern Reference (0.7)')
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='x')
    ax1.set_xlim(0, 1)
    
    # 子图2: δ¹⁵N 范围对比
    ax2 = fig.add_subplot(gs[0, 1])
    
    delta_ranges = [get_scenario_info(s[0])['delta15N_sed_range'] for s in scenarios]
    delta_means = [(r[0] + r[1])/2 for r in delta_ranges]
    delta_errors = [[(r[0]+r[1])/2 - r[0], r[1] - (r[0]+r[1])/2] for r in delta_ranges]
    
    ax2.barh(y_pos, delta_means, xerr=np.array(delta_errors).T,
            color=colors_list, alpha=0.7, edgecolor='black', capsize=5)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([s[1] for s in scenarios], fontsize=9)
    ax2.set_xlabel('Sediment δ¹⁵N (‰)', fontsize=11)
    ax2.set_title('Expected δ¹⁵N$_{sed}$ Range\n(Model Prediction)', fontsize=12)
    ax2.grid(True, alpha=0.3, axis='x')
    
    # 子图3: 时间序列示意图
    ax3 = fig.add_subplot(gs[1, :])
    
    # 模拟时间序列
    np.random.seed(42)
    ages = np.linspace(252, 247, 100)  # 早三叠世
    
    # 创建模拟趋势
    trend = np.piecewise(ages, 
        [ages > 250.5, (ages <= 250.5) & (ages > 248.8), ages <= 248.8],
        [lambda x: 1.0 + 0.5*np.sin((252-x)*2),  # Stage I
         lambda x: 3.0 + 1.0*np.sin((250.5-x)*3), # Stage II
         lambda x: 2.0 + 0.8*np.sin((248.8-x)*2)]) # Stage III
    
    noise = np.random.normal(0, 0.3, len(ages))
    delta15N_data = trend + noise
    
    ax3.fill_between(ages, delta15N_data - 0.5, delta15N_data + 0.5, 
                     alpha=0.3, color='gray', label='Uncertainty')
    ax3.plot(ages, delta15N_data, 'ko-', markersize=3, alpha=0.6, label='Observed δ¹⁵N')
    
    # 标记阶段
    ax3.axvspan(251.9, 250.6, alpha=0.2, color='red', label='Stage I (Griesbachian-Smithian)')
    ax3.axvspan(250.6, 248.8, alpha=0.2, color='green', label='Stage II (Spathian early)')
    ax3.axvspan(248.8, 247.2, alpha=0.2, color='orange', label='Stage III (Spathian late)')
    
    # 标记关键界线
    ax3.axvline(x=251.902, color='black', linestyle='--', linewidth=1.5)
    ax3.text(251.902, ax3.get_ylim()[1]*0.9, 'PTB', rotation=90, 
            verticalalignment='top', fontsize=9)
    ax3.axvline(x=250.6, color='black', linestyle='--', linewidth=1.5)
    ax3.text(250.6, ax3.get_ylim()[1]*0.9, 'SSB', rotation=90,
            verticalalignment='top', fontsize=9)
    
    ax3.set_xlabel('Age (Ma)', fontsize=11)
    ax3.set_ylabel('Sediment δ¹⁵N (‰)', fontsize=11)
    ax3.set_title('Simulated Early Triassic δ¹⁵N$_{sed}$ Record\n' + 
                  '(Based on Ma et al. 2025)', fontsize=12)
    ax3.legend(loc='lower right', fontsize=8)
    ax3.grid(True, alpha=0.3)
    ax3.invert_xaxis()
    
    plt.suptitle('Nitrogen Cycle Evolution Through Geological Time', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"  已保存: {save_path}")
    
    plt.close()


def plot_4_inverse_model_demo(save_path=None):
    """
    图4: 反向模型演示
    展示如何从观测的 δ¹⁵N 反演 f_assimilator
    """
    print("\n正在生成图4: 反向模型演示...")
    
    n_system = NIsotopeSystem(scenario='early_triassic')
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 左图: 正向曲线 + 反演示例
    ax1 = axes[0]
    
    f_range = np.linspace(0, 0.48, 100)  # 缺氧环境范围
    delta_forward = [n_system.forward_model(f) for f in f_range]
    
    ax1.plot(f_range, delta_forward, 'b-', linewidth=2.5, label='Forward Model')
    ax1.fill_between(f_range, np.array(delta_forward)-0.5, 
                    np.array(delta_forward)+0.5, alpha=0.2, color='blue')
    
    # 示例反演点
    example_observations = [0.5, 1.5, 3.0, 4.5, 2.0]
    example_colors = ['red', 'orange', 'green', 'blue', 'purple']
    
    for obs, color in zip(example_observations, example_colors):
        result = n_system.inverse_model(obs, f_range=(0, 0.48))
        f_inv = result['f_assimilator']
        
        # 绘制水平线
        ax1.axhline(y=obs, color=color, linestyle='--', alpha=0.5)
        # 绘制垂直线
        ax1.axvline(x=f_inv, color=color, linestyle='--', alpha=0.5)
        # 标记交点
        ax1.plot(f_inv, obs, 'o', color=color, markersize=10,
                label=f'δ¹⁵N={obs}‰ → f={f_inv:.3f}')
    
    ax1.set_xlabel('Nitrate Assimilator Fraction (f$_{assimilator}$)', fontsize=12)
    ax1.set_ylabel('Sediment δ¹⁵N (‰)', fontsize=12)
    ax1.set_title('Inverse Model Demonstration\n(Oxygen-Deficient Ocean Scenario)', fontsize=12)
    ax1.legend(loc='lower right', fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 0.5)
    ax1.set_ylim(-1, 6)
    
    # 右图: 反演结果汇总
    ax2 = axes[1]
    
    # 生成一系列观测值并反演
    obs_range = np.linspace(0, 5, 50)
    f_inverted = []
    
    for obs in obs_range:
        result = n_system.inverse_model(obs, f_range=(0, 0.48))
        f_inverted.append(result['f_assimilator'])
    
    ax2.plot(obs_range, f_inverted, 'g-', linewidth=2.5)
    ax2.fill_between(obs_range, 0, f_inverted, alpha=0.3, color='green')
    
    # 标记阶段区域
    ax2.axhspan(0, 0.1, alpha=0.2, color='red', label='Stage I (Anoxic)')
    ax2.axhspan(0.15, 0.25, alpha=0.2, color='green', label='Stage II (Oxic)')
    ax2.axhspan(0.05, 0.15, alpha=0.2, color='orange', label='Stage III (Re-anoxic)')
    
    ax2.set_xlabel('Observed δ¹⁵N$_{sed}$ (‰)', fontsize=12)
    ax2.set_ylabel('Inferred f$_{assimilator}$', fontsize=12)
    ax2.set_title('Inversion Result: δ¹⁵N → f$_{assimilator}$\n(Early Triassic Application)', fontsize=12)
    ax2.legend(loc='upper left', fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 5)
    ax2.set_ylim(0, 0.5)
    
    plt.suptitle('Nitrogen Isotope Inverse Model', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"  已保存: {save_path}")
    
    plt.close()


def plot_5_box_model_diagram(save_path=None):
    """
    图5: 箱式模型示意图
    """
    print("\n正在生成图5: 箱式模型示意图...")
    
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    
    # 标题
    ax.text(5, 9.5, 'Nitrogen Cycle Box Model', fontsize=18, 
           ha='center', fontweight='bold')
    ax.text(5, 9.0, 'Based on Kang et al. (2023) and Ma et al. (2025)', 
           fontsize=11, ha='center', style='italic')
    
    # 大气框
    atm_box = FancyBboxPatch((3.5, 7), 3, 1.2, boxstyle="round,pad=0.1",
                            edgecolor='black', facecolor='lightblue', linewidth=2)
    ax.add_patch(atm_box)
    ax.text(5, 7.6, 'Atmosphere N₂', fontsize=12, ha='center', fontweight='bold')
    ax.text(5, 7.2, 'δ¹⁵N = 0‰', fontsize=10, ha='center')
    
    # 固氮生物-铵储库
    fixer_box = FancyBboxPatch((0.5, 4), 3.5, 1.8, boxstyle="round,pad=0.1",
                              edgecolor='black', facecolor='lightyellow', linewidth=2)
    ax.add_patch(fixer_box)
    ax.text(2.25, 5.4, 'N$_{fixer}$/ammonium', fontsize=11, ha='center', fontweight='bold')
    ax.text(2.25, 5.0, '(Nitrogen Fixers)', fontsize=9, ha='center', style='italic')
    ax.text(2.25, 4.5, 'δ¹⁵N$_{ammonium}$ = 0‰ + ε$_{fix}$', fontsize=9, ha='center')
    ax.text(2.25, 4.2, 'ε$_{fix}$ = -2‰ to +1‰', fontsize=8, ha='center', color='blue')
    
    # 硝酸盐同化生物-硝酸盐储库
    assim_box = FancyBboxPatch((6, 4), 3.5, 1.8, boxstyle="round,pad=0.1",
                              edgecolor='black', facecolor='lightgreen', linewidth=2)
    ax.add_patch(assim_box)
    ax.text(7.75, 5.4, 'N$_{assimilator}$/nitrate', fontsize=11, ha='center', fontweight='bold')
    ax.text(7.75, 5.0, '(Nitrate Assimilators)', fontsize=9, ha='center', style='italic')
    ax.text(7.75, 4.5, 'δ¹⁵N$_{nitrate}$ = δ¹⁵N$_{ammonium}$ - Δ', fontsize=9, ha='center')
    ax.text(7.75, 4.2, 'Δ = (F$_{wcd}$×ε$_{wcd}$ + F$_{sd}$×ε$_{sd}$)/F$_{remin}$', 
           fontsize=7, ha='center', color='blue')
    
    # 沉积物
    sed_box = FancyBboxPatch((3.5, 1), 3, 1.5, boxstyle="round,pad=0.1",
                            edgecolor='black', facecolor='wheat', linewidth=2)
    ax.add_patch(sed_box)
    ax.text(5, 2.1, 'Sediment', fontsize=11, ha='center', fontweight='bold')
    ax.text(5, 1.7, 'δ¹⁵N$_{sed}$ = (1-f)×δ¹⁵N$_{ammonium}$', fontsize=9, ha='center')
    ax.text(5, 1.35, '+ f×δ¹⁵N$_{nitrate}$', fontsize=9, ha='center')
    
    # 箭头 - 固氮作用
    arrow1 = FancyArrowPatch((4.5, 7), (2.5, 5.8),
                            arrowstyle='->', mutation_scale=25, 
                            linewidth=2.5, color='darkgreen')
    ax.add_patch(arrow1)
    ax.text(3.2, 6.6, 'F$_{fix}$\nN fixation', fontsize=9, ha='center', color='darkgreen')
    ax.text(3.2, 6.2, 'ε$_{fix}$', fontsize=8, ha='center', color='darkgreen')
    
    # 箭头 - 再矿化
    arrow2 = FancyArrowPatch((4, 4.9), (6, 4.9),
                            arrowstyle='->', mutation_scale=20,
                            linewidth=2, color='purple')
    ax.add_patch(arrow2)
    ax.text(5, 5.3, 'F$_{remin}$', fontsize=9, ha='center', color='purple')
    
    # 箭头 - 水柱反硝化
    arrow3 = FancyArrowPatch((7.5, 4), (7.5, 2.5),
                            arrowstyle='->', mutation_scale=20,
                            linewidth=2, color='red')
    ax.add_patch(arrow3)
    ax.text(8.3, 3.2, 'F$_{wcd}$\nWater column\ndenitrification', 
           fontsize=8, ha='center', color='red')
    ax.text(8.3, 2.7, 'ε$_{wcd}$=-30‰~-22‰', fontsize=7, ha='center', color='red')
    
    # 箭头 - 沉积反硝化
    arrow4 = FancyArrowPatch((6.5, 4), (5, 2.5),
                            arrowstyle='->', mutation_scale=20,
                            linewidth=2, color='orange', linestyle='--')
    ax.add_patch(arrow4)
    ax.text(5.5, 3.4, 'F$_{sd}$', fontsize=8, ha='center', color='orange')
    ax.text(5.5, 3.1, 'ε$_{sd}$=0‰', fontsize=7, ha='center', color='orange')
    
    # 箭头 - 埋藏
    arrow5 = FancyArrowPatch((2.5, 4), (4, 2.5),
                            arrowstyle='->', mutation_scale=20,
                            linewidth=2, color='brown')
    ax.add_patch(arrow5)
    ax.text(2.8, 3.2, 'F$_{fixer}$', fontsize=8, ha='center', color='brown')
    
    arrow6 = FancyArrowPatch((7.5, 4), (6, 2.5),
                            arrowstyle='->', mutation_scale=20,
                            linewidth=2, color='brown')
    ax.add_patch(arrow6)
    ax.text(7.2, 3.2, 'F$_{assimilator}$', fontsize=8, ha='center', color='brown')
    
    # 说明框
    info_box = FancyBboxPatch((0.2, 0.1), 9.6, 0.7, boxstyle="round,pad=0.05",
                             edgecolor='black', facecolor='white', linewidth=1.5)
    ax.add_patch(info_box)
    ax.text(5, 0.6, 'f$_{assimilator}$ = F$_{assimilator}$/(F$_{assimilator}$+F$_{fixer}$) = Nitrate availability indicator',
           fontsize=10, ha='center', fontweight='bold')
    ax.text(5, 0.25, 'Key relationship: f$_{assimilator}$↑ → F$_{wcd}$↓ → δ¹⁵N$_{nitrate}$↓ → δ¹⁵N$_{sed}$ peak at f ≈ 0.48',
           fontsize=9, ha='center', style='italic')
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"  已保存: {save_path}")
    
    plt.close()


def main():
    """生成所有图表"""
    print("=" * 70)
    print("氮同位素体系可视化示例")
    print("=" * 70)
    
    # 创建输出目录
    output_dir = "examples/output"
    os.makedirs(output_dir, exist_ok=True)
    
    # 生成所有图表
    plot_1_forward_model_curve(f"{output_dir}/01_forward_model_curve.png")
    plot_2_monte_carlo_analysis(f"{output_dir}/02_monte_carlo_analysis.png")
    plot_3_geological_comparison(f"{output_dir}/03_geological_comparison.png")
    plot_4_inverse_model_demo(f"{output_dir}/04_inverse_model_demo.png")
    plot_5_box_model_diagram(f"{output_dir}/05_box_model_diagram.png")
    
    print("\n" + "=" * 70)
    print("所有图表生成完成!")
    print(f"输出目录: {output_dir}/")
    print("=" * 70)
    print("\n生成的文件:")
    for f in sorted(os.listdir(output_dir)):
        if f.endswith('.png'):
            print(f"  - {f}")


if __name__ == '__main__':
    main()
