"""
硅酸盐风化通量模拟 - 基于碎屑组分Mg同位素
基于 Hu et al. (2023) 框架

核心功能：
1. 从碎屑沉积物 δ²⁶Mg 反演风化程度
2. 计算硅酸盐风化通量 (F_sili)
3. 评估硅酸盐风化强度指数 (SWI)
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple
import matplotlib.pyplot as plt


@dataclass
class SilicateWeatheringParams:
    """硅酸盐风化模拟参数"""
    # 同位素端元 (‰)
    d26Mg_UCC: float = -0.25           # 上地壳基准值
    d26Mg_carbonate: float = -2.0       # 碳酸盐风化端元
    d26Mg_river_water: float = -1.14    # 河水 δ²⁶Mg (实测)
    
    # 分馏因子 (‰)
    Delta_fluid_protolith: float = -0.50  # 流体与原岩分馏 (范围: -0.60 ~ -0.40)
    Delta_release_clay: float = -0.45     # 释放Mg与黏土的分馏
    
    # 通量参数 (mol/yr)
    F_river_total: float = 30e10        # 河流总Mg通量 (长江: 30×10¹⁰ mol/yr)
    
    # 不确定性
    uncertainty_clay: float = 0.07      # 碎屑 δ²⁶Mg 分析误差


class SilicateWeatheringModel:
    """
    硅酸盐风化通量模型
    
    从碎屑沉积物（黏土）Mg同位素反演硅酸盐风化贡献
    """
    
    def __init__(self, params: Optional[SilicateWeatheringParams] = None):
        self.params = params or SilicateWeatheringParams()
    
    def calculate_weathering_degree(self, d26Mg_clay: float) -> float:
        """
        从碎屑 δ²⁶Mg 计算风化程度 (f_Mg)
        
        使用 Rayleigh 分馏模型:
        δ²⁶Mg_clay = δ²⁶Mg_UCC - Δ·ln(f_Mg)
        
        Parameters:
            d26Mg_clay: 碎屑黏土 δ²⁶Mg (‰)
            
        Returns:
            f_Mg: 保留在碎屑中的Mg比例 (0-1)
                  越小表示风化程度越高
        """
        p = self.params
        
        # 检查输入合理性
        if d26Mg_clay < p.d26Mg_UCC:
            raise ValueError(f"碎屑δ²⁶Mg ({d26Mg_clay}‰) 不应低于UCC ({p.d26Mg_UCC}‰)")
        
        # Rayleigh 分馏反演
        numerator = p.d26Mg_UCC - d26Mg_clay
        denominator = p.Delta_fluid_protolith
        
        f_Mg = np.exp(numerator / denominator)
        
        return f_Mg
    
    def calculate_silicate_endmember(self, d26Mg_clay: float) -> float:
        """
        计算硅酸盐风化端元 δ²⁶Mg_sili
        
        假设：释放到流体中的Mg与残余黏土存在分馏
        δ²⁶Mg_sili = δ²⁶Mg_clay + Δ_release_clay
        
        Parameters:
            d26Mg_clay: 碎屑黏土 δ²⁶Mg (‰)
            
        Returns:
            δ²⁶Mg_sili: 硅酸盐风化端元 (‰)
        """
        return d26Mg_clay + self.params.Delta_release_clay
    
    def calculate_silicate_flux(self, d26Mg_clay: float) -> Tuple[float, float, float]:
        """
        计算硅酸盐风化通量
        
        使用双端元质量平衡:
        F_riv = F_sili + F_car
        F_riv·δ²⁶Mg_riv = F_sili·δ²⁶Mg_sili + F_car·δ²⁶Mg_car
        
        Parameters:
            d26Mg_clay: 碎屑黏土 δ²⁶Mg (‰)
            
        Returns:
            (F_sili, F_car, SWI): 
                - 硅酸盐风化通量 (mol/yr)
                - 碳酸盐风化通量 (mol/yr)
                - 硅酸盐风化强度指数 (%)
        """
        p = self.params
        
        # 计算硅酸盐风化端元
        d26Mg_sili = self.calculate_silicate_endmember(d26Mg_clay)
        
        # 检查端元合理性（硅酸盐应重于碳酸盐）
        if d26Mg_sili <= p.d26Mg_carbonate:
            raise ValueError(
                f"硅酸盐端元({d26Mg_sili}‰)必须重于碳酸盐({p.d26Mg_carbonate}‰)"
            )
        
        if d26Mg_sili <= p.d26Mg_river_water:
            raise ValueError(
                f"硅酸盐端元({d26Mg_sili}‰)必须重于河水({p.d26Mg_river_water}‰)"
            )
        
        # 求解质量平衡方程组
        # F_sili = F_riv * (δ²⁶Mg_riv - δ²⁶Mg_car) / (δ²⁶Mg_sili - δ²⁶Mg_car)
        numerator = p.d26Mg_river_water - p.d26Mg_carbonate
        denominator = d26Mg_sili - p.d26Mg_carbonate
        
        ratio = numerator / denominator
        F_sili = p.F_river_total * ratio
        F_car = p.F_river_total - F_sili
        
        # 硅酸盐风化强度指数
        SWI = ratio * 100  # 百分比
        
        return F_sili, F_car, SWI
    
    def forward_model(self, d26Mg_clay_range: np.ndarray) -> dict:
        """
        正演模型：给定碎屑 δ²⁶Mg 范围，计算所有输出
        
        Parameters:
            d26Mg_clay_range: 碎屑 δ²⁶Mg 数组 (‰)
            
        Returns:
            包含所有计算结果的字典
        """
        results = {
            'd26Mg_clay': d26Mg_clay_range,
            'f_Mg': [],
            'd26Mg_sili': [],
            'F_sili': [],
            'F_car': [],
            'SWI': []
        }
        
        for d26Mg in d26Mg_clay_range:
            try:
                f_Mg = self.calculate_weathering_degree(d26Mg)
                d26Mg_sili = self.calculate_silicate_endmember(d26Mg)
                F_sili, F_car, SWI = self.calculate_silicate_flux(d26Mg)
                
                results['f_Mg'].append(f_Mg)
                results['d26Mg_sili'].append(d26Mg_sili)
                results['F_sili'].append(F_sili)
                results['F_car'].append(F_car)
                results['SWI'].append(SWI)
            except ValueError as e:
                # 超出合理范围，填充NaN
                results['f_Mg'].append(np.nan)
                results['d26Mg_sili'].append(np.nan)
                results['F_sili'].append(np.nan)
                results['F_car'].append(np.nan)
                results['SWI'].append(np.nan)
        
        # 转换为数组
        for key in ['f_Mg', 'd26Mg_sili', 'F_sili', 'F_car', 'SWI']:
            results[key] = np.array(results[key])
        
        return results
    
    def monte_carlo_uncertainty(
        self, 
        d26Mg_clay: float,
        n_iterations: int = 10000
    ) -> dict:
        """
        Monte Carlo 不确定性分析
        
        Parameters:
            d26Mg_clay: 观测的碎屑 δ²⁶Mg
            n_iterations: 迭代次数
            
        Returns:
            包含统计结果的字典
        """
        p = self.params
        
        # 生成随机参数
        np.random.seed(42)
        d26Mg_clay_samples = np.random.normal(d26Mg_clay, p.uncertainty_clay, n_iterations)
        Delta_samples = np.random.uniform(-0.60, -0.40, n_iterations)  # 分馏因子范围
        d26Mg_river_samples = np.random.normal(p.d26Mg_river_water, 0.15, n_iterations)
        
        F_sili_samples = []
        SWI_samples = []
        
        for i in range(n_iterations):
            # 临时修改参数
            original_Delta = p.Delta_fluid_protolith
            original_d26Mg_river = p.d26Mg_river_water
            
            p.Delta_fluid_protolith = Delta_samples[i]
            p.d26Mg_river_water = d26Mg_river_samples[i]
            
            try:
                F_sili, _, SWI = self.calculate_silicate_flux(d26Mg_clay_samples[i])
                F_sili_samples.append(F_sili)
                SWI_samples.append(SWI)
            except ValueError:
                pass
            
            # 恢复参数
            p.Delta_fluid_protolith = original_Delta
            p.d26Mg_river_water = original_d26Mg_river
        
        F_sili_samples = np.array(F_sili_samples)
        SWI_samples = np.array(SWI_samples)
        
        return {
            'F_sili_mean': np.mean(F_sili_samples),
            'F_sili_std': np.std(F_sili_samples),
            'F_sili_95ci': np.percentile(F_sili_samples, [2.5, 97.5]),
            'SWI_mean': np.mean(SWI_samples),
            'SWI_std': np.std(SWI_samples),
            'SWI_95ci': np.percentile(SWI_samples, [2.5, 97.5]),
            'F_sili_distribution': F_sili_samples,
            'SWI_distribution': SWI_samples
        }


def plot_weathering_model_results(results: dict, save_path: Optional[str] = None):
    """
    可视化风化模型结果
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. 碎屑 δ²⁶Mg vs 风化程度
    ax = axes[0, 0]
    ax.plot(results['d26Mg_clay'], results['f_Mg'] * 100, 'b-', linewidth=2)
    ax.set_xlabel('碎屑黏土 δ²⁶Mg (‰)', fontsize=12)
    ax.set_ylabel('保留Mg比例 f_Mg (%)', fontsize=12)
    ax.set_title('风化程度 vs 碎屑同位素组成', fontsize=13)
    ax.axhline(y=100, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=-0.25, color='gray', linestyle='--', alpha=0.5, label='UCC')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # 2. 碎屑 δ²⁶Mg vs 硅酸盐风化通量
    ax = axes[0, 1]
    ax.plot(results['d26Mg_clay'], np.array(results['F_sili']) / 1e10, 'r-', linewidth=2, label='$F_{sili}$')
    ax.plot(results['d26Mg_clay'], np.array(results['F_car']) / 1e10, 'g--', linewidth=2, label='$F_{car}$')
    ax.set_xlabel('碎屑黏土 δ²⁶Mg (‰)', fontsize=12)
    ax.set_ylabel('Mg 通量 (×10¹⁰ mol/yr)', fontsize=12)
    ax.set_title('风化通量分解', fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 3. 硅酸盐风化强度指数 (SWI)
    ax = axes[1, 0]
    ax.plot(results['d26Mg_clay'], results['SWI'], 'purple', linewidth=2)
    ax.fill_between(results['d26Mg_clay'], 0, results['SWI'], alpha=0.3, color='purple')
    ax.set_xlabel('碎屑黏土 δ²⁶Mg (‰)', fontsize=12)
    ax.set_ylabel('SWI (%)', fontsize=12)
    ax.set_title('硅酸盐风化强度指数', fontsize=13)
    ax.axhline(y=50, color='gray', linestyle='--', alpha=0.5)
    ax.grid(True, alpha=0.3)
    
    # 4. 端元混合图
    ax = axes[1, 1]
    # 混合线
    x_mix = np.linspace(0, 100, 100)
    d26Mg_mix = (
        x_mix/100 * results['d26Mg_sili'][len(results['d26Mg_sili'])//2] + 
        (1-x_mix/100) * (-2.0)
    )
    ax.plot(x_mix, d26Mg_mix, 'k-', linewidth=1, alpha=0.5)
    
    # 端元
    ax.scatter([0, 100], [-2.0, results['d26Mg_sili'][len(results['d26Mg_sili'])//2]], 
               c=['green', 'red'], s=100, zorder=5)
    ax.text(0, -2.1, '碳酸盐\n风化', ha='center', fontsize=10)
    ax.text(100, results['d26Mg_sili'][len(results['d26Mg_sili'])//2]-0.1, '硅酸盐\n风化', 
            ha='center', fontsize=10)
    
    # 河水观测
    ax.scatter([50], [-1.14], c='blue', s=150, marker='*', zorder=5, label='河水实测')
    ax.set_xlabel('硅酸盐贡献 (%)', fontsize=12)
    ax.set_ylabel('混合 δ²⁶Mg (‰)', fontsize=12)
    ax.set_title('双端元混合模型', fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.show()


def example_usage():
    """
    示例：使用 Hu et al. (2023) 长江数据
    """
    print("=" * 60)
    print("碎屑组分 Mg 同位素 - 硅酸盐风化通量模拟")
    print("基于 Hu et al. (2023) 长江数据")
    print("=" * 60)
    
    # 初始化模型
    params = SilicateWeatheringParams(
        d26Mg_UCC=-0.25,
        d26Mg_carbonate=-2.0,
        d26Mg_river_water=-1.14,
        Delta_fluid_protolith=-0.50,
        F_river_total=30e10
    )
    
    model = SilicateWeatheringModel(params)
    
    # 观测的碎屑黏土 δ²⁶Mg 范围 (Hu et al.)
    d26Mg_clay_observed = np.array([-0.15, -0.10, -0.05, 0.00])
    
    print("\n【观测数据】")
    print(f"碎屑黏土 δ²⁶Mg 范围: {d26Mg_clay_observed.min():.2f}‰ ~ {d26Mg_clay_observed.max():.2f}‰")
    
    print("\n【计算结果】")
    print("-" * 60)
    print(f"{'δ²⁶Mg_黏土':<12} {'f_Mg':<10} {'δ²⁶Mg_sili':<12} {'F_sili':<15} {'SWI':<8}")
    print("-" * 60)
    
    for d26Mg in d26Mg_clay_observed:
        try:
            f_Mg = model.calculate_weathering_degree(d26Mg)
            d26Mg_sili = model.calculate_silicate_endmember(d26Mg)
            F_sili, F_car, SWI = model.calculate_silicate_flux(d26Mg)
            
            print(f"{d26Mg:<12.2f} {f_Mg*100:<10.1f} {d26Mg_sili:<12.2f} "
                  f"{F_sili/1e10:<15.1f} {SWI:<8.1f}")
        except ValueError as e:
            print(f"{d26Mg:<12.2f} 错误: {e}")
    
    # 正演模拟
    print("\n【正演模拟】")
    d26Mg_range = np.linspace(-0.20, 0.05, 50)
    results = model.forward_model(d26Mg_range)
    
    # Monte Carlo 不确定性
    print("\n【Monte Carlo 不确定性分析】")
    print("-" * 60)
    mc_results = model.monte_carlo_uncertainty(d26Mg_clay=-0.10, n_iterations=5000)
    
    print(f"输入碎屑 δ²⁶Mg: -0.10 ± 0.07‰")
    print(f"硅酸盐风化通量 F_sili:")
    print(f"  均值: {mc_results['F_sili_mean']/1e10:.2f} × 10¹⁰ mol/yr")
    print(f"  标准差: {mc_results['F_sili_std']/1e10:.2f} × 10¹⁰ mol/yr")
    print(f"  95% CI: [{mc_results['F_sili_95ci'][0]/1e10:.2f}, {mc_results['F_sili_95ci'][1]/1e10:.2f}] × 10¹⁰ mol/yr")
    print(f"\n硅酸盐风化强度指数 SWI:")
    print(f"  均值: {mc_results['SWI_mean']:.1f}%")
    print(f"  标准差: {mc_results['SWI_std']:.1f}%")
    print(f"  95% CI: [{mc_results['SWI_95ci'][0]:.1f}, {mc_results['SWI_95ci'][1]:.1f}]%")
    
    # 可视化
    print("\n【生成可视化】")
    plot_weathering_model_results(results)
    
    return model, results, mc_results


if __name__ == "__main__":
    example_usage()
