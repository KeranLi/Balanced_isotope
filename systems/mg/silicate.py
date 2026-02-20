"""
Mg同位素碎屑岩（硅酸盐）体系模型
基于 Hu et al. (2023) 框架

核心功能：
1. 从碎屑沉积物 δ²⁶Mg 反演风化程度
2. 计算硅酸盐风化通量 (F_sili)
3. 评估硅酸盐风化强度指数 (SWI)

适用场景：
- 分析碎屑岩组分（黏土矿物）的 Mg 同位素
- 区分硅酸盐 vs 碳酸盐风化贡献
- 约束大陆尺度硅酸盐风化通量

与碳酸盐体系的区别：
- 碳酸盐体系：基于海相碳酸盐岩 δ²⁶Mg 反演风化
- 碎屑岩体系：基于陆源碎屑沉积物 δ²⁶Mg 反演风化
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple, Dict, Any


@dataclass
class SilicateWeatheringParams:
    """硅酸盐风化模拟参数 - 基于 Hu et al. (2023)"""
    
    # 同位素端元 (‰ vs DSM3)
    d26Mg_UCC: float = -0.25           # 上地壳基准值
    d26Mg_carbonate: float = -2.0       # 碳酸盐风化端元
    d26Mg_river_water: float = -1.14    # 河水 δ²⁶Mg (长江实测)
    
    # 分馏因子 (‰)
    Delta_fluid_protolith: float = -0.50  # 流体与原岩分馏 (范围: -0.60 ~ -0.40)
    Delta_release_clay: float = -0.45     # 释放Mg与黏土的分馏
    
    # 通量参数 (mol/yr)
    F_river_total: float = 30e10        # 河流总Mg通量 (长江: 30×10¹⁰ mol/yr)
    
    # 不确定性
    uncertainty_clay: float = 0.07      # 碎屑 δ²⁶Mg 分析误差
    uncertainty_river: float = 0.15     # 河水 δ²⁶Mg 变异
    
    # 模拟参数
    weathering_stages: Dict[str, Tuple[float, float]] = field(default_factory=lambda: {
        'incipient': (0.8, 1.0),        # 初始风化: 保留80-100% Mg
        'intermediate': (0.5, 0.8),     # 中等风化: 保留50-80% Mg
        'advanced': (0.2, 0.5),         # 高级风化: 保留20-50% Mg
        'extreme': (0.0, 0.2),          # 极端风化: 保留0-20% Mg
    })


@dataclass
class SilicateModelResult:
    """硅酸盐体系模型结果"""
    # 基础字段
    success: bool = True
    message: str = ""
    data: Dict[str, Any] = field(default_factory=dict)
    
    # 风化参数
    f_Mg: float = None                    # 保留Mg比例
    weathering_stage: str = None          # 风化阶段
    
    # 同位素端元
    d26Mg_silicate: float = None          # 硅酸盐风化端元
    d26Mg_carbonate: float = -2.0         # 碳酸盐风化端元
    
    # 通量结果
    F_silicate: float = None              # 硅酸盐风化通量 (mol/yr)
    F_carbonate: float = None             # 碳酸盐风化通量 (mol/yr)
    SWI: float = None                     # 硅酸盐风化强度指数 (%)
    
    # 质量平衡
    mass_balance_check: float = None      # 质量平衡检验
    
    def get(self, key: str, default=None):
        """获取数据字段"""
        return self.data.get(key, default)


class SilicateMgSystem:
    """
    Mg同位素碎屑岩（硅酸盐）体系
    
    基于 Hu et al. (2023) 实现，从碎屑沉积物（黏土）Mg同位素
    反演硅酸盐风化通量。
    
    核心假设：
    1. 碎屑沉积物中的Mg主要来自硅酸盐风化残余
    2. 碳酸盐Mg在搬运过程中已基本溶解
    3. 风化过程遵循Rayleigh分馏
    
    与碳酸盐体系的关键区别：
    - 输入数据：碎屑岩 δ²⁶Mg vs 碳酸盐岩 δ²⁶Mg
    - 模型原理：风化残余分馏 vs 海水沉淀分馏
    - 应用对象：陆源碎屑沉积物 vs 海相碳酸盐岩
    """
    
    COMPONENT_TYPE = 'siliciclastic'  # 体系类型标识
    ELEMENT = 'mg'
    NAME = 'Magnesium (Siliciclastic)'
    
    def __init__(self, params: Optional[SilicateWeatheringParams] = None, 
                 basin: str = 'changjiang'):
        """
        初始化硅酸盐体系
        
        Parameters:
            params: 自定义参数，默认使用 Hu et al. (2023) 长江参数
            basin: 流域名称 ('changjiang', 'global', 'custom')
        """
        self.params = params or self._get_default_params(basin)
        self.basin = basin
        
        # 内部状态
        self._last_result: Optional[SilicateModelResult] = None
    
    def _get_default_params(self, basin: str) -> SilicateWeatheringParams:
        """获取默认参数集"""
        if basin == 'changjiang':
            # Hu et al. (2023) 长江参数
            return SilicateWeatheringParams()
        elif basin == 'global':
            # 全球平均参数 (需根据文献调整)
            return SilicateWeatheringParams(
                d26Mg_river_water=-1.0,
                F_river_total=100e10,  # 全球河流总通量
            )
        else:
            return SilicateWeatheringParams()
    
    def calculate_weathering_degree(self, d26Mg_clay: float) -> Tuple[float, str]:
        """
        从碎屑 δ²⁶Mg 计算风化程度
        
        使用 Rayleigh 分馏模型:
        δ²⁶Mg_clay = δ²⁶Mg_UCC - Δ·ln(f_Mg)
        
        Parameters:
            d26Mg_clay: 碎屑黏土 δ²⁶Mg (‰)
            
        Returns:
            (f_Mg, stage): 保留Mg比例 (0-1) 和 风化阶段描述
        """
        p = self.params
        
        # 检查输入合理性
        if d26Mg_clay < p.d26Mg_UCC - 0.5:
            # 允许一定范围内的偏差
            pass
        
        # Rayleigh 分馏反演
        numerator = p.d26Mg_UCC - d26Mg_clay
        denominator = p.Delta_fluid_protolith
        
        # 计算 f_Mg
        f_Mg = np.exp(numerator / denominator)
        
        # 限制在合理范围
        f_Mg = np.clip(f_Mg, 0.01, 2.0)
        
        # 判断风化阶段
        stage = self._classify_weathering_stage(f_Mg)
        
        return f_Mg, stage
    
    def _classify_weathering_stage(self, f_Mg: float) -> str:
        """分类风化阶段"""
        stages = self.params.weathering_stages
        for stage_name, (f_min, f_max) in stages.items():
            if f_min <= f_Mg <= f_max:
                return stage_name
        # 边界情况
        if f_Mg > 1.0:
            return 'unweathered'
        return 'extreme'
    
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
    
    def calculate_weathering_flux(self, d26Mg_clay: float) -> SilicateModelResult:
        """
        计算硅酸盐风化通量 - 主计算函数
        
        使用双端元质量平衡:
        F_riv = F_sili + F_car
        F_riv·δ²⁶Mg_riv = F_sili·δ²⁶Mg_sili + F_car·δ²⁶Mg_car
        
        Parameters:
            d26Mg_clay: 碎屑黏土 δ²⁶Mg (‰)
            
        Returns:
            SilicateModelResult: 完整计算结果
        """
        p = self.params
        
        # 步骤1: 计算风化程度
        f_Mg, stage = self.calculate_weathering_degree(d26Mg_clay)
        
        # 步骤2: 计算硅酸盐风化端元
        d26Mg_sili = self.calculate_silicate_endmember(d26Mg_clay)
        
        # 步骤3: 求解质量平衡
        # F_sili = F_riv * (δ²⁶Mg_riv - δ²⁶Mg_car) / (δ²⁶Mg_sili - δ²⁶Mg_car)
        numerator = p.d26Mg_river_water - p.d26Mg_carbonate
        denominator = d26Mg_sili - p.d26Mg_carbonate
        
        if abs(denominator) < 1e-10:
            # 避免除零
            return SilicateModelResult(
                success=False,
                message=f"端元值过于接近: δ²⁶Mg_sili = {d26Mg_sili:.3f}‰, "
                        f"δ²⁶Mg_car = {p.d26Mg_carbonate:.3f}‰"
            )
        
        ratio = numerator / denominator
        F_sili = p.F_river_total * ratio
        F_car = p.F_river_total - F_sili
        
        # 硅酸盐风化强度指数
        SWI = ratio * 100  # 百分比
        
        # 质量平衡检验
        check_lhs = p.d26Mg_river_water
        check_rhs = (F_sili * d26Mg_sili + F_car * p.d26Mg_carbonate) / p.F_river_total
        mass_balance = abs(check_lhs - check_rhs)
        
        # 构建结果
        result = SilicateModelResult(
            success=True,
            data={
                'd26Mg_clay': d26Mg_clay,
                'f_Mg': f_Mg,
                'weathering_stage': stage,
                'd26Mg_silicate': d26Mg_sili,
                'd26Mg_carbonate': p.d26Mg_carbonate,
                'F_silicate': F_sili,
                'F_carbonate': F_car,
                'SWI': SWI,
                'F_river_total': p.F_river_total,
            },
            f_Mg=f_Mg,
            weathering_stage=stage,
            d26Mg_silicate=d26Mg_sili,
            d26Mg_carbonate=p.d26Mg_carbonate,
            F_silicate=F_sili,
            F_carbonate=F_car,
            SWI=SWI,
            mass_balance_check=mass_balance
        )
        
        self._last_result = result
        return result
    
    def calculate_from_row(self, row_data: Dict[str, Any]) -> SilicateModelResult:
        """
        从数据行计算（用于批量处理）
        
        Parameters:
            row_data: 包含 'delta_26_Mg_iso' 的字典或Series
            
        Returns:
            SilicateModelResult
        """
        # 支持多种列名
        possible_cols = ['delta_26_Mg_iso', 'd26Mg', 'δ26Mg', 'delta_26Mg']
        d26Mg_clay = None
        
        for col in possible_cols:
            if col in row_data:
                d26Mg_clay = float(row_data[col])
                break
        
        if d26Mg_clay is None:
            return SilicateModelResult(
                success=False,
                message=f"未找到有效的Mg同位素列，支持的列名: {possible_cols}"
            )
        
        return self.calculate_weathering_flux(d26Mg_clay)
    
    def monte_carlo_analysis(self, d26Mg_clay: float, 
                            d26Mg_clay_std: float = None,
                            n_iterations: int = 10000) -> Dict[str, Any]:
        """
        Monte Carlo 不确定性分析
        
        Parameters:
            d26Mg_clay: 观测的碎屑 δ²⁶Mg
            d26Mg_clay_std: 测量标准差，默认使用参数中的 uncertainty_clay
            n_iterations: 迭代次数
            
        Returns:
            统计结果字典
        """
        p = self.params
        d26Mg_clay_std = d26Mg_clay_std or p.uncertainty_clay
        
        # 生成随机参数
        np.random.seed(42)
        d26Mg_clay_samples = np.random.normal(d26Mg_clay, d26Mg_clay_std, n_iterations)
        Delta_samples = np.random.uniform(-0.60, -0.40, n_iterations)
        d26Mg_river_samples = np.random.normal(p.d26Mg_river_water, p.uncertainty_river, n_iterations)
        
        F_sili_samples = []
        SWI_samples = []
        d26Mg_sili_samples = []
        
        for i in range(n_iterations):
            # 临时修改参数
            orig_Delta = p.Delta_fluid_protolith
            orig_d26Mg_river = p.d26Mg_river_water
            
            p.Delta_fluid_protolith = Delta_samples[i]
            p.d26Mg_river_water = d26Mg_river_samples[i]
            
            try:
                result = self.calculate_weathering_flux(d26Mg_clay_samples[i])
                if result.success:
                    F_sili_samples.append(result.F_silicate)
                    SWI_samples.append(result.SWI)
                    d26Mg_sili_samples.append(result.d26Mg_silicate)
            except (ValueError, ZeroDivisionError):
                pass
            
            # 恢复参数
            p.Delta_fluid_protolith = orig_Delta
            p.d26Mg_river_water = orig_d26Mg_river
        
        F_sili_samples = np.array(F_sili_samples)
        SWI_samples = np.array(SWI_samples)
        d26Mg_sili_samples = np.array(d26Mg_sili_samples)
        
        return {
            'F_silicate': {
                'mean': np.mean(F_sili_samples),
                'std': np.std(F_sili_samples),
                'median': np.median(F_sili_samples),
                'ci_95': np.percentile(F_sili_samples, [2.5, 97.5]),
                'ci_68': np.percentile(F_sili_samples, [16, 84]),
                'samples': F_sili_samples,
            },
            'SWI': {
                'mean': np.mean(SWI_samples),
                'std': np.std(SWI_samples),
                'median': np.median(SWI_samples),
                'ci_95': np.percentile(SWI_samples, [2.5, 97.5]),
                'ci_68': np.percentile(SWI_samples, [16, 84]),
                'samples': SWI_samples,
            },
            'd26Mg_silicate': {
                'mean': np.mean(d26Mg_sili_samples),
                'std': np.std(d26Mg_sili_samples),
                'ci_95': np.percentile(d26Mg_sili_samples, [2.5, 97.5]),
            },
            'n_valid': len(F_sili_samples),
            'n_total': n_iterations,
        }
    
    def simulate_transient(self, 
                          d26Mg_clay_timeseries: np.ndarray,
                          times: np.ndarray) -> list:
        """
        模拟时间序列的碎屑 δ²⁶Mg 变化
        
        Parameters:
            d26Mg_clay_timeseries: 碎屑 δ²⁶Mg 时间序列
            times: 对应的时间点 (Myr)
            
        Returns:
            SilicateModelResult 列表
        """
        results = []
        for d26Mg in d26Mg_clay_timeseries:
            result = self.calculate_weathering_flux(d26Mg)
            results.append(result)
        return results
    
    def get_info(self) -> Dict[str, Any]:
        """获取体系信息"""
        p = self.params
        return {
            'component_type': self.COMPONENT_TYPE,
            'basin': self.basin,
            'reference': 'Hu et al. (2023)',
            'description': '碎屑岩（硅酸盐）Mg同位素体系 - 基于黏土矿物δ²⁶Mg反演风化通量',
            'end_members': {
                'UCC': {'d26Mg': p.d26Mg_UCC, 'description': '上地壳'},
                'carbonate': {'d26Mg': p.d26Mg_carbonate, 'description': '碳酸盐风化'},
                'river_water': {'d26Mg': p.d26Mg_river_water, 'description': '河水'},
            },
            'fractionation_factors': {
                'Delta_fluid_protolith': p.Delta_fluid_protolith,
                'Delta_release_clay': p.Delta_release_clay,
            },
            'fluxes': {
                'F_river_total': p.F_river_total,
            },
            'applicability': [
                '陆源碎屑沉积物（黏土粒级）',
                '河流沉积物',
                '大陆边缘海相泥质岩',
            ],
        }


# 别名，方便导入
SilicateWeatheringModel = SilicateMgSystem
