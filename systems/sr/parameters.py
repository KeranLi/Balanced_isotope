"""
Sr同位素参数模块
基于Wang等(2021)及其他文献的标准参数

参考：
- Wang et al., 2021 - Earth-Science Reviews 218: 103679
- Peucker-Ehrenbrink and Fiske, 2019 - Continental weathering
- Coogan and Dosso, 2015 - Low-temperature hydrothermal flux
- Li and Elderfield, 2013 - High-temperature hydrothermal flux
"""

from systems.base.isotope_system import IsotopeParameters
from typing import Dict, Tuple
import numpy as np


# =============================================================================
# 现代Sr通量参数（mol/yr）
# =============================================================================

MODERN_FLUXES = {
    # 河流输入
    'river': 47.6e9,           # Peucker-Ehrenbrink & Fiske, 2019
    
    # 高温热液输入（洋中脊）
    'hydrothermal_highT': 8.04e9,  # Li & Elderfield, 2013
    
    # 低温热液输入（洋壳蚀变）
    'hydrothermal_lowT': 10e9,     # Coogan & Dosso, 2015
    
    # 成岩作用输入
    'diagenetic': 3.4e9,           # 估算值
}

# =============================================================================
# 端元⁸⁷Sr/⁸⁶Sr值
# =============================================================================

END_MEMBERS = {
    'river': {
        'ratio': 0.71107,          # 现代河流平均
        'uncertainty': 0.0005,
        'range': (0.7025, 0.75),   # 可变范围（不同岩性）
        'description': '河流输入（大陆风化）'
    },
    'hydrothermal_highT': {
        'ratio': 0.7037,           # 地幔值
        'uncertainty': 0.0001,
        'range': (0.703, 0.704),   # Antonelli et al., 2017
        'description': '高温热液（洋中脊）'
    },
    'hydrothermal_lowT': {
        'ratio': 0.7084,           # Coogan & Dosso, 2015
        'uncertainty': 0.0005,
        'range': (0.7025, 0.7084), # 温度依赖
        'description': '低温热液（洋壳蚀变）'
    },
    'diagenetic': {
        'ratio': 0.7084,           # 与海水平衡
        'uncertainty': 0.0005,
        'description': '成岩作用'
    },
    'seawater_modern': {
        'ratio': 0.70917,          # 现代海水
        'uncertainty': 0.00001,
        'description': '现代海水'
    },
    'mantle': {
        'ratio': 0.7037,           # 地幔值
        'uncertainty': 0.0001,
        'description': '地幔/玄武岩'
    }
}

# =============================================================================
# 海洋储库参数
# =============================================================================

RESERVOIR = {
    'mass': 1.2e17,                # 海水Sr总量 (mol)
    'residence_time': 2.6e6,       # 停留时间 (yr) - Wang et al., 2021
}

# =============================================================================
# 二叠纪特定参数
# =============================================================================

PERMIAN_PARAMETERS = {
    # 时间范围 (Ma)
    'time_range': (299, 252),
    
    # 观测到的Sr同位素范围（基于brachiopod数据）
    'sr_ratio_range': (0.706832, 0.708271),
    
    # 变化率 (每Myr)
    'decreasing_rate': -0.000053,  # Asselian-Kungurian下降
    'increasing_rate': 0.000076,   # Changhsingian上升
}


# =============================================================================
# 模型参数范围（用于随机模拟）
# =============================================================================

# 基于Wang et al. (2021) Table 1
STOCHASTIC_RANGES = {
    'F_river': (10e9, 190e9),              # 河流通量范围 (mol/yr)
    'R_river': (0.7025, 0.75),             # 河流同位素范围
    'F_hydrothermal_lowT': (2.5e9, 40e9),  # 低温热液通量
    'R_hydrothermal_lowT': (0.7025, 0.7084),
    'F_hydrothermal_highT': (2e9, 35e9),   # 高温热液通量
    'R_hydrothermal_highT': (0.703, 0.704),
    'F_diagenetic': (3.4e9, 3.4e9),        # 保持恒定
    'R_diagenetic': (0.7084, 0.7084),      # 保持恒定
}

# 参数不确定性（1σ）
PARAMETER_UNCERTAINTIES = {
    'F_river': 0.10,              # 10%
    'R_river': 0.0005,
    'F_hydrothermal_lowT': 0.10,
    'R_hydrothermal_lowT': 0.0005,
    'F_hydrothermal_highT': 0.10,
    'R_hydrothermal_highT': 0.0001,
}


def get_sr_parameters() -> IsotopeParameters:
    """
    获取Sr同位素体系的标准参数
    
    Returns
    -------
    IsotopeParameters
        Sr同位素参数对象
    """
    return IsotopeParameters(
        element='sr',
        name='Strontium',
        
        # 参考标准：NBS-987
        reference_standard='NBS-987',
        reference_ratios={
            '87/86': 0.710247,     # NBS-987标准值
            '86/88': 0.1194,
        },
        
        # 分馏系数（Sr同位素无显著质量分馏，主要考虑混合）
        fractionation_factors={},
        
        # 端元值
        end_members=END_MEMBERS.copy(),
        
        # 海洋储库
        reservoir_mass=RESERVOIR['mass'],
        
        # 输入通量 (mol/yr)
        input_fluxes=MODERN_FLUXES.copy(),
        
        # 输出通量 (mol/yr) - 稳态时等于输入总和
        output_fluxes={
            'total': sum(MODERN_FLUXES.values()),
        }
    )


def get_modern_parameters() -> IsotopeParameters:
    """
    获取现代Sr同位素体系参数
    
    Returns
    -------
    IsotopeParameters
        现代参数设置
    """
    params = get_sr_parameters()
    
    # 使用现代观测值
    params.end_members['seawater'] = END_MEMBERS['seawater_modern'].copy()
    
    return params


def get_permian_parameters() -> IsotopeParameters:
    """
    获取二叠纪Sr同位素体系参数
    基于Wang et al. (2021)论文设置
    
    Returns
    -------
    IsotopeParameters
        二叠纪参数设置
    """
    params = get_sr_parameters()
    
    # 二叠纪特定参数
    params.end_members['seawater'] = {
        'ratio': 0.7070,           # 二叠纪平均值（近似）
        'range': PERMIAN_PARAMETERS['sr_ratio_range'],
        'description': '二叠纪海水'
    }
    
    # 添加随机模拟范围
    params.fractionation_factors['stochastic_ranges'] = STOCHASTIC_RANGES.copy()
    params.fractionation_factors['uncertainties'] = PARAMETER_UNCERTAINTIES.copy()
    
    return params


# =============================================================================
# 工具函数
# =============================================================================

def calculate_seawater_sr(
    F_river: float,
    R_river: float,
    F_highT: float,
    R_highT: float,
    F_lowT: float,
    R_lowT: float,
    F_dia: float = 3.4e9,
    R_dia: float = 0.7084
) -> float:
    """
    计算海水⁸⁷Sr/⁸⁶Sr值
    
    基于质量平衡方程：
    R_sw = (F_riv×R_riv + F_highT×R_highT + F_lowT×R_lowT + F_dia×R_dia) / 
           (F_riv + F_highT + F_lowT + F_dia)
    
    Parameters
    ----------
    F_river, R_river : float
        河流输入通量(mol/yr)和同位素比值
    F_highT, R_highT : float
        高温热液通量和同位素比值
    F_lowT, R_lowT : float
        低温热液通量和同位素比值
    F_dia, R_dia : float
        成岩作用通量和同位素比值
        
    Returns
    -------
    float
        海水⁸⁷Sr/⁸⁶Sr比值
    """
    total_input = F_river + F_highT + F_lowT + F_dia
    
    if total_input == 0:
        return 0.708  # 默认值
    
    weighted_sum = (
        F_river * R_river +
        F_highT * R_highT +
        F_lowT * R_lowT +
        F_dia * R_dia
    )
    
    return weighted_sum / total_input


def generate_random_parameters(
    ranges: Dict[str, Tuple[float, float]] = None,
    distribution: str = 'uniform'
) -> Dict[str, float]:
    """
    生成随机参数组合
    
    Parameters
    ----------
    ranges : dict, optional
        参数范围字典，默认使用STOCHASTIC_RANGES
    distribution : str
        分布类型 ('uniform' 或 'loguniform')
        
    Returns
    -------
    dict
        随机参数组合
    """
    if ranges is None:
        ranges = STOCHASTIC_RANGES
    
    params = {}
    
    for key, (min_val, max_val) in ranges.items():
        if distribution == 'uniform':
            params[key] = np.random.uniform(min_val, max_val)
        elif distribution == 'loguniform':
            # 对数均匀分布（适用于通量）
            log_min, log_max = np.log(min_val), np.log(max_val)
            params[key] = np.exp(np.random.uniform(log_min, log_max))
        else:
            raise ValueError(f"Unknown distribution: {distribution}")
    
    return params


def filter_by_observations(
    model_ratios: np.ndarray,
    observed_ratios: np.ndarray,
    observed_ages: np.ndarray,
    model_ages: np.ndarray,
    confidence_interval: float = 0.95
) -> bool:
    """
    检查模型结果是否在观测数据的置信区间内
    
    Parameters
    ----------
    model_ratios : array_like
        模型计算的Sr同位素比值
    observed_ratios : array_like
        观测的Sr同位素比值
    observed_ages, model_ages : array_like
        对应的时间点（Ma）
    confidence_interval : float
        置信区间（默认95%）
        
    Returns
    -------
    bool
        是否在置信区间内
    """
    # 简化的过滤逻辑：检查是否在给定误差范围内匹配
    # 实际应用中可能需要更复杂的统计检验
    
    tolerance = (1 - confidence_interval) * 0.1  # 简化的容差
    
    # 插值观测数据到模型时间点
    observed_interp = np.interp(model_ages, observed_ages, observed_ratios)
    
    # 计算差异
    diff = np.abs(model_ratios - observed_interp)
    max_diff = np.max(diff)
    
    return max_diff < tolerance
