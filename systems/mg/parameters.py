"""
Mg同位素参数模块
基于Kasemann等(2014)及其他文献的标准参数

参考标准说明：
- DSM3 = D3MS (Dead Sea Magnesium Standard)
- 这是同一个标准的两种写法，数值完全等价

数据来源：
- Kasemann et al., 2014 - Continental weathering following a Cryogenian glaciation
- Tipper et al., 2006a, 2006b - Mg isotope systematics
- Li et al., 2010 - Mg isotopes in modern environments
"""

from systems.base.isotope_system import IsotopeParameters
from typing import Dict, Optional
import numpy as np


# =============================================================================
# 论文标准参数（Kasemann et al., 2014）
# =============================================================================

# 端元δ²⁶Mg值（‰ vs DSM3）
END_MEMBERS = {
    'silicate': {
        'delta26': -0.3,           # 硅酸盐岩石（上地壳）
        'uncertainty': 0.2,
        'range': (-0.5, 0.0),
        'description': '硅酸盐风化端元（大陆地壳）'
    },
    'carbonate': {
        'delta26': -2.5,           # 碳酸盐岩
        'uncertainty': 0.3,
        'range': (-2.8, -2.2),
        'description': '碳酸盐风化端元（Cryogenian碳酸盐）'
    },
    'seawater': {
        'delta26': -0.83,          # 现代海水
        'uncertainty': 0.09,
        'description': '现代海水'
    },
    'mantle': {
        'delta26': -0.25,          # 地幔
        'uncertainty': 0.04,
        'description': '地幔/原始地球'
    },
    'chondrite': {
        'delta26': -0.30,
        'uncertainty': 0.06,
        'description': '球粒陨石'
    }
}

# 现代海洋通量参数（mol/yr）
MODERN_FLUXES = {
    'river_total': 5.6e12,         # 河流总Mg通量
    'hydrothermal': 2.0e12,        # 海底热液通量（Ca模型中的值，Mg类似）
    'carbonate_burial': 4.0e12,    # 碳酸盐埋藏
}

# 海洋储库参数
RESERVOIR = {
    'mass': 7.3e19,                # 海水Mg总量 (mol)
    'residence_time': 13e6,        # 停留时间 (yr)
}

# 分馏系数（‰）
FRACTIONATION = {
    'calcite_precipitation': -2.7,     # 方解石沉淀（vs 海水）
    'dolomite_precipitation': -1.8,    # 白云石沉淀（vs 海水）
    'silicate_weathering': -0.3,       # 硅酸盐风化
}


def get_mg_parameters() -> IsotopeParameters:
    """
    获取Mg同位素体系的标准参数
    基于Kasemann等(2014)论文参数
    
    Returns
    -------
    IsotopeParameters
        Mg同位素参数对象
    """
    return IsotopeParameters(
        element='mg',
        name='Magnesium',
        
        # 参考标准：DSM3
        reference_standard='DSM3',
        reference_ratios={
            '25/24': 0.12663,
            '26/24': 0.13932
        },
        
        # 分馏系数（ε值，‰）
        fractionation_factors={
            'calcite_precipitation': -2.7,
            'dolomite_precipitation': -1.8,
            'generic_carbonate': -2.7,      # 通用碳酸盐分馏
            'silicate_weathering': -0.3,
            'hydrothermal': 0.0,            # 热液无分馏
        },
        
        # 端元值
        end_members=END_MEMBERS.copy(),
        
        # 海洋储库
        reservoir_mass=RESERVOIR['mass'],
        
        # 输入通量 (mol/yr)
        input_fluxes={
            'river_total': MODERN_FLUXES['river_total'],
            'rivers_silicate': MODERN_FLUXES['river_total'] * 0.5,  # 假设一半来自硅酸盐
            'rivers_carbonate': MODERN_FLUXES['river_total'] * 0.5, # 一半来自碳酸盐
            'hydrothermal': MODERN_FLUXES['hydrothermal'],
        },
        
        # 输出通量 (mol/yr)
        output_fluxes={
            'carbonate_burial': MODERN_FLUXES['carbonate_burial'],
            'hydrothermal_sink': MODERN_FLUXES['hydrothermal'],
        }
    )


def get_cryogenian_parameters() -> IsotopeParameters:
    """
    获取Cryogenian时期（新元古代冰期）的参数
    根据Kasemann等(2014)论文设置
    
    Returns
    -------
    IsotopeParameters
        Cryogenian时期参数
    """
    params = get_mg_parameters()
    
    # Cryogenian碳酸盐同位素特征
    params.end_members['carbonate']['delta26'] = -1.6  # Cryogenian碳酸盐平均值
    params.end_members['carbonate']['range'] = (-2.2, -1.1)
    
    # 估计的高风化通量场景
    params.input_fluxes['high_weathering_silicate'] = 6 * MODERN_FLUXES['river_total']
    params.input_fluxes['high_weathering_carbonate'] = 9 * MODERN_FLUXES['river_total']
    
    return params


# =============================================================================
# 风化比例计算工具函数
# =============================================================================

def calculate_river_delta26(
    f_silicate: float,
    delta_silicate: Optional[float] = None,
    delta_carbonate: Optional[float] = None
) -> float:
    """
    计算河流Mg同位素组成（双端元混合）
    
    δ²⁶Mg_riv = f_sil × δ²⁶Mg_sil + (1-f_sil) × δ²⁶Mg_carb
    
    Parameters
    ----------
    f_silicate : float
        硅酸盐风化比例 (0-1)
    delta_silicate : float, optional
        硅酸盐端元δ²⁶Mg，默认-0.3‰
    delta_carbonate : float, optional
        碳酸盐端元δ²⁶Mg，默认-2.5‰
        
    Returns
    -------
    float
        河流δ²⁶Mg值
    """
    if delta_silicate is None:
        delta_silicate = END_MEMBERS['silicate']['delta26']
    if delta_carbonate is None:
        delta_carbonate = END_MEMBERS['carbonate']['delta26']
    
    f_sil = np.clip(f_silicate, 0, 1)
    return f_sil * delta_silicate + (1 - f_sil) * delta_carbonate


def solve_f_silicate(
    delta_riv: float,
    delta_silicate: Optional[float] = None,
    delta_carbonate: Optional[float] = None
) -> float:
    """
    从河流同位素组成反演硅酸盐风化比例
    
    f_sil = (δ²⁶Mg_riv - δ²⁶Mg_carb) / (δ²⁶Mg_sil - δ²⁶Mg_carb)
    
    Parameters
    ----------
    delta_riv : float
        河流δ²⁶Mg值
    delta_silicate : float, optional
        硅酸盐端元δ²⁶Mg
    delta_carbonate : float, optional
        碳酸盐端元δ²⁶Mg
        
    Returns
    -------
    float
        硅酸盐风化比例 (0-1)
    """
    if delta_silicate is None:
        delta_silicate = END_MEMBERS['silicate']['delta26']
    if delta_carbonate is None:
        delta_carbonate = END_MEMBERS['carbonate']['delta26']
    
    denom = delta_silicate - delta_carbonate
    if abs(denom) < 1e-10:
        return 0.5  # 避免除零
    
    f_sil = (delta_riv - delta_carbonate) / denom
    return np.clip(f_sil, 0, 1)


# =============================================================================
# 时间依赖的风化比例函数
# =============================================================================

def weathering_transition_linear(
    t: float,
    t_start: float,
    t_end: float,
    f_initial: float,
    f_final: float
) -> float:
    """
    线性风化比例转变
    
    Parameters
    ----------
    t : float
        当前时间
    t_start, t_end : float
        转变开始和结束时间
    f_initial, f_final : float
        初始和最终硅酸盐风化比例
        
    Returns
    -------
    float
        当前时间的硅酸盐风化比例
    """
    if t <= t_start:
        return f_initial
    elif t >= t_end:
        return f_final
    else:
        ratio = (t - t_start) / (t_end - t_start)
        return f_initial + ratio * (f_final - f_initial)


def weathering_transition_exponential(
    t: float,
    tau: float,
    f_initial: float,
    f_final: float
) -> float:
    """
    指数风化比例转变
    
    f(t) = f_final + (f_initial - f_final) × exp(-t/τ)
    
    Parameters
    ----------
    t : float
        时间
    tau : float
        时间常数
    f_initial, f_final : float
        初始和最终硅酸盐风化比例
        
    Returns
    -------
    float
        当前时间的硅酸盐风化比例
    """
    return f_final + (f_initial - f_final) * np.exp(-t / tau)


def weathering_pulse(
    t: float,
    t_peak: float,
    duration: float,
    f_baseline: float,
    f_peak: float
) -> float:
    """
    高斯型风化脉冲
    
    Parameters
    ----------
    t : float
        时间
    t_peak : float
        脉冲峰值时间
    duration : float
        脉冲持续时间（半高宽）
    f_baseline : float
        基线硅酸盐风化比例
    f_peak : float
        峰值硅酸盐风化比例
        
    Returns
    -------
    float
        当前时间的硅酸盐风化比例
    """
    sigma = duration / 2.355  # FWHM to sigma
    pulse = np.exp(-0.5 * ((t - t_peak) / sigma) ** 2)
    return f_baseline + (f_peak - f_baseline) * pulse
