"""
Mg同位素参数模块
包含Mg同位素体系的所有物理参数和端元值
数据来源：文献综合（Tipper et al., 2006; Teng et al., 2007等）
"""

from systems.base.isotope_system import IsotopeParameters
from typing import Dict


def get_mg_parameters() -> IsotopeParameters:
    """
    获取Mg同位素体系的标准参数
    
    Returns
    -------
    IsotopeParameters
        Mg同位素参数对象
    """
    return IsotopeParameters(
        element='mg',
        name='Magnesium',
        
        # 参考标准：DSM3 (Dead Sea Magnesium)
        reference_standard='DSM3',
        reference_ratios={
            '25/24': 0.12663,   # ²⁵Mg/²⁴Mg
            '26/24': 0.13932    # ²⁶Mg/²⁴Mg
        },
        
        # 分馏系数（ε值，‰）
        # 碳酸盐-海水分馏（室温）
        fractionation_factors={
            # 平衡分馏
            'carb_sw_equilibrium': -3.5,      # 碳酸盐vs海水（‰）
            'silicate_sw_equilibrium': -0.5,  # 硅酸盐vs海水（‰）
            
            # 动力学分馏
            'carb_precipitation': -2.5,       # 碳酸盐沉淀
            'silicate_weathering': -0.3,      # 硅酸盐风化
            
            # 高温分馏
            'mantle_basalt': -0.2,            # 地幔/玄武岩
            'chondrite': -0.3                 # 球粒陨石
        },
        
        # 端元值（δ²⁶Mg，‰ vs DSM3）
        end_members={
            'seawater': {
                'delta26': -0.83,
                'uncertainty': 0.09,
                'description': '现代海水'
            },
            'carbonate': {
                'delta26': -4.3,              # 典型碳酸盐岩
                'uncertainty': 0.5,
                'range': (-5.0, -3.5),
                'description': '碳酸盐风化端元'
            },
            'silicate': {
                'delta26': -0.3,              # 上地壳平均值
                'uncertainty': 0.3,
                'range': (-0.6, 0.0),
                'description': '硅酸盐风化端元'
            },
            'basalt': {
                'delta26': -0.2,
                'uncertainty': 0.1,
                'description': '玄武岩（地幔来源）'
            },
            'granite': {
                'delta26': -0.4,
                'uncertainty': 0.2,
                'description': '花岗岩（长英质）'
            },
            'mantle': {
                'delta26': -0.25,
                'uncertainty': 0.04,
                'description': '地幔（全球平均）'
            },
            'chondrite': {
                'delta26': -0.30,
                'uncertainty': 0.06,
                'description': '球粒陨石（太阳系原始组成）'
            },
            'rivers': {
                'delta26': -1.2,
                'uncertainty': 0.4,
                'range': (-2.0, -0.5),
                'description': '全球河流平均值'
            }
        },
        
        # 现代海洋储库参数
        reservoir_mass=5.1e19,  # 海水Mg总量 (mol)
        
        # 现代输入通量 (mol/Ma)
        input_fluxes={
            'rivers_carbonate': 2.5e18,   # 碳酸盐风化
            'rivers_silicate': 2.5e18,    # 硅酸盐风化
            'hydrothermal': 0.5e18,       # 海底热液
            'dolomitization': 1.0e18      # 白云岩化
        },
        
        # 现代输出通量 (mol/Ma)
        output_fluxes={
            'carbonate_burial': 4.0e18,   # 碳酸盐埋藏
            'clay_formation': 1.0e18,     # 粘土形成
            'hydrothermal_sink': 0.5e18   # 热液汇
        }
    )


def get_ancient_parameters(age_ma: float) -> IsotopeParameters:
    """
    获取特定地质时期的Mg同位素参数
    考虑地质演化过程中的参数变化
    
    Parameters
    ----------
    age_ma : float
        年龄（Ma）
        
    Returns
    -------
    IsotopeParameters
        调整后的参数
    """
    params = get_mg_parameters()
    
    # 示例：根据地质时期调整参数
    # 这里可以根据文献数据添加更多条件
    
    if age_ma > 1000:  # 前寒武纪
        # 假设前寒武纪海水Mg同位素与现代不同
        params.end_members['seawater']['delta26'] = -0.5
        # 硅酸盐风化可能更强
        params.input_fluxes['rivers_silicate'] *= 1.5
    
    return params


# 分馏系数温度依赖公式参数
# 1000·ln(α) = A/T² + B/T （T单位为K）
FRACTIONATION_TEMP_PARAMS = {
    'calcite_water': {
        'A': 0.0,      # 待补充
        'B': -0.8,     # 简化值
        'reference': 'Li et al., 2012'
    },
    'dolomite_water': {
        'A': 0.0,
        'B': -1.2,
        'reference': 'Li et al., 2015'
    }
}


def calculate_temperature_dependent_fractionation(
    temperature_c: float,
    mineral: str = 'calcite'
) -> float:
    """
    计算温度依赖的分馏系数
    
    Parameters
    ----------
    temperature_c : float
        温度（摄氏度）
    mineral : str
        矿物类型（calcite/dolomite）
        
    Returns
    -------
    float
        分馏系数ε（‰）
    """
    T_k = temperature_c + 273.15
    
    if mineral == 'calcite':
        params = FRACTIONATION_TEMP_PARAMS['calcite_water']
    elif mineral == 'dolomite':
        params = FRACTIONATION_TEMP_PARAMS['dolomite_water']
    else:
        raise ValueError(f"Unknown mineral: {mineral}")
    
    # 1000·ln(α) = A/T² + B/T
    ln_alpha = params['A'] / T_k**2 + params['B'] / T_k
    
    # 转换为ε（‰）
    epsilon = (np.exp(ln_alpha / 1000) - 1) * 1000
    
    return epsilon


import numpy as np
