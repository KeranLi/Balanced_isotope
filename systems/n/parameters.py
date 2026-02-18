"""
氮同位素体系参数模块
基于 Kang et al. (2023) 和 Ma et al. (2025) 的氮循环箱式模型

参考资料:
- Kang et al. (2023) Nitrate limitation in early Neoproterozoic oceans delayed 
  the ecological rise of eukaryotes. Science Advances.
- Ma et al. (2025) Prolonged nitrate depletion delayed marine ecosystem recovery 
  after the end-Permian mass extinction. Science China Earth Sciences.
"""

from systems.base.isotope_system import IsotopeParameters
from typing import Dict, Tuple


# 标准物质参考值
NITROGEN_STANDARDS = {
    'air': {
        'name': 'Atmospheric N2',
        'delta15': 0.0,  # 定义标准
        'description': '大气氮气，δ¹⁵N定义零点'
    }
}


# 分馏系数范围 (‰)
FRACTIONATION_RANGES = {
    'fixation': {'min': -2.0, 'max': 1.0, 'description': '固氮作用分馏'},
    'water_column_denitrification': {'min': -30.0, 'max': -22.0, 'description': '水柱反硝化分馏'},
    'sedimentary_denitrification': {'min': 0.0, 'max': 0.0, 'description': '沉积反硝化分馏'},
}


# 氮循环通量参数 (Tg N/a)
FLUX_PARAMETERS = {
    'modern': {
        'F_fix': 205.0,              # 固氮通量
        'F_total_burial': 25.0,       # 总埋藏通量
        'F_wcd': 140.0,              # 水柱反硝化通量 (估算)
        'F_sd': 40.0,                # 沉积反硝化通量 (估算)
        'description': '现代海洋氮循环参数 (Kang et al. 2023)'
    },
    'early_triassic': {
        'F_fix': 205.0,              # 固氮通量 (假设与现代相同)
        'F_total_burial': 25.0,       # 总埋藏通量 (假设与现代相同)
        'description': '早三叠世氮循环参数 (Ma et al. 2025)'
    },
    'neoproterozoic': {
        'F_fix': 205.0,              # 固氮通量
        'F_total_burial': 25.0,       # 总埋藏通量
        'description': '新元古代氮循环参数 (Kang et al. 2023)'
    }
}


# 典型情景参数
SCENARIO_PARAMETERS = {
    'modern_oxic': {
        'f_assimilator_range': (0.3, 0.7),  # 现代海洋硝酸盐占比范围
        'delta15N_sed_range': (4.0, 6.0),   # 现代沉积物δ¹⁵N范围
        'description': '现代氧化海洋'
    },
    'early_triassic_stage_I': {
        'f_assimilator_range': (0.0, 0.1),  # 第I阶段 (格里斯巴赫-斯密斯)
        'delta15N_sed_range': (0.0, 2.0),
        'description': '早三叠世第I阶段 (缺氧高温)'
    },
    'early_triassic_stage_II': {
        'f_assimilator_range': (0.15, 0.25),  # 第II阶段 (斯帕斯早期)
        'delta15N_sed_range': (3.0, 5.0),
        'description': '早三叠世第II阶段 (氧化降温)'
    },
    'early_triassic_stage_III': {
        'f_assimilator_range': (0.05, 0.15),  # 第III阶段 (斯帕斯晚期)
        'delta15N_sed_range': (1.0, 3.0),
        'description': '早三叠世第III阶段 (再缺氧)'
    },
    'neoproterozoic_pre_800Ma': {
        'f_assimilator_range': (0.05, 0.20),  # 800Ma前
        'delta15N_sed_range': (0.5, 3.0),
        'description': '新元古代早期 (<800Ma)'
    },
    'neoproterozoic_post_800Ma': {
        'f_assimilator_range': (0.15, 0.35),  # 800Ma后
        'delta15N_sed_range': (3.0, 6.0),
        'description': '新元古代晚期 (>800Ma)'
    },
    'anoxic_nitrate_depleted': {
        'f_assimilator_range': (0.0, 0.1),
        'delta15N_sed_range': (-1.0, 2.0),
        'description': '缺氧硝酸盐匮乏环境'
    }
}


def get_n_parameters(scenario: str = 'modern') -> IsotopeParameters:
    """
    获取氮同位素体系参数
    
    Parameters
    ----------
    scenario : str
        情景名称：'modern', 'early_triassic', 'neoproterozoic'
        
    Returns
    -------
    IsotopeParameters
        氮同位素参数对象
    """
    
    # 基础参数
    base_params = {
        'element': 'n',
        'name': 'Nitrogen',
        'reference_standard': 'Air-N2',
        'reference_ratios': {'15/14': 0.003676},  # 大气氮同位素比值
    }
    
    # 获取通量参数
    flux_params = FLUX_PARAMETERS.get(scenario, FLUX_PARAMETERS['modern'])
    
    return IsotopeParameters(
        **base_params,
        
        # 储库参数 (海洋溶解氮总量，mol N)
        reservoir_mass=5.7e16,  # 约 5700 Tg N
        
        # 输入通量 (mol/Ma 或 Tg/a)
        input_fluxes={
            'fixation': flux_params['F_fix'],           # 固氮作用
            'atmospheric_deposition': 5.0,               # 大气沉降 (次要)
            'river_input': 20.0,                         # 河流输入 (次要)
        },
        
        # 输出通量
        output_fluxes={
            'water_column_denitrification': flux_params.get('F_wcd', 140.0),
            'sedimentary_denitrification': flux_params.get('F_sd', 40.0),
            'burial': flux_params['F_total_burial'],
        },
        
        # 端元同位素值 (‰)
        end_members={
            'atmosphere': {
                'delta15': 0.0,
                'description': '大气N₂ (定义标准)'
            },
            'nitrogen_fixer': {
                'delta15': -1.0,  # 固氮生物典型值
                'range': (-2.0, 1.0),
                'description': '固氮生物 (蓝细菌等)'
            },
            'nitrate_assimilator': {
                'delta15': 5.0,   # 硝酸盐同化生物典型值
                'range': (3.0, 8.0),
                'description': '硝酸盐同化生物 (真核藻类等)'
            },
            'ammonium_dominant_ocean': {
                'delta15': -0.5,
                'range': (-1.0, 1.0),
                'description': '铵主导海洋 (缺氧环境)'
            },
            'nitrate_dominant_ocean': {
                'delta15': 5.0,
                'range': (3.0, 7.0),
                'description': '硝酸盐主导海洋 (氧化环境)'
            }
        },
        
        # 分馏参数 (‰)
        fractionation_factors={
            'epsilon_fixation': -0.5,           # 固氮分馏 (典型值)
            'epsilon_fixation_min': -2.0,       # 固氮分馏最小值
            'epsilon_fixation_max': 1.0,        # 固氮分馏最大值
            'epsilon_wcd': -26.0,               # 水柱反硝化分馏 (典型值)
            'epsilon_wcd_min': -30.0,           # 水柱反硝化分馏最小值
            'epsilon_wcd_max': -22.0,           # 水柱反硝化分馏最大值
            'epsilon_sd': 0.0,                  # 沉积反硝化分馏
        }
    )


def get_scenario_info(scenario_name: str) -> Dict:
    """
    获取情景信息
    
    Parameters
    ----------
    scenario_name : str
        情景名称
        
    Returns
    -------
    dict
        情景参数
    """
    return SCENARIO_PARAMETERS.get(scenario_name, {
        'description': '未知情景',
        'f_assimilator_range': (0.0, 1.0),
        'delta15N_sed_range': (-5.0, 10.0)
    })


def get_fractionation_ranges() -> Dict[str, Tuple[float, float]]:
    """
    获取分馏系数范围
    
    Returns
    -------
    dict
        各过程的分馏系数范围 (min, max)
    """
    return {
        'fixation': (FRACTIONATION_RANGES['fixation']['min'], 
                     FRACTIONATION_RANGES['fixation']['max']),
        'water_column_denitrification': (
            FRACTIONATION_RANGES['water_column_denitrification']['min'],
            FRACTIONATION_RANGES['water_column_denitrification']['max']
        ),
        'sedimentary_denitrification': (
            FRACTIONATION_RANGES['sedimentary_denitrification']['min'],
            FRACTIONATION_RANGES['sedimentary_denitrification']['max']
        )
    }
