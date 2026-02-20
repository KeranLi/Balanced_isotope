"""
铀同位素体系参数模块
基于海洋铀循环模型 (Algeo et al., 2023; Lau et al., 2016; Elrick et al., 2017)

参考资料:
- Algeo et al. (2023) Marine uranium cycle modeling
- Montoya-Pino et al. (2010) Geology
- Brennecka et al. (2011) PNAS
- Lau et al. (2016) PNAS
- Elrick et al. (2017) Geology
- Tissot & Dauphas (2015) GCA
- Weyer et al. (2008) GCA
- Romaniello et al. (2013) Chem. Geol.
- Stirling et al. (2015) GCA

标准物质:
- CRM-112a (NBL): 重铀酸铵标准，定义 δ²³⁸U = 0‰
- 现代海水 δ²³⁸U ≈ -0.39‰ (Tissot & Dauphas, 2015)
"""

from systems.base.isotope_system import IsotopeParameters
from typing import Dict, Tuple


# 标准物质参考值
URANIUM_STANDARDS = {
    'nbl': {
        'name': 'CRM-112a (NBL)',
        'delta238': 0.0,  # 定义标准
        'description': 'NIST SRM 950a / CRM-112a，重铀酸铵标准'
    },
    'seawater': {
        'name': 'Modern Seawater',
        'delta238': -0.392,  # Tissot & Dauphas (2015)
        'range': (-0.40, -0.35),
        'description': '现代海水'
    }
}


# 分馏系数 (‰)
FRACTIONATION_PARAMETERS = {
    'delta_sw_ox': {
        'value': 0.0,
        'range': (0.0, 0.0),
        'description': '氧化/亚氧化汇分馏 (Chen et al., 2016)'
    },
    'delta_sw_anox': {
        'value': 0.77,
        'range': (0.73, 0.81),
        'description': '缺氧/硫化汇分馏 (Stirling et al., 2015)'
    },
    'delta_diag': {
        'value': 0.40,
        'range': (0.30, 0.50),
        'description': '成岩校正因子 (Romaniello et al., 2013; Elrick et al., 2017)'
    },
    'delta_river': {
        'value': -0.29,
        'range': (-0.72, 0.06),
        'description': '河流输入分馏 (Andersen et al., 2016)'
    }
}


# 海洋铀储库参数
RESERVOIR_PARAMETERS = {
    'modern': {
        'mass': 3.7e16,           # mol U (现代海洋总铀储量)
        'concentration': 13.6e-9,  # mol/kg (13.6 nM)
        'delta238_seawater': -0.392,  # ‰
        'residence_time': 0.5,     # Myr (停留时间)
        'description': '现代海洋铀储库'
    },
    'frasnian_famennian': {
        'mass': 3.7e16,           # mol U (假设与现代类似)
        'delta238_seawater': None,  # 待计算
        'residence_time': 0.5,
        'description': 'Frasnian-Famennian过渡期'
    }
}


# 通量参数 (mol/Ma)
FLUX_PARAMETERS = {
    'modern': {
        'F_river': 7.4e16,        # 河流输入通量
        'F_total': 7.4e16,        # 总输出通量 (稳态)
        'description': '现代海洋铀循环通量'
    },
    'high_anoxia': {
        'F_river': 7.4e16,        # 河流输入通量
        'F_total': 7.4e16,        # 总输出通量
        'description': '高缺氧环境'
    }
}


# 典型地质情景参数
SCENARIO_PARAMETERS = {
    'modern_oxic': {
        'f_anox_range': (0.05, 0.15),  # 现代海洋缺氧汇比例
        'delta238_carb_range': (-0.60, -0.40),  # 碳酸盐δ²³⁸U范围
        'anoxic_area': 0.5,  # % 海底缺氧面积
        'description': '现代氧化海洋'
    },
    'oceanic_anoxic_event': {
        'f_anox_range': (0.50, 0.90),
        'delta238_carb_range': (-1.10, -0.70),
        'anoxic_area': 5.0,
        'description': '大洋缺氧事件 (OAE)'
    },
    'end_permain': {
        'f_anox_range': (0.60, 0.95),
        'delta238_carb_range': (-1.20, -0.80),
        'anoxic_area': 10.0,
        'description': '二叠纪末大灭绝'
    },
    'frasnian_famennian': {
        'f_anox_range': (0.30, 0.85),
        'delta238_carb_range': (-0.95, -0.55),
        'anoxic_area': 1.5,
        'description': 'Frasnian-Famennian过渡 (F-F事件)'
    }
}


def get_u_parameters(scenario: str = 'modern') -> IsotopeParameters:
    """
    获取铀同位素体系参数
    
    Parameters
    ----------
    scenario : str
        情景名称：'modern', 'frasnian_famennian', 'oceanic_anoxic_event'
        
    Returns
    -------
    IsotopeParameters
        铀同位素参数对象
    """
    
    # 基础参数
    base_params = {
        'element': 'u',
        'name': 'Uranium',
        'reference_standard': 'CRM-112a (NBL)',
        'reference_ratios': {'238/235': 137.818},  # ²³⁸U/²³⁵U 自然丰度比
    }
    
    # 获取储库和通量参数
    reservoir = RESERVOIR_PARAMETERS.get(scenario, RESERVOIR_PARAMETERS['modern'])
    flux = FLUX_PARAMETERS.get(scenario, FLUX_PARAMETERS['modern'])
    
    return IsotopeParameters(
        **base_params,
        
        # 储库参数 (mol U)
        reservoir_mass=reservoir['mass'],
        
        # 输入通量 (mol/Ma)
        input_fluxes={
            'river': flux['F_river'],
            'hydrothermal': 1.0e15,  # 热液输入 (次要)
        },
        
        # 输出通量 (mol/Ma)
        output_fluxes={
            'oxic_sink': flux['F_total'] * 0.8,   # 氧化汇 (初始估计)
            'anoxic_sink': flux['F_total'] * 0.2,  # 缺氧汇 (初始估计)
            'total': flux['F_total'],
        },
        
        # 端元同位素值 (‰)
        end_members={
            'river': {
                'delta238': -0.29,
                'range': (-0.72, 0.06),
                'description': '河流输入 (平均地壳值)'
            },
            'seawater': {
                'delta238': -0.392,
                'range': (-0.40, -0.35),
                'description': '现代海水'
            },
            'oxic_sink': {
                'delta238': -0.392,  # 与海水相同 (Δ=0)
                'range': (-0.45, -0.35),
                'description': '氧化/亚氧化汇 (碳酸盐等)'
            },
            'anoxic_sink': {
                'delta238': -1.16,  # -0.392 - 0.77
                'range': (-1.20, -1.10),
                'description': '缺氧/硫化汇 (富有机质沉积)'
            },
            'hydrothermal': {
                'delta238': -0.29,
                'range': (-0.35, -0.25),
                'description': '热液输入'
            }
        },
        
        # 分馏参数 (‰)
        fractionation_factors={
            'delta_sw_ox': 0.0,      # 氧化汇分馏
            'delta_sw_anox': 0.77,   # 缺氧汇分馏
            'delta_diag': 0.40,      # 成岩校正
            'delta_river': -0.29,    # 河流输入
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
        'f_anox_range': (0.0, 1.0),
        'delta238_carb_range': (-1.5, 0.0),
        'anoxic_area': 0.0
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
        'delta_sw_ox': FRACTIONATION_PARAMETERS['delta_sw_ox']['range'],
        'delta_sw_anox': FRACTIONATION_PARAMETERS['delta_sw_anox']['range'],
        'delta_diag': FRACTIONATION_PARAMETERS['delta_diag']['range'],
        'delta_river': FRACTIONATION_PARAMETERS['delta_river']['range'],
    }
