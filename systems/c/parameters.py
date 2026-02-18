"""
C同位素参数模块
包含碳同位素体系的所有物理参数
数据来源：Li et al. 2020, Precambrian Research等
"""

from systems.base.isotope_system import IsotopeParameters


def get_c_parameters(scenario: str = 'dice') -> IsotopeParameters:
    """
    获取C同位素体系参数
    
    Parameters
    ----------
    scenario : str
        'dice' - DICE事件参数
        'modern' - 现代碳循环参数
        'perturbation' - 一般扰动情景
        
    Returns
    -------
    IsotopeParameters
        C同位素参数对象
    """
    
    # 基础参数（所有情景共用）
    base_params = {
        'element': 'c',
        'name': 'Carbon',
        'reference_standard': 'VPDB',
        'reference_ratios': {'13/12': 0.0112372},
    }
    
    if scenario == 'dice':
        # DICE (Doushantuo Carbon Isotope Excursion) 参数
        # 来自 Li et al. 2020, Precambrian Research
        return IsotopeParameters(
            **base_params,
            
            # 海洋DIC库储量 (mol)
            reservoir_mass=4.0e18,
            
            # 输入通量 (mol/Ma)
            input_fluxes={
                'weathering': 25.0e18,      # 风化输入
                'doc_remineralization': 0,  # DOC再矿化（可变）
                'mantle': 0.1e18            # 地幔脱气
            },
            
            # 输出通量 (mol/Ma)
            output_fluxes={
                'organic_burial': 3.5e18,   # 有机碳埋藏
                'carbonate_burial': 21.5e18 # 碳酸盐埋藏
            },
            
            # 端元同位素值（‰）
            end_members={
                'weathering': {
                    'delta13': -4.0,
                    'description': '风化输入（平均）'
                },
                'mantle': {
                    'delta13': -6.0,
                    'description': '地幔CO₂'
                },
                'marine_doc': {
                    'delta13': -30.0,
                    'description': '海洋DOC（溶解有机碳）'
                },
                'organic_carbon': {
                    'delta13': -30.0,
                    'description': '有机碳（与DOC相同）'
                },
                'seawater_dic': {
                    'delta13': 0.0,
                    'range': (-6, 2),
                    'description': '海水DIC（现代≈0）'
                },
                'carbonate': {
                    'delta13': 0.0,
                    'description': '碳酸盐岩'
                }
            },
            
            # 分馏参数
            fractionation_factors={
                'delta_carb_org': 30.0,      # 碳酸盐-有机碳分馏Δ¹³C
                'organic_burial_fraction': 0.14,  # f_org
                'event_duration_ma': 1.0     # 事件持续时间
            }
        )
    
    elif scenario == 'modern':
        # 现代碳循环参数
        return IsotopeParameters(
            **base_params,
            reservoir_mass=3.8e18,  # 现代DIC略低
            
            input_fluxes={
                'weathering': 20e18,
                'mantle': 0.1e18
            },
            
            output_fluxes={
                'organic_burial': 2.8e18,
                'carbonate_burial': 17.2e18
            },
            
            end_members={
                'weathering': {'delta13': -4.0},
                'mantle': {'delta13': -6.0},
                'seawater_dic': {'delta13': 0.5},
                'organic_carbon': {'delta13': -25.0}
            },
            
            fractionation_factors={
                'delta_carb_org': 25.0,
                'organic_burial_fraction': 0.14
            }
        )
    
    else:
        raise ValueError(f"Unknown scenario: {scenario}")


# 氧化剂消耗情景（Li et al. 2020）
OXIDANT_SCENARIOS = {
    'modern_o2_high_sulfate': {
        'name': 'Morden_O2+12mM_Sul',
        'k1': 0.52,  # 被O₂氧化的DOC比例
        'k2': 0.48,  # 被硫酸盐氧化的DOC比例
        'description': '现代氧气水平 + 12mM硫酸盐'
    },
    'modern_o2_low_sulfate': {
        'name': 'Morden_O2+2mM_Sul',
        'k1': 0.87,
        'k2': 0.13,
        'description': '现代氧气水平 + 2mM硫酸盐'
    },
    'low_o2_high_sulfate': {
        'name': '10per_Morden_O2+12mM_Sul',
        'k1': 0.10,
        'k2': 0.90,
        'description': '10%现代氧气 + 12mM硫酸盐'
    },
    'low_o2_low_sulfate': {
        'name': '10per_Morden_O2+2mM_sul',
        'k1': 0.39,
        'k2': 0.61,
        'description': '10%现代氧气 + 2mM硫酸盐'
    },
    'complete_oxidation': {
        'name': 'DOC_complete_Oxi',
        'k1': 1.0,
        'k2': 0.0,
        'description': 'DOC完全氧化'
    }
}


# 硫酸盐还原化学计量系数
# 1 mol有机碳被硫酸盐氧化需要 15/8 mol硫酸盐
SULFATE_REDUCTION_STOICHIOMETRY = {
    'sulfate_per_carbon': 15.0 / 8.0,
    'carbon_per_sulfate': 8.0 / 15.0
}
