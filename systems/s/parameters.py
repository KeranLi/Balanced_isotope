"""
S同位素参数模块（模板）
多硫同位素体系（³²S, ³³S, ³⁴S, ³⁶S）
"""

from systems.base.isotope_system import IsotopeParameters


def get_s_parameters() -> IsotopeParameters:
    """
    获取S同位素体系参数
    
    TODO: 补充完整参数
    """
    return IsotopeParameters(
        element='s',
        name='Sulfur',
        
        # 参考标准
        reference_standard='VCDT',
        reference_ratios={
            '33/32': 0.007877,
            '34/32': 0.0441626,
            '36/32': 0.0001547
        },
        
        # 分馏系数
        fractionation_factors={
            'sulfate_reduction': 25.0,    # 硫酸盐还原
            'sulfide_oxidation': -5.0     # 硫化物氧化
        },
        
        # 端元值
        end_members={
            'seawater_sulfate': {
                'delta34': 20.0,
                'delta33': 0.0,  # Δ³³S ≈ 0（现代海水）
                'description': '现代海水硫酸盐'
            },
            'mantle': {
                'delta34': 0.0,
                'description': '地幔硫'
            }
        },
        
        # 储库参数
        reservoir_mass=2.9e19,  # 海水硫酸盐 (mol)
        
        input_fluxes={
            'rivers': 1.0e18,
            'hydrothermal': 0.5e18
        },
        
        output_fluxes={
            'pyrite_burial': 1.0e18,
            'gypsum_burial': 0.5e18
        }
    )
