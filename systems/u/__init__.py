"""
铀同位素体系 (Uranium Isotope System)

实现海洋铀循环的质量平衡模型，用于研究古海洋氧化还原条件变化。

主要功能:
- 稳态模型: 从碳酸盐δ²³⁸U计算缺氧汇比例 f_anox
- 非稳态模型: 模拟铀同位素的时间演化 (ODE积分)
- 反演模型: 从观测数据推断古环境参数

示例用法:
    from systems.u import UIsotopeSystem
    
    # 创建体系实例
    u_system = UIsotopeSystem(scenario='modern')
    
    # 稳态计算: 从碳酸盐δ²³⁸U计算 f_anox
    result = u_system.calculate_f_anox_steady_state(delta238_carb=-0.65)
    print(f"缺氧汇比例: {result['f_anox']:.1%}")
    
    # 非稳态模拟: 模拟缺氧事件
    transient_result = u_system.simulate_anoxic_event(
        event_duration=1.0,  # 1 Myr
        peak_f_anox=0.8
    )
"""

from systems.u.model import UIsotopeSystem, UraniumCycleFluxes, FractionationFactors
from systems.u.parameters import (
    get_u_parameters,
    get_scenario_info,
    get_fractionation_ranges,
    URANIUM_STANDARDS,
    SCENARIO_PARAMETERS
)
from systems.u.uncertainty import (
    UncertaintyAnalyzer,
    UncertaintyConfig,
    ParameterUncertainty,
    analyze_likelihood
)

__all__ = [
    'UIsotopeSystem',
    'UraniumCycleFluxes',
    'FractionationFactors',
    'UncertaintyAnalyzer',
    'UncertaintyConfig',
    'ParameterUncertainty',
    'analyze_likelihood',
    'get_u_parameters',
    'get_scenario_info',
    'get_fractionation_ranges',
    'URANIUM_STANDARDS',
    'SCENARIO_PARAMETERS',
]
