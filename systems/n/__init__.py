"""
氮同位素体系 (Nitrogen Isotope System)

基于 Kang et al. (2023) 和 Ma et al. (2025) 的双箱稳态氮循环模型

主要功能:
1. 正向模型: 从硝酸盐占比(f_assimilator)计算沉积物氮同位素(δ¹⁵N_sed)
2. 反向模型: 从沉积物氮同位素反演硝酸盐占比
3. 蒙特卡洛模拟: 评估参数不确定性

使用示例:
    >>> from systems.n import NIsotopeSystem
    >>> 
    >>> # 创建模型实例
    >>> n_system = NIsotopeSystem(scenario='modern')
    >>> 
    >>> # 正向计算
    >>> delta15N = n_system.forward_model(f_assimilator=0.5)
    >>> print(f"δ¹⁵N_sed = {delta15N:.2f}‰")
    >>> 
    >>> # 反向反演
    >>> result = n_system.inverse_model(delta15N_sed=3.0)
    >>> print(f"f_assimilator = {result['f_assimilator']:.3f}")
    >>> 
    >>> # 蒙特卡洛不确定性分析
    >>> mc_result = n_system.monte_carlo_simulation(f_assimilator=0.2)
    >>> print(f"δ¹⁵N_sed = {mc_result['delta15N_sed_mean']:.2f} ± {mc_result['delta15N_sed_std']:.2f}‰")

应用场景:
- 新元古代真核生物崛起研究 (Kang et al. 2023)
- 早三叠世生态系统复苏研究 (Ma et al. 2025)
- 现代海洋氮循环研究
"""

from systems.n.model import NIsotopeSystem, NitrogenCycleFluxes, FractionationFactors
from systems.n.parameters import (
    get_n_parameters,
    get_scenario_info,
    get_fractionation_ranges,
    SCENARIO_PARAMETERS,
    FLUX_PARAMETERS
)

__all__ = [
    'NIsotopeSystem',
    'NitrogenCycleFluxes',
    'FractionationFactors',
    'get_n_parameters',
    'get_scenario_info',
    'get_fractionation_ranges',
    'SCENARIO_PARAMETERS',
    'FLUX_PARAMETERS'
]
