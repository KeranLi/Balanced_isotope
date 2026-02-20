"""
Mg同位素体系模块
基于Kasemann等(2014)论文的风化通量模型
"""

from systems.mg.model import MgIsotopeSystem, WeatheringFluxConfig
from systems.mg.parameters import (
    get_mg_parameters, 
    get_cryogenian_parameters,
    calculate_river_delta26,
    solve_f_silicate
)

__all__ = [
    'MgIsotopeSystem', 
    'get_mg_parameters', 
    'get_cryogenian_parameters',
    'calculate_river_delta26',
    'solve_f_silicate',
    'WeatheringFluxConfig'
]
