"""
Mg同位素体系模块
"""

from systems.mg.model import MgIsotopeSystem
from systems.mg.parameters import get_mg_parameters, get_ancient_parameters

__all__ = ['MgIsotopeSystem', 'get_mg_parameters', 'get_ancient_parameters']
