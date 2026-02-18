"""
同位素体系层
各元素同位素体系的实现
"""

from systems.base.isotope_system import (
    IsotopeSystem,
    MultiIsotopeSystem,
    RadiogenicSystem,
    IsotopeParameters,
    ModelResult
)

# 延迟导入具体体系（避免循环导入）
def get_mg_system():
    from systems.mg.model import MgIsotopeSystem
    return MgIsotopeSystem

def get_c_system():
    from systems.c.model import CIsotopeSystem
    return CIsotopeSystem

def get_n_system():
    from systems.n.model import NIsotopeSystem
    return NIsotopeSystem

__all__ = [
    'IsotopeSystem',
    'MultiIsotopeSystem', 
    'RadiogenicSystem',
    'IsotopeParameters',
    'ModelResult',
    'get_mg_system',
    'get_c_system',
    'get_n_system'
]
