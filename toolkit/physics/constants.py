"""
物理化学常数模块
包含同位素地球化学中常用的物理常数和基本参数
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict


@dataclass(frozen=True)
class PhysicalConstants:
    """物理常数"""
    
    # 基本常数
    AVOGADRO: float = 6.02214076e23  # mol^-1
    BOLTZMANN: float = 1.380649e-23  # J/K
    GAS_CONSTANT: float = 8.314462618  # J/(mol·K)
    
    # 标准状况
    STANDARD_TEMPERATURE: float = 273.15  # K
    STANDARD_PRESSURE: float = 101325  # Pa
    
    # 时间单位
    YEAR_TO_SECOND: float = 365.25 * 24 * 3600
    MA_TO_YEAR: float = 1e6  # 百万年到年


@dataclass
class IsotopeConstants:
    """
    同位素相关常数基类
    每个具体同位素体系继承并扩展此类
    """
    
    # 参考标准（delta值计算基准）
    reference_standard: str = ""
    reference_ratio: float = 0.0
    
    # 自然丰度（近似值）
    natural_abundance: Dict[str, float] = None
    
    # 原子量
    atomic_weights: Dict[str, float] = None
    
    def __post_init__(self):
        if self.natural_abundance is None:
            self.natural_abundance = {}
        if self.atomic_weights is None:
            self.atomic_weights = {}


class FractionationTheory:
    """
    分馏理论基础公式
    """
    
    @staticmethod
    def equilibrium_fractionation(T: float, A: float, B: float) -> float:
        """
        平衡分馏系数与温度关系
        通用形式: 1000·ln(α) = A/T^2 + B/T
        
        Parameters
        ----------
        T : float
            温度（K）
        A, B : float
            经验系数
            
        Returns
        -------
        float
            1000·ln(α)
        """
        return A / T**2 + B / T
    
    @staticmethod
    def alpha_to_epsilon(alpha: float) -> float:
        """
        分馏系数α转换为ε值（‰）
        ε = (α - 1) × 1000
        """
        return (alpha - 1) * 1000
    
    @staticmethod
    def epsilon_to_alpha(epsilon: float) -> float:
        """
        ε值（‰）转换为分馏系数α
        α = ε/1000 + 1
        """
        return epsilon / 1000 + 1
    
    @staticmethod
    def kinetic_fractionation(D: float, beta: float = 0.5) -> float:
        """
        动力学分馏（简化模型）
        ε ∝ D^(-β)，其中D为扩散系数比
        
        Parameters
        ----------
        D : float
            扩散系数比（轻/重同位素）
        beta : float
            指数（通常0.5左右）
            
        Returns
        -------
        float
            预期的分馏值
        """
        return (D**(-beta) - 1) * 1000
    
    @staticmethod
    def mass_dependent_fractionation(delta_major: float, 
                                     mass_ratio: float,
                                     exponent: float = 0.5) -> float:
        """
        基于质量依赖关系计算次要同位素的分馏
        例如：从δ³⁰Si推算δ²⁹Si
        
        δ_minor = δ_major × (1/m_ratio)^exponent
        
        Parameters
        ----------
        delta_major : float
            主要同位素的delta值
        mass_ratio : float
            质量比（m_minor/m_major）
        exponent : float
            指数（理论值0.5）
            
        Returns
        -------
        float
            预测的同位素delta值
        """
        return delta_major * (mass_ratio ** exponent)


class ReservoirConstants:
    """
    地球化学储库常数
    用于箱模型计算
    """
    
    # 海洋相关 (mol)
    SEAWATER_MG_MOL = 5.1e19  # 海水Mg储量
    SEAWATER_CA_MOL = 1.0e19  # 海水Ca储量  
    SEAWATER_SR_MOL = 1.2e17  # 海水Sr储量
    SEAWATER_DIC_MOL = 4.0e18  # 海水DIC储量
    SEAWATER_SULFATE_MOL = 2.9e19  # 海水硫酸盐储量
    
    # 河流输入通量 (mol/Ma)
    RIVER_MG_FLUX = 5.0e18  # Mg输入通量
    RIVER_CA_FLUX = 1.3e19  # Ca输入通量
    RIVER_SR_FLUX = 3.3e16  # Sr输入通量
    RIVER_DIC_FLUX = 25e18  # DIC输入通量
    RIVER_SULFATE_FLUX = 1.0e18  # 硫酸盐输入通量
    
    # 海水停留时间 (Ma)
    RESIDENCE_TIME_MG = 10  # Mg停留时间
    RESIDENCE_TIME_CA = 1  # Ca停留时间
    RESIDENCE_TIME_SR = 5  # Sr停留时间
    RESIDENCE_TIME_C = 0.1  # C停留时间


def get_element_info(element: str) -> Dict:
    """
    获取元素的基本信息
    
    Parameters
    ----------
    element : str
        元素符号（小写）
        
    Returns
    -------
    dict
        元素信息字典
    """
    info = {
        'c': {
            'name': 'Carbon',
            'atomic_number': 6,
            'isotopes': ['12C', '13C', '14C'],
            'major_isotope': '12C',
            'minor_isotope': '13C',
            'reference': 'VPDB',
            'reference_ratio': 0.0112372
        },
        'mg': {
            'name': 'Magnesium',
            'atomic_number': 12,
            'isotopes': ['24Mg', '25Mg', '26Mg'],
            'major_isotope': '24Mg',
            'minor_isotopes': ['25Mg', '26Mg'],
            'reference': 'DSM3',
            'reference_ratios': {'25/24': 0.12663, '26/24': 0.13932}
        },
        's': {
            'name': 'Sulfur',
            'atomic_number': 16,
            'isotopes': ['32S', '33S', '34S', '36S'],
            'major_isotope': '32S',
            'minor_isotopes': ['33S', '34S', '36S'],
            'reference': 'VCDT',
            'reference_ratio': 0.0441626
        },
        'sr': {
            'name': 'Strontium',
            'atomic_number': 38,
            'isotopes': ['84Sr', '86Sr', '87Sr', '88Sr'],
            'major_isotope': '88Sr',
            'minor_isotope': '87Sr',
            'reference': 'NIST SRM 987',
            'reference_ratio': 0.71034
        },
        'nd': {
            'name': 'Neodymium',
            'atomic_number': 60,
            'isotopes': ['142Nd', '143Nd', '144Nd', '145Nd', '146Nd', '148Nd', '150Nd'],
            'major_isotope': '144Nd',
            'minor_isotope': '143Nd',
            'reference': 'CHUR',
            'reference_ratio': 0.512638
        },
        'os': {
            'name': 'Osmium',
            'atomic_number': 76,
            'isotopes': ['184Os', '186Os', '187Os', '188Os', '189Os', '190Os', '192Os'],
            'major_isotope': '188Os',
            'minor_isotope': '187Os',
            'reference': 'UMd',
            'reference_ratio': 1.0
        },
        'li': {
            'name': 'Lithium',
            'atomic_number': 3,
            'isotopes': ['6Li', '7Li'],
            'major_isotope': '7Li',
            'minor_isotope': '6Li',
            'reference': 'L-SVEC',
            'reference_ratio': 0.083062
        },
        'ca': {
            'name': 'Calcium',
            'atomic_number': 20,
            'isotopes': ['40Ca', '42Ca', '43Ca', '44Ca', '46Ca', '48Ca'],
            'major_isotope': '40Ca',
            'minor_isotopes': ['42Ca', '43Ca', '44Ca'],
            'reference': 'SRM 915a',
            'reference_ratios': {'44/40': 0.021212}
        },
        'fe': {
            'name': 'Iron',
            'atomic_number': 26,
            'isotopes': ['54Fe', '56Fe', '57Fe', '58Fe'],
            'major_isotope': '56Fe',
            'minor_isotope': '57Fe',
            'reference': 'IRMM-014',
            'reference_ratio': 0.003978
        },
        'mo': {
            'name': 'Molybdenum',
            'atomic_number': 42,
            'isotopes': ['92Mo', '94Mo', '95Mo', '96Mo', '97Mo', '98Mo', '100Mo'],
            'major_isotope': '98Mo',
            'minor_isotopes': ['95Mo', '97Mo'],
            'reference': 'NIST SRM 3134',
            'reference_ratios': {'97/95': 0.6029, '98/95': 1.5308}
        },
        'u': {
            'name': 'Uranium',
            'atomic_number': 92,
            'isotopes': ['234U', '235U', '238U'],
            'major_isotope': '238U',
            'minor_isotope': '235U',
            'reference': 'NIST SRM 960',
            'reference_ratio': 0.007253
        }
    }
    
    return info.get(element.lower(), {})
