"""
同位素数学公式模块
提供同位素地球化学中通用的数学计算
包括delta值转换、质量平衡、瑞利分馏等
"""

import numpy as np
from typing import Union, Tuple, Optional
from dataclasses import dataclass


@dataclass
class DeltaValue:
    """
    Delta值数据类
    包含值和误差信息
    """
    value: float
    error: Optional[float] = None
    unit: str = "‰"
    
    def __post_init__(self):
        if self.unit not in ["‰", "permil", "epsilon"]:
            raise ValueError(f"Unknown unit: {self.unit}")
    
    def __repr__(self):
        if self.error is not None:
            return f"δ = {self.value:.3f} ± {self.error:.3f} {self.unit}"
        return f"δ = {self.value:.3f} {self.unit}"


class DeltaCalculator:
    """
    Delta值计算工具
    """
    
    @staticmethod
    def ratio_to_delta(measured_ratio: float, 
                       standard_ratio: float) -> float:
        """
        将同位素比值转换为delta值
        
        δ = ((R_sample / R_standard) - 1) × 1000
        
        Parameters
        ----------
        measured_ratio : float
            样品测量比值
        standard_ratio : float
            标准物质比值
            
        Returns
        -------
        float
            delta值（‰）
        """
        return ((measured_ratio / standard_ratio) - 1) * 1000
    
    @staticmethod
    def delta_to_ratio(delta: float, 
                       standard_ratio: float) -> float:
        """
        将delta值转换为同位素比值
        
        R_sample = R_standard × (δ/1000 + 1)
        
        Parameters
        ----------
        delta : float
            delta值（‰）
        standard_ratio : float
            标准物质比值
            
        Returns
        -------
        float
            同位素比值
        """
        return standard_ratio * (delta / 1000 + 1)
    
    @staticmethod
    def combine_deltas(deltas: np.ndarray,
                      fractions: np.ndarray) -> float:
        """
        按组分比例混合多个delta值（质量加权平均）
        
        注意：delta值本身不能直接加权平均，需要转换回比值
        
        Parameters
        ----------
        deltas : array_like
            各组分的delta值
        fractions : array_like
            各组分的质量分数（和为1）
            
        Returns
        -------
        float
            混合后的delta值
        """
        # 简化为线性近似（当delta值较小时足够准确）
        return np.sum(deltas * fractions)


class MassBalance:
    """
    同位素质量平衡计算
    """
    
    @staticmethod
    def two_component_mixing(
        delta_A: float,
        delta_B: float,
        f_A: float,
        standard_ratio: float = 1.0
    ) -> float:
        """
        两元混合模型
        
        Parameters
        ----------
        delta_A, delta_B : float
            端元A和B的delta值
        f_A : float
            端元A的质量分数（0-1）
        standard_ratio : float
            标准物质比值（用于精确计算）
            
        Returns
        -------
        float
            混合物的delta值
        """
        if standard_ratio == 1.0:
            # 线性近似（delta值较小时）
            return f_A * delta_A + (1 - f_A) * delta_B
        else:
            # 精确计算
            R_A = DeltaCalculator.delta_to_ratio(delta_A, standard_ratio)
            R_B = DeltaCalculator.delta_to_ratio(delta_B, standard_ratio)
            R_mix = f_A * R_A + (1 - f_A) * R_B
            return DeltaCalculator.ratio_to_delta(R_mix, standard_ratio)
    
    @staticmethod
    def multi_component_mixing(
        deltas: np.ndarray,
        fractions: np.ndarray,
        standard_ratio: float = 1.0
    ) -> float:
        """
        多元混合模型
        
        Parameters
        ----------
        deltas : array_like
            各端元的delta值
        fractions : array_like
            各端元的质量分数（和为1）
        standard_ratio : float
            标准物质比值
            
        Returns
        -------
        float
            混合物的delta值
        """
        if len(deltas) != len(fractions):
            raise ValueError("deltas和fractions长度必须相同")
        
        if not np.isclose(np.sum(fractions), 1.0):
            raise ValueError("fractions之和必须等于1")
        
        if standard_ratio == 1.0:
            return np.sum(deltas * fractions)
        else:
            ratios = np.array([
                DeltaCalculator.delta_to_ratio(d, standard_ratio) 
                for d in deltas
            ])
            R_mix = np.sum(ratios * fractions)
            return DeltaCalculator.ratio_to_delta(R_mix, standard_ratio)
    
    @staticmethod
    def flux_balance(
        delta_inputs: np.ndarray,
        flux_inputs: np.ndarray,
        delta_outputs: np.ndarray,
        flux_outputs: np.ndarray,
        M: float,
        d_delta_dt: float = 0,
        standard_ratio: float = 1.0
    ) -> bool:
        """
        验证稳态质量平衡
        
        d(M×δ)/dt = Σ(F_in×δ_in) - Σ(F_out×δ_out)
        
        Parameters
        ----------
        delta_inputs, flux_inputs : array_like
            输入同位素值和通量
        delta_outputs, flux_outputs : array_like
            输出同位素值和通量
        M : float
            储库质量
        d_delta_dt : float
            储库delta值变化率（稳态时为0）
        standard_ratio : float
            标准物质比值
            
        Returns
        -------
        bool
            是否满足质量平衡
        """
        if standard_ratio == 1.0:
            input_sum = np.sum(delta_inputs * flux_inputs)
            output_sum = np.sum(delta_outputs * flux_outputs)
        else:
            R_inputs = [DeltaCalculator.delta_to_ratio(d, standard_ratio) 
                       for d in delta_inputs]
            R_outputs = [DeltaCalculator.delta_to_ratio(d, standard_ratio) 
                        for d in delta_outputs]
            input_sum = np.sum(np.array(R_inputs) * flux_inputs)
            output_sum = np.sum(np.array(R_outputs) * flux_outputs)
        
        balance = input_sum - output_sum - M * d_delta_dt
        return np.isclose(balance, 0, atol=1e-10)


class RayleighFractionation:
    """
    瑞利分馏模型
    描述单步不可逆过程（如蒸发、光合作用等）
    """
    
    @staticmethod
    def residual_fraction(
        f: float,
        alpha: float,
        delta_0: float = 0.0
    ) -> float:
        """
        瑞利分馏 - 残余相delta值
        
        δ_residual = (δ_0 + 1000) × f^(α-1) - 1000
        
        Parameters
        ----------
        f : float
            残余分数（0-1）
        alpha : float
            分馏系数
        delta_0 : float
            初始delta值
            
        Returns
        -------
        float
            残余相delta值
        """
        return (delta_0 + 1000) * f**(alpha - 1) - 1000
    
    @staticmethod
    def instantaneous_product(
        f: float,
        alpha: float,
        delta_0: float = 0.0
    ) -> float:
        """
        瑞利分馏 - 瞬时产物delta值
        
        δ_product = (δ_0 + 1000) × α × f^(α-1) - 1000
        
        Parameters
        ----------
        f : float
            残余分数（0-1）
        alpha : float
            分馏系数
        delta_0 : float
            初始delta值
            
        Returns
        -------
        float
            瞬时产物delta值
        """
        return (delta_0 + 1000) * alpha * f**(alpha - 1) - 1000
    
    @staticmethod
    def accumulated_product(
        f: float,
        alpha: float,
        delta_0: float = 0.0
    ) -> float:
        """
        瑞利分馏 - 累积产物delta值
        
        δ_accumulated = (δ_0 + 1000) × (1 - f^α)/(1 - f) - 1000
        
        Parameters
        ----------
        f : float
            残余分数（0-1）
        alpha : float
            分馏系数
        delta_0 : float
            初始delta值
            
        Returns
        -------
        float
            累积产物delta值
        """
        if np.isclose(f, 1.0):
            return delta_0
        return (delta_0 + 1000) * (1 - f**alpha) / (1 - f) - 1000


class NonTraditionalFractionation:
    """
    非传统稳定同位素分馏（多同位素体系）
    如：Δ³³S, Δ³⁶S, Δ¹⁷O等
    """
    
    @staticmethod
    def capital_delta(delta_minor: float,
                     delta_major: float,
                     theta: float) -> float:
        """
        计算Capital Delta（质量无关分馏指标）
        
        Δ' = δ_minor - θ × δ_major
        
        Parameters
        ----------
        delta_minor : float
            次要同位素delta值
        delta_major : float
            主要同位素delta值
        theta : float
            质量依赖斜率
            
        Returns
        -------
        float
            Capital Delta值
        """
        return delta_minor - theta * delta_major
    
    @staticmethod
    def theoretical_theta(mass_minor: float,
                         mass_major: float,
                         exponent: float = 0.5) -> float:
        """
        计算理论质量依赖斜率
        
        θ = (1/m1 - 1/m2) / (1/m1 - 1/m3)
        
        Parameters
        ----------
        mass_minor : float
            次要同位素质量
        mass_major : float
            主要同位素质量
        exponent : float
            分馏指数
            
        Returns
        -------
        float
            理论theta值
        """
        # 简化计算，实际需要考虑三个同位素
        return (mass_major / mass_minor) ** exponent


class EvolutionEquations:
    """
    同位素演化方程
    用于箱模型和时间演化模拟
    """
    
    @staticmethod
    def box_model_evolution(
        t: float,
        delta: float,
        F_in: float,
        delta_in: float,
        F_out: float,
        M: float,
        alpha: float = 1.0
    ) -> float:
        """
        单箱模型同位素演化微分方程
        
        dδ/dt = [F_in × (δ_in - δ) + F_out × δ × (α - 1)] / M
        
        Parameters
        ----------
        t : float
            时间（仅用于接口兼容）
        delta : float
            当前储库delta值
        F_in : float
            输入通量
        delta_in : float
            输入同位素值
        F_out : float
            输出通量
        M : float
            储库质量
        alpha : float
            输出分馏系数
            
        Returns
        -------
        float
            dδ/dt
        """
        return (F_in * (delta_in - delta) + F_out * delta * (alpha - 1)) / M
    
    @staticmethod
    def radioactive_decay(
        t: float,
        parent_ratio: float,
        decay_constant: float
    ) -> float:
        """
        放射性衰变方程
        
        dP/dt = -λ × P
        
        Parameters
        ----------
        t : float
            时间
        parent_ratio : float
            母同位素比值
        decay_constant : float
            衰变常数
            
        Returns
        -------
        float
            dP/dt
        """
        return -decay_constant * parent_ratio
    
    @staticmethod
    def radiogenic_growth(
        t: float,
        daughter_ratio: float,
        parent_ratio: float,
        decay_constant: float
    ) -> float:
        """
        放射成因同位素增长方程
        
        dD/dt = λ × P
        
        Parameters
        ----------
        t : float
            时间
        daughter_ratio : float
            子同位素比值
        parent_ratio : float
            母同位素比值
        decay_constant : float
            衰变常数
            
        Returns
        -------
        float
            dD/dt
        """
        return decay_constant * parent_ratio


def epsilon_to_ratio(epsilon: float, reference_ratio: float) -> float:
    """
    ε表示法转换为同位素比值
    
    用于Nd, Hf等同位素体系
    
    Parameters
    ----------
    epsilon : float
        ε值
    reference_ratio : float
        参考比值
        
    Returns
    -------
    float
        同位素比值
    """
    return reference_ratio * (epsilon / 10000 + 1)


def ratio_to_epsilon(ratio: float, reference_ratio: float) -> float:
    """
    同位素比值转换为ε表示法
    
    Parameters
    ----------
    ratio : float
        同位素比值
    reference_ratio : float
        参考比值
        
    Returns
    -------
    float
        ε值
    """
    return ((ratio / reference_ratio) - 1) * 10000
