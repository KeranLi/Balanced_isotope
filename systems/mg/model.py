"""
Mg同位素体系模型
实现风化-沉积体系的质量平衡和分馏计算
"""

import numpy as np
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass

from systems.base.isotope_system import IsotopeSystem, ModelResult, IsotopeParameters
from systems.mg.parameters import get_mg_parameters
from toolkit.isotope.formulas import MassBalance, EvolutionEquations
from toolkit.math.numerical import ODESolver, Interpolator


@dataclass
class MgWeatheringFluxes:
    """Mg风化通量数据类"""
    F_carb_in: float      # 碳酸盐风化输入 (mol/Ma)
    F_sil_in: float       # 硅酸盐风化输入 (mol/Ma)
    F_carb_out: float     # 碳酸盐沉积输出 (mol/Ma)
    F_sil_out: float      # 硅酸盐相关输出 (mol/Ma)
    
    # 同位素值（δ²⁶Mg，‰）
    delta_carb_in: float = -4.3   # 碳酸盐风化输入
    delta_sil_in: float = -0.3    # 硅酸盐风化输入
    delta_carb_out: float = -4.0  # 碳酸盐沉积
    delta_sil_out: float = -0.5   # 硅酸盐相关输出


class MgIsotopeSystem(IsotopeSystem):
    """
    Mg同位素体系
    
    主要应用：
    1. 海水Mg同位素演化模拟
    2. 风化通量反演
    3. 碳酸盐-硅酸盐风化比例估算
    """
    
    ELEMENT = 'mg'
    NAME = 'Magnesium'
    ISOTOPES = ['24Mg', '25Mg', '26Mg']
    
    def __init__(self, parameters: Optional[IsotopeParameters] = None):
        super().__init__(parameters or get_mg_parameters())
        
        # Mg专用缓存
        self._swpre_cache: Optional[np.ndarray] = None
    
    def _default_parameters(self) -> IsotopeParameters:
        return get_mg_parameters()
    
    # ============== 核心计算接口实现 ==============
    
    def mass_balance_equation(self,
                             state: np.ndarray,
                             fluxes: Dict[str, float],
                             time: Optional[float] = None) -> np.ndarray:
        """
        Mg同位素质量平衡微分方程
        
        state[0] = M_sw  (海水Mg浓度)
        state[1] = δ²⁶Mg_sw  (海水Mg同位素)
        
        Parameters
        ----------
        state : array_like, shape (2,)
            [M_sw, delta_sw]
        fluxes : dict
            包含 F_in, F_out, delta_in, delta_out
            
        Returns
        -------
        array_like, shape (2,)
            [dM/dt, dδ/dt]
        """
        M_sw, delta_sw = state
        
        # 总输入输出通量
        F_in_total = fluxes.get('F_in_total', 0)
        F_out_total = fluxes.get('F_out_total', 0)
        
        # 输入同位素（加权平均）
        delta_in = fluxes.get('delta_in', 0)
        
        # 输出分馏系数
        alpha_out = fluxes.get('alpha_out', 1.0)
        
        # Mg总量变化
        dM_dt = F_in_total - F_out_total
        
        # 同位素演化（简化模型）
        # dδ/dt = [F_in(δ_in - δ) + F_out·δ(α-1)] / M
        d_delta_dt = (
            F_in_total * (delta_in - delta_sw) +
            F_out_total * delta_sw * (alpha_out - 1)
        ) / M_sw if M_sw > 0 else 0
        
        return np.array([dM_dt, d_delta_dt])
    
    def fractionation_factor(self,
                            process: str,
                            temperature: Optional[float] = None,
                            **kwargs) -> float:
        """
        获取Mg分馏系数
        
        Parameters
        ----------
        process : str
            'carb_sw', 'sil_sw', 'precipitation', 'weathering'
        temperature : float, optional
            温度（K）
            
        Returns
        -------
        float
            分馏系数 α
        """
        epsilon = self.params.fractionation_factors.get(process, 0)
        
        # 温度校正（如果有）
        if temperature is not None and process in ['carb_sw', 'precipitation']:
            # 简化温度依赖：每升高100K，分馏减小约0.1‰
            temp_correction = -0.001 * (temperature - 298) / 100
            epsilon += temp_correction
        
        # ε转换为α
        return 1 + epsilon / 1000
    
    def mixing_model(self,
                    end_member_values: np.ndarray,
                    proportions: np.ndarray) -> float:
        """
        端元混合模型
        
        Parameters
        ----------
        end_member_values : array_like
            各端元的δ²⁶Mg值
        proportions : array_like
            各端元的比例（和为1）
            
        Returns
        -------
        float
            混合后的δ²⁶Mg值
        """
        return MassBalance.multi_component_mixing(
            end_member_values, proportions, 
            standard_ratio=self.params.reference_ratios['26/24']
        )
    
    # ============== Mg专用方法 ==============
    
    def calculate_weathering_ratio(self,
                                  delta_sample: float,
                                  delta_seawater: float,
                                  method: str = 'two_endmember') -> Dict[str, float]:
        """
        计算风化端元比例
        
        Parameters
        ----------
        delta_sample : float
            样品δ²⁶Mg
        delta_seawater : float
            同期海水δ²⁶Mg
        method : str
            'two_endmember' 或 'three_endmember'
            
        Returns
        -------
        dict
            各端元比例
        """
        end_members = self.params.end_members
        
        if method == 'two_endmember':
            # 简化的二元混合：碳酸盐 vs 硅酸盐
            delta_carb = end_members['carbonate']['delta26']
            delta_sil = end_members['silicate']['delta26']
            
            # 质量平衡：delta_sample = f_carb * delta_carb + (1-f_carb) * delta_sil
            # 求解 f_carb
            if abs(delta_carb - delta_sil) < 0.001:
                f_carb = 0.5  # 避免除零
            else:
                f_carb = (delta_sample - delta_sil) / (delta_carb - delta_sil)
            
            # 限制在合理范围
            f_carb = np.clip(f_carb, 0, 1)
            
            return {
                'f_carbonate': float(f_carb),
                'f_silicate': float(1 - f_carb)
            }
        
        elif method == 'three_endmember':
            # 三元混合：碳酸盐 + 硅酸盐 + 海水
            # 需要额外的约束条件
            raise NotImplementedError("Three end-member mixing not yet implemented")
        
        else:
            raise ValueError(f"Unknown method: {method}")
    
    def calculate_swpre(self,
                       rm_data: np.ndarray,
                       decay_constant: float = 0.001) -> np.ndarray:
        """
        计算风化前信号（swpre）
        这是原mass_balance_model.py中的核心算法
        
        Parameters
        ----------
        rm_data : array_like
            沉积物Mg同位素数据
        decay_constant : float
            衰减常数
            
        Returns
        -------
        array_like
            swpre值
        """
        n = len(rm_data)
        swpre = np.zeros(n)
        
        for i in range(n):
            decay = np.arange(i, n)
            Rdecay = np.exp(-decay_constant * (decay - i) / 100)
            swpre[i] = np.sum(rm_data[i:] * Rdecay) / np.sum(Rdecay)
        
        self._swpre_cache = swpre
        return swpre
    
    def seawater_evolution(self,
                          time_span: Tuple[float, float],
                          initial_delta: float = -0.5,
                          flux_scenario: str = 'modern',
                          n_points: int = 1000) -> ModelResult:
        """
        模拟海水Mg同位素演化
        
        Parameters
        ----------
        time_span : tuple
            (start_age, end_age) in Ma
        initial_delta : float
            初始海水δ²⁶Mg
        flux_scenario : str
            'modern', 'high_weathering', 'low_weathering'
        n_points : int
            时间点数
            
        Returns
        -------
        ModelResult
        """
        # 根据情景设置通量
        if flux_scenario == 'modern':
            F_carb = self.params.input_fluxes['rivers_carbonate']
            F_sil = self.params.input_fluxes['rivers_silicate']
        elif flux_scenario == 'high_weathering':
            F_carb = self.params.input_fluxes['rivers_carbonate'] * 1.5
            F_sil = self.params.input_fluxes['rivers_silicate'] * 2.0
        elif flux_scenario == 'low_weathering':
            F_carb = self.params.input_fluxes['rivers_carbonate'] * 0.5
            F_sil = self.params.input_fluxes['rivers_silicate'] * 0.5
        else:
            raise ValueError(f"Unknown scenario: {flux_scenario}")
        
        # 计算输入同位素（加权平均）
        delta_carb = self.params.end_members['carbonate']['delta26']
        delta_sil = self.params.end_members['silicate']['delta26']
        
        F_total = F_carb + F_sil
        delta_in = (F_carb * delta_carb + F_sil * delta_sil) / F_total
        
        # 输出通量（简化为稳态）
        F_out = F_total
        
        # 分馏系数
        alpha_out = self.fractionation_factor('carb_sw_equilibrium')
        
        fluxes = {
            'F_in_total': F_total,
            'F_out_total': F_out,
            'delta_in': delta_in,
            'alpha_out': alpha_out
        }
        
        # 初始状态 [M_sw, delta_sw]
        M_sw = self.params.reservoir_mass
        initial_state = np.array([M_sw, initial_delta])
        
        return self.time_evolution(initial_state, time_span, fluxes, n_points)
    
    def inverse_weathering_flux(self,
                               age_data: np.ndarray,
                               delta_data: np.ndarray,
                               delta_uncertainty: Optional[np.ndarray] = None) -> ModelResult:
        """
        从观测数据反演风化通量历史
        
        Parameters
        ----------
        age_data : array_like
            年龄数据（Ma）
        delta_data : array_like
            δ²⁶Mg观测值
        delta_uncertainty : array_like, optional
            观测误差
            
        Returns
        -------
        ModelResult
            包含反演的风化通量历史
        """
        # 简化的反演：假设稳态，求解通量比
        n = len(age_data)
        f_carb_history = np.zeros(n)
        
        delta_carb = self.params.end_members['carbonate']['delta26']
        delta_sil = self.params.end_members['silicate']['delta26']
        
        for i in range(n):
            # 二元混合反演
            delta = delta_data[i]
            if abs(delta_carb - delta_sil) > 0.001:
                f_carb = (delta - delta_sil) / (delta_carb - delta_sil)
                f_carb_history[i] = np.clip(f_carb, 0, 1)
        
        result = ModelResult(success=True)
        result.add('age', age_data)
        result.add('f_carbonate', f_carb_history)
        result.add('f_silicate', 1 - f_carb_history)
        
        return result
    
    # ============== 工具方法实现 ==============
    
    def state_dimension(self) -> int:
        """状态变量维度：M_sw, δ²⁶Mg_sw"""
        return 2
    
    def validate_data(self, data: Dict[str, np.ndarray]) -> Tuple[bool, str]:
        """
        验证Mg同位素数据
        """
        required_fields = ['delta_26_mg']
        
        for field in required_fields:
            if field not in data and 'delta_26_Mg' not in data:
                return False, f"Missing required field: {field}"
        
        # 检查数值范围
        delta_values = data.get('delta_26_mg', data.get('delta_26_Mg'))
        if np.any(delta_values < -10) or np.any(delta_values > 5):
            return False, "Delta values out of reasonable range (-10 to 5)"
        
        return True, "Validation passed"
