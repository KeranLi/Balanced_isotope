"""
C同位素体系模型
实现碳循环和DOC氧化模型
基于 Li et al. 2020, Precambrian Research
"""

import numpy as np
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass

from systems.base.isotope_system import IsotopeSystem, ModelResult, IsotopeParameters
from systems.c.parameters import get_c_parameters, OXIDANT_SCENARIOS, SULFATE_REDUCTION_STOICHIOMETRY
from toolkit.isotope.formulas import MassBalance, EvolutionEquations
from toolkit.math.numerical import ODESolver


@dataclass
class OxidantConsumption:
    """氧化剂消耗结果"""
    o2_consumption: float      # O₂消耗通量 (mol/Ma)
    sulfate_consumption: float # 硫酸盐消耗通量 (mol/Ma)
    total_consumption: float   # 总氧化剂消耗
    
    def to_dict(self) -> Dict:
        return {
            'o2': self.o2_consumption,
            'sulfate': self.sulfate_consumption,
            'total': self.total_consumption
        }


class CIsotopeSystem(IsotopeSystem):
    """
    C同位素体系
    
    主要应用：
    1. 碳循环质量平衡
    2. DOC氧化与碳同位素负漂
    3. 氧化剂消耗估算
    """
    
    ELEMENT = 'c'
    NAME = 'Carbon'
    ISOTOPES = ['12C', '13C', '14C']
    
    def __init__(self, parameters: Optional[IsotopeParameters] = None,
                 scenario: str = 'dice'):
        self.scenario = scenario
        super().__init__(parameters or get_c_parameters(scenario))
    
    def _default_parameters(self) -> IsotopeParameters:
        return get_c_parameters(self.scenario)
    
    # ============== 核心计算接口实现 ==============
    
    def mass_balance_equation(self,
                             state: np.ndarray,
                             fluxes: Dict[str, float],
                             time: Optional[float] = None) -> np.ndarray:
        """
        C同位素质量平衡微分方程
        
        state[0] = δ¹³C_org (有机碳同位素)
        基于 Li et al. 2020 公式3
        
        d(δ¹³C_org)/dt = [δ¹³C_w·F_w + δ¹³C_odoc·F_odoc
                        - (δ¹³C_org + Δ¹³C)·F_carb - δ¹³C_org·F_org] / M_DIC
        
        Parameters
        ----------
        state : array_like, shape (1,)
            [delta13C_org]
        fluxes : dict
            包含 F_w, F_odoc, delta_w, delta_odoc等
            
        Returns
        -------
        array_like
            dδ¹³C_org/dt
        """
        delta_org = state[0]
        
        # 参数
        M_DIC = self.params.reservoir_mass
        Delta_13C = self.params.fractionation_factors['delta_carb_org']
        f_org = self.params.fractionation_factors['organic_burial_fraction']
        
        # 通量
        F_w = fluxes.get('F_weathering', self.params.input_fluxes.get('weathering', 0))
        F_odoc = fluxes.get('F_odoc', 0)  # DOC再矿化通量（可变）
        
        delta_w = fluxes.get('delta_weathering', 
                           self.params.end_members['weathering']['delta13'])
        delta_odoc = fluxes.get('delta_doc', 
                               self.params.end_members['marine_doc']['delta13'])
        
        # 输出通量（稳态）
        F_out = F_w + F_odoc
        F_org = f_org * F_out
        F_carb = (1 - f_org) * F_out
        
        # 计算导数
        d_delta_org_dt = (
            delta_w * F_w + delta_odoc * F_odoc
            - (delta_org + Delta_13C) * F_carb
            - delta_org * F_org
        ) / M_DIC
        
        return np.array([d_delta_org_dt])
    
    def fractionation_factor(self,
                            process: str,
                            temperature: Optional[float] = None,
                            **kwargs) -> float:
        """
        获取C分馏系数
        
        Parameters
        ----------
        process : str
            'carb_org', 'organic_burial'
        temperature : float
            温度（不使用）
            
        Returns
        -------
        float
            分馏系数 α
        """
        if process == 'carb_org':
            epsilon = self.params.fractionation_factors.get('delta_carb_org', 30.0)
        elif process == 'organic_burial':
            # 有机碳埋藏分馏（简化）
            epsilon = -1.0  # 略富集轻同位素
        else:
            epsilon = 0
        
        return 1 + epsilon / 1000
    
    def mixing_model(self,
                    end_member_values: np.ndarray,
                    proportions: np.ndarray) -> float:
        """
        端元混合模型
        
        Parameters
        ----------
        end_member_values : array_like
            各端元的δ¹³C值
        proportions : array_like
            各端元的比例
            
        Returns
        -------
        float
            混合后的δ¹³C值
        """
        return MassBalance.multi_component_mixing(
            end_member_values, proportions,
            standard_ratio=self.params.reference_ratios['13/12']
        )
    
    # ============== C同位素专用方法 ==============
    
    def solve_steady_state(self,
                          F_odoc: float,
                          t_span: float = 10.0,
                          n_points: int = 1000) -> ModelResult:
        """
        求解特定DOC通量下的稳态碳同位素
        
        Parameters
        ----------
        F_odoc : float
            DOC再矿化通量 (mol/Ma)
        t_span : float
            积分时间范围
        n_points : int
            点数
            
        Returns
        -------
        ModelResult
            包含稳态δ¹³C_carb和δ¹³C_org
        """
        # 初始条件
        Delta_13C = self.params.fractionation_factors['delta_carb_org']
        delta_carb_initial = self.params.end_members['seawater_dic']['delta13']
        delta_org_0 = delta_carb_initial - Delta_13C
        
        # 设置通量
        fluxes = {'F_odoc': F_odoc}
        
        # 求解ODE
        from scipy.integrate import odeint
        
        def ode_func(y, t):
            result = self.mass_balance_equation(np.array(y).reshape(-1), fluxes, t)
            return float(result.item())  # 返回标量
        
        t = np.linspace(0, t_span, n_points)
        solution = odeint(ode_func, delta_org_0, t)
        
        delta_org_final = solution[-1, 0]
        delta_carb_final = delta_org_final + Delta_13C
        
        result = ModelResult(success=True)
        result.add('delta13C_carb', float(delta_carb_final))
        result.add('delta13C_org', float(delta_org_final))
        result.add('solution', solution)
        result.add('time', t)
        
        return result
    
    def calculate_oxidant_consumption(self,
                                     F_odoc: float,
                                     scenario_name: str) -> OxidantConsumption:
        """
        计算DOC氧化所需的氧化剂消耗
        
        基于 Li et al. 2020 公式：
        d(oxidants)/dt = k1·F_odoc + k2·(8/15)·F_odoc
        
        Parameters
        ----------
        F_odoc : float
            DOC通量 (mol/Ma)
        scenario_name : str
            氧化剂情景名称
            
        Returns
        -------
        OxidantConsumption
        """
        scenario = OXIDANT_SCENARIOS.get(scenario_name)
        if scenario is None:
            raise ValueError(f"Unknown scenario: {scenario_name}")
        
        k1 = scenario['k1']
        k2 = scenario['k2']
        
        # O₂消耗
        o2_consume = k1 * F_odoc
        
        # 硫酸盐消耗（考虑化学计量）
        sulfate_factor = SULFATE_REDUCTION_STOICHIOMETRY['carbon_per_sulfate']
        sulfate_consume = k2 * sulfate_factor * F_odoc
        
        return OxidantConsumption(
            o2_consumption=o2_consume,
            sulfate_consumption=sulfate_consume,
            total_consumption=o2_consume + sulfate_consume
        )
    
    def doc_excursion_model(self,
                           F_odoc_range: Tuple[float, float] = (0, 10e18),
                           n_points: int = 300) -> ModelResult:
        """
        计算不同DOC通量下的碳同位素偏移
        
        Parameters
        ----------
        F_odoc_range : tuple
            DOC通量范围 (mol/Ma)
        n_points : int
            计算点数
            
        Returns
        -------
        ModelResult
            包含完整的模型结果
        """
        F_odoc_values = np.linspace(F_odoc_range[0], F_odoc_range[1], n_points)
        
        # 计算碳同位素偏移
        delta_delta13C = np.zeros(n_points)
        
        initial_delta_carb = self.params.end_members['seawater_dic']['delta13']
        
        for i, F_odoc in enumerate(F_odoc_values):
            steady_result = self.solve_steady_state(F_odoc, t_span=10.0, n_points=500)
            if steady_result.success:
                delta_carb = steady_result.get('delta13C_carb')
                delta_delta13C[i] = delta_carb - initial_delta_carb
        
        # 计算各情景的氧化剂消耗
        oxidant_results = {}
        for scenario_key, scenario in OXIDANT_SCENARIOS.items():
            oxidant_consumption = []
            for F_odoc in F_odoc_values:
                ox = self.calculate_oxidant_consumption(F_odoc, scenario_key)
                oxidant_consumption.append(ox.total_consumption)
            oxidant_results[scenario_key] = {
                'name': scenario['name'],
                'consumption': np.array(oxidant_consumption),
                'color': self._get_scenario_color(scenario_key),
                'linestyle': self._get_scenario_linestyle(scenario_key)
            }
        
        result = ModelResult(success=True)
        result.add('F_odoc', F_odoc_values)
        result.add('delta_delta13C', delta_delta13C)
        result.add('oxidant_scenarios', oxidant_results)
        
        return result
    
    def _get_scenario_color(self, scenario_key: str) -> str:
        """获取情景绘图颜色"""
        colors = {
            'modern_o2_high_sulfate': 'black',
            'modern_o2_low_sulfate': 'black',
            'low_o2_high_sulfate': 'green',
            'low_o2_low_sulfate': 'green',
            'complete_oxidation': 'blue'
        }
        return colors.get(scenario_key, 'gray')
    
    def _get_scenario_linestyle(self, scenario_key: str) -> str:
        """获取情景绘图线型"""
        linestyles = {
            'modern_o2_high_sulfate': '-',
            'modern_o2_low_sulfate': '--',
            'low_o2_high_sulfate': '-',
            'low_o2_low_sulfate': '--',
            'complete_oxidation': '-'
        }
        return linestyles.get(scenario_key, '-')
    
    def find_doc_for_excursion(self,
                              target_excursion: float,
                              tolerance: float = 0.01) -> Dict:
        """
        查找产生特定碳同位素偏移所需的DOC通量
        
        Parameters
        ----------
        target_excursion : float
            目标偏移值（负值，如-4表示4‰负漂）
        tolerance : float
            容差
            
        Returns
        -------
        dict
            结果字典
        """
        # 使用二分查找
        F_low, F_high = 0, 20e18
        
        for _ in range(50):  # 最大迭代次数
            F_mid = (F_low + F_high) / 2
            
            result = self.solve_steady_state(F_mid, t_span=10.0, n_points=500)
            if not result.success:
                break
            
            delta_carb = result.get('delta13C_carb')
            initial_delta = self.params.end_members['seawater_dic']['delta13']
            excursion = delta_carb - initial_delta
            
            if abs(excursion - target_excursion) < tolerance:
                return {
                    'F_odoc': F_mid,
                    'excursion': excursion,
                    'iterations': _
                }
            
            if excursion > target_excursion:  # 偏移不够负，增加DOC
                F_low = F_mid
            else:
                F_high = F_mid
        
        return {
            'F_odoc': F_mid,
            'excursion': excursion,
            'converged': False
        }
    
    # ============== 工具方法实现 ==============
    
    def state_dimension(self) -> int:
        """状态变量维度：δ¹³C_org"""
        return 1
    
    def validate_data(self, data: Dict[str, np.ndarray]) -> Tuple[bool, str]:
        """验证C同位素数据"""
        if 'delta_13_c' not in data and 'delta13C' not in data:
            return False, "Missing carbon isotope data"
        
        delta_values = data.get('delta_13_c', data.get('delta13C'))
        if np.any(delta_values < -50) or np.any(delta_values > 20):
            return False, "Delta13C values out of reasonable range (-50 to 20)"
        
        return True, "Validation passed"
