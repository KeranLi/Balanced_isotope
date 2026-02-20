"""
Mg同位素风化通量模型
基于Kasemann等(2014)论文的海洋箱模型实现

核心功能：
1. 双端元风化模型（硅酸盐 vs 碳酸盐）
2. 海水Mg同位素演化模拟
3. 风化通量反演
4. 风化比例时间演化
"""

import numpy as np
from typing import Dict, Tuple, Optional, List, Callable, Union
from dataclasses import dataclass
from scipy.optimize import minimize

from systems.base.isotope_system import IsotopeSystem, ModelResult, IsotopeParameters
from systems.mg.parameters import (
    get_mg_parameters, get_cryogenian_parameters,
    calculate_river_delta26, solve_f_silicate,
    weathering_transition_linear, weathering_transition_exponential,
    weathering_pulse
)
from toolkit.math.numerical import ODESolver


@dataclass
class WeatheringFluxConfig:
    """风化通量配置"""
    f_silicate: float = 0.5           # 硅酸盐风化比例
    F_riv_multiplier: float = 1.0     # 河流总通量放大系数（相对于现代）
    
    # 端元值
    delta_silicate: float = -0.3      # 硅酸盐端元δ²⁶Mg
    delta_carbonate: float = -2.5     # 碳酸盐端元δ²⁶Mg
    
    # 分馏系数
    Delta_carb: float = -2.7          # 碳酸盐沉淀分馏（‰）


class MgWeatheringModel:
    """
    Mg同位素风化模型（碳酸盐体系）
    实现论文中的双端元混合和质量平衡计算
    
    适用场景：
    - 分析海相碳酸盐岩的 Mg 同位素
    - 从碳酸盐岩 δ²⁶Mg 反演古海水组成
    - 模拟海洋 Mg 循环演化
    
    与碎屑岩体系的区别：
    - 输入数据：碳酸盐岩 δ²⁶Mg vs 碎屑岩 δ²⁶Mg
    - 模型原理：海水沉淀分馏 vs 风化残余分馏
    - 应用对象：海相碳酸盐岩 vs 陆源碎屑沉积物
    """
    
    def __init__(self, params: Optional[IsotopeParameters] = None):
        self.params = params or get_mg_parameters()
        self.M_sw = self.params.reservoir_mass  # 海水Mg储库 (mol)
        self.F_riv_0 = self.params.input_fluxes['river_total']  # 现代河流通量
    
    def river_input(self, config: WeatheringFluxConfig) -> Tuple[float, float]:
        """
        计算河流输入通量和同位素组成
        
        Returns
        -------
        (F_riv, delta_riv) : tuple
            河流通量 (mol/yr) 和 δ²⁶Mg值
        """
        F_riv = config.F_riv_multiplier * self.F_riv_0
        delta_riv = calculate_river_delta26(
            f_silicate=config.f_silicate,
            delta_silicate=config.delta_silicate,
            delta_carbonate=config.delta_carbonate
        )
        return F_riv, delta_riv
    
    def steady_state_seawater(self, config: WeatheringFluxConfig) -> float:
        """
        计算稳态海水δ²⁶Mg值
        
        δ_sw = δ_riv - (F_carb/F_riv) × Δ_carb
        
        假设：F_carb ≈ F_riv - F_hydro，且F_hydro定量移除Mg无分馏
        
        Parameters
        ----------
        config : WeatheringFluxConfig
            风化配置
            
        Returns
        -------
        float
            稳态海水δ²⁶Mg
        """
        F_riv, delta_riv = self.river_input(config)
        F_hydro = self.params.input_fluxes.get('hydrothermal', 0.2 * self.F_riv_0)
        F_carb = F_riv - F_hydro  # 碳酸盐沉淀通量
        
        # 稳态方程
        delta_sw = delta_riv - (F_carb / F_riv) * config.Delta_carb
        return delta_sw
    
    def derivative(self, state: np.ndarray, t: float, 
                   config_func: Callable[[float], WeatheringFluxConfig]) -> np.ndarray:
        """
        计算状态变量的时间导数
        
        dM_sw/dt = F_riv - F_carb - F_hydro
        dδ_sw/dt = [F_riv(δ_riv - δ_sw) - F_carb·Δ_carb] / M_sw
        
        Parameters
        ----------
        state : array_like, shape (2,)
            [M_sw, delta_sw]
        t : float
            时间
        config_func : callable
            返回当前时间WeatheringFluxConfig的函数
            
        Returns
        -------
        array_like, shape (2,)
            [dM/dt, dδ/dt]
        """
        M_sw, delta_sw = state
        
        # 获取当前配置
        config = config_func(t)
        
        # 河流输入
        F_riv, delta_riv = self.river_input(config)
        
        # 热液输出（定量移除，无分馏）
        F_hydro = self.params.input_fluxes.get('hydrothermal', 0.2 * self.F_riv_0)
        
        # 碳酸盐沉淀（与海水平衡，产生分馏）
        F_carb = max(0, F_riv - F_hydro)  # 确保非负
        
        # 质量平衡
        dM_dt = F_riv - F_carb - F_hydro
        
        # 同位素演化
        if M_sw > 0:
            d_delta_dt = (
                F_riv * (delta_riv - delta_sw) 
                - F_carb * config.Delta_carb
            ) / M_sw
        else:
            d_delta_dt = 0
        
        return np.array([dM_dt, d_delta_dt])


class MgIsotopeSystem(IsotopeSystem):
    """
    Mg同位素体系 - 风化通量模型（碳酸盐体系）
    
    基于Kasemann等(2014)论文实现：
    - 双端元风化模型（硅酸盐 vs 碳酸盐）
    - 海洋箱模型演化
    - 风化通量反演
    
    适用场景：
    - 分析海相碳酸盐岩的 Mg 同位素
    - 从碳酸盐岩 δ²⁶Mg 反演古海水组成
    - 模拟海洋 Mg 循环演化
    
    与碎屑岩体系的区别：
    - 输入数据：碳酸盐岩 δ²⁶Mg vs 碎屑岩 δ²⁶Mg
    - 模型原理：海水沉淀分馏 vs 风化残余分馏
    - 应用对象：海相碳酸盐岩 vs 陆源碎屑沉积物
    """
    
    ELEMENT = 'mg'
    NAME = 'Magnesium'
    ISOTOPES = ['24Mg', '25Mg', '26Mg']
    COMPONENT_TYPE = 'carbonate'  # 体系类型标识
    
    def __init__(self, parameters: Optional[IsotopeParameters] = None,
                 scenario: str = 'modern'):
        """
        初始化Mg同位素体系
        
        Parameters
        ----------
        parameters : IsotopeParameters, optional
            自定义参数
        scenario : str
            'modern' - 现代参数
            'cryogenian' - Cryogenian时期参数
        """
        if parameters is None:
            if scenario == 'cryogenian':
                parameters = get_cryogenian_parameters()
            else:
                parameters = get_mg_parameters()
        
        super().__init__(parameters)
        
        # 初始化风化模型
        self.weathering_model = MgWeatheringModel(self.params)
        
        # 默认端元值
        self._delta_sil = self.params.end_members['silicate']['delta26']
        self._delta_carb = self.params.end_members['carbonate']['delta26']
    
    def _default_parameters(self) -> IsotopeParameters:
        return get_mg_parameters()
    
    # ============== 核心接口实现 ==============
    
    def mass_balance_equation(self,
                             state: np.ndarray,
                             fluxes: Dict[str, float],
                             time: Optional[float] = None) -> np.ndarray:
        """
        质量平衡微分方程（基类接口）
        
        注意：推荐使用evolve_seawater()进行演化计算
        """
        M_sw, delta_sw = state
        
        F_riv = fluxes.get('F_riv', self.params.input_fluxes['river_total'])
        delta_riv = fluxes.get('delta_riv', -1.2)
        F_carb = fluxes.get('F_carb', F_riv * 0.8)
        Delta_carb = fluxes.get('Delta_carb', -2.7)
        
        dM_dt = fluxes.get('dM_dt', 0)  # 通常假设储库恒定
        
        if M_sw > 0:
            d_delta_dt = (F_riv * (delta_riv - delta_sw) - F_carb * Delta_carb) / M_sw
        else:
            d_delta_dt = 0
        
        return np.array([dM_dt, d_delta_dt])
    
    def fractionation_factor(self,
                            process: str,
                            temperature: Optional[float] = None,
                            **kwargs) -> float:
        """获取分馏系数"""
        epsilon = self.params.fractionation_factors.get(process, 0)
        return 1 + epsilon / 1000
    
    def mixing_model(self,
                    end_member_values: np.ndarray,
                    proportions: np.ndarray) -> float:
        """端元混合模型"""
        proportions = np.array(proportions)
        proportions = proportions / np.sum(proportions)  # 归一化
        return float(np.sum(end_member_values * proportions))
    
    def state_dimension(self) -> int:
        """状态维度"""
        return 2
    
    # ============== Mg专用方法 ==============
    
    def calculate_weathering_ratio(self,
                                  delta_sample: float,
                                  delta_seawater: Optional[float] = None,
                                  method: str = 'two_endmember') -> Dict[str, float]:
        """
        计算风化端元比例
        
        Parameters
        ----------
        delta_sample : float
            样品δ²⁶Mg（沉积碳酸盐）
        delta_seawater : float, optional
            同期海水δ²⁶Mg，默认使用稳态值
        method : str
            'two_endmember' - 碳酸盐 vs 硅酸盐
            'from_river' - 从河流组成计算
            
        Returns
        -------
        dict
            {
                'f_silicate': 硅酸盐风化比例,
                'f_carbonate': 碳酸盐风化比例,
                'delta_river': 估算的河流δ²⁶Mg
            }
        """
        if delta_seawater is None:
            delta_seawater = self.params.end_members['seawater']['delta26']
        
        # 从碳酸盐δ²⁶Mg反演海水δ²⁶Mg
        # δ_carb = δ_sw + Δ_carb  →  δ_sw = δ_carb - Δ_carb
        delta_riv_est = delta_sample - (-2.7)  # 简化：假设分馏系数-2.7‰
        
        if method == 'two_endmember':
            f_sil = solve_f_silicate(delta_riv_est, self._delta_sil, self._delta_carb)
            
            return {
                'f_silicate': float(f_sil),
                'f_carbonate': float(1 - f_sil),
                'delta_river': float(delta_riv_est)
            }
        
        elif method == 'from_river':
            # 直接给定河流组成
            f_sil = solve_f_silicate(delta_sample, self._delta_sil, self._delta_carb)
            return {
                'f_silicate': float(f_sil),
                'f_carbonate': float(1 - f_sil),
                'delta_river': float(delta_sample)
            }
        
        else:
            raise ValueError(f"Unknown method: {method}")
    
    def evolve_seawater(self,
                       time_span: Tuple[float, float],
                       f_silicate_func: Callable[[float], float],
                       flux_multiplier_func: Optional[Callable[[float], float]] = None,
                       initial_delta: float = -0.83,
                       n_points: int = 1000) -> ModelResult:
        """
        模拟海水Mg同位素演化
        
        核心方程：
        dδ_sw/dt = [F_riv(δ_riv - δ_sw) - F_carb·Δ_carb] / M_sw
        
        其中：
        - δ_riv = f_sil × δ_sil + (1-f_sil) × δ_carb
        - F_riv = φ(t) × F_riv_0 （φ为通量放大系数）
        
        Parameters
        ----------
        time_span : tuple
            (t_start, t_end)，单位：年
        f_silicate_func : callable
            硅酸盐风化比例函数 f_sil(t)
        flux_multiplier_func : callable, optional
            通量放大系数函数 φ(t)，默认恒为1
        initial_delta : float
            初始海水δ²⁶Mg
        n_points : int
            时间点数
            
        Returns
        -------
        ModelResult
            包含时间序列的结果
        """
        if flux_multiplier_func is None:
            flux_multiplier_func = lambda t: 1.0
        
        # 创建配置函数
        def config_func(t: float) -> WeatheringFluxConfig:
            return WeatheringFluxConfig(
                f_silicate=f_silicate_func(t),
                F_riv_multiplier=flux_multiplier_func(t),
                delta_silicate=self._delta_sil,
                delta_carbonate=self._delta_carb,
                Delta_carb=-2.7
            )
        
        # 初始状态 [M_sw, delta_sw]
        initial_state = np.array([self.params.reservoir_mass, initial_delta])
        
        # 定义ODE
        def ode_func(t: float, y: np.ndarray) -> np.ndarray:
            return self.weathering_model.derivative(y, t, config_func)
        
        # 求解
        solver = ODESolver()
        result = solver.solve(ode_func, initial_state, time_span, n_points=n_points)
        
        if result.success:
            # 计算伴随的通量信息
            times = result.t
            f_sil_history = np.array([f_silicate_func(t) for t in times])
            flux_mult_history = np.array([flux_multiplier_func(t) for t in times])
            delta_riv_history = np.array([
                calculate_river_delta26(f, self._delta_sil, self._delta_carb)
                for f in f_sil_history
            ])
            
            return ModelResult(
                success=True,
                time=times,
                values=result.y,
                data={
                    'M_sw': result.y[:, 0],
                    'delta_sw': result.y[:, 1],
                    'f_silicate': f_sil_history,
                    'f_carbonate': 1 - f_sil_history,
                    'flux_multiplier': flux_mult_history,
                    'delta_river': delta_riv_history
                }
            )
        else:
            return ModelResult(success=False, message=result.message)
    
    def simulate_weathering_transition(self,
                                      time_span: Tuple[float, float],
                                      transition: Dict[str, float],
                                      initial_delta: float = -0.83,
                                      n_points: int = 1000) -> ModelResult:
        """
        模拟风化转变情景
        
        Parameters
        ----------
        time_span : tuple
            (t_start, t_end)，单位：Ma
        transition : dict
            {
                't_start': 转变开始时间 (Ma),
                't_end': 转变结束时间 (Ma),
                'f_initial': 初始硅酸盐比例,
                'f_final': 最终硅酸盐比例,
                'mode': 'linear' or 'exponential'
            }
        initial_delta : float
            初始海水δ²⁶Mg
        n_points : int
            时间点数
            
        Returns
        -------
        ModelResult
        """
        t_start = transition.get('t_start', time_span[0])
        t_end = transition.get('t_end', time_span[1])
        f_initial = transition.get('f_initial', 0.2)  # 初始以碳酸盐为主
        f_final = transition.get('f_final', 0.8)      # 最终以硅酸盐为主
        mode = transition.get('mode', 'linear')
        
        # 转换为年
        t_start_yr = t_start * 1e6
        t_end_yr = t_end * 1e6
        time_span_yr = (time_span[0] * 1e6, time_span[1] * 1e6)
        
        if mode == 'linear':
            f_func = lambda t: weathering_transition_linear(
                t, t_start_yr, t_end_yr, f_initial, f_final
            )
        elif mode == 'exponential':
            tau = (t_end_yr - t_start_yr) / 3  # 3τ达到95%
            f_func = lambda t: weathering_transition_exponential(
                t - t_start_yr, tau, f_initial, f_final
            )
        else:
            raise ValueError(f"Unknown mode: {mode}")
        
        return self.evolve_seawater(
            time_span=time_span_yr,
            f_silicate_func=f_func,
            initial_delta=initial_delta,
            n_points=n_points
        )
    
    def simulate_cryogenian_scenario(self,
                                    duration_ma: float = 3.0,
                                    n_points: int = 1000) -> ModelResult:
        """
        模拟Cryogenian冰期后情景（基于Kasemann等2014论文）
        
        情景设定：
        - 阶段1（0-0.5 Myr）：混合风化，高总通量（9×现代）
        - 阶段2（0.5-1.5 Myr）：硅酸盐主导风化（6×现代）
        
        Parameters
        ----------
        duration_ma : float
            总模拟时长（Myr）
        n_points : int
            时间点数
            
        Returns
        -------
        ModelResult
        """
        # 定义风化比例随时间变化
        def f_silicate_func(t: float) -> float:
            t_ma = t / 1e6
            if t_ma < 0.5:
                return 0.5  # 混合风化
            else:
                return 0.8  # 硅酸盐主导
        
        # 定义总通量放大系数
        def flux_multiplier_func(t: float) -> float:
            t_ma = t / 1e6
            if t_ma < 0.5:
                return 9.0   # 高风化通量
            else:
                return 6.0   # 中等硅酸盐风化
        
        # 初始条件：假设冰期后以碳酸盐风化为主
        initial_delta = calculate_river_delta26(0.2, self._delta_sil, self._delta_carb)
        initial_delta -= 2.7  # 减去分馏得到海水值
        
        return self.evolve_seawater(
            time_span=(0, duration_ma * 1e6),
            f_silicate_func=f_silicate_func,
            flux_multiplier_func=flux_multiplier_func,
            initial_delta=initial_delta,
            n_points=n_points
        )
    
    def inverse_weathering_flux(self,
                               age_data: np.ndarray,
                               delta_carb_data: np.ndarray,
                               delta_uncertainty: Optional[np.ndarray] = None,
                               assume_steady_state: bool = False) -> ModelResult:
        """
        从观测数据反演风化通量历史
        
        Parameters
        ----------
        age_data : array_like
            年龄数据（Ma）
        delta_carb_data : array_like
            碳酸盐δ²⁶Mg观测值
        delta_uncertainty : array_like, optional
            观测误差
        assume_steady_state : bool
            是否假设稳态（简化计算）
            
        Returns
        -------
        ModelResult
            包含反演结果
        """
        n = len(age_data)
        
        # 反演硅酸盐风化比例
        f_sil_history = np.zeros(n)
        delta_riv_history = np.zeros(n)
        
        for i in range(n):
            # δ_carb = δ_sw + Δ_carb
            # 假设沉积碳酸盐与海水平衡
            delta_sw = delta_carb_data[i] - (-2.7)  # +2.7‰分馏
            
            if assume_steady_state:
                # 稳态：δ_sw = δ_riv - (F_carb/F_riv) × Δ_carb
                # 简化：假设 δ_riv ≈ δ_sw + 0.5 × Δ_carb
                delta_riv = delta_sw + 0.5 * (-2.7)
            else:
                # 非稳态：直接使用δ_sw作为河流估算
                delta_riv = delta_sw
            
            delta_riv_history[i] = delta_riv
            f_sil_history[i] = solve_f_silicate(delta_riv, self._delta_sil, self._delta_carb)
        
        result = ModelResult(success=True)
        result.add('age_ma', age_data)
        result.add('f_silicate', f_sil_history)
        result.add('f_carbonate', 1 - f_sil_history)
        result.add('delta_river_inferred', delta_riv_history)
        
        if delta_uncertainty is not None:
            # 误差传播（简化）
            result.add('f_silicate_error', delta_uncertainty / abs(self._delta_sil - self._delta_carb))
        
        return result
    
    def validate_data(self, data: Dict[str, np.ndarray]) -> Tuple[bool, str]:
        """验证Mg同位素数据"""
        if 'delta_26_mg' not in data and 'delta26Mg' not in data:
            return False, "Missing required field: delta_26_mg or delta26Mg"
        
        delta_values = data.get('delta_26_mg', data.get('delta26Mg'))
        
        if np.any(delta_values < -10) or np.any(delta_values > 5):
            return False, "Delta values out of reasonable range (-10 to 5 ‰)"
        
        return True, "Validation passed"
