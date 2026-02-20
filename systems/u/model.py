"""
铀同位素体系模型
实现海洋铀循环的稳态和非稳态模型

基于:
- Algeo et al. (2023) Marine uranium cycle modeling
- Montoya-Pino et al. (2010) Geology
- Brennecka et al. (2011) PNAS
- Lau et al. (2016) PNAS
- Elrick et al. (2017) Geology

模型特点:
1. 稳态模型: 适用于长期稳定或缓慢变化的系统
2. 非稳态模型: 适用于短期快速变化事件（积分方法）
"""

import numpy as np
from typing import Dict, Tuple, Optional, List, Union, Callable
from dataclasses import dataclass
from scipy.optimize import fsolve, minimize_scalar

from systems.base.isotope_system import IsotopeSystem, ModelResult, IsotopeParameters
from systems.u.parameters import get_u_parameters, get_fractionation_ranges


@dataclass
class UraniumCycleFluxes:
    """
    铀循环通量数据类 (单位: mol/Ma)
    """
    F_river: float = 7.4e16       # 河流输入通量
    F_total: float = 7.4e16       # 总输出通量
    
    @property
    def F_oxic(self) -> float:
        """氧化汇通量"""
        return self.F_total * (1 - self.f_anox)
    
    @property
    def F_anox(self) -> float:
        """缺氧汇通量"""
        return self.F_total * self.f_anox
    
    @property
    def f_anox(self) -> float:
        """缺氧汇比例 (默认20%)"""
        return 0.2


@dataclass
class FractionationFactors:
    """分馏系数数据类 (‰)"""
    delta_sw_ox: float = 0.0      # 氧化汇分馏
    delta_sw_anox: float = 0.77   # 缺氧汇分馏
    delta_diag: float = 0.40      # 成岩校正因子
    delta_river: float = -0.29    # 河流输入
    
    def get_alpha(self, process: str) -> float:
        """
        将δ值转换为α分馏因子
        
        α = 1 + δ/1000
        """
        delta_map = {
            'oxic': self.delta_sw_ox,
            'anox': self.delta_sw_anox,
            'anoxic': self.delta_sw_anox,
            'river': self.delta_river,
            'diag': self.delta_diag,
        }
        delta = delta_map.get(process, 0.0)
        return 1 + delta / 1000


class UIsotopeSystem(IsotopeSystem):
    """
    铀同位素体系
    
    实现海洋铀循环的质量平衡模型，包括稳态和非稳态模拟。
    
    核心方程:
    
    稳态 (S1):
        δ²³⁸U_river = (1-f_anox) × δ²³⁸U_ox + f_anox × δ²³⁸U_anox
    
    稳态 f_anox 计算 (S2):
        f_anox = (δ²³⁸U_source - δ²³⁸U_ox) / (δ²³⁸U_anox - δ²³⁸U_ox)
    
    成岩校正 (S3):
        δ²³⁸U_ox = δ²³⁸U_meas + Δ_diag
    
    非稳态质量平衡:
        dM/dt = F_river - F_oxic - F_anox
        d(M×δ_sw)/dt = F_river×δ_river - F_oxic×δ_ox - F_anox×δ_anox
    
    Attributes
    ----------
    ELEMENT : str
        元素符号 'u'
    NAME : str
        元素名称 'Uranium'
    ISOTOPES : list
        涉及的同位素 ['238U', '235U']
    """
    
    ELEMENT = 'u'
    NAME = 'Uranium'
    ISOTOPES = ['238U', '235U']
    
    # 参考值
    DELTA238U_RIVER = -0.29      # ‰，河流/地壳平均值
    DELTA238U_SEAWATER_MODERN = -0.392  # ‰，现代海水
    
    def __init__(self, 
                 parameters: Optional[IsotopeParameters] = None,
                 scenario: str = 'modern',
                 model_type: str = 'steady'):
        """
        初始化铀同位素体系
        
        Parameters
        ----------
        parameters : IsotopeParameters, optional
            自定义参数，如果为None则使用默认参数
        scenario : str
            情景名称：'modern', 'oceanic_anoxic_event', 'end_permain', 
                     'frasnian_famennian'
        model_type : str
            模型类型：'steady' (稳态) 或 'transient' (非稳态)
        """
        self.scenario = scenario
        self.model_type = model_type
        
        super().__init__(parameters or get_u_parameters(scenario))
        
        # 初始化通量和分馏系数
        self.fluxes = self._init_fluxes()
        self.fractionation = self._init_fractionation()
    
    def _default_parameters(self) -> IsotopeParameters:
        """返回默认参数"""
        return get_u_parameters(self.scenario)
    
    def _init_fluxes(self) -> UraniumCycleFluxes:
        """初始化铀循环通量"""
        params = self.params
        return UraniumCycleFluxes(
            F_river=params.input_fluxes.get('river', 7.4e16),
            F_total=params.output_fluxes.get('total', 7.4e16)
        )
    
    def _init_fractionation(self) -> FractionationFactors:
        """初始化分馏系数"""
        ff = self.params.fractionation_factors
        return FractionationFactors(
            delta_sw_ox=ff.get('delta_sw_ox', 0.0),
            delta_sw_anox=ff.get('delta_sw_anox', 0.77),
            delta_diag=ff.get('delta_diag', 0.40),
            delta_river=ff.get('delta_river', -0.29)
        )
    
    # ============== 核心计算接口实现 ==============
    
    def mass_balance_equation(self,
                             state: np.ndarray,
                             fluxes: Dict[str, float],
                             time: Optional[float] = None) -> np.ndarray:
        """
        铀同位素质量平衡微分方程 (用于非稳态模型)
        
        状态变量: [M, M×δ_sw]
        
        dM/dt = F_river - F_oxic - F_anox
        d(M×δ_sw)/dt = F_river×δ_river - F_oxic×δ_ox - F_anox×δ_anox
        
        Parameters
        ----------
        state : array_like, shape (2,)
            [M, M_delta] 其中 M_delta = M × δ_sw
        fluxes : dict
            通量字典，必须包含:
            - 'F_river': 河流输入通量
            - 'f_anox': 缺氧汇比例 (函数或常数)
            - 'delta_river': 河流同位素值
        time : float, optional
            当前时间 (Myr)
            
        Returns
        -------
        array_like, shape (2,)
            d(state)/dt
        """
        M = state[0]
        M_delta = state[1]
        
        if M <= 0:
            return np.array([0.0, 0.0])
        
        delta_sw = M_delta / M  # 当前海水δ值
        
        # 获取通量参数
        F_river = fluxes.get('F_river', self.fluxes.F_river)
        delta_river = fluxes.get('delta_river', self.fractionation.delta_river)
        
        # f_anox 可以是常数或时间函数
        f_anox_func = fluxes.get('f_anox', lambda t: 0.2)
        if callable(f_anox_func):
            f_anox = f_anox_func(time) if time is not None else 0.2
        else:
            f_anox = f_anox_func
        
        # 计算汇的同位素值
        delta_ox = delta_sw - self.fractionation.delta_sw_ox
        delta_anox = delta_sw + self.fractionation.delta_sw_anox
        
        # 计算输出通量 (假设总通量与储库质量成正比以保持稳态)
        F_total = F_river  # 稳态时输入=输出
        F_oxic = F_total * (1 - f_anox)
        F_anox = F_total * f_anox
        
        # 质量平衡方程
        dM_dt = F_river - F_oxic - F_anox
        
        # 同位素质量平衡
        dM_delta_dt = (F_river * delta_river - 
                       F_oxic * delta_ox - 
                       F_anox * delta_anox)
        
        return np.array([dM_dt, dM_delta_dt])
    
    def fractionation_factor(self,
                            process: str,
                            temperature: Optional[float] = None,
                            **kwargs) -> float:
        """
        获取铀同位素分馏系数 α
        
        Parameters
        ----------
        process : str
            过程名称：
            - 'oxic', 'ox': 氧化/亚氧化汇
            - 'anoxic', 'anox': 缺氧/硫化汇
            - 'river': 河流输入
            - 'diagenesis', 'diag': 成岩过程
        temperature : float, optional
            温度 (K)，目前不使用
        **kwargs
            其他参数
            
        Returns
        -------
        float
            分馏系数 α = 1 + δ/1000
        """
        return self.fractionation.get_alpha(process)
    
    def mixing_model(self,
                    end_member_values: np.ndarray,
                    proportions: np.ndarray) -> float:
        """
        端元混合模型
        
        Parameters
        ----------
        end_member_values : array_like
            各端元的δ²³⁸U值
        proportions : array_like
            各端元的比例（和为1）
            
        Returns
        -------
        float
            混合后的δ²³⁸U值
        """
        from toolkit.isotope.formulas import MassBalance
        
        return MassBalance.multi_component_mixing(
            end_member_values, proportions,
            standard_ratio=self.params.reference_ratios['238/235']
        )
    
    # ============== 稳态模型方法 ==============
    
    def calculate_f_anox_steady_state(self,
                                     delta238_carb: float,
                                     apply_diagenetic_correction: bool = True,
                                     delta_diag: Optional[float] = None) -> Dict:
        """
        稳态模型: 从碳酸盐δ²³⁸U计算缺氧汇比例 f_anox
        
        基于方程 (S2):
        f_anox = (δ_source - δ_ox) / (δ_anox - δ_ox)
        
        其中:
        δ_source = δ_river
        δ_ox = δ_carb_corrected - Δ_sw_ox
        δ_anox = δ_sw + Δ_sw_anox
        
        Parameters
        ----------
        delta238_carb : float
            实测碳酸盐δ²³⁸U值 (‰)
        apply_diagenetic_correction : bool
            是否应用成岩校正
        delta_diag : float, optional
            成岩校正因子，默认使用参数中的值
            
        Returns
        -------
        dict
            {
                'f_anox': 缺氧汇比例 (0-1),
                'f_oxic': 氧化汇比例 (0-1),
                'delta238_seawater': 计算的海水δ²³⁸U值,
                'delta238_carb_corrected': 校正后的碳酸盐δ²³⁸U值,
                'delta238_anoxic_sink': 缺氧汇同位素值,
                'model': 'steady_state'
            }
        """
        # 成岩校正
        # 注意：成岩作用使碳酸盐富集²³⁸U，因此测量值比原始值重
        # 校正公式：δ_corrected = δ_meas - Δ_diag
        if delta_diag is None:
            delta_diag = self.fractionation.delta_diag
        
        if apply_diagenetic_correction:
            delta238_carb_corrected = delta238_carb - delta_diag
        else:
            delta238_carb_corrected = delta238_carb
        
        # 计算海水同位素值 (假设氧化汇与海水平衡)
        # δ_ox = δ_sw - Δ_sw_ox, 当Δ_sw_ox=0时, δ_ox = δ_sw
        delta238_sw = delta238_carb_corrected + self.fractionation.delta_sw_ox
        
        # 计算各汇的同位素值
        delta238_ox = delta238_sw - self.fractionation.delta_sw_ox  # 氧化汇
        delta238_anox = delta238_sw + self.fractionation.delta_sw_anox  # 缺氧汇
        
        # 计算 f_anox (方程 S2)
        delta238_source = self.fractionation.delta_river
        
        numerator = delta238_source - delta238_ox
        denominator = delta238_anox - delta238_ox
        
        if abs(denominator) < 1e-10:
            f_anox = 0.0
        else:
            f_anox = numerator / denominator
        
        # 限制在合理范围内
        f_anox = np.clip(f_anox, 0.0, 1.0)
        f_oxic = 1 - f_anox
        
        return {
            'f_anox': float(f_anox),
            'f_oxic': float(f_oxic),
            'delta238_seawater': float(delta238_sw),
            'delta238_carb_measured': float(delta238_carb),
            'delta238_carb_corrected': float(delta238_carb_corrected),
            'delta238_oxic_sink': float(delta238_ox),
            'delta238_anoxic_sink': float(delta238_anox),
            'delta_diag': float(delta_diag),
            'model': 'steady_state'
        }
    
    def calculate_seawater_delta_steady_state(self,
                                             f_anox: float,
                                             delta_river: Optional[float] = None) -> Dict:
        """
        稳态模型: 从缺氧汇比例计算海水δ²³⁸U
        
        基于方程 (S1):
        δ_source = (1-f_anox) × δ_ox + f_anox × δ_anox
        
        Parameters
        ----------
        f_anox : float
            缺氧汇比例 (0-1)
        delta_river : float, optional
            河流输入δ值，默认使用参数中的值
            
        Returns
        -------
        dict
            {
                'delta238_seawater': 海水δ²³⁸U值,
                'delta238_carbonate': 碳酸盐δ²³⁸U值 (已校正),
                'f_anox': 缺氧汇比例,
                'model': 'steady_state'
            }
        """
        if delta_river is None:
            delta_river = self.fractionation.delta_river
        
        f_anox = np.clip(f_anox, 0.0, 1.0)
        f_oxic = 1 - f_anox
        
        # 需要迭代求解，因为δ_anox依赖于δ_sw
        # δ_source = (1-f_anox) × (δ_sw - Δ_ox) + f_anox × (δ_sw + Δ_anox)
        # δ_source = δ_sw - (1-f_anox)×Δ_ox + f_anox×Δ_anox
        # δ_sw = δ_source + (1-f_anox)×Δ_ox - f_anox×Δ_anox
        
        delta_sw = (delta_river + 
                   f_oxic * self.fractionation.delta_sw_ox - 
                   f_anox * self.fractionation.delta_sw_anox)
        
        # 碳酸盐值 (未校正)
        delta238_carb = delta_sw - self.fractionation.delta_diag
        
        return {
            'delta238_seawater': float(delta_sw),
            'delta238_carbonate': float(delta238_carb),
            'delta238_carb_corrected': float(delta_sw),
            'f_anox': float(f_anox),
            'f_oxic': float(f_oxic),
            'model': 'steady_state'
        }
    
    # ============== 非稳态模型方法 ==============
    
    def solve_transient(self,
                       initial_delta238: float,
                       time_span: Tuple[float, float],
                       f_anox_function: Union[float, Callable[[float], float]],
                       n_points: int = 1000,
                       F_river: Optional[float] = None) -> ModelResult:
        """
        非稳态模型: 求解铀同位素的时间演化
        
        使用ODE积分求解质量平衡方程。
        
        Parameters
        ----------
        initial_delta238 : float
            初始海水δ²³⁸U值 (‰)
        time_span : tuple
            (t_start, t_end)，单位 Myr
        f_anox_function : float or callable
            缺氧汇比例，可以是常数或时间函数 f(t)
        n_points : int
            输出点数
        F_river : float, optional
            河流输入通量，默认使用参数中的值
            
        Returns
        -------
        ModelResult
            包含时间序列的结果:
            - time: 时间点数组 (Myr)
            - delta_seawater: 海水δ²³⁸U演化
            - f_anox: 缺氧汇比例
            - reservoir_mass: 储库质量演化
        """
        from toolkit.math.numerical import ODESolver
        
        # 准备初始条件
        M0 = self.params.reservoir_mass
        initial_state = np.array([M0, M0 * initial_delta238])
        
        # 准备通量参数
        if F_river is None:
            F_river = self.fluxes.F_river
        
        # 确保 f_anox 是函数
        if not callable(f_anox_function):
            f_anox_val = float(f_anox_function)
            f_anox_func = lambda t: f_anox_val
        else:
            f_anox_func = f_anox_function
        
        fluxes = {
            'F_river': F_river,
            'f_anox': f_anox_func,
            'delta_river': self.fractionation.delta_river
        }
        
        # 定义ODE函数
        def ode_func(t, y):
            return self.mass_balance_equation(y, fluxes, t)
        
        # 求解ODE
        result = ODESolver.solve(
            ode_func, initial_state, time_span, 
            method='RK45', n_points=n_points
        )
        
        if not result.success:
            return ModelResult(
                success=False,
                message=f"Transient solve failed: {result.message}"
            )
        
        # 提取结果
        M_t = result.y[:, 0]
        M_delta_t = result.y[:, 1]
        delta_sw_t = M_delta_t / M_t  # 海水同位素值
        
        # 计算 f_anox 时间序列
        f_anox_t = np.array([f_anox_func(t) for t in result.t])
        
        return ModelResult(
            success=True,
            time=result.t,
            data={
                'time_myr': result.t,
                'reservoir_mass': M_t,
                'delta_seawater': delta_sw_t,
                'f_anox': f_anox_t,
                'model': 'transient'
            }
        )
    
    def solve_equilibration_time(self,
                                target_f_anox: float,
                                initial_f_anox: float = 0.2,
                                tolerance: float = 0.01) -> Dict:
        """
        计算达到新的稳态所需的平衡时间
        
        当 f_anox 发生突变时，估算系统达到新平衡态的特征时间。
        
        Parameters
        ----------
        target_f_anox : float
            目标缺氧汇比例
        initial_f_anox : float
            初始缺氧汇比例
        tolerance : float
            达到平衡的判断容差
            
        Returns
        -------
        dict
            {
                'equilibration_time': 平衡特征时间 (Myr),
                'residence_time': 铀停留时间 (Myr),
                'initial_delta': 初始δ值,
                'final_delta': 最终平衡δ值
            }
        """
        # 铀停留时间
        tau = self.params.reservoir_mass / self.fluxes.F_river
        
        # 计算初始和最终平衡态
        result_initial = self.calculate_seawater_delta_steady_state(initial_f_anox)
        result_final = self.calculate_seawater_delta_steady_state(target_f_anox)
        
        # 达到 95% 平衡的时间 ≈ 3τ
        equilibration_time = 3 * tau
        
        return {
            'equilibration_time': float(equilibration_time),
            'residence_time': float(tau),
            'initial_delta': result_initial['delta238_seawater'],
            'final_delta': result_final['delta238_seawater'],
            'delta_change': float(result_final['delta238_seawater'] - 
                                 result_initial['delta238_seawater'])
        }
    
    def simulate_anoxic_event(self,
                             event_duration: float,
                             peak_f_anox: float,
                             background_f_anox: float = 0.2,
                             n_points: int = 1000) -> ModelResult:
        """
        模拟一次缺氧事件的非稳态响应
        
        使用高斯型脉冲函数描述 f_anox 的变化。
        
        Parameters
        ----------
        event_duration : float
            事件持续时间 (Myr)
        peak_f_anox : float
            峰值缺氧汇比例
        background_f_anox : float
            背景缺氧汇比例
        n_points : int
            输出点数
            
        Returns
        -------
        ModelResult
            包含完整演化过程的结果
        """
        # 定义时间跨度 (事件前后各2倍持续时间)
        t_start = -2 * event_duration
        t_end = 3 * event_duration
        
        # 定义 f_anox 随时间变化 (高斯脉冲)
        event_center = 0.5 * event_duration
        sigma = event_duration / 4  # 脉冲宽度
        
        def f_anox_func(t):
            # 背景值 + 脉冲
            pulse = (peak_f_anox - background_f_anox) * np.exp(
                -0.5 * ((t - event_center) / sigma) ** 2
            )
            return background_f_anox + pulse
        
        # 获取背景稳态作为初始条件
        bg_result = self.calculate_seawater_delta_steady_state(background_f_anox)
        initial_delta = bg_result['delta238_seawater']
        
        # 求解
        return self.solve_transient(
            initial_delta238=initial_delta,
            time_span=(t_start, t_end),
            f_anox_function=f_anox_func,
            n_points=n_points
        )
    
    # ============== 批量处理和反演方法 ==============
    
    def batch_steady_state_calculation(self,
                                      delta238_carb_array: np.ndarray,
                                      apply_diagenetic_correction: bool = True,
                                      delta_diag: Optional[float] = None) -> Dict:
        """
        批量稳态计算
        
        Parameters
        ----------
        delta238_carb_array : array_like
            碳酸盐δ²³⁸U值数组
        apply_diagenetic_correction : bool
            是否应用成岩校正
        delta_diag : float, optional
            成岩校正因子
            
        Returns
        -------
        dict
            {
                'f_anox': 缺氧汇比例数组,
                'delta_seawater': 海水δ²³⁸U数组,
                ...其他字段
            }
        """
        results = {
            'f_anox': [],
            'f_oxic': [],
            'delta_seawater': [],
            'delta_carb_corrected': []
        }
        
        for delta_carb in delta238_carb_array:
            result = self.calculate_f_anox_steady_state(
                delta_carb,
                apply_diagenetic_correction=apply_diagenetic_correction,
                delta_diag=delta_diag
            )
            results['f_anox'].append(result['f_anox'])
            results['f_oxic'].append(result['f_oxic'])
            results['delta_seawater'].append(result['delta238_seawater'])
            results['delta_carb_corrected'].append(result['delta238_carb_corrected'])
        
        # 转换为数组
        for key in results:
            results[key] = np.array(results[key])
        
        return results
    
    def inverse_model_steady_state(self,
                                  observed_delta238_carb: float,
                                  uncertainty: float = 0.1,
                                  n_monte_carlo: int = 10000) -> Dict:
        """
        稳态反演模型: 从观测数据推断 f_anox 及不确定性
        
        考虑分馏系数的不确定性。
        
        Parameters
        ----------
        observed_delta238_carb : float
            观测的碳酸盐δ²³⁸U值
        uncertainty : float
            观测不确定性 (1σ, ‰)
        n_monte_carlo : int
            蒙特卡洛采样次数
            
        Returns
        -------
        dict
            {
                'f_anox_mean': 平均f_anox,
                'f_anox_std': 标准差,
                'f_anox_ci95': 95%置信区间,
                'samples': 所有样本结果
            }
        """
        # 获取分馏系数范围
        frac_ranges = get_fractionation_ranges()
        
        # 蒙特卡洛采样
        samples = []
        for _ in range(n_monte_carlo):
            # 采样观测值 (考虑测量误差)
            delta_carb_sample = np.random.normal(
                observed_delta238_carb, uncertainty
            )
            
            # 采样分馏系数
            delta_diag_sample = np.random.uniform(*frac_ranges['delta_diag'])
            delta_anox_sample = np.random.uniform(*frac_ranges['delta_sw_anox'])
            
            # 临时修改分馏系数
            original_diag = self.fractionation.delta_diag
            original_anox = self.fractionation.delta_sw_anox
            
            self.fractionation.delta_diag = delta_diag_sample
            self.fractionation.delta_sw_anox = delta_anox_sample
            
            # 计算 f_anox
            result = self.calculate_f_anox_steady_state(delta_carb_sample)
            samples.append(result['f_anox'])
            
            # 恢复原值
            self.fractionation.delta_diag = original_diag
            self.fractionation.delta_sw_anox = original_anox
        
        samples = np.array(samples)
        
        return {
            'f_anox_mean': float(np.mean(samples)),
            'f_anox_std': float(np.std(samples)),
            'f_anox_median': float(np.median(samples)),
            'f_anox_ci95': (float(np.percentile(samples, 2.5)),
                           float(np.percentile(samples, 97.5))),
            'samples': samples
        }
    
    # ============== 工具方法实现 ==============
    
    def state_dimension(self) -> int:
        """状态变量维度: [M, M×δ_sw]"""
        return 2
    
    def validate_data(self, data: Dict[str, np.ndarray]) -> Tuple[bool, str]:
        """验证铀同位素数据"""
        if 'delta238U' not in data and 'delta_238_u' not in data:
            return False, "Missing uranium isotope data (delta238U or delta_238_u)"
        
        delta_values = data.get('delta238U', data.get('delta_238_u'))
        if np.any(delta_values < -2.0) or np.any(delta_values > 1.0):
            return False, "Delta238U values out of reasonable range (-2.0 to 1.0)"
        
        return True, "Validation passed"
    
    def get_model_info(self) -> Dict:
        """获取模型信息"""
        return {
            'element': self.ELEMENT,
            'name': self.NAME,
            'model_type': self.model_type,
            'reservoir_mass': self.params.reservoir_mass,
            'residence_time': self.params.reservoir_mass / self.fluxes.F_river,
            'fluxes': {
                'F_river': self.fluxes.F_river,
                'F_total': self.fluxes.F_total
            },
            'fractionation': {
                'delta_sw_ox': self.fractionation.delta_sw_ox,
                'delta_sw_anox': self.fractionation.delta_sw_anox,
                'delta_diag': self.fractionation.delta_diag,
                'delta_river': self.fractionation.delta_river
            }
        }
    
    def estimate_anoxic_area(self,
                            f_anox: float,
                            area_modern: float = 0.5) -> float:
        """
        从 f_anox 估算海底缺氧面积
        
        假设缺氧面积与 f_anox 成正比。
        
        Parameters
        ----------
        f_anox : float
            缺氧汇比例
        area_modern : float
            现代缺氧面积比例 (%)
            
        Returns
        -------
        float
            估算的缺氧面积比例 (%)
        """
        # 现代 f_anox ≈ 0.2，缺氧面积 ≈ 0.5%
        f_anox_modern = 0.2
        return area_modern * (f_anox / f_anox_modern)
