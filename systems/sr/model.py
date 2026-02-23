"""
Sr同位素海洋箱模型
基于Wang等(2021)论文实现

核心功能：
1. 质量平衡计算（四端元混合）
2. 随机蒙特卡洛模拟
3. 时间演化模拟
4. 情景分析
"""

import numpy as np
from typing import Dict, Tuple, Optional, List, Callable, Union
from dataclasses import dataclass, field
from scipy import stats
from scipy.interpolate import interp1d

from systems.base.isotope_system import (
    IsotopeSystem, 
    ModelResult, 
    IsotopeParameters,
    RadiogenicSystem
)
from systems.sr.parameters import (
    get_sr_parameters,
    get_permian_parameters,
    calculate_seawater_sr,
    generate_random_parameters,
    STOCHASTIC_RANGES,
    PARAMETER_UNCERTAINTIES,
    PERMIAN_PARAMETERS,
)


@dataclass
class SrFluxConfig:
    """
    Sr通量配置
    四端元输入模型
    """
    # 河流输入
    F_river: float = 47.6e9           # mol/yr
    R_river: float = 0.71107          # ⁸⁷Sr/⁸⁶Sr
    
    # 高温热液输入（洋中脊）
    F_hydrothermal_highT: float = 8.04e9
    R_hydrothermal_highT: float = 0.7037
    
    # 低温热液输入（洋壳蚀变）
    F_hydrothermal_lowT: float = 10e9
    R_hydrothermal_lowT: float = 0.7084
    
    # 成岩作用输入
    F_diagenetic: float = 3.4e9
    R_diagenetic: float = 0.7084
    
    def to_dict(self) -> Dict[str, float]:
        """转换为字典"""
        return {
            'F_river': self.F_river,
            'R_river': self.R_river,
            'F_hydrothermal_highT': self.F_hydrothermal_highT,
            'R_hydrothermal_highT': self.R_hydrothermal_highT,
            'F_hydrothermal_lowT': self.F_hydrothermal_lowT,
            'R_hydrothermal_lowT': self.R_hydrothermal_lowT,
            'F_diagenetic': self.F_diagenetic,
            'R_diagenetic': self.R_diagenetic,
        }


@dataclass
class SrModelResult(ModelResult):
    """
    Sr模型结果数据类
    扩展基类以包含Sr特定的结果
    """
    # 成功运行的参数
    successful_params: List[Dict[str, float]] = field(default_factory=list)
    
    # 参数密度分布
    param_density: Dict[str, np.ndarray] = field(default_factory=dict)
    
    # 统计信息
    statistics: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    def add_successful_run(self, params: Dict[str, float]):
        """添加一次成功的运行"""
        self.successful_params.append(params)
    
    def get_parameter_array(self, param_name: str) -> np.ndarray:
        """获取特定参数的所有成功值"""
        return np.array([p[param_name] for p in self.successful_params])


class StochasticSrModel:
    """
    随机Sr同位素模型
    实现Wang et al. (2021)的蒙特卡洛模拟方法
    """
    
    def __init__(
        self,
        time_points: np.ndarray,
        observed_data: Optional[Tuple[np.ndarray, np.ndarray]] = None,
        confidence_interval: float = 0.95
    ):
        """
        初始化随机模型
        
        Parameters
        ----------
        time_points : array_like
            模拟时间点（Ma）
        observed_data : tuple, optional
            (ages, sr_ratios) 观测数据
            如果为None，则不进行过滤
        confidence_interval : float
            过滤用的置信区间
        """
        self.time_points = time_points
        self.observed_ages = None
        self.observed_ratios = None
        self.confidence_interval = confidence_interval
        
        # 如果提供了观测数据，使用它
        if observed_data is not None:
            self.observed_ages = observed_data[0]
            self.observed_ratios = observed_data[1]
            self._setup_interpolation()
    
    def _setup_interpolation(self):
        """设置观测数据插值"""
        if self.observed_ages is not None and self.observed_ratios is not None:
            self._interp_func = interp1d(
                self.observed_ages, self.observed_ratios, kind='linear', 
                bounds_error=False, fill_value='extrapolate'
            )
    
    def run_single(
        self,
        params: Dict[str, float],
        time_dependent: bool = False
    ) -> Tuple[np.ndarray, bool]:
        """
        运行单次模拟
        
        Parameters
        ----------
        params : dict
            参数字典
        time_dependent : bool
            是否考虑时间依赖的参数变化
            
        Returns
        -------
        (sr_ratios, is_valid) : tuple
            Sr同位素比值数组和是否通过过滤
        """
        sr_ratios = np.zeros(len(self.time_points))
        
        # 默认值（用于固定参数的情景）
        defaults = {
            'F_river': 47.6e9,
            'R_river': 0.71107,
            'F_hydrothermal_highT': 8.04e9,
            'R_hydrothermal_highT': 0.7037,
            'F_hydrothermal_lowT': 10e9,
            'R_hydrothermal_lowT': 0.7084,
            'F_diagenetic': 3.4e9,
            'R_diagenetic': 0.7084,
        }
        
        for i, t in enumerate(self.time_points):
            if time_dependent:
                # 时间依赖的参数（可扩展）
                p = self._get_time_dependent_params(params, t)
            else:
                p = params
            
            # 使用传入参数或默认值
            sr_ratios[i] = calculate_seawater_sr(
                F_river=p.get('F_river', defaults['F_river']),
                R_river=p.get('R_river', defaults['R_river']),
                F_highT=p.get('F_hydrothermal_highT', defaults['F_hydrothermal_highT']),
                R_highT=p.get('R_hydrothermal_highT', defaults['R_hydrothermal_highT']),
                F_lowT=p.get('F_hydrothermal_lowT', defaults['F_hydrothermal_lowT']),
                R_lowT=p.get('R_hydrothermal_lowT', defaults['R_hydrothermal_lowT']),
                F_dia=p.get('F_diagenetic', defaults['F_diagenetic']),
                R_dia=p.get('R_diagenetic', defaults['R_diagenetic'])
            )
        
        # 检查是否通过过滤
        is_valid = self._filter_result(sr_ratios)
        
        return sr_ratios, is_valid
    
    def _get_time_dependent_params(
        self, 
        base_params: Dict[str, float], 
        time: float
    ) -> Dict[str, float]:
        """
        获取时间依赖的参数
        可根据需要扩展为更复杂的函数
        """
        params = base_params.copy()
        
        # 示例：线性变化的热液通量（基于Wang et al.假设）
        # 二叠纪期间热液活动增强
        t_normalized = (299 - time) / (299 - 252)  # 0 to 1
        
        # 高温热液通量可能增加
        if 'F_hydrothermal_highT' in params:
            params['F_hydrothermal_highT'] *= (1 + 0.5 * t_normalized)
        
        return params
    
    def _filter_result(self, sr_ratios: np.ndarray) -> bool:
        """
        过滤模型结果
        检查是否匹配观测数据
        """
        # 如果没有观测数据，不进行过滤
        if self.observed_ratios is None:
            return True
        
        # 插值观测数据到模型时间点
        if hasattr(self, '_interp_func'):
            observed_interp = self._interp_func(self.time_points)
        else:
            observed_interp = np.interp(
                self.time_points, 
                self.observed_ages, 
                self.observed_ratios
            )
        
        # 计算95%置信区间
        diff = np.abs(sr_ratios - observed_interp)
        
        # 简化：检查最大偏差是否在合理范围内
        # 实际应用中可能需要更严格的统计检验
        max_allowed_diff = 0.0005  # 根据Wang et al.的不确定性
        
        return np.all(diff < max_allowed_diff)
    
    def run_monte_carlo(
        self,
        n_runs: int = 5000,
        param_ranges: Optional[Dict[str, Tuple[float, float]]] = None,
        distribution: str = 'uniform',
        time_dependent: bool = False,
        verbose: bool = True
    ) -> SrModelResult:
        """
        运行蒙特卡洛模拟
        
        Parameters
        ----------
        n_runs : int
            运行次数
        param_ranges : dict, optional
            参数范围
        distribution : str
            分布类型
        time_dependent : bool
            是否考虑时间依赖
        verbose : bool
            是否输出进度
            
        Returns
        -------
        SrModelResult
            模拟结果
        """
        if param_ranges is None:
            param_ranges = STOCHASTIC_RANGES
        
        result = SrModelResult(success=True)
        result.time = self.time_points
        
        # 存储所有成功的结果
        successful_ratios = []
        
        for i in range(n_runs):
            if verbose and (i + 1) % 1000 == 0:
                print(f"  Run {i + 1}/{n_runs}...")
            
            # 生成随机参数
            params = generate_random_parameters(param_ranges, distribution)
            
            # 运行单次模拟
            sr_ratios, is_valid = self.run_single(params, time_dependent)
            
            if is_valid:
                successful_ratios.append(sr_ratios)
                result.add_successful_run(params)
        
        if len(successful_ratios) == 0:
            result.success = False
            result.message = "No successful runs found"
            return result
        
        # 计算统计信息
        successful_array = np.array(successful_ratios)
        
        result.values = successful_array.mean(axis=0)
        result.uncertainties = successful_array.std(axis=0)
        
        result.data = {
            'mean_ratio': successful_array.mean(axis=0),
            'std_ratio': successful_array.std(axis=0),
            'percentile_2.5': np.percentile(successful_array, 2.5, axis=0),
            'percentile_97.5': np.percentile(successful_array, 97.5, axis=0),
            'n_successful': len(successful_ratios),
            'success_rate': len(successful_ratios) / n_runs,
        }
        
        # 计算参数统计
        for param_name in param_ranges.keys():
            param_values = result.get_parameter_array(param_name)
            result.statistics[param_name] = {
                'mean': float(np.mean(param_values)),
                'std': float(np.std(param_values)),
                'median': float(np.median(param_values)),
                'min': float(np.min(param_values)),
                'max': float(np.max(param_values)),
            }
        
        return result


class SrIsotopeSystem(RadiogenicSystem):
    """
    Sr同位素体系 - 海洋箱模型
    
    基于Wang et al. (2021)实现：
    - 四端元混合模型
    - 随机蒙特卡洛模拟
    - 多种情景分析
    
    继承自RadiogenicSystem因为Sr包含放射成因的⁸⁷Sr
    （来自⁸⁷Rb衰变，但在海洋时间尺度上可忽略）
    """
    
    ELEMENT = 'sr'
    NAME = 'Strontium'
    ISOTOPES = ['84Sr', '86Sr', '87Sr', '88Sr']
    
    # ⁸⁷Rb衰变常数（用于地质年代学，但海洋模型中通常忽略）
    DECAY_CONSTANT = 1.42e-5  # 1/Ma
    
    def __init__(
        self, 
        parameters: Optional[IsotopeParameters] = None,
        scenario: str = 'modern'
    ):
        """
        初始化Sr同位素体系
        
        Parameters
        ----------
        parameters : IsotopeParameters, optional
            自定义参数
        scenario : str
            'modern' - 现代参数
            'permian' - 二叠纪参数
        """
        if parameters is None:
            if scenario == 'permian':
                parameters = get_permian_parameters()
            else:
                parameters = get_sr_parameters()
        
        self.scenario = scenario
        super().__init__(parameters)
        
        # 初始化随机模型
        self._stochastic_model = None
    
    def _default_parameters(self) -> IsotopeParameters:
        return get_sr_parameters()
    
    # ============== 核心接口实现 ==============
    
    def mass_balance_equation(
        self,
        state: np.ndarray,
        fluxes: Dict[str, float],
        time: Optional[float] = None
    ) -> np.ndarray:
        """
        质量平衡微分方程
        
        dR_sw/dt = [Σ(F_in×(R_in - R_sw))] / M_sw
        
        稳态时 dR_sw/dt = 0
        """
        R_sw = state[0]  # 当前海水Sr同位素比值
        M_sw = self.params.reservoir_mass
        
        # 计算总输入
        total_flux = 0
        weighted_input = 0
        
        for key, flux in fluxes.items():
            if key.startswith('F_'):
                param_base = key[2:]  # 去掉'F_'前缀
                R_key = f'R_{param_base}'
                R_in = fluxes.get(R_key, R_sw)
                
                total_flux += flux
                weighted_input += flux * (R_in - R_sw)
        
        dR_dt = weighted_input / M_sw if M_sw > 0 else 0
        
        return np.array([dR_dt])
    
    def fractionation_factor(
        self,
        process: str,
        temperature: Optional[float] = None,
        **kwargs
    ) -> float:
        """
        获取分馏系数
        
        Sr同位素分馏很小，主要考虑混合
        """
        # Sr同位素质量分馏很小，通常忽略
        return 1.0
    
    def mixing_model(
        self,
        end_member_values: np.ndarray,
        proportions: np.ndarray
    ) -> float:
        """端元混合模型"""
        proportions = np.array(proportions)
        proportions = proportions / np.sum(proportions)
        return float(np.sum(end_member_values * proportions))
    
    def state_dimension(self) -> int:
        """状态维度"""
        return 1
    
    def calculate_model_age(
        self,
        parent_daughter_ratio: float,
        radiogenic_ratio: float,
        initial_ratio: float
    ) -> float:
        """
        计算Rb-Sr模式年龄
        
        t = 1/λ × ln[1 + (R_sample - R_initial) / (⁸⁷Rb/⁸⁶Sr)]
        
        注意：在海洋模型中通常不使用
        """
        return (1 / self.DECAY_CONSTANT) * np.log(
            1 + (radiogenic_ratio - initial_ratio) / parent_daughter_ratio
        )
    
    def epsilon_value(self, ratio: float) -> float:
        """
        计算ε值
        
        对于Sr通常不使用ε表示法，这里提供接口
        """
        # Sr通常不使用epsilon表示法
        return 0.0
    
    # ============== Sr专用方法 ==============
    
    def calculate_seawater_ratio(
        self,
        config: Optional[SrFluxConfig] = None,
        **kwargs
    ) -> float:
        """
        计算海水⁸⁷Sr/⁸⁶Sr值
        
        Parameters
        ----------
        config : SrFluxConfig, optional
            通量配置
        **kwargs : dict
            可直接传入参数
            
        Returns
        -------
        float
            海水⁸⁷Sr/⁸⁶Sr比值
        """
        if config is not None:
            params = config.to_dict()
        else:
            params = kwargs
        
        return calculate_seawater_sr(
            F_river=params.get('F_river', 47.6e9),
            R_river=params.get('R_river', 0.71107),
            F_highT=params.get('F_hydrothermal_highT', 8.04e9),
            R_highT=params.get('R_hydrothermal_highT', 0.7037),
            F_lowT=params.get('F_hydrothermal_lowT', 10e9),
            R_lowT=params.get('R_hydrothermal_lowT', 0.7084),
            F_dia=params.get('F_diagenetic', 3.4e9),
            R_dia=params.get('R_diagenetic', 0.7084)
        )
    
    def monte_carlo_simulation(
        self,
        time_span: Tuple[float, float] = (299, 252),
        n_runs: int = 5000,
        n_time_points: int = 50,
        param_ranges: Optional[Dict[str, Tuple[float, float]]] = None,
        observed_data: Optional[Tuple[np.ndarray, np.ndarray]] = None,
        filter_by_data: bool = True,
        verbose: bool = True
    ) -> SrModelResult:
        """
        运行蒙特卡洛模拟
        
        实现Wang et al. (2021)的随机模拟方法
        
        Parameters
        ----------
        time_span : tuple
            时间范围 (Ma)，默认二叠纪 (299-252 Ma)
        n_runs : int
            蒙特卡洛运行次数
        n_time_points : int
            时间点数量
        param_ranges : dict, optional
            参数范围，默认使用论文中的范围
        observed_data : tuple, optional
            (ages, sr_ratios) 观测数据用于过滤
        filter_by_data : bool
            是否根据观测数据过滤结果
        verbose : bool
            是否输出进度
            
        Returns
        -------
        SrModelResult
            包含成功运行的参数和统计信息
            
        Examples
        --------
        >>> sr = SrIsotopeSystem(scenario='permian')
        >>> result = sr.monte_carlo_simulation(n_runs=5000)
        >>> print(f"成功率: {result.data['success_rate']:.2%}")
        >>> print(f"平均河流通量: {result.statistics['F_river']['mean']:.2e} mol/yr")
        """
        if verbose:
            print("\n" + "="*60)
            print("Sr同位素随机海洋箱模型模拟")
            print(f"基于 Wang et al. (2021) Earth-Science Reviews")
            print("="*60)
            print(f"模拟时间范围: {time_span[0]} - {time_span[1]} Ma")
            print(f"蒙特卡洛运行次数: {n_runs}")
            print(f"时间点数量: {n_time_points}")
        
        # 生成时间点
        time_points = np.linspace(time_span[0], time_span[1], n_time_points)
        
        # 创建随机模型
        stochastic = StochasticSrModel(
            time_points=time_points,
            observed_data=observed_data if filter_by_data else None
        )
        
        # 运行模拟
        result = stochastic.run_monte_carlo(
            n_runs=n_runs,
            param_ranges=param_ranges,
            verbose=verbose
        )
        
        if verbose and result.success:
            print(f"\n模拟完成!")
            print(f"成功率: {result.data['success_rate']:.2%}")
            print(f"成功运行次数: {result.data['n_successful']}/{n_runs}")
        
        return result
    
    def scenario_analysis(
        self,
        scenario_name: str,
        time_span: Tuple[float, float] = (299, 252),
        n_runs: int = 5000,
        **kwargs
    ) -> SrModelResult:
        """
        情景分析
        
        测试Wang et al. (2021)中的不同约束情景
        
        Parameters
        ----------
        scenario_name : str
            情景名称:
            - 'scenario1': 所有参数可变（除成岩）
            - 'scenario2': 固定河流同位素
            - 'scenario3': 固定河流通量
            - 'scenario4': 固定高温热液通量
            - 'scenario5': 固定低温热液通量
            - 'scenario6': 固定所有热液通量
            - 'scenario7': 固定河流参数
        time_span : tuple
            时间范围
        n_runs : int
            运行次数
        **kwargs : dict
            其他参数
            
        Returns
        -------
        SrModelResult
            情景分析结果
        """
        # 定义各情景的固定参数
        fixed_params = {
            'scenario1': ['F_diagenetic', 'R_diagenetic'],
            'scenario2': ['R_river', 'F_diagenetic', 'R_diagenetic'],
            'scenario3': ['F_river', 'F_diagenetic', 'R_diagenetic'],
            'scenario4': ['F_hydrothermal_highT', 'F_diagenetic', 'R_diagenetic'],
            'scenario5': ['F_hydrothermal_lowT', 'F_diagenetic', 'R_diagenetic'],
            'scenario6': ['F_hydrothermal_highT', 'F_hydrothermal_lowT', 
                         'F_diagenetic', 'R_diagenetic'],
            'scenario7': ['F_river', 'R_river', 'F_diagenetic', 'R_diagenetic'],
        }
        
        if scenario_name not in fixed_params:
            raise ValueError(
                f"未知情景: {scenario_name}\n"
                f"支持的情景: {list(fixed_params.keys())}"
            )
        
        # 获取参数范围并移除固定参数
        param_ranges = STOCHASTIC_RANGES.copy()
        for param in fixed_params[scenario_name]:
            if param in param_ranges:
                del param_ranges[param]
        
        print(f"\n运行情景: {scenario_name}")
        print(f"可变参数: {list(param_ranges.keys())}")
        
        return self.monte_carlo_simulation(
            time_span=time_span,
            n_runs=n_runs,
            param_ranges=param_ranges,
            **kwargs
        )
    
    def sensitivity_analysis(
        self,
        param_name: str,
        param_range: Tuple[float, float],
        n_steps: int = 50,
        fixed_params: Optional[Dict[str, float]] = None
    ) -> ModelResult:
        """
        参数敏感性分析
        
        Parameters
        ----------
        param_name : str
            要分析的参数名
        param_range : tuple
            参数范围
        n_steps : int
            步数
        fixed_params : dict, optional
            其他固定参数
            
        Returns
        -------
        ModelResult
            敏感性分析结果
        """
        param_values = np.linspace(param_range[0], param_range[1], n_steps)
        sr_ratios = np.zeros(n_steps)
        
        if fixed_params is None:
            fixed_params = {}
        
        # 参数名映射 (CLI名称 -> 函数参数名)
        param_mapping = {
            'F_river': 'F_river',
            'R_river': 'R_river',
            'F_highT': 'F_hydrothermal_highT',
            'R_highT': 'R_hydrothermal_highT',
            'F_lowT': 'F_hydrothermal_lowT',
            'R_lowT': 'R_hydrothermal_lowT',
            'F_dia': 'F_diagenetic',
            'R_dia': 'R_diagenetic',
        }
        
        func_param_name = param_mapping.get(param_name, param_name)
        
        for i, val in enumerate(param_values):
            params = fixed_params.copy()
            params[func_param_name] = val
            
            sr_ratios[i] = self.calculate_seawater_ratio(**params)
        
        result = ModelResult(success=True)
        result.data = {
            'param_name': param_name,
            'param_values': param_values,
            'sr_ratios': sr_ratios,
            'sensitivity': np.gradient(sr_ratios, param_values),
        }
        
        return result
    
    def validate_data(self, data: Dict[str, np.ndarray]) -> Tuple[bool, str]:
        """
        验证Sr同位素数据
        
        Parameters
        ----------
        data : dict
            输入数据字典
            
        Returns
        -------
        (is_valid, message) : tuple
        """
        if 'sr_ratio' not in data and '87Sr/86Sr' not in data:
            return False, "Missing required field: sr_ratio or 87Sr/86Sr"
        
        sr_values = data.get('sr_ratio', data.get('87Sr/86Sr'))
        
        # Sr同位素合理范围（显生宙）
        if np.any(sr_values < 0.706) or np.any(sr_values > 0.730):
            return False, "Sr values out of reasonable range (0.706 to 0.730)"
        
        return True, "Validation passed"
