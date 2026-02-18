"""
同位素体系基类模块
定义所有同位素体系必须实现的接口
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Optional, Callable, Any
from dataclasses import dataclass, field
import numpy as np
from pathlib import Path


@dataclass
class IsotopeParameters:
    """
    同位素参数数据类
    各体系继承并扩展此类
    """
    # 基本标识
    element: str = ""                          # 元素符号
    name: str = ""                             # 元素名称
    
    # 标准物质
    reference_standard: str = ""               # 参考标准名称
    reference_ratios: Dict[str, float] = field(default_factory=dict)  # 标准比值
    
    # 分馏系数
    fractionation_factors: Dict[str, float] = field(default_factory=dict)
    
    # 端元值（用于混合模型）
    end_members: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    # 储库参数（用于箱模型）
    reservoir_mass: float = 0.0                # 储库质量 (mol)
    input_fluxes: Dict[str, float] = field(default_factory=dict)
    output_fluxes: Dict[str, float] = field(default_factory=dict)
    
    def to_dict(self) -> Dict:
        """转换为字典"""
        return {
            'element': self.element,
            'name': self.name,
            'reference_standard': self.reference_standard,
            'reference_ratios': self.reference_ratios,
            'fractionation_factors': self.fractionation_factors,
            'end_members': self.end_members,
            'reservoir_mass': self.reservoir_mass,
            'input_fluxes': self.input_fluxes,
            'output_fluxes': self.output_fluxes
        }


@dataclass
class ModelResult:
    """
    模型结果数据类
    """
    success: bool = True
    message: str = ""
    data: Dict[str, Any] = field(default_factory=dict)
    
    # 时间序列数据
    time: Optional[np.ndarray] = None
    values: Optional[np.ndarray] = None
    uncertainties: Optional[np.ndarray] = None
    
    # 参数估计结果
    fitted_parameters: Optional[Dict[str, float]] = None
    parameter_errors: Optional[Dict[str, float]] = None
    
    def get(self, key: str, default=None):
        """获取数据字段"""
        return self.data.get(key, default)
    
    def add(self, key: str, value: Any):
        """添加数据字段"""
        self.data[key] = value


class IsotopeSystem(ABC):
    """
    同位素体系抽象基类
    所有具体同位素体系必须继承此类
    
    设计原则：
    1. 与数据I/O解耦 - 只负责计算逻辑
    2. 与可视化解耦 - 返回数据而非绘图
    3. 可配置化 - 通过parameters配置而非硬编码
    """
    
    # 类属性：体系标识
    ELEMENT: str = ""           # 元素符号（必须小写）
    NAME: str = ""              # 元素名称
    ISOTOPES: List[str] = []    # 涉及的同位素列表
    
    def __init__(self, parameters: Optional[IsotopeParameters] = None):
        """
        初始化同位素体系
        
        Parameters
        ----------
        parameters : IsotopeParameters, optional
            体系参数，如果为None则使用默认参数
        """
        self.params = parameters or self._default_parameters()
        self._validate_parameters()
        
        # 内部状态
        self._cache: Dict[str, Any] = {}  # 计算缓存
    
    @abstractmethod
    def _default_parameters(self) -> IsotopeParameters:
        """
        返回默认参数
        子类必须实现此方法
        """
        pass
    
    def _validate_parameters(self):
        """
        验证参数有效性
        子类可重写此方法添加特定验证
        """
        if not self.params.element:
            self.params.element = self.ELEMENT
        if not self.params.name:
            self.params.name = self.NAME
    
    # ============== 核心计算接口（必须实现） ==============
    
    @abstractmethod
    def mass_balance_equation(self, 
                             state: np.ndarray,
                             fluxes: Dict[str, float],
                             time: Optional[float] = None) -> np.ndarray:
        """
        质量平衡微分方程
        
        dC/dt = f(state, fluxes, t)
        
        Parameters
        ----------
        state : array_like
            当前状态变量（浓度、同位素值等）
        fluxes : dict
            通量字典
        time : float, optional
            当前时间
            
        Returns
        -------
        array_like
            d(state)/dt
        """
        pass
    
    @abstractmethod
    def fractionation_factor(self,
                            process: str,
                            temperature: Optional[float] = None,
                            **kwargs) -> float:
        """
        获取分馏系数
        
        Parameters
        ----------
        process : str
            过程名称（如"equilibrium", "kinetic"）
        temperature : float, optional
            温度（K），用于温度依赖的分馏
        **kwargs
            其他参数
            
        Returns
        -------
        float
            分馏系数α
        """
        pass
    
    @abstractmethod
    def mixing_model(self,
                    end_member_values: np.ndarray,
                    proportions: np.ndarray) -> float:
        """
        端元混合模型
        
        Parameters
        ----------
        end_member_values : array_like
            各端元的同位素值
        proportions : array_like
            各端元的比例（和为1）
            
        Returns
        -------
        float
            混合后的同位素值
        """
        pass
    
    # ============== 高级计算接口（可选实现） ==============
    
    def solve_steady_state(self,
                          fluxes: Dict[str, float],
                          initial_guess: Optional[np.ndarray] = None) -> ModelResult:
        """
        求解稳态（dC/dt = 0）
        
        Parameters
        ----------
        fluxes : dict
            通量条件
        initial_guess : array_like, optional
            初始猜测值
            
        Returns
        -------
        ModelResult
            稳态解
        """
        # 默认实现：使用根求解
        from scipy.optimize import fsolve
        
        def steady_eq(state):
            return self.mass_balance_equation(state, fluxes)
        
        if initial_guess is None:
            initial_guess = np.zeros(self.state_dimension())
        
        try:
            solution = fsolve(steady_eq, initial_guess)
            return ModelResult(
                success=True,
                data={'steady_state': solution}
            )
        except Exception as e:
            return ModelResult(
                success=False,
                message=f"Steady state solve failed: {e}"
            )
    
    def time_evolution(self,
                      initial_state: np.ndarray,
                      time_span: Tuple[float, float],
                      fluxes: Dict[str, float],
                      n_points: int = 1000) -> ModelResult:
        """
        计算时间演化
        
        Parameters
        ----------
        initial_state : array_like
            初始状态
        time_span : tuple
            (t_start, t_end)
        fluxes : dict
            通量条件（可以是常数或函数）
        n_points : int
            时间点数量
            
        Returns
        -------
        ModelResult
            包含时间序列的结果
        """
        from toolkit.math.numerical import ODESolver
        
        def ode_func(t, y):
            return self.mass_balance_equation(y, fluxes, t)
        
        result = ODESolver.solve(
            ode_func, initial_state, time_span, n_points=n_points
        )
        
        if result.success:
            return ModelResult(
                success=True,
                time=result.t,
                values=result.y
            )
        else:
            return ModelResult(
                success=False,
                message=result.message
            )
    
    def inverse_model(self,
                     observations: np.ndarray,
                     observation_times: np.ndarray,
                     parameters_to_fit: List[str],
                     initial_guess: Dict[str, float]) -> ModelResult:
        """
        反演模型（从观测数据推断参数）
        
        Parameters
        ----------
        observations : array_like
            观测数据
        observation_times : array_like
            观测时间点
        parameters_to_fit : list
            需要拟合的参数名称列表
        initial_guess : dict
            参数初始猜测值
            
        Returns
        -------
        ModelResult
            拟合结果
        """
        # 默认实现：使用最小二乘拟合
        # 子类应该根据具体模型重写此方法
        raise NotImplementedError("Inverse model not implemented for this system")
    
    def uncertainty_propagation(self,
                               model_result: ModelResult,
                               parameter_errors: Dict[str, float],
                               n_monte_carlo: int = 10000) -> ModelResult:
        """
        误差传播计算
        
        Parameters
        ----------
        model_result : ModelResult
            模型结果
        parameter_errors : dict
            参数误差
        n_monte_carlo : int
            蒙特卡洛采样次数
            
        Returns
        -------
        ModelResult
            包含不确定性的结果
        """
        # 默认实现：简单的高斯误差传播
        # 子类可根据需要重写
        return model_result
    
    # ============== 工具方法 ==============
    
    @abstractmethod
    def state_dimension(self) -> int:
        """返回状态变量的维度"""
        pass
    
    def clear_cache(self):
        """清除计算缓存"""
        self._cache.clear()
    
    def get_info(self) -> Dict:
        """获取体系信息"""
        return {
            'element': self.ELEMENT,
            'name': self.NAME,
            'isotopes': self.ISOTOPES,
            'parameters': self.params.to_dict()
        }
    
    def validate_data(self, data: Dict[str, np.ndarray]) -> Tuple[bool, str]:
        """
        验证输入数据
        
        Parameters
        ----------
        data : dict
            输入数据字典
            
        Returns
        -------
        (is_valid, message) : tuple
        """
        # 基础验证：检查必要字段是否存在
        return True, "Data validation passed"


class MultiIsotopeSystem(IsotopeSystem):
    """
    多同位素体系基类（如Mg-25/Mg-24/Mg-26, S-33/S-34/S-36等）
    """
    
    @abstractmethod
    def mass_independent_fractionation(self,
                                      deltas: Dict[str, float]) -> Dict[str, float]:
        """
        计算质量无关分馏（如Δ³³S, Δ³⁶S）
        
        Parameters
        ----------
        deltas : dict
            各同位素的delta值
            
        Returns
        -------
        dict
            质量无关分馏指标
        """
        pass
    
    @abstractmethod
    def triple_isotope_plot(self,
                           data: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        准备三同位素图解数据
        
        Parameters
        ----------
        data : dict
            同位素数据
            
        Returns
        -------
        dict
            绘图所需数据
        """
        pass


class RadiogenicSystem(IsotopeSystem):
    """
    放射成因同位素体系基类（如Sr, Nd, Os等）
    """
    
    DECAY_CONSTANT: float = 0.0  # 衰变常数 (1/Ma)
    
    @abstractmethod
    def calculate_model_age(self,
                           parent_daughter_ratio: float,
                           radiogenic_ratio: float,
                           initial_ratio: float) -> float:
        """
        计算模式年龄
        
        t = 1/λ × ln[1 + (R_sample - R_initial) / (P/D)]
        
        Parameters
        ----------
        parent_daughter_ratio : float
            母/子同位素比值
        radiogenic_ratio : float
            放射成因同位素比值
        initial_ratio : float
            初始比值
            
        Returns
        -------
        float
            模式年龄（Ma）
        """
        pass
    
    @abstractmethod
    def epsilon_value(self, ratio: float) -> float:
        """
        计算ε值
        
        Parameters
        ----------
        ratio : float
            同位素比值
            
        Returns
        -------
        float
            ε值
        """
        pass
