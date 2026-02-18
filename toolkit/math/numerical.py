"""
核心数学工具模块
提供数值计算、ODE求解、插值、优化等通用数学功能
与具体同位素体系无关
"""

import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.interpolate import interp1d, CubicSpline
from scipy.optimize import minimize, curve_fit
from typing import Callable, Tuple, Optional, Union, List
from dataclasses import dataclass


@dataclass
class ODEResult:
    """ODE求解结果容器"""
    t: np.ndarray
    y: np.ndarray
    success: bool
    message: str = ""


class ODESolver:
    """
    常微分方程求解器
    支持多种求解方法，适用于各类同位素演化模型
    """
    
    METHODS = ['odeint', 'RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA']
    
    @staticmethod
    def solve(
        func: Callable,
        y0: Union[float, np.ndarray],
        t_span: Tuple[float, float],
        args: tuple = (),
        method: str = 'rk45',
        n_points: int = 1000,
        **kwargs
    ) -> ODEResult:
        """
        求解ODE
        
        Parameters
        ----------
        func : Callable
            dy/dt = func(t, y, *args)
        y0 : array_like
            初始条件
        t_span : tuple
            (t_start, t_end)
        args : tuple
            传递给func的额外参数
        method : str
            求解方法 ('odeint', 'RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA')
        n_points : int
            输出点数
            
        Returns
        -------
        ODEResult
        """
        try:
            if method == 'odeint':
                # 使用odeint（旧API，但某些情况更稳定）
                t = np.linspace(t_span[0], t_span[1], n_points)
                
                # 调整函数签名以匹配odeint要求（y, t）而非（t, y）
                # odeint要求func(y, t, ...)，而solve_ivp要求func(t, y, ...)
                def func_odeint(y, t):
                    return func(t, y, *args)
                
                solution = odeint(func_odeint, y0, t, **kwargs)
                return ODEResult(t=t, y=solution.flatten() if solution.ndim > 1 else solution, success=True)
            else:
                # 使用solve_ivp（新API，支持更多方法）
                t_eval = np.linspace(t_span[0], t_span[1], n_points)
                
                def func_ivp(t, y):
                    return func(t, y, *args)
                
                # 转换method名称以匹配scipy要求
                method_map = {
                    'rk45': 'RK45', 'rk23': 'RK23', 'dop853': 'DOP853',
                    'radau': 'Radau', 'bdf': 'BDF', 'lsoda': 'LSODA'
                }
                method_normalized = method_map.get(method.lower(), method)
                
                sol = solve_ivp(
                    func_ivp, t_span, [y0] if np.isscalar(y0) else y0,
                    method=method_normalized, t_eval=t_eval, **kwargs
                )
                
                return ODEResult(
                    t=sol.t, 
                    y=sol.y.T if sol.y.ndim > 1 else sol.y,
                    success=sol.success,
                    message=sol.message if hasattr(sol, 'message') else ""
                )
        except Exception as e:
            return ODEResult(
                t=np.array([]), 
                y=np.array([]), 
                success=False, 
                message=str(e)
            )


class Interpolator:
    """
    插值工具类
    提供多种插值方法
    """
    
    METHODS = ['linear', 'cubic', 'spline', 'nearest']
    
    @staticmethod
    def interpolate(
        x: np.ndarray,
        y: np.ndarray,
        x_new: np.ndarray,
        method: str = 'linear',
        fill_value: str = 'extrapolate',
        **kwargs
    ) -> np.ndarray:
        """
        一维插值
        
        Parameters
        ----------
        x, y : array_like
            原始数据点
        x_new : array_like
            需要插值的点
        method : str
            插值方法
        fill_value : str or float
            外推处理方式
            
        Returns
        -------
        array_like
            插值结果
        """
        if method == 'cubic':
            cs = CubicSpline(x, y, **kwargs)
            return cs(x_new)
        else:
            kind = 'linear' if method == 'spline' else method
            f = interp1d(x, y, kind=kind, fill_value=fill_value, bounds_error=False)
            return f(x_new)
    
    @staticmethod
    def resample(
        data: np.ndarray,
        original_scale: np.ndarray,
        new_scale: np.ndarray,
        axis: int = 0
    ) -> np.ndarray:
        """
        按新的尺度重采样数据
        
        Parameters
        ----------
        data : array_like
            原始数据
        original_scale : array_like
            原始尺度（如年龄、深度）
        new_scale : array_like
            新尺度
        axis : int
            重采样轴
            
        Returns
        -------
        array_like
            重采样后的数据
        """
        f = interp1d(original_scale, data, axis=axis, 
                     fill_value='extrapolate', bounds_error=False)
        return f(new_scale)


class Optimizer:
    """
    优化工具类
    用于参数反演、拟合等
    """
    
    @staticmethod
    def curve_fit_wrapper(
        f: Callable,
        xdata: np.ndarray,
        ydata: np.ndarray,
        p0: Optional[np.ndarray] = None,
        bounds: Optional[Tuple] = None,
        **kwargs
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        曲线拟合包装器
        
        Returns
        -------
        popt : array
            最优参数
        pcov : 2D array
            参数协方差矩阵
        """
        try:
            popt, pcov = curve_fit(f, xdata, ydata, p0=p0, bounds=bounds, **kwargs)
            return popt, pcov
        except Exception as e:
            raise RuntimeError(f"Curve fitting failed: {e}")
    
    @staticmethod
    def minimize_wrapper(
        func: Callable,
        x0: np.ndarray,
        method: str = 'L-BFGS-B',
        bounds: Optional[List[Tuple]] = None,
        **kwargs
    ):
        """
        最小化优化包装器
        """
        result = minimize(func, x0, method=method, bounds=bounds, **kwargs)
        return result


class StatisticalTools:
    """
    统计工具类
    """
    
    @staticmethod
    def weighted_average(values: np.ndarray, errors: np.ndarray) -> Tuple[float, float]:
        """
        计算加权平均值和误差
        
        Parameters
        ----------
        values : array_like
            测量值
        errors : array_like
            误差（1σ）
            
        Returns
        -------
        mean, error : float
            加权平均值和误差
        """
        weights = 1.0 / errors**2
        mean = np.sum(values * weights) / np.sum(weights)
        error = np.sqrt(1.0 / np.sum(weights))
        return mean, error
    
    @staticmethod
    def monte_carlo_uncertainty(
        model_func: Callable,
        param_means: np.ndarray,
        param_stds: np.ndarray,
        x_data: np.ndarray,
        n_samples: int = 10000
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        使用蒙特卡洛方法计算模型输出的不确定性
        
        Parameters
        ----------
        model_func : Callable
            模型函数 f(params, x) -> y
        param_means : array_like
            参数均值
        param_stds : array_like
            参数标准差
        x_data : array_like
            输入数据点
        n_samples : int
            采样次数
            
        Returns
        -------
        mean, std : array_like
            预测的均值和标准差
        """
        n_params = len(param_means)
        results = np.zeros((n_samples, len(x_data)))
        
        for i in range(n_samples):
            # 从正态分布中采样参数
            params_sampled = np.random.normal(param_means, param_stds)
            results[i] = model_func(params_sampled, x_data)
        
        return np.mean(results, axis=0), np.std(results, axis=0)


# 常用数学函数
class MathUtils:
    """通用数学工具"""
    
    @staticmethod
    def moving_average(data: np.ndarray, window: int) -> np.ndarray:
        """移动平均"""
        return np.convolve(data, np.ones(window)/window, mode='same')
    
    @staticmethod
    def exponential_decay(t: np.ndarray, tau: float, A: float = 1.0) -> np.ndarray:
        """指数衰减 A * exp(-t/tau)"""
        return A * np.exp(-t / tau)
    
    @staticmethod
    def box_model_residence_time(M: float, F_in: float, F_out: float) -> float:
        """
        计算箱模型停留时间
        
        Parameters
        ----------
        M : float
            库储量
        F_in, F_out : float
            输入/输出通量
            
        Returns
        -------
        float
            停留时间
        """
        if F_out == 0:
            return np.inf
        return M / F_out
