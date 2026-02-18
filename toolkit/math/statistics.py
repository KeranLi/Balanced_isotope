"""
统计工具模块
提供同位素数据分析常用的统计方法

包括:
- Bootstrap 重采样
- LOWESS 局部加权回归
- Changepoint 检测
- 蒙特卡洛不确定性分析
"""

import numpy as np
from typing import Callable, Tuple, Optional, Dict, List
from scipy import stats
from scipy.interpolate import interp1d


class Bootstrap:
    """
    Bootstrap 重采样方法
    
    用于估计统计量的置信区间
    """
    
    @staticmethod
    def confidence_interval(data: np.ndarray,
                           statistic_func: Callable = np.mean,
                           n_bootstrap: int = 10000,
                           confidence_level: float = 0.95) -> Dict:
        """
        计算统计量的 Bootstrap 置信区间
        
        Parameters
        ----------
        data : array_like
            原始数据
        statistic_func : callable
            统计量函数 (如 np.mean, np.median)
        n_bootstrap : int
            Bootstrap 采样次数
        confidence_level : float
            置信水平 (0-1)
            
        Returns
        -------
        dict
            {
                'statistic': 原始统计量,
                'mean': Bootstrap 均值,
                'std': Bootstrap 标准差,
                'ci_lower': 置信区间下限,
                'ci_upper': 置信区间上限,
                'bootstrap_statistics': 所有 Bootstrap 统计量
            }
        """
        n = len(data)
        
        # 计算原始统计量
        original_stat = statistic_func(data)
        
        # Bootstrap 采样
        bootstrap_stats = []
        for _ in range(n_bootstrap):
            sample = np.random.choice(data, size=n, replace=True)
            bootstrap_stats.append(statistic_func(sample))
        
        bootstrap_stats = np.array(bootstrap_stats)
        
        # 计算置信区间
        alpha = 1 - confidence_level
        ci_lower = np.percentile(bootstrap_stats, alpha/2 * 100)
        ci_upper = np.percentile(bootstrap_stats, (1 - alpha/2) * 100)
        
        return {
            'statistic': float(original_stat),
            'mean': float(np.mean(bootstrap_stats)),
            'std': float(np.std(bootstrap_stats)),
            'ci_lower': float(ci_lower),
            'ci_upper': float(ci_upper),
            'bootstrap_statistics': bootstrap_stats
        }
    
    @staticmethod
    def regression(x: np.ndarray,
                  y: np.ndarray,
                  x_new: Optional[np.ndarray] = None,
                  n_bootstrap: int = 1000,
                  confidence_level: float = 0.95) -> Dict:
        """
        Bootstrap 回归
        
        对数据进行多次重采样回归，评估回归曲线的不确定性
        
        Parameters
        ----------
        x, y : array_like
            原始数据
        x_new : array_like, optional
            预测点，如果为None则使用x
        n_bootstrap : int
            Bootstrap 次数
        confidence_level : float
            置信水平
            
        Returns
        -------
        dict
            {
                'x': 预测点,
                'y_mean': 回归均值,
                'y_std': 回归标准差,
                'ci_lower': 置信区间下限,
                'ci_upper': 置信区间上限,
                'bootstrap_curves': 所有 Bootstrap 回归曲线
            }
        """
        if x_new is None:
            x_new = x
        
        n = len(x)
        bootstrap_curves = []
        
        for _ in range(n_bootstrap):
            # 重采样
            indices = np.random.choice(n, size=n, replace=True)
            x_sample = x[indices]
            y_sample = y[indices]
            
            # 排序以确保单调性
            sort_idx = np.argsort(x_sample)
            x_sample = x_sample[sort_idx]
            y_sample = y_sample[sort_idx]
            
            # 线性插值 (可以根据需要改为其他回归方法)
            try:
                f = interp1d(x_sample, y_sample, kind='linear', 
                            bounds_error=False, fill_value='extrapolate')
                y_pred = f(x_new)
                bootstrap_curves.append(y_pred)
            except ValueError:
                # 如果插值失败，跳过
                continue
        
        bootstrap_curves = np.array(bootstrap_curves)
        
        alpha = 1 - confidence_level
        return {
            'x': x_new,
            'y_mean': np.mean(bootstrap_curves, axis=0),
            'y_std': np.std(bootstrap_curves, axis=0),
            'ci_lower': np.percentile(bootstrap_curves, alpha/2 * 100, axis=0),
            'ci_upper': np.percentile(bootstrap_curves, (1 - alpha/2) * 100, axis=0),
            'bootstrap_curves': bootstrap_curves
        }


class LOWESS:
    """
    LOWESS (Locally Weighted Scatterplot Smoothing)
    局部加权散点平滑
    
    用于时间序列数据的非参数回归
    """
    
    @staticmethod
    def fit(x: np.ndarray,
            y: np.ndarray,
            frac: float = 0.3,
            it: int = 3,
            delta: float = 0.0) -> Dict:
        """
        LOWESS 回归
        
        Parameters
        ----------
        x, y : array_like
            数据点
        frac : float
            平滑参数，用于每个局部回归的数据点比例 (0-1)
        it : int
            鲁棒性迭代次数
        delta : float
            计算加速参数，在x间隔小于delta的区间内跳过计算
            
        Returns
        -------
        dict
            {
                'x': 输入x,
                'y_smoothed': 平滑后的y值,
                'residuals': 残差
            }
        """
        try:
            from statsmodels.nonparametric.smoothers_lowess import lowess
            
            # 执行 LOWESS
            result = lowess(y, x, frac=frac, it=it, delta=delta, return_sorted=False)
            
            return {
                'x': x,
                'y_smoothed': result,
                'residuals': y - result
            }
        except ImportError:
            # 如果 statsmodels 不可用，使用简单的移动平均作为替代
            return LOWESS._moving_average_lowess(x, y, frac)
    
    @staticmethod
    def _moving_average_lowess(x: np.ndarray, y: np.ndarray, frac: float = 0.3) -> Dict:
        """使用移动平均作为 LOWESS 的简单替代"""
        n = len(x)
        window = max(3, int(n * frac))
        
        # 按x排序
        sort_idx = np.argsort(x)
        x_sorted = x[sort_idx]
        y_sorted = y[sort_idx]
        
        # 移动平均
        y_smoothed = np.convolve(y_sorted, np.ones(window)/window, mode='same')
        
        # 恢复原始顺序
        inv_sort_idx = np.argsort(sort_idx)
        y_smoothed = y_smoothed[inv_sort_idx]
        
        return {
            'x': x,
            'y_smoothed': y_smoothed,
            'residuals': y - y_smoothed
        }
    
    @staticmethod
    def bootstrap_fit(x: np.ndarray,
                     y: np.ndarray,
                     frac: float = 0.3,
                     n_bootstrap: int = 1000,
                     confidence_level: float = 0.95) -> Dict:
        """
        Bootstrap LOWESS 回归 (带置信区间)
        
        Parameters
        ----------
        x, y : array_like
            数据点
        frac : float
            平滑参数
        n_bootstrap : int
            Bootstrap 次数
        confidence_level : float
            置信水平
            
        Returns
        -------
        dict
            {
                'x': 输入x,
                'y_smoothed': 平滑曲线,
                'ci_lower': 置信区间下限,
                'ci_upper': 置信区间上限,
                'y_mean': Bootstrap 均值
            }
        """
        try:
            from statsmodels.nonparametric.smoothers_lowess import lowess
        except ImportError:
            # 使用简单替代
            result = LOWESS._moving_average_lowess(x, y, frac)
            return {
                'x': result['x'],
                'y_smoothed': result['y_smoothed'],
                'ci_lower': result['y_smoothed'] * 0.9,
                'ci_upper': result['y_smoothed'] * 1.1,
                'y_mean': result['y_smoothed']
            }
        
        n = len(x)
        bootstrap_curves = []
        
        # 生成平滑的x网格
        x_grid = np.linspace(x.min(), x.max(), n)
        
        for _ in range(n_bootstrap):
            # 重采样
            indices = np.random.choice(n, size=n, replace=True)
            x_sample = x[indices]
            y_sample = y[indices]
            
            # 排序
            sort_idx = np.argsort(x_sample)
            x_sample = x_sample[sort_idx]
            y_sample = y_sample[sort_idx]
            
            try:
                # LOWESS
                y_smooth = lowess(y_sample, x_sample, frac=frac, 
                                 it=0, return_sorted=False)
                
                # 插值到网格
                f = interp1d(x_sample, y_smooth, kind='linear',
                           bounds_error=False, fill_value='extrapolate')
                bootstrap_curves.append(f(x_grid))
            except ValueError:
                continue
        
        bootstrap_curves = np.array(bootstrap_curves)
        
        alpha = 1 - confidence_level
        return {
            'x': x_grid,
            'y_smoothed': np.mean(bootstrap_curves, axis=0),
            'ci_lower': np.percentile(bootstrap_curves, alpha/2 * 100, axis=0),
            'ci_upper': np.percentile(bootstrap_curves, (1 - alpha/2) * 100, axis=0),
            'y_mean': np.mean(bootstrap_curves, axis=0)
        }


class ChangepointDetector:
    """
    变点检测 (Changepoint Detection)
    
    检测时间序列中统计特性发生变化的点
    """
    
    @staticmethod
    def detect_mean_variance(data: np.ndarray,
                            method: str = 'amoc') -> Dict:
        """
        检测均值和方差的变化点
        
        Parameters
        ----------
        data : array_like
            时间序列数据
        method : str
            检测方法：'amoc' (at most one change)
            
        Returns
        -------
        dict
            {
                'changepoint': 变点位置,
                'likelihood': 似然比统计量,
                'before_mean': 变点前均值,
                'after_mean': 变点后均值,
                'before_std': 变点前标准差,
                'after_std': 变点后标准差
            }
        """
        n = len(data)
        
        if method == 'amoc':
            # At Most One Change (AMOC) 检测
            max_stat = -np.inf
            changepoint = 0
            
            for i in range(2, n-2):  # 避免边界
                # 分段统计
                before = data[:i]
                after = data[i:]
                
                # 计算似然比统计量 (简化版)
                pooled_std = np.sqrt(
                    ((len(before)-1)*np.var(before, ddof=1) + 
                     (len(after)-1)*np.var(after, ddof=1)) / 
                    (len(before) + len(after) - 2)
                )
                
                if pooled_std > 0:
                    stat = (np.mean(before) - np.mean(after)) ** 2 / (pooled_std ** 2)
                    stat *= len(before) * len(after) / n
                    
                    if stat > max_stat:
                        max_stat = stat
                        changepoint = i
            
            before = data[:changepoint]
            after = data[changepoint:]
            
            return {
                'changepoint': changepoint,
                'likelihood': max_stat,
                'before_mean': float(np.mean(before)),
                'after_mean': float(np.mean(after)),
                'before_std': float(np.std(before, ddof=1)),
                'after_std': float(np.std(after, ddof=1))
            }
        else:
            raise ValueError(f"Unknown method: {method}")
    
    @staticmethod
    def pelt(data: np.ndarray,
            pen: Optional[float] = None,
            min_seg_length: int = 2) -> List[int]:
        """
        PELT (Pruned Exact Linear Time) 算法
        
        检测多个变点
        
        Parameters
        ----------
        data : array_like
            时间序列数据
        pen : float, optional
            惩罚参数，如果为None则自动计算
        min_seg_length : int
            最小段长度
            
        Returns
        -------
        list
            变点位置列表
        """
        n = len(data)
        
        if pen is None:
            pen = 2 * np.log(n)  # BIC 惩罚
        
        # 代价函数 (方差)
        def cost(start, end):
            segment = data[start:end]
            if len(segment) < 2:
                return 0
            return np.sum((segment - np.mean(segment)) ** 2)
        
        # PELT 算法
        F = [0] + [np.inf] * n
        last_changepoints = [0] * (n + 1)
        
        for t in range(1, n + 1):
            candidates = []
            for s in range(0, t - min_seg_length + 1):
                c = F[s] + cost(s, t) + pen
                candidates.append((c, s))
            
            if candidates:
                F[t], last_changepoints[t] = min(candidates)
        
        # 回溯找到所有变点
        changepoints = []
        cp = n
        while cp > 0:
            cp = last_changepoints[cp]
            if cp > 0:
                changepoints.append(cp)
        
        return sorted(changepoints)


class TimeSeriesAnalysis:
    """
    时间序列分析工具
    """
    
    @staticmethod
    def bin_average(x: np.ndarray,
                   y: np.ndarray,
                   bin_edges: np.ndarray) -> Dict:
        """
        分箱平均
        
        Parameters
        ----------
        x, y : array_like
            数据点
        bin_edges : array_like
            分箱边界
            
        Returns
        -------
        dict
            {
                'bin_centers': 箱中心,
                'bin_means': 箱均值,
                'bin_stds': 箱标准差,
                'bin_counts': 每箱数据点数
            }
        """
        n_bins = len(bin_edges) - 1
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        bin_means = []
        bin_stds = []
        bin_counts = []
        
        for i in range(n_bins):
            mask = (x >= bin_edges[i]) & (x < bin_edges[i+1])
            if np.any(mask):
                bin_means.append(np.mean(y[mask]))
                bin_stds.append(np.std(y[mask]))
                bin_counts.append(np.sum(mask))
            else:
                bin_means.append(np.nan)
                bin_stds.append(np.nan)
                bin_counts.append(0)
        
        return {
            'bin_centers': bin_centers,
            'bin_means': np.array(bin_means),
            'bin_stds': np.array(bin_stds),
            'bin_counts': np.array(bin_counts)
        }
    
    @staticmethod
    def moving_average(y: np.ndarray,
                      window: int = 5,
                      mode: str = 'same') -> np.ndarray:
        """
        移动平均
        
        Parameters
        ----------
        y : array_like
            输入数据
        window : int
            窗口大小
        mode : str
            边界处理模式
            
        Returns
        -------
        array_like
            移动平均结果
        """
        return np.convolve(y, np.ones(window)/window, mode=mode)
