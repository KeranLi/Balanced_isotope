"""
铀同位素模型不确定度分析模块

提供蒙特卡洛模拟、Bootstrap和敏感性分析等方法
用于评估模型参数不确定度对结果的影响

基于文献中的参数不确定度范围:
- 海水 δ²³⁸U: -0.392 ± 0.005 ‰ (Tissot & Dauphas, 2015)
- 河流 δ²³⁸U: -0.27 ± 0.16 ‰ (Andersen et al., 2016)
- Δ_sw-anox: 0.77 ± 0.04 ‰ (Stirling et al., 2015)
- Δ_diag: 0.30-0.50 ‰ (Elrick et al., 2017)
"""

import numpy as np
from typing import Dict, Tuple, Optional, List, Callable, Union
from dataclasses import dataclass
from scipy import stats


@dataclass
class ParameterUncertainty:
    """参数不确定度数据类"""
    mean: float
    std: Optional[float] = None  # 标准差 (正态分布)
    lower: Optional[float] = None  # 下限 (均匀分布)
    upper: Optional[float] = None  # 上限 (均匀分布)
    distribution: str = 'normal'  # 'normal' 或 'uniform'
    
    def sample(self, n: int = 1) -> np.ndarray:
        """从分布中采样"""
        if self.distribution == 'normal':
            if self.std is None:
                raise ValueError("Standard deviation required for normal distribution")
            return np.random.normal(self.mean, self.std, n)
        elif self.distribution == 'uniform':
            if self.lower is None or self.upper is None:
                raise ValueError("Lower and upper bounds required for uniform distribution")
            return np.random.uniform(self.lower, self.upper, n)
        else:
            raise ValueError(f"Unknown distribution: {self.distribution}")


@dataclass
class UncertaintyConfig:
    """不确定度配置"""
    # 分馏系数不确定度
    delta_sw_ox: ParameterUncertainty = None
    delta_sw_anox: ParameterUncertainty = None
    delta_diag: ParameterUncertainty = None
    delta_river: ParameterUncertainty = None
    
    # 测量值不确定度
    delta_measurement: ParameterUncertainty = None
    
    def __post_init__(self):
        # 设置默认值
        if self.delta_sw_ox is None:
            self.delta_sw_ox = ParameterUncertainty(
                mean=0.0, std=0.0, distribution='normal'
            )
        if self.delta_sw_anox is None:
            self.delta_sw_anox = ParameterUncertainty(
                mean=0.77, std=0.04, distribution='normal'
            )
        if self.delta_diag is None:
            self.delta_diag = ParameterUncertainty(
                mean=0.40, lower=0.30, upper=0.50, distribution='uniform'
            )
        if self.delta_river is None:
            self.delta_river = ParameterUncertainty(
                mean=-0.29, std=0.16, distribution='normal'
            )
        if self.delta_measurement is None:
            self.delta_measurement = ParameterUncertainty(
                mean=0.0, std=0.05, distribution='normal'
            )


class UncertaintyAnalyzer:
    """
    不确定度分析器
    
    为U同位素模型提供蒙特卡洛、Bootstrap和敏感性分析
    """
    
    def __init__(self, u_system):
        """
        初始化不确定度分析器
        
        Parameters
        ----------
        u_system : UIsotopeSystem
            铀同位素体系实例
        """
        self.u_system = u_system
        self.config = UncertaintyConfig()
    
    def monte_carlo_steady_state(
        self,
        delta238_carb: float,
        n_samples: int = 10000,
        confidence_level: float = 0.95,
        apply_diagenetic_correction: bool = True,
        random_seed: Optional[int] = None
    ) -> Dict:
        """
        蒙特卡洛模拟稳态模型不确定度
        
        通过随机采样参数，评估f_anox的不确定度范围
        
        Parameters
        ----------
        delta238_carb : float
            实测碳酸盐δ²³⁸U值
        n_samples : int
            蒙特卡洛采样次数
        confidence_level : float
            置信水平 (默认95%)
        apply_diagenetic_correction : bool
            是否应用成岩校正
        random_seed : int, optional
            随机种子，用于结果复现
            
        Returns
        -------
        dict
            {
                'f_anox_mean': f_anox均值,
                'f_anox_std': 标准差,
                'f_anox_ci': (下限, 上限) 置信区间,
                'f_anox_median': 中位数,
                'f_anox_samples': 所有样本,
                'parameter_samples': 参数采样记录,
                'convergence': 收敛性指标
            }
        """
        if random_seed is not None:
            np.random.seed(random_seed)
        
        # 保存原始参数
        original_fractionation = {
            'delta_sw_ox': self.u_system.fractionation.delta_sw_ox,
            'delta_sw_anox': self.u_system.fractionation.delta_sw_anox,
            'delta_diag': self.u_system.fractionation.delta_diag,
            'delta_river': self.u_system.fractionation.delta_river,
        }
        
        # 采样并计算
        f_anox_samples = []
        param_samples = {
            'delta_sw_anox': [],
            'delta_diag': [],
            'delta_river': [],
            'delta_meas': []
        }
        
        for _ in range(n_samples):
            # 采样参数
            delta_sw_anox = self.config.delta_sw_anox.sample(1)[0]
            delta_diag = self.config.delta_diag.sample(1)[0]
            delta_river = self.config.delta_river.sample(1)[0]
            delta_meas_noise = self.config.delta_measurement.sample(1)[0]
            
            delta_meas_sampled = delta238_carb + delta_meas_noise
            
            # 临时修改参数
            self.u_system.fractionation.delta_sw_anox = delta_sw_anox
            self.u_system.fractionation.delta_diag = delta_diag
            self.u_system.fractionation.delta_river = delta_river
            
            # 计算f_anox
            try:
                result = self.u_system.calculate_f_anox_steady_state(
                    delta238_carb=delta_meas_sampled,
                    apply_diagenetic_correction=apply_diagenetic_correction,
                    delta_diag=delta_diag
                )
                f_anox_samples.append(result['f_anox'])
                
                # 记录参数
                param_samples['delta_sw_anox'].append(delta_sw_anox)
                param_samples['delta_diag'].append(delta_diag)
                param_samples['delta_river'].append(delta_river)
                param_samples['delta_meas'].append(delta_meas_sampled)
            except Exception:
                continue
        
        # 恢复原始参数
        for key, value in original_fractionation.items():
            setattr(self.u_system.fractionation, key, value)
        
        f_anox_samples = np.array(f_anox_samples)
        
        # 计算统计量
        alpha = 1 - confidence_level
        ci_lower = np.percentile(f_anox_samples, alpha/2 * 100)
        ci_upper = np.percentile(f_anox_samples, (1 - alpha/2) * 100)
        
        # 收敛性检查 (Gelman-Rubin-like statistic)
        n_chains = 4
        chain_size = len(f_anox_samples) // n_chains
        if chain_size > 100:
            chains = f_anox_samples[:n_chains*chain_size].reshape(n_chains, chain_size)
            chain_means = np.mean(chains, axis=1)
            chain_vars = np.var(chains, axis=1, ddof=1)
            B = chain_size * np.var(chain_means, ddof=1)
            W = np.mean(chain_vars)
            V_hat = (1 - 1/chain_size) * W + B/chain_size
            r_hat = np.sqrt(V_hat / W) if W > 0 else 1.0
        else:
            r_hat = 1.0
        
        return {
            'f_anox_mean': float(np.mean(f_anox_samples)),
            'f_anox_std': float(np.std(f_anox_samples)),
            'f_anox_median': float(np.median(f_anox_samples)),
            'f_anox_ci': (float(ci_lower), float(ci_upper)),
            'f_anox_samples': f_anox_samples,
            'parameter_samples': {k: np.array(v) for k, v in param_samples.items()},
            'convergence': {
                'r_hat': float(r_hat),
                'n_effective': len(f_anox_samples),
                'converged': r_hat < 1.1
            }
        }
    
    def bootstrap_analysis(
        self,
        delta238_measurements: np.ndarray,
        n_bootstrap: int = 10000,
        confidence_level: float = 0.95,
        random_seed: Optional[int] = None
    ) -> Dict:
        """
        Bootstrap重采样分析
        
        适用于有多个测量样本的情况，评估统计不确定度
        
        Parameters
        ----------
        delta238_measurements : array_like
            多个δ²³⁸U测量值数组
        n_bootstrap : int
            Bootstrap重采样次数
        confidence_level : float
            置信水平
        random_seed : int, optional
            随机种子
            
        Returns
        -------
        dict
            {
                'f_anox_mean': 均值,
                'f_anox_std': 标准差,
                'f_anox_ci': 置信区间,
                'f_anox_bootstrap': Bootstrap样本,
                'original_estimate': 原始样本估计值
            }
        """
        if random_seed is not None:
            np.random.seed(random_seed)
        
        n_measurements = len(delta238_measurements)
        
        # 原始样本估计
        original_mean = np.mean(delta238_measurements)
        original_result = self.u_system.calculate_f_anox_steady_state(original_mean)
        original_f_anox = original_result['f_anox']
        
        # Bootstrap重采样
        bootstrap_f_anox = []
        
        for _ in range(n_bootstrap):
            # 有放回地重采样
            resampled = np.random.choice(
                delta238_measurements, 
                size=n_measurements, 
                replace=True
            )
            resampled_mean = np.mean(resampled)
            
            try:
                result = self.u_system.calculate_f_anox_steady_state(resampled_mean)
                bootstrap_f_anox.append(result['f_anox'])
            except Exception:
                continue
        
        bootstrap_f_anox = np.array(bootstrap_f_anox)
        
        # 偏差校正 (BCa方法简化版)
        bias = np.mean(bootstrap_f_anox) - original_f_anox
        
        # 百分位置信区间
        alpha = 1 - confidence_level
        ci_lower = np.percentile(bootstrap_f_anox, alpha/2 * 100)
        ci_upper = np.percentile(bootstrap_f_anox, (1 - alpha/2) * 100)
        
        return {
            'f_anox_mean': float(np.mean(bootstrap_f_anox)),
            'f_anox_std': float(np.std(bootstrap_f_anox)),
            'f_anox_ci': (float(ci_lower), float(ci_upper)),
            'f_anox_bootstrap': bootstrap_f_anox,
            'original_estimate': float(original_f_anox),
            'bias': float(bias)
        }
    
    def sensitivity_analysis(
        self,
        delta238_carb: float,
        parameter_ranges: Optional[Dict[str, Tuple[float, float]]] = None,
        n_points: int = 50
    ) -> Dict:
        """
        敏感性分析
        
        评估各参数变化对 f_anox 结果的影响程度
        
        Parameters
        ----------
        delta238_carb : float
            实测碳酸盐δ²³⁸U值
        parameter_ranges : dict, optional
            参数变化范围，如 {'delta_sw_anox': (0.6, 0.9)}
        n_points : int
            每个参数的分析点数
            
        Returns
        -------
        dict
            {
                'sensitivities': 各参数的敏感性系数,
                'tornado_data': 龙卷风图数据
            }
        """
        if parameter_ranges is None:
            # 默认参数范围（基于文献不确定度）
            parameter_ranges = {
                'delta_sw_anox': (0.69, 0.85),  # 0.77 ± 0.08
                'delta_diag': (0.30, 0.50),      # 0.40 ± 0.10
                'delta_river': (-0.45, -0.13),   # -0.29 ± 0.16
            }
        
        # 基准计算
        baseline = self.u_system.calculate_f_anox_steady_state(delta238_carb)
        baseline_f_anox = baseline['f_anox']
        
        sensitivities = {}
        tornado_data = []
        
        for param_name, (p_min, p_max) in parameter_ranges.items():
            # 保存原始值
            original_value = getattr(self.u_system.fractionation, param_name)
            
            # 计算参数变化的影响
            values = np.linspace(p_min, p_max, n_points)
            f_anox_values = []
            
            for val in values:
                setattr(self.u_system.fractionation, param_name, val)
                try:
                    result = self.u_system.calculate_f_anox_steady_state(delta238_carb)
                    f_anox_values.append(result['f_anox'])
                except Exception:
                    f_anox_values.append(np.nan)
            
            # 恢复原始值
            setattr(self.u_system.fractionation, param_name, original_value)
            
            f_anox_values = np.array(f_anox_values)
            valid_mask = ~np.isnan(f_anox_values)
            
            if np.sum(valid_mask) > 1:
                # 计算敏感性系数 (Δf_anox / Δparam)
                sensitivity = np.polyfit(
                    values[valid_mask], 
                    f_anox_values[valid_mask], 
                    1
                )[0]
                
                # 归一化敏感性 (考虑参数典型范围)
                param_range = p_max - p_min
                normalized_sensitivity = abs(sensitivity) * param_range
                
                sensitivities[param_name] = {
                    'slope': float(sensitivity),
                    'normalized': float(normalized_sensitivity),
                    'f_anox_range': (float(np.min(f_anox_values[valid_mask])),
                                    float(np.max(f_anox_values[valid_mask])))
                }
                
                tornado_data.append({
                    'parameter': param_name,
                    'min_effect': float(np.min(f_anox_values[valid_mask]) - baseline_f_anox),
                    'max_effect': float(np.max(f_anox_values[valid_mask]) - baseline_f_anox),
                    'range': float(normalized_sensitivity)
                })
        
        # 按影响大小排序
        tornado_data.sort(key=lambda x: abs(x['range']), reverse=True)
        
        return {
            'baseline_f_anox': float(baseline_f_anox),
            'sensitivities': sensitivities,
            'tornado_data': tornado_data,
            'most_important': tornado_data[0]['parameter'] if tornado_data else None
        }
    
    def measurement_uncertainty_only(
        self,
        delta238_carb: float,
        measurement_std: float,
        n_samples: int = 10000
    ) -> Dict:
        """
        仅考虑测量值不确定度
        
        当参数确定，仅测量值有误差时使用
        
        Parameters
        ----------
        delta238_carb : float
            实测值
        measurement_std : float
            测量标准差
        n_samples : int
            采样次数
            
        Returns
        -------
        dict
            不确定度分析结果
        """
        # 临时保存配置
        original_config = self.config
        
        # 仅设置测量不确定度，其他参数确定
        self.config = UncertaintyConfig(
            delta_sw_ox=ParameterUncertainty(mean=0.0, std=0.0),
            delta_sw_anox=ParameterUncertainty(mean=self.u_system.fractionation.delta_sw_anox, std=0.0),
            delta_diag=ParameterUncertainty(mean=self.u_system.fractionation.delta_diag, std=0.0),
            delta_river=ParameterUncertainty(mean=self.u_system.fractionation.delta_river, std=0.0),
            delta_measurement=ParameterUncertainty(mean=0.0, std=measurement_std)
        )
        
        result = self.monte_carlo_steady_state(
            delta238_carb=delta238_carb,
            n_samples=n_samples
        )
        
        # 恢复配置
        self.config = original_config
        
        return result


def analyze_likelihood(
    u_system,
    observed_delta238: float,
    f_anox_grid: np.ndarray,
    measurement_std: float = 0.05
) -> Dict:
    """
    基于似然的f_anox估计
    
    计算在给定观测值下，不同f_anox值的似然度
    
    Parameters
    ----------
    u_system : UIsotopeSystem
        铀同位素体系
    observed_delta238 : float
        观测的δ²³⁸U值
    f_anox_grid : array_like
        f_anox搜索网格
    measurement_std : float
        测量标准差
        
    Returns
    -------
    dict
        {
            'f_anox_mle': 最大似然估计,
            'log_likelihood': 对数似然值,
            'likelihood': 似然值 (归一化),
            'credible_interval': 可信区间
        }
    """
    log_likelihoods = []
    
    for f_anox in f_anox_grid:
        # 正向计算预测的δ值
        result = u_system.calculate_seawater_delta_steady_state(f_anox=f_anox)
        predicted_delta = result['delta238_carbonate']
        
        # 计算对数似然 (高斯似然)
        log_lik = -0.5 * ((observed_delta238 - predicted_delta) / measurement_std) ** 2
        log_likelihoods.append(log_lik)
    
    log_likelihoods = np.array(log_likelihoods)
    likelihoods = np.exp(log_likelihoods - np.max(log_likelihoods))
    likelihoods = likelihoods / np.sum(likelihoods)  # 归一化
    
    # 最大似然估计
    mle_idx = np.argmax(likelihoods)
    f_anox_mle = f_anox_grid[mle_idx]
    
    # 计算95%可信区间
    cumulative = np.cumsum(likelihoods)
    ci_lower_idx = np.searchsorted(cumulative, 0.025)
    ci_upper_idx = np.searchsorted(cumulative, 0.975)
    
    return {
        'f_anox_mle': float(f_anox_mle),
        'f_anox_grid': f_anox_grid,
        'log_likelihood': log_likelihoods,
        'likelihood': likelihoods,
        'credible_interval': (float(f_anox_grid[ci_lower_idx]), 
                             float(f_anox_grid[ci_upper_idx]))
    }
