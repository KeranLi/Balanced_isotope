"""
氮同位素体系模型
实现基于 Kang et al. (2023) 和 Ma et al. (2025) 的双箱稳态氮循环模型

该模型用于：
1. 从硝酸盐占比(f_assimilator)计算沉积物氮同位素(δ¹⁵N_sed)
2. 从沉积物氮同位素反演硝酸盐占比
3. 蒙特卡洛不确定性分析
"""

import numpy as np
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass

from systems.base.isotope_system import IsotopeSystem, ModelResult, IsotopeParameters
from systems.n.parameters import get_n_parameters, get_fractionation_ranges


@dataclass
class NitrogenCycleFluxes:
    """氮循环通量数据类 (单位: Tg N/a)"""
    F_fix: float = 205.0           # 固氮通量
    F_total_burial: float = 25.0    # 总埋藏通量
    F_wcd: float = 140.0           # 水柱反硝化通量
    F_sd: float = 40.0             # 沉积反硝化通量
    
    @property
    def F_remin(self) -> float:
        """再矿化通量 = 固氮通量 - 固氮生物埋藏"""
        # 注意：这里假设 F_fixer_burial 很小，近似 F_fix ≈ F_remin
        return self.F_fix
    
    @property
    def F_total_denit(self) -> float:
        """总反硝化通量"""
        return self.F_wcd + self.F_sd


@dataclass
class FractionationFactors:
    """分馏系数数据类 (‰)"""
    epsilon_fix: float = -0.5       # 固氮分馏
    epsilon_wcd: float = -26.0      # 水柱反硝化分馏
    epsilon_sd: float = 0.0         # 沉积反硝化分馏
    
    def to_alpha(self, epsilon_name: str) -> float:
        """将ε转换为α"""
        epsilon = getattr(self, epsilon_name)
        return 1 + epsilon / 1000


class NIsotopeSystem(IsotopeSystem):
    """
    氮同位素体系
    
    基于双箱稳态模型，模拟海洋氮循环和同位素分馏
    
    主要应用：
    1. 计算沉积物氮同位素组成
    2. 估算硝酸盐可利用性 (f_assimilator)
    3. 评估氮循环演化
    
    Attributes
    ----------
    ELEMENT : str
        元素符号 'n'
    NAME : str
        元素名称 'Nitrogen'
    ISOTOPES : list
        涉及的同位素 ['14N', '15N']
    """
    
    ELEMENT = 'n'
    NAME = 'Nitrogen'
    ISOTOPES = ['14N', '15N']
    
    # 大气氮同位素标准值
    DELTA15N_ATMOSPHERE = 0.0
    
    def __init__(self, parameters: Optional[IsotopeParameters] = None,
                 scenario: str = 'modern'):
        """
        初始化氮同位素体系
        
        Parameters
        ----------
        parameters : IsotopeParameters, optional
            自定义参数，如果为None则使用默认参数
        scenario : str
            情景名称：'modern', 'early_triassic', 'neoproterozoic'
        """
        self.scenario = scenario
        super().__init__(parameters or get_n_parameters(scenario))
        
        # 初始化通量和分馏系数
        self.fluxes = self._init_fluxes()
        self.fractionation = self._init_fractionation()
    
    def _default_parameters(self) -> IsotopeParameters:
        """返回默认参数"""
        return get_n_parameters(self.scenario)
    
    def _init_fluxes(self) -> NitrogenCycleFluxes:
        """初始化氮循环通量"""
        params = self.params
        return NitrogenCycleFluxes(
            F_fix=params.input_fluxes.get('fixation', 205.0),
            F_total_burial=params.output_fluxes.get('burial', 25.0),
            F_wcd=params.output_fluxes.get('water_column_denitrification', 140.0),
            F_sd=params.output_fluxes.get('sedimentary_denitrification', 40.0)
        )
    
    def _init_fractionation(self) -> FractionationFactors:
        """初始化分馏系数"""
        ff = self.params.fractionation_factors
        return FractionationFactors(
            epsilon_fix=ff.get('epsilon_fixation', -0.5),
            epsilon_wcd=ff.get('epsilon_wcd', -26.0),
            epsilon_sd=ff.get('epsilon_sd', 0.0)
        )
    
    # ============== 核心计算接口实现 ==============
    
    def mass_balance_equation(self,
                             state: np.ndarray,
                             fluxes: Dict[str, float],
                             time: Optional[float] = None) -> np.ndarray:
        """
        氮同位素质量平衡微分方程 (稳态时 d/dt = 0)
        
        注意：本模型主要使用稳态解，此方法主要用于接口兼容
        
        Parameters
        ----------
        state : array_like, shape (2,)
            [delta15N_ammonium, delta15N_nitrate]
        fluxes : dict
            通量字典
        time : float, optional
            时间 (不适用)
            
        Returns
        -------
        array_like, shape (2,)
            d(state)/dt (稳态时应该接近0)
        """
        # 在稳态模型中，状态不随时间变化
        return np.array([0.0, 0.0])
    
    def fractionation_factor(self,
                            process: str,
                            temperature: Optional[float] = None,
                            **kwargs) -> float:
        """
        获取氮同位素分馏系数 α
        
        Parameters
        ----------
        process : str
            过程名称：
            - 'fixation': 固氮作用
            - 'water_column_denitrification': 水柱反硝化
            - 'sedimentary_denitrification': 沉积反硝化
            - 'nitrification': 硝化作用
            - 'assimilation': 同化作用
        temperature : float, optional
            温度 (K)，目前不使用
        **kwargs
            其他参数
            
        Returns
        -------
        float
            分馏系数 α
        """
        epsilon_map = {
            'fixation': self.fractionation.epsilon_fix,
            'nitrogen_fixation': self.fractionation.epsilon_fix,
            'water_column_denitrification': self.fractionation.epsilon_wcd,
            'wcd': self.fractionation.epsilon_wcd,
            'sedimentary_denitrification': self.fractionation.epsilon_sd,
            'sd': self.fractionation.epsilon_sd,
        }
        
        epsilon = epsilon_map.get(process, 0.0)
        return 1 + epsilon / 1000
    
    def mixing_model(self,
                    end_member_values: np.ndarray,
                    proportions: np.ndarray) -> float:
        """
        端元混合模型
        
        Parameters
        ----------
        end_member_values : array_like
            各端元的δ¹⁵N值
        proportions : array_like
            各端元的比例（和为1）
            
        Returns
        -------
        float
            混合后的δ¹⁵N值
        """
        from toolkit.isotope.formulas import MassBalance
        
        return MassBalance.multi_component_mixing(
            end_member_values, proportions,
            standard_ratio=self.params.reference_ratios['15/14']
        )
    
    # ============== 氮同位素专用方法 ==============
    
    def calculate_reservoir_isotopes(self,
                                     f_assimilator: float,
                                     epsilon_fix: Optional[float] = None,
                                     epsilon_wcd: Optional[float] = None,
                                     epsilon_sd: Optional[float] = None) -> Dict[str, float]:
        """
        计算储库同位素组成
        
        基于双箱模型方程 (3) 和 (4)
        
        关键关系：
        - f_assimilator 增加 → 硝酸盐更充足 → 反硝化减少 → δ¹⁵N_nitrate 降低
        - 这导致 δ¹⁵N_sed 与 f_assimilator 呈非线性关系，在 f ≈ 0.48 时达到峰值
        
        Parameters
        ----------
        f_assimilator : float
            硝酸盐同化埋藏通量比例 (0-1)
        epsilon_fix : float, optional
            固氮分馏系数 (‰)，默认使用模型设定值
        epsilon_wcd : float, optional
            水柱反硝化分馏系数 (‰)
        epsilon_sd : float, optional
            沉积反硝化分馏系数 (‰)
            
        Returns
        -------
        dict
            {
                'delta15N_ammonium': 铵储库δ¹⁵N,
                'delta15N_nitrate': 硝酸盐储库δ¹⁵N
            }
        """
        # 使用默认或自定义分馏系数
        eps_fix = epsilon_fix if epsilon_fix is not None else self.fractionation.epsilon_fix
        eps_wcd = epsilon_wcd if epsilon_wcd is not None else self.fractionation.epsilon_wcd
        eps_sd = epsilon_sd if epsilon_sd is not None else self.fractionation.epsilon_sd
        
        # 计算铵储库同位素 (方程 3)
        delta15N_ammonium = self.DELTA15N_ATMOSPHERE + eps_fix
        
        # 根据 Kang et al. (2023) 模型:
        # F_assimilator_burial 与总反硝化通量成反比
        # 即 f_assimilator ↑ → (F_wcd + F_sd) ↓
        # 
        # 简化的线性关系：当 f 从 0 增加到 1 时，
        # 水柱反硝化从最大值线性减少到 0
        # 
        # 现代海洋参考值：f ≈ 0.7 时，F_wcd ≈ 140, F_sd ≈ 40
        # 完全缺氧时 (f → 0)：F_wcd 最大，F_sd 保持相对稳定
        # 完全氧化时 (f → 1)：F_wcd → 0，F_sd 保持相对稳定
        
        # 假设最大水柱反硝化通量 (完全缺氧)
        F_wcd_max = self.fluxes.F_wcd * 1.5  # 约 210 Tg/a
        F_wcd_min = 0.0  # 完全氧化时
        
        # 水柱反硝化通量随 f 增加而减少 (线性近似)
        # 注意：这里使用 (1-f) 的关系，f 越大，反硝化越少
        F_wcd_effective = F_wcd_max * (1 - f_assimilator)
        
        # 沉积反硝化通量相对稳定
        F_sd_effective = self.fluxes.F_sd
        
        # 再矿化通量 (近似等于固氮通量)
        F_remin = self.fluxes.F_fix
        
        # 计算硝酸盐储库同位素 (方程 4)
        # δ¹⁵N_nitrate = δ¹⁵N_ammonium - (F_wcd*ε_wcd + F_sd*ε_sd) / F_remin
        weighted_epsilon = (F_wcd_effective * eps_wcd + 
                           F_sd_effective * eps_sd) / F_remin
        delta15N_nitrate = delta15N_ammonium - weighted_epsilon
        
        return {
            'delta15N_ammonium': delta15N_ammonium,
            'delta15N_nitrate': delta15N_nitrate
        }
    
    def forward_model(self,
                     f_assimilator: float,
                     epsilon_fix: Optional[float] = None,
                     epsilon_wcd: Optional[float] = None,
                     epsilon_sd: Optional[float] = None) -> float:
        """
        正向模型：从硝酸盐占比计算沉积物δ¹⁵N
        
        基于方程 (6)
        
        Parameters
        ----------
        f_assimilator : float
            硝酸盐同化埋藏通量比例 (0-1)
        epsilon_fix : float, optional
            固氮分馏系数 (‰)
        epsilon_wcd : float, optional
            水柱反硝化分馏系数 (‰)
        epsilon_sd : float, optional
            沉积反硝化分馏系数 (‰)
            
        Returns
        -------
        float
            沉积物δ¹⁵N (‰)
        """
        # 获取储库同位素值
        reservoirs = self.calculate_reservoir_isotopes(
            f_assimilator, epsilon_fix, epsilon_wcd, epsilon_sd
        )
        
        delta15N_ammonium = reservoirs['delta15N_ammonium']
        delta15N_nitrate = reservoirs['delta15N_nitrate']
        
        # 计算沉积物同位素 (方程 6)
        delta15N_sed = ((1 - f_assimilator) * delta15N_ammonium + 
                       f_assimilator * delta15N_nitrate)
        
        return delta15N_sed
    
    def inverse_model(self,
                     delta15N_sed: float,
                     epsilon_fix: Optional[float] = None,
                     epsilon_wcd: Optional[float] = None,
                     epsilon_sd: Optional[float] = None,
                     f_range: Tuple[float, float] = (0.0, 0.48)) -> Dict:
        """
        反向模型（反演）：从沉积物δ¹⁵N计算硝酸盐占比
        
        注意：由于非线性关系，相同的δ¹⁵N_sed可能对应两个f_assimilator值
        在缺氧环境下，通常选择较小的值
        
        Parameters
        ----------
        delta15N_sed : float
            沉积物δ¹⁵N (‰)
        epsilon_fix : float, optional
            固氮分馏系数 (‰)
        epsilon_wcd : float, optional
            水柱反硝化分馏系数 (‰)
        epsilon_sd : float, optional
            沉积反硝化分馏系数 (‰)
        f_range : tuple
            合理的f_assimilator范围，默认(0.0, 0.48)对应缺氧环境
            
        Returns
        -------
        dict
            {
                'f_assimilator': 硝酸盐占比,
                'delta15N_sed_calculated': 计算的δ¹⁵N_sed,
                'residual': 残差
            }
        """
        from scipy.optimize import minimize_scalar
        
        def objective(f):
            calculated = self.forward_model(f, epsilon_fix, epsilon_wcd, epsilon_sd)
            return (calculated - delta15N_sed) ** 2
        
        # 在合理范围内寻找最优解
        result = minimize_scalar(objective, bounds=f_range, method='bounded')
        
        f_optimal = result.x
        delta_calc = self.forward_model(f_optimal, epsilon_fix, epsilon_wcd, epsilon_sd)
        
        return {
            'f_assimilator': float(f_optimal),
            'delta15N_sed_calculated': float(delta_calc),
            'residual': float(delta_calc - delta15N_sed)
        }
    
    def monte_carlo_simulation(self,
                              f_assimilator: float,
                              n_samples: int = 10000,
                              epsilon_fix_range: Tuple[float, float] = (-2.0, 1.0),
                              epsilon_wcd_range: Tuple[float, float] = (-30.0, -22.0)) -> Dict:
        """
        蒙特卡洛模拟评估不确定性
        
        Parameters
        ----------
        f_assimilator : float
            硝酸盐同化埋藏通量比例
        n_samples : int
            采样次数
        epsilon_fix_range : tuple
            固氮分馏系数范围 (min, max)
        epsilon_wcd_range : tuple
            水柱反硝化分馏系数范围 (min, max)
            
        Returns
        -------
        dict
            {
                'delta15N_sed_mean': 均值,
                'delta15N_sed_std': 标准差,
                'delta15N_sed_median': 中位数,
                'delta15N_sed_ci68': 68%置信区间,
                'delta15N_sed_ci95': 95%置信区间,
                'samples': 所有样本结果
            }
        """
        # 从均匀分布中采样分馏系数
        epsilon_fix_samples = np.random.uniform(
            epsilon_fix_range[0], epsilon_fix_range[1], n_samples
        )
        epsilon_wcd_samples = np.random.uniform(
            epsilon_wcd_range[0], epsilon_wcd_range[1], n_samples
        )
        
        # 计算沉积物同位素
        delta15N_sed_samples = np.array([
            self.forward_model(f_assimilator, eps_fix, eps_wcd)
            for eps_fix, eps_wcd in zip(epsilon_fix_samples, epsilon_wcd_samples)
        ])
        
        # 统计结果
        mean_val = np.mean(delta15N_sed_samples)
        std_val = np.std(delta15N_sed_samples)
        median_val = np.median(delta15N_sed_samples)
        ci68 = (np.percentile(delta15N_sed_samples, 16), 
                np.percentile(delta15N_sed_samples, 84))
        ci95 = (np.percentile(delta15N_sed_samples, 2.5), 
                np.percentile(delta15N_sed_samples, 97.5))
        
        return {
            'delta15N_sed_mean': float(mean_val),
            'delta15N_sed_std': float(std_val),
            'delta15N_sed_median': float(median_val),
            'delta15N_sed_ci68': (float(ci68[0]), float(ci68[1])),
            'delta15N_sed_ci95': (float(ci95[0]), float(ci95[1])),
            'samples': delta15N_sed_samples
        }
    
    def calculate_f_assimilator_curve(self,
                                     f_range: Tuple[float, float] = (0.0, 1.0),
                                     n_points: int = 100,
                                     n_monte_carlo: int = 1000) -> Dict:
        """
        计算 f_assimilator 与 δ¹⁵N_sed 的关系曲线
        
        Parameters
        ----------
        f_range : tuple
            f_assimilator范围
        n_points : int
            曲线点数
        n_monte_carlo : int
            每个点的蒙特卡洛采样次数
            
        Returns
        -------
        dict
            {
                'f_assimilator': f值数组,
                'delta15N_sed_mean': 均值曲线,
                'delta15N_sed_ci68': 68%置信区间,
                'delta15N_sed_ci95': 95%置信区间
            }
        """
        f_values = np.linspace(f_range[0], f_range[1], n_points)
        
        mean_curve = []
        ci68_lower = []
        ci68_upper = []
        ci95_lower = []
        ci95_upper = []
        
        for f in f_values:
            mc_result = self.monte_carlo_simulation(f, n_samples=n_monte_carlo)
            mean_curve.append(mc_result['delta15N_sed_mean'])
            ci68_lower.append(mc_result['delta15N_sed_ci68'][0])
            ci68_upper.append(mc_result['delta15N_sed_ci68'][1])
            ci95_lower.append(mc_result['delta15N_sed_ci95'][0])
            ci95_upper.append(mc_result['delta15N_sed_ci95'][1])
        
        return {
            'f_assimilator': f_values,
            'delta15N_sed_mean': np.array(mean_curve),
            'delta15N_sed_ci68_lower': np.array(ci68_lower),
            'delta15N_sed_ci68_upper': np.array(ci68_upper),
            'delta15N_sed_ci95_lower': np.array(ci95_lower),
            'delta15N_sed_ci95_upper': np.array(ci95_upper)
        }
    
    # ============== 工具方法实现 ==============
    
    def state_dimension(self) -> int:
        """状态变量维度：δ¹⁵N_ammonium, δ¹⁵N_nitrate"""
        return 2
    
    def validate_data(self, data: Dict[str, np.ndarray]) -> Tuple[bool, str]:
        """验证氮同位素数据"""
        if 'delta15N' not in data and 'delta_15_n' not in data:
            return False, "Missing nitrogen isotope data (delta15N or delta_15_n)"
        
        delta_values = data.get('delta15N', data.get('delta_15_n'))
        if np.any(delta_values < -10) or np.any(delta_values > 15):
            return False, "Delta15N values out of reasonable range (-10 to 15)"
        
        return True, "Validation passed"
    
    def get_model_info(self) -> Dict:
        """获取模型信息"""
        return {
            'element': self.ELEMENT,
            'name': self.NAME,
            'fluxes': {
                'F_fix': self.fluxes.F_fix,
                'F_total_burial': self.fluxes.F_total_burial,
                'F_wcd': self.fluxes.F_wcd,
                'F_sd': self.fluxes.F_sd
            },
            'fractionation': {
                'epsilon_fix': self.fractionation.epsilon_fix,
                'epsilon_wcd': self.fractionation.epsilon_wcd,
                'epsilon_sd': self.fractionation.epsilon_sd
            }
        }
