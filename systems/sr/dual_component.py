"""
Sr同位素双组分河流模型

用于模拟峨眉山LIP期间的玄武岩风化脉冲

关键创新:
- 将河流输入拆分为地壳和玄武岩两个端元
- 玄武岩风化使用高斯脉冲函数模拟
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, Callable
from scipy.interpolate import interp1d


@dataclass
class BasaltPulseConfig:
    """
    玄武岩风化脉冲配置

    高斯脉冲函数:
    F_basalt(t) = F_max * exp(-(t - t_peak)^2 / (2 * sigma^2)) * H(t - t_onset)
    """
    F_max: float = 20e9          # 峰值通量 (mol/yr)
    t_onset: float = 265.5       # 脉冲开始时间 (Ma)
    t_peak: float = 263.5        # 峰值时间 (Ma)
    sigma: float = 1.5           # 脉冲宽度 (Myr)
    R_basalt: float = 0.7040     # 玄武岩同位素比值

    def pulse_function(self, time: float) -> float:
        """
        计算特定时间的玄武岩风化通量

        注意: 地质时间是从老到新递减的 (e.g., 300 Ma -> 200 Ma)
        所以 t < t_onset 表示在脉冲开始之后

        Parameters
        ----------
        time : float
            时间 (Ma)

        Returns
        -------
        float
            玄武岩风化通量 (mol/yr)
        """
        # Heaviside阶跃函数 (注意地质时间方向)
        # t_onset = 265.5 Ma, 时间小于265.5表示在脉冲开始之后
        if time <= self.t_onset:
            # 高斯脉冲
            exponent = -((time - self.t_peak) ** 2) / (2 * self.sigma ** 2)
            return self.F_max * np.exp(exponent)
        else:
            return 0.0

    def pulse_array(self, times: np.ndarray) -> np.ndarray:
        """数组版本: 计算多个时间点的通量"""
        return np.array([self.pulse_function(t) for t in times])


@dataclass
class DualComponentRiverConfig:
    """
    双组分河流配置

    F_riv * R_riv = F_crust * R_crust + F_basalt(t) * R_basalt
    """
    # 地壳组分 (背景)
    F_crust: float = 35e9          # 地壳风化通量 (mol/yr)
    R_crust: float = 0.71107       # 地壳同位素比值

    # 玄武岩组分 (脉冲)
    basalt_config: Optional[BasaltPulseConfig] = None

    # 气候增强系数 (Scenario C)
    climate_enhancement: float = 1.0  # 温室效应增强风化 (1.2-1.8)

    def __post_init__(self):
        if self.basalt_config is None:
            self.basalt_config = BasaltPulseConfig()

    def get_total_river_flux(self, time: float) -> float:
        """获取总河流风化通量"""
        F_basalt = self.basalt_config.pulse_function(time) if self.basalt_config else 0.0
        return (self.F_crust + F_basalt) * self.climate_enhancement

    def get_effective_river_ratio(self, time: float) -> float:
        """
        计算有效河流同位素比值

        R_riv_eff = (F_crust * R_crust + F_basalt * R_basalt) / (F_crust + F_basalt)
        """
        F_basalt = self.basalt_config.pulse_function(time) if self.basalt_config else 0.0
        F_total = self.F_crust + F_basalt

        if F_total == 0:
            return self.R_crust

        weighted_sum = (
            self.F_crust * self.R_crust +
            F_basalt * self.basalt_config.R_basalt
        )

        return weighted_sum / F_total

    def get_fluxes_at_time(self, time: float) -> dict:
        """获取特定时间的详细通量信息"""
        F_basalt = self.basalt_config.pulse_function(time) if self.basalt_config else 0.0
        F_total = (self.F_crust + F_basalt) * self.climate_enhancement
        R_eff = self.get_effective_river_ratio(time)

        return {
            'time_ma': time,
            'F_crust': self.F_crust * self.climate_enhancement,
            'F_basalt': F_basalt * self.climate_enhancement,
            'F_total': F_total,
            'R_crust': self.R_crust,
            'R_basalt': self.basalt_config.R_basalt if self.basalt_config else None,
            'R_river_effective': R_eff,
            'f_basalt_fraction': F_basalt / (self.F_crust + F_basalt) if (self.F_crust + F_basalt) > 0 else 0,
            'climate_enhancement': self.climate_enhancement
        }


def calculate_seawater_sr_dual(
    river_config: DualComponentRiverConfig,
    time: float,
    F_highT: float = 8.04e9,
    R_highT: float = 0.7037,
    F_lowT: float = 10e9,
    R_lowT: float = 0.7084,
    F_dia: float = 3.4e9,
    R_dia: float = 0.7084
) -> float:
    """
    使用双组分河流模型计算海水Sr同位素

    Parameters
    ----------
    river_config : DualComponentRiverConfig
        双组分河流配置
    time : float
        时间 (Ma)
    F_highT, R_highT : float
        高温热液通量和同位素
    F_lowT, R_lowT : float
        低温热液通量和同位素
    F_dia, R_dia : float
        成岩通量和同位素

    Returns
    -------
    float
        海水87Sr/86Sr比值
    """
    # 获取河流输入
    F_river = river_config.get_total_river_flux(time)
    R_river = river_config.get_effective_river_ratio(time)

    # 质量平衡计算
    total_flux = F_river + F_highT + F_lowT + F_dia

    if total_flux == 0:
        return 0.708  # 默认值

    weighted_sum = (
        F_river * R_river +
        F_highT * R_highT +
        F_lowT * R_lowT +
        F_dia * R_dia
    )

    return weighted_sum / total_flux


class TimeDependentSrModel:
    """
    时间依赖的Sr同位素模型

    支持随时间变化的参数
    """

    def __init__(
        self,
        river_config: Optional[DualComponentRiverConfig] = None,
        time_span: tuple = (299, 252)
    ):
        """
        初始化时间依赖模型

        Parameters
        ----------
        river_config : DualComponentRiverConfig
            河流配置
        time_span : tuple
            模拟时间范围 (Ma)
        """
        self.river_config = river_config or DualComponentRiverConfig()
        self.time_span = time_span

    def simulate_time_series(
        self,
        times: np.ndarray,
        F_highT_func: Optional[Callable[[float], float]] = None,
        F_lowT_func: Optional[Callable[[float], float]] = None,
    ) -> np.ndarray:
        """
        模拟时间序列

        Parameters
        ----------
        times : array_like
            时间点数组 (Ma)
        F_highT_func : callable, optional
            高温热液通量随时间变化的函数
        F_lowT_func : callable, optional
            低温热液通量随时间变化的函数

        Returns
        -------
        array
            Sr同位素比值数组
        """
        sr_ratios = np.zeros_like(times)

        for i, t in enumerate(times):
            # 动态热液通量
            F_highT = F_highT_func(t) if F_highT_func else 8.04e9
            F_lowT = F_lowT_func(t) if F_lowT_func else 10e9

            sr_ratios[i] = calculate_seawater_sr_dual(
                self.river_config,
                t,
                F_highT=F_highT,
                F_lowT=F_lowT
            )

        return sr_ratios

    def get_pulse_visualization(self, times: np.ndarray) -> dict:
        """
        获取玄武岩脉冲可视化数据

        Returns
        -------
        dict
            包含时间序列的各种通量和比值
        """
        data = {
            'times': times,
            'F_basalt': [],
            'F_crust': [],
            'F_total_river': [],
            'R_river': [],
            'R_sw': []
        }

        for t in times:
            fluxes = self.river_config.get_fluxes_at_time(t)
            data['F_basalt'].append(fluxes['F_basalt'])
            data['F_crust'].append(fluxes['F_crust'])
            data['F_total_river'].append(fluxes['F_total'])
            data['R_river'].append(fluxes['R_river_effective'])

            # 计算海水Sr
            R_sw = calculate_seawater_sr_dual(self.river_config, t)
            data['R_sw'].append(R_sw)

        # 转换为数组
        for key in ['F_basalt', 'F_crust', 'F_total_river', 'R_river', 'R_sw']:
            data[key] = np.array(data[key])

        return data


if __name__ == '__main__':
    # 测试双组分模型
    print("="*60)
    print("双组分河流模型测试")
    print("="*60)

    # 测试1: 无玄武岩脉冲 (背景状态)
    print("\n[1] 背景状态 (无玄武岩风化):")
    config_bg = DualComponentRiverConfig(
        F_crust=35e9,
        R_crust=0.71107,
        basalt_config=None,
        climate_enhancement=1.0
    )

    for t in [268, 265, 263, 261, 259]:
        R_sw = calculate_seawater_sr_dual(config_bg, t)
        print(f"  t={t} Ma: R_sw = {R_sw:.5f}")

    # 测试2: 有玄武岩脉冲 (Scenario C)
    print("\n[2] Scenario C (峨眉山玄武岩风化脉冲):")
    basalt_config = BasaltPulseConfig(
        F_max=25e9,
        t_onset=265.5,
        t_peak=263.5,
        sigma=1.5,
        R_basalt=0.7040
    )

    config_scenC = DualComponentRiverConfig(
        F_crust=35e9,
        R_crust=0.71107,
        basalt_config=basalt_config,
        climate_enhancement=1.5
    )

    print("\n  时间序列:")
    times = np.arange(266, 258, -0.5)
    for t in times:
        fluxes = config_scenC.get_fluxes_at_time(t)
        R_sw = calculate_seawater_sr_dual(config_scenC, t)
        print(f"  t={t:5.1f} Ma: F_basalt={fluxes['F_basalt']/1e9:5.1f}×10⁹, "
              f"R_riv={fluxes['R_river_effective']:.5f}, R_sw={R_sw:.5f}")

    # 测试3: 时间序列模拟
    print("\n[3] 时间序列模拟:")
    model = TimeDependentSrModel(config_scenC)
    times = np.linspace(266, 258, 100)
    viz_data = model.get_pulse_visualization(times)

    print(f"  模拟时间点: {len(times)}个")
    print(f"  R_sw范围: {viz_data['R_sw'].min():.5f} - {viz_data['R_sw'].max():.5f}")
    print(f"  F_basalt峰值: {viz_data['F_basalt'].max()/1e9:.1f}×10⁹ mol/yr")

    # 找出最低Sr值
    min_idx = np.argmin(viz_data['R_sw'])
    print(f"  最低R_sw: {viz_data['R_sw'][min_idx]:.5f} at t={times[min_idx]:.1f} Ma")
