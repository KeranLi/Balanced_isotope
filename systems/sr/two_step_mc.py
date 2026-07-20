"""
Sr同位素两步蒙特卡洛模拟模块

基于Wang et al. (2021) Earth-Science Reviews的方法:
1. Step 1: 5000次均匀随机采样，筛选匹配目标曲线
2. Step 2: 基于Step 1成功参数的高斯分布，2000次细化采样

参考:
Wang et al. (2021) "Revisiting the Permian seawater Sr/86Sr record:
New perspectives from brachiopod proxy data and stochastic oceanic box models"
Earth-Science Reviews 218: 103679
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Tuple, List, Optional, Callable
from scipy import stats
from scipy.interpolate import interp1d

from systems.sr.scenarios import BaseScenario


@dataclass
class Step1Result:
    """Step 1 结果"""
    successful_params: List[Dict] = field(default_factory=list)
    successful_curves: List[np.ndarray] = field(default_factory=list)
    ages: np.ndarray = None
    target_sr: np.ndarray = None
    n_total_runs: int = 0

    def get_param_statistics(self) -> Dict:
        """获取成功参数的统计信息"""
        if not self.successful_params:
            return {}

        stats_dict = {}
        param_names = self.successful_params[0].keys()

        for param_name in param_names:
            values = np.array([p[param_name] for p in self.successful_params])
            stats_dict[param_name] = {
                'mean': float(np.mean(values)),
                'std': float(np.std(values)),
                'median': float(np.median(values)),
                'min': float(np.min(values)),
                'max': float(np.max(values)),
                'percentile_16': float(np.percentile(values, 16)),
                'percentile_84': float(np.percentile(values, 84)),
            }

        return stats_dict

    def get_success_rate(self) -> float:
        """获取成功率"""
        if self.n_total_runs == 0:
            return 0.0
        return len(self.successful_params) / self.n_total_runs


@dataclass
class Step2Result:
    """Step 2 结果"""
    successful_params: List[Dict] = field(default_factory=list)
    successful_curves: List[np.ndarray] = field(default_factory=list)
    ages: np.ndarray = None
    target_sr: np.ndarray = None
    step1_statistics: Dict = None
    n_total_runs: int = 0

    def get_refined_statistics(self) -> Dict:
        """获取细化后的统计信息"""
        if not self.successful_params:
            return {}

        stats_dict = {}
        param_names = self.successful_params[0].keys()

        for param_name in param_names:
            values = np.array([p[param_name] for p in self.successful_params])
            stats_dict[param_name] = {
                'mean': float(np.mean(values)),
                'std': float(np.std(values)),
                'median': float(np.median(values)),
                'ci_95_lower': float(np.percentile(values, 2.5)),
                'ci_95_upper': float(np.percentile(values, 97.5)),
            }

        return stats_dict


@dataclass
class TwoStepMCResult:
    """完整的两步蒙特卡洛结果"""
    step1: Step1Result = None
    step2: Step2Result = None

    # 最终预测
    final_mean_sr: np.ndarray = None
    final_std_sr: np.ndarray = None
    final_ci_lower: np.ndarray = None
    final_ci_upper: np.ndarray = None

    def summary(self) -> str:
        """生成结果摘要"""
        lines = []
        lines.append("="*70)
        lines.append("两步蒙特卡洛模拟结果摘要")
        lines.append("="*70)

        if self.step1:
            lines.append(f"\nStep 1 (均匀采样):")
            lines.append(f"  总运行: {self.step1.n_total_runs}")
            lines.append(f"  成功: {len(self.step1.successful_params)}")
            lines.append(f"  成功率: {self.step1.get_success_rate()*100:.2f}%")

        if self.step2:
            lines.append(f"\nStep 2 (高斯细化):")
            lines.append(f"  总运行: {self.step2.n_total_runs}")
            lines.append(f"  成功: {len(self.step2.successful_params)}")
            lines.append(f"  成功率: {len(self.step2.successful_params)/self.step2.n_total_runs*100:.2f}%")

        if self.final_mean_sr is not None:
            lines.append(f"\n最终预测:")
            lines.append(f"  Sr范围: {self.final_mean_sr.min():.5f} - {self.final_mean_sr.max():.5f}")

        return "\n".join(lines)


class TwoStepMonteCarlo:
    """
    两步蒙特卡洛模拟器

    实现Wang et al. (2021)的两步方法
    """

    def __init__(
        self,
        scenario: BaseScenario,
        target_ages: np.ndarray,
        target_sr: np.ndarray,
        target_ci_lower: Optional[np.ndarray] = None,
        target_ci_upper: Optional[np.ndarray] = None,
        tolerance: float = 0.0002,
        adaptive_tolerance: bool = True,
    ):
        """
        初始化两步蒙特卡洛模拟器

        Parameters
        ----------
        scenario : BaseScenario
            情景对象 (A, B, or C)
        target_ages : array_like
            目标曲线的年龄点
        target_sr : array_like
            目标曲线的Sr值
        target_ci_lower : array_like, optional
            目标曲线95% CI下限
        target_ci_upper : array_like, optional
            目标曲线95% CI上限
        tolerance : float
            基础匹配容差 (默认0.0002)
        adaptive_tolerance : bool
            是否使用自适应容差
        """
        self.scenario = scenario
        self.target_ages = target_ages
        self.target_sr = target_sr
        self.base_tolerance = tolerance
        self.adaptive_tolerance = adaptive_tolerance

        # 计算目标Sr的变化幅度
        self.target_sr_range = (float(target_sr.min()), float(target_sr.max()))
        self.target_sr_amplitude = abs(self.target_sr_range[1] - self.target_sr_range[0])

        # 如果没有提供CI，使用基于数据变化的动态容差
        if target_ci_lower is None:
            # 根据数据变化幅度设置默认容差
            dynamic_tol = max(tolerance, self.target_sr_amplitude * 0.1)
            target_ci_lower = target_sr - dynamic_tol
        if target_ci_upper is None:
            dynamic_tol = max(tolerance, self.target_sr_amplitude * 0.1)
            target_ci_upper = target_sr + dynamic_tol

        self.target_ci_lower = target_ci_lower
        self.target_ci_upper = target_ci_upper

        # 创建目标插值函数
        self._target_interp = interp1d(
            target_ages, target_sr,
            kind='linear',
            bounds_error=False,
            fill_value='extrapolate'
        )
        self._ci_lower_interp = interp1d(
            target_ages, target_ci_lower,
            kind='linear',
            bounds_error=False,
            fill_value='extrapolate'
        )
        self._ci_upper_interp = interp1d(
            target_ages, target_ci_upper,
            kind='linear',
            bounds_error=False,
            fill_value='extrapolate'
        )

    def run(
        self,
        n_runs_step1: int = 5000,
        n_runs_step2: int = 2000,
        tolerance: float = 0.0002,
        verbose: bool = True
    ) -> TwoStepMCResult:
        """
        运行完整的两步蒙特卡洛模拟

        Parameters
        ----------
        n_runs_step1 : int
            Step 1运行次数 (默认5000)
        n_runs_step2 : int
            Step 2运行次数 (默认2000)
        tolerance : float
            匹配容差
        verbose : bool
            是否输出进度

        Returns
        -------
        TwoStepMCResult
            完整结果
        """
        result = TwoStepMCResult()

        # Step 1: 均匀随机采样
        if verbose:
            print("="*70)
            print("Step 1: 均匀随机采样")
            print("="*70)

        result.step1 = self._run_step1(
            n_runs=n_runs_step1,
            tolerance=tolerance,
            verbose=verbose
        )

        if not result.step1.successful_params:
            if verbose:
                print("Step 1未找到成功参数，调整容差重试...")
            # 扩大容差重试
            result.step1 = self._run_step1(
                n_runs=n_runs_step1,
                tolerance=tolerance*2,
                verbose=verbose
            )

        if not result.step1.successful_params:
            if verbose:
                print("警告: Step 1未找到任何成功参数")
            return result

        # Step 2: 高斯细化
        if verbose:
            print("\n" + "="*70)
            print("Step 2: 高斯细化采样")
            print("="*70)

        result.step2 = self._run_step2(
            step1_stats=result.step1.get_param_statistics(),
            n_runs=n_runs_step2,
            tolerance=tolerance,
            verbose=verbose
        )

        # 计算最终结果
        self._compute_final_results(result)

        if verbose:
            print("\n" + result.summary())

        return result

    def _run_step1(
        self,
        n_runs: int,
        tolerance: float,
        verbose: bool
    ) -> Step1Result:
        """运行Step 1: 均匀随机采样"""
        result = Step1Result(
            ages=self.target_ages,
            target_sr=self.target_sr,
            n_total_runs=n_runs
        )

        param_ranges = self.scenario.get_param_ranges()

        if verbose:
            print(f"参数范围:")
            for name, (min_v, max_v) in param_ranges.items():
                print(f"  {name}: {min_v} - {max_v}")
            print(f"\n运行 {n_runs} 次模拟...")

        for i in range(n_runs):
            if verbose and (i + 1) % 1000 == 0:
                print(f"  完成 {i+1}/{n_runs}... (成功: {len(result.successful_params)})")

            # 生成随机参数
            params = self._generate_uniform_params(param_ranges)

            # 计算Sr曲线
            sr_curve = np.array([
                self.scenario.calculate_sr(params, age)
                for age in self.target_ages
            ])

            # 检查是否匹配目标
            if self._check_match(sr_curve, tolerance):
                result.successful_params.append(params)
                result.successful_curves.append(sr_curve)

        if verbose:
            print(f"\nStep 1完成: {len(result.successful_params)}/{n_runs} 成功 "
                  f"({len(result.successful_params)/n_runs*100:.2f}%)")

        return result

    def _run_step2(
        self,
        step1_stats: Dict,
        n_runs: int,
        tolerance: float,
        verbose: bool
    ) -> Step2Result:
        """运行Step 2: 高斯细化采样"""
        result = Step2Result(
            ages=self.target_ages,
            target_sr=self.target_sr,
            step1_statistics=step1_stats,
            n_total_runs=n_runs
        )

        if not step1_stats:
            if verbose:
                print("Step 1统计信息为空，跳过Step 2")
            return result

        if verbose:
            print(f"基于Step 1的高斯分布采样...")
            print(f"运行 {n_runs} 次模拟...")

        for i in range(n_runs):
            if verbose and (i + 1) % 500 == 0:
                print(f"  完成 {i+1}/{n_runs}... (成功: {len(result.successful_params)})")

            # 从高斯分布生成参数
            params = self._generate_gaussian_params(step1_stats)

            # 计算Sr曲线
            sr_curve = np.array([
                self.scenario.calculate_sr(params, age)
                for age in self.target_ages
            ])

            # 检查是否匹配目标
            if self._check_match(sr_curve, tolerance):
                result.successful_params.append(params)
                result.successful_curves.append(sr_curve)

        if verbose:
            success_rate = len(result.successful_params) / n_runs * 100
            print(f"\nStep 2完成: {len(result.successful_params)}/{n_runs} 成功 "
                  f"({success_rate:.2f}%)")

        return result

    def _generate_uniform_params(
        self,
        param_ranges: Dict[str, Tuple[float, float]]
    ) -> Dict:
        """生成均匀分布的随机参数"""
        params = {}
        for param_name, (min_val, max_val) in param_ranges.items():
            if param_name.startswith('R_'):
                # 同位素比值用均匀分布
                params[param_name] = np.random.uniform(min_val, max_val)
            else:
                # 通量用对数均匀分布
                log_min, log_max = np.log(min_val), np.log(max_val)
                params[param_name] = np.exp(np.random.uniform(log_min, log_max))
        return params

    def _generate_gaussian_params(self, step1_stats: Dict) -> Dict:
        """基于Step 1统计信息生成高斯分布的参数"""
        params = {}
        for param_name, stats in step1_stats.items():
            mean = stats['mean']
            std = stats['std']

            # 从高斯分布采样
            value = np.random.normal(mean, std)

            # 确保在合理范围内
            if param_name.startswith('F_'):
                value = max(value, 0)  # 通量不能为负

            params[param_name] = value

        return params

    def _check_match(
        self,
        sr_curve: np.ndarray,
        tolerance: float
    ) -> bool:
        """
        检查Sr曲线是否匹配观测数据目标

        使用多种指标综合判断:
        1. 是否在CI范围内（主要标准）
        2. MSE是否足够小
        3. 曲线范围是否重叠

        注：文献约束仅用于参考和验证，不直接参与匹配
        """
        # 自适应容差
        if self.adaptive_tolerance:
            dynamic_tolerance = max(tolerance, self.target_sr_amplitude * 0.15)
        else:
            dynamic_tolerance = tolerance

        # ===== 方法1: 检查是否在CI范围内 =====
        in_ci = (
            (sr_curve >= self.target_ci_lower - dynamic_tolerance) &
            (sr_curve <= self.target_ci_upper + dynamic_tolerance)
        )

        ci_match_ratio = np.mean(in_ci)
        if ci_match_ratio >= 0.6:
            return True

        # ===== 方法2: 检查MSE =====
        mse = np.mean((sr_curve - self.target_sr) ** 2)
        max_mse = (dynamic_tolerance * 2) ** 2
        if mse < max_mse:
            return True

        # ===== 方法3: 检查曲线范围 =====
        sr_min, sr_max = sr_curve.min(), sr_curve.max()
        target_min, target_max = self.target_sr_range

        range_overlap = (
            sr_min <= target_max + dynamic_tolerance and
            sr_max >= target_min - dynamic_tolerance
        )

        if range_overlap and ci_match_ratio >= 0.4:
            return True

        return False

    def _compute_final_results(self, result: TwoStepMCResult):
        """计算最终结果统计"""
        # 优先使用Step 2的结果
        if result.step2 and result.step2.successful_curves:
            curves = result.step2.successful_curves
        elif result.step1 and result.step1.successful_curves:
            curves = result.step1.successful_curves
        else:
            return

        curves_array = np.array(curves)

        result.final_mean_sr = np.mean(curves_array, axis=0)
        result.final_std_sr = np.std(curves_array, axis=0)
        result.final_ci_lower = np.percentile(curves_array, 2.5, axis=0)
        result.final_ci_upper = np.percentile(curves_array, 97.5, axis=0)


def calculate_optimal_tolerance(target_sr: np.ndarray) -> float:
    """
    基于目标数据特征计算最优容差

    Parameters
    ----------
    target_sr : array_like
        目标Sr值

    Returns
    -------
    float
        建议的容差值
    """
    # 计算数据特征
    sr_std = np.std(target_sr)
    sr_range = np.max(target_sr) - np.min(target_sr)
    n_points = len(target_sr)

    # 基于标准差和范围的容差
    tolerance_from_std = sr_std * 2.5  # 2.5倍标准差覆盖约99%的数据
    tolerance_from_range = sr_range * 0.5  # 范围的一半

    # 取较大值，但设置上下限
    optimal_tolerance = np.clip(
        max(tolerance_from_std, tolerance_from_range * 0.3),
        0.0001,  # 最小容差: 1 ppm
        0.0010   # 最大容差: 10 ppm
    )

    # 如果数据点少，适当增加容差
    if n_points < 10:
        optimal_tolerance *= 1.5

    return float(optimal_tolerance)


def run_two_step_mc_for_scenario(
    scenario_name: str,
    target_ages: np.ndarray,
    target_sr: np.ndarray,
    n_runs_step1: int = 5000,
    n_runs_step2: int = 2000,
    tolerance: Optional[float] = None,
    verbose: bool = True
) -> TwoStepMCResult:
    """
    便捷函数: 为特定情景运行两步蒙特卡洛

    Parameters
    ----------
    scenario_name : str
        'A', 'B', 或 'C'
    target_ages : array_like
        目标年龄点
    target_sr : array_like
        目标Sr值
    n_runs_step1 : int
    n_runs_step2 : int
    tolerance : float, optional
        匹配容差，默认自动计算
    verbose : bool

    Returns
    -------
    TwoStepMCResult
    """
    from systems.sr.scenarios import ScenarioManager

    manager = ScenarioManager()
    scenario = manager.get_scenario(scenario_name)

    # 自动计算最优容差
    if tolerance is None:
        tolerance = calculate_optimal_tolerance(target_sr)
        if verbose:
            print(f"  自动计算容差: {tolerance:.5f} ({tolerance*10000:.1f} ppm)")

    # 创建模拟器
    mc = TwoStepMonteCarlo(
        scenario=scenario,
        target_ages=target_ages,
        target_sr=target_sr,
        tolerance=tolerance,
        adaptive_tolerance=True
    )

    # 运行
    return mc.run(
        n_runs_step1=n_runs_step1,
        n_runs_step2=n_runs_step2,
        verbose=verbose
    )


if __name__ == '__main__':
    # 测试两步蒙特卡洛
    import sys
    sys.path.insert(0, '.')

    print("="*70)
    print("两步蒙特卡洛模拟测试")
    print("="*70)

    from systems.sr.scenarios import ScenarioManager

    # 使用Scenario A进行测试
    manager = ScenarioManager()
    scenario = manager.get_scenario('A')

    # 创建简单的目标曲线
    ages = np.linspace(268, 273, 20)
    # Roadian下降: 0.70762 -> 0.70701
    target = np.linspace(0.70701, 0.70762, 20)

    print(f"\n测试情景: {scenario.name}")
    print(f"年龄范围: {ages.min():.1f} - {ages.max():.1f} Ma")
    print(f"目标Sr范围: {target.min():.5f} - {target.max():.5f}")

    # 运行两步MC (减少运行次数用于测试)
    mc = TwoStepMonteCarlo(scenario, ages, target)
    result = mc.run(n_runs_step1=500, n_runs_step2=200, verbose=True)

    print("\n" + result.summary())
