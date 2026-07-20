"""
Sr同位素可视化模块

提供以下图表:
1. Sr演化曲线图 (Seawater Sr evolution)
2. 参数密度热图 (Parameter density heatmaps)
3. 多情景对比图 (Multi-scenario comparison)
4. 蒙特卡洛结果展示
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde
from typing import Dict, List, Optional, Tuple
import warnings

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


class SrVisualization:
    """Sr同位素可视化类"""

    def __init__(self, style: str = 'default'):
        """
        初始化可视化器

        Parameters
        ----------
        style : str
            matplotlib样式
        """
        self.style = style
        plt.style.use(style)

        # 定义颜色方案
        self.colors = {
            'scenario_a': '#1f77b4',  # 蓝色
            'scenario_b': '#ff7f0e',  # 橙色
            'scenario_c': '#2ca02c',  # 绿色
            'observed': '#d62728',     # 红色
            'target': '#9467bd',       # 紫色
            'model_mean': '#8c564b',   # 棕色
            'model_ci': '#e377c2',     # 粉色
        }

    def plot_sr_evolution(
        self,
        ages: np.ndarray,
        observed_sr: np.ndarray,
        model_mean: Optional[np.ndarray] = None,
        model_std: Optional[np.ndarray] = None,
        model_ci_lower: Optional[np.ndarray] = None,
        model_ci_upper: Optional[np.ndarray] = None,
        title: str = "Seawater ⁸⁷Sr/⁸⁶Sr Evolution",
        xlabel: str = "Age (Ma)",
        ylabel: str = "⁸⁷Sr/⁸⁶Sr",
        figsize: Tuple[int, int] = (12, 6),
        save_path: Optional[str] = None,
        show: bool = True
    ) -> plt.Figure:
        """
        绘制Sr同位素演化曲线

        Parameters
        ----------
        ages : array_like
            年龄数组
        observed_sr : array_like
            观测的Sr值
        model_mean : array_like, optional
            模型预测的均值
        model_std : array_like, optional
            模型预测的标准差
        model_ci_lower, model_ci_upper : array_like, optional
            置信区间
        title : str
            图表标题
        save_path : str, optional
            保存路径
        show : bool
            是否显示图表

        Returns
        -------
        matplotlib.figure.Figure
        """
        fig, ax = plt.subplots(figsize=figsize)

        # 绘制观测数据
        ax.scatter(ages, observed_sr, c=self.colors['observed'],
                  s=50, alpha=0.6, label='Observed', zorder=3)

        # 绘制模型结果
        if model_mean is not None:
            # 绘制均值曲线
            ax.plot(ages, model_mean, color=self.colors['model_mean'],
                   linewidth=2, label='Model Mean', zorder=2)

            # 绘制置信区间
            if model_ci_lower is not None and model_ci_upper is not None:
                ax.fill_between(ages, model_ci_lower, model_ci_upper,
                              alpha=0.3, color=self.colors['model_ci'],
                              label='95% CI', zorder=1)
            elif model_std is not None:
                ax.fill_between(ages, model_mean - 2*model_std,
                              model_mean + 2*model_std,
                              alpha=0.3, color=self.colors['model_ci'],
                              label='±2σ', zorder=1)

        # 标注三个下降期
        self._annotate_intervals(ax)

        # 设置标签和标题
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')

        # 反转x轴（地质时间方向）
        ax.invert_xaxis()

        # 添加图例
        ax.legend(loc='best', fontsize=10)

        # 添加网格
        ax.grid(True, alpha=0.3, linestyle='--')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"图表已保存: {save_path}")

        if show:
            plt.show()

        return fig

    def plot_scenario_comparison(
        self,
        scenario_results: Dict[str, Dict],
        figsize: Tuple[int, int] = (14, 8),
        save_path: Optional[str] = None,
        show: bool = True
    ) -> plt.Figure:
        """
        绘制三个情景的对比图

        Parameters
        ----------
        scenario_results : dict
            {scenario_name: {'ages': [], 'mean_sr': [], 'ci_lower': [], 'ci_upper': []}}
        figsize : tuple
        save_path : str, optional
        show : bool

        Returns
        -------
        matplotlib.figure.Figure
        """
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle('Middle Permian Seawater ⁸⁷Sr/⁸⁶Sr: Three Declines',
                    fontsize=14, fontweight='bold')

        scenario_names = ['A', 'B', 'C']
        scenario_colors = [self.colors['scenario_a'],
                          self.colors['scenario_b'],
                          self.colors['scenario_c']]

        # 绘制每个情景
        for idx, (name, color) in enumerate(zip(scenario_names, scenario_colors)):
            if name in scenario_results:
                ax = axes[idx // 2, idx % 2]
                result = scenario_results[name]

                ages = result['ages']
                mean_sr = result['mean_sr']
                ci_lower = result.get('ci_lower')
                ci_upper = result.get('ci_upper')

                # 绘制模型结果
                ax.plot(ages, mean_sr, color=color, linewidth=2,
                       label=f'Scenario {name} Model')

                if ci_lower is not None and ci_upper is not None:
                    ax.fill_between(ages, ci_lower, ci_upper,
                                  alpha=0.3, color=color, label='95% CI')

                # 绘制原始实测数据点（实心圆点）
                if 'observed_ages_raw' in result and 'observed_sr_raw' in result:
                    ax.scatter(result['observed_ages_raw'], result['observed_sr_raw'],
                             c=self.colors['observed'], s=50, alpha=0.9,
                             edgecolors='darkred', linewidths=0.5,
                             label='Observed (measured)', zorder=5)
                # 备用：使用插值后的观测数据
                elif 'observed_sr' in result:
                    ax.scatter(ages, result['observed_sr'],
                             c=self.colors['observed'], s=40, alpha=0.7,
                             label='Observed', zorder=5)

                ax.set_xlabel('Age (Ma)', fontsize=10)
                ax.set_ylabel('⁸⁷Sr/⁸⁶Sr', fontsize=10)
                ax.set_title(f'Scenario {name}: {self._get_scenario_title(name)}',
                           fontsize=11)
                ax.invert_xaxis()
                ax.legend(fontsize=8, loc='best')
                ax.grid(True, alpha=0.3)

        # 第四个图：综合对比（包含实测点）
        ax = axes[1, 1]
        for name, color in zip(scenario_names, scenario_colors):
            if name in scenario_results:
                result = scenario_results[name]
                # 绘制模型曲线（实线）
                ax.plot(result['ages'], result['mean_sr'],
                       color=color, linewidth=2.5, label=f'Scenario {name} Model')
                # 绘制原始实测点（散点）
                if 'observed_ages_raw' in result and 'observed_sr_raw' in result:
                    ax.scatter(result['observed_ages_raw'], result['observed_sr_raw'],
                             c=color, s=30, alpha=0.6, edgecolors='black', linewidths=0.3,
                             label=f'Scenario {name} Observed', zorder=5)

        ax.set_xlabel('Age (Ma)', fontsize=10)
        ax.set_ylabel('⁸⁷Sr/⁸⁶Sr', fontsize=10)
        ax.set_title('All Scenarios Comparison (Model vs Observed)', fontsize=11)
        ax.invert_xaxis()
        ax.legend(fontsize=8, loc='best')
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"图表已保存: {save_path}")

        if show:
            plt.show()

        return fig

    def plot_parameter_density(
        self,
        params_list: List[Dict],
        param_names: Optional[List[str]] = None,
        figsize: Tuple[int, int] = (14, 10),
        save_path: Optional[str] = None,
        show: bool = True
    ) -> plt.Figure:
        """
        绘制参数密度热图（Wang et al. 2021风格）

        Parameters
        ----------
        params_list : list of dict
            成功的参数组合列表
        param_names : list, optional
            要绘制的参数名，默认自动选择
        figsize : tuple
        save_path : str, optional
        show : bool

        Returns
        -------
        matplotlib.figure.Figure
        """
        if not params_list:
            print("警告: 没有参数数据可供绘图")
            return None

        # 自动选择关键参数
        if param_names is None:
            all_params = list(params_list[0].keys())
            # 优先选择通量参数
            param_names = [p for p in all_params if p.startswith('F_')][:4]
            # 添加同位素参数
            param_names += [p for p in all_params if p.startswith('R_')][:2]

        n_params = len(param_names)

        # 创建子图
        fig, axes = plt.subplots(n_params, n_params, figsize=figsize)
        fig.suptitle('Parameter Density Distributions (Wang et al. 2021 style)',
                    fontsize=14, fontweight='bold')

        # 准备数据
        data = {name: np.array([p[name] for p in params_list])
                for name in param_names}

        for i, name_i in enumerate(param_names):
            for j, name_j in enumerate(param_names):
                ax = axes[i, j]

                if i == j:
                    # 对角线：绘制1D密度图
                    self._plot_1d_density(ax, data[name_i], name_i)
                elif i > j:
                    # 下三角：绘制2D密度热图
                    self._plot_2d_density(ax, data[name_j], data[name_i],
                                        name_j, name_i)
                else:
                    # 上三角：隐藏
                    ax.axis('off')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"图表已保存: {save_path}")

        if show:
            plt.show()

        return fig

    def plot_monte_carlo_results(
        self,
        ages: np.ndarray,
        all_curves: List[np.ndarray],
        target_sr: Optional[np.ndarray] = None,
        title: str = "Monte Carlo Simulation Results",
        figsize: Tuple[int, int] = (12, 7),
        save_path: Optional[str] = None,
        show: bool = True
    ) -> plt.Figure:
        """
        绘制蒙特卡洛结果展示图

        Parameters
        ----------
        ages : array_like
            年龄数组
        all_curves : list of array
            所有成功的Sr曲线
        target_sr : array_like, optional
            目标曲线
        title : str
        figsize : tuple
        save_path : str, optional
        show : bool

        Returns
        -------
        matplotlib.figure.Figure
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize,
                                       gridspec_kw={'height_ratios': [3, 1]})

        # 转换为数组
        curves_array = np.array(all_curves)

        # 计算统计量
        mean_curve = np.mean(curves_array, axis=0)
        std_curve = np.std(curves_array, axis=0)
        p2_5 = np.percentile(curves_array, 2.5, axis=0)
        p97_5 = np.percentile(curves_array, 97.5, axis=0)

        # 上图：所有曲线和统计
        # 绘制所有成功的曲线（透明度低）
        for curve in all_curves[:min(100, len(all_curves))]:  # 限制显示数量
            ax1.plot(ages, curve, color='gray', alpha=0.1, linewidth=0.5)

        # 绘制统计曲线
        ax1.plot(ages, mean_curve, color=self.colors['model_mean'],
                linewidth=2.5, label='Mean')
        ax1.fill_between(ages, p2_5, p97_5, alpha=0.3,
                        color=self.colors['model_ci'], label='95% CI')

        # 绘制目标曲线
        if target_sr is not None:
            ax1.plot(ages, target_sr, color=self.colors['target'],
                    linewidth=2, linestyle='--', label='Target', zorder=5)

        ax1.set_ylabel('⁸⁷Sr/⁸⁶Sr', fontsize=11)
        ax1.set_title(title, fontsize=13, fontweight='bold')
        ax1.legend(loc='best')
        ax1.invert_xaxis()
        ax1.grid(True, alpha=0.3)

        # 下图：标准差
        ax2.fill_between(ages, 0, std_curve, alpha=0.5,
                        color=self.colors['model_ci'])
        ax2.plot(ages, std_curve, color=self.colors['model_mean'], linewidth=1.5)
        ax2.set_xlabel('Age (Ma)', fontsize=11)
        ax2.set_ylabel('Std Dev', fontsize=11)
        ax2.invert_xaxis()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"图表已保存: {save_path}")

        if show:
            plt.show()

        return fig

    def _plot_1d_density(self, ax, data: np.ndarray, name: str):
        """绘制1D密度图"""
        # 核密度估计
        kde = gaussian_kde(data)
        x_range = np.linspace(data.min(), data.max(), 100)
        density = kde(x_range)

        ax.fill_between(x_range, density, alpha=0.6)
        ax.plot(x_range, density, linewidth=1.5)
        ax.set_ylabel('Density')

        # 格式化标签
        label = self._format_param_name(name)
        ax.set_xlabel(label, fontsize=9)

        # 添加统计信息
        mean_val = np.mean(data)
        std_val = np.std(data)
        ax.axvline(mean_val, color='red', linestyle='--', linewidth=1.5,
                  label=f'μ={mean_val:.3f}')
        ax.legend(fontsize=7)

    def _plot_2d_density(self, ax, x_data: np.ndarray, y_data: np.ndarray,
                         x_name: str, y_name: str):
        """绘制2D密度热图"""
        # 创建2D直方图
        hist, x_edges, y_edges = np.histogram2d(x_data, y_data, bins=30)

        # 绘制热图
        extent = [x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]]
        im = ax.imshow(hist.T, origin='lower', extent=extent,
                      aspect='auto', cmap='YlOrRd', interpolation='bilinear')

        # 设置标签
        if ax.get_subplotspec().is_last_row():
            ax.set_xlabel(self._format_param_name(x_name), fontsize=9)
        if ax.get_subplotspec().is_first_col():
            ax.set_ylabel(self._format_param_name(y_name), fontsize=9)

    def _annotate_intervals(self, ax):
        """标注三个下降期"""
        intervals = [
            (273, 268, 'Decline 1\n(Roadian)'),
            (268, 265, 'Decline 2\n(Wordian)'),
            (265, 259.1, 'Decline 3\n(Capitanian)'),
        ]

        colors_intervals = ['#1f77b4', '#ff7f0e', '#2ca02c']

        for (age_start, age_end, label), color in zip(intervals, colors_intervals):
            ax.axvspan(age_end, age_start, alpha=0.1, color=color)
            ax.text((age_start + age_end) / 2, ax.get_ylim()[1], label,
                   ha='center', va='bottom', fontsize=8, color=color,
                   fontweight='bold', alpha=0.8)

    def _get_scenario_title(self, name: str) -> str:
        """获取情景标题"""
        titles = {
            'A': 'Roadian Glaciation',
            'B': 'Wordian Rifting',
            'C': 'Capitanian LIP',
        }
        return titles.get(name, name)

    def _format_param_name(self, name: str) -> str:
        """格式化参数名"""
        # 替换下划线并格式化
        name_map = {
            'F_river': 'F_river (10⁹ mol/yr)',
            'R_river': 'R_river',
            'F_highT': 'F_highT (10⁹ mol/yr)',
            'F_lowT': 'F_lowT (10⁹ mol/yr)',
            'F_crust': 'F_crust (10⁹ mol/yr)',
            'F_basalt_max': 'F_basalt^max (10⁹ mol/yr)',
            'R_basalt': 'R_basalt',
            'climate_enhancement': 'Climate Enhancement',
        }
        return name_map.get(name, name.replace('_', ' '))


def plot_dukou_data(
    dukou_data,
    show_target_curve: bool = True,
    figsize: Tuple[int, int] = (14, 8),
    save_path: Optional[str] = None,
    show: bool = True
):
    """
    便捷函数: 绘制Dukou剖面数据

    Parameters
    ----------
    dukou_data : DukouData
        Dukou数据对象
    show_target_curve : bool
        是否显示平滑曲线
    figsize : tuple
    save_path : str, optional
    show : bool
    """
    viz = SrVisualization()

    # 准备数据 - 按年龄排序以确保一致性
    sort_idx = np.argsort(dukou_data.ages)
    ages = dukou_data.ages[sort_idx]
    sr_ratios = dukou_data.sr_ratios[sort_idx]

    # 平滑曲线数据（已经是排序的）
    target_sr = dukou_data.target_sr if show_target_curve else None
    target_ages = dukou_data.target_ages if show_target_curve else None

    ci_lower = dukou_data.target_ci_lower if show_target_curve else None
    ci_upper = dukou_data.target_ci_upper if show_target_curve else None

    return viz.plot_sr_evolution(
        ages=ages,  # 排序后的原始数据年龄
        observed_sr=sr_ratios,  # 排序后的观测值
        model_mean=target_sr,  # 平滑曲线
        model_ci_lower=ci_lower,
        model_ci_upper=ci_upper,
        title="Dukou Section: Seawater ⁸⁷Sr/⁸⁶Sr Evolution",
        figsize=figsize,
        save_path=save_path,
        show=show
    )


if __name__ == '__main__':
    # 测试可视化
    import sys
    sys.path.insert(0, '.')

    print("="*60)
    print("Sr同位素可视化模块测试")
    print("="*60)

    # 测试Dukou数据可视化
    from toolkit.io.dukou_data import load_dukou_data

    print("\n1. 加载Dukou数据...")
    data = load_dukou_data()

    print("\n2. 绘制Sr演化曲线...")
    fig = plot_dukou_data(data, show=False,
                         save_path='results/dukou_sr_evolution.png')

    print("\n3. 测试完成！")
    print("   输出文件: results/dukou_sr_evolution.png")
