#!/usr/bin/env python3
"""
中二叠世Sr同位素三步情景分析完整流程

运行流程:
1. 加载Dukou剖面数据
2. 运行三步情景蒙特卡洛模拟
3. 生成可视化图表
4. 输出结果报告

使用方法:
    python run_sr_three_declines.py --all
    python run_sr_three_declines.py --scenario A --visualize
    python run_sr_three_declines.py --help
"""

import sys
import argparse
import numpy as np
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))


def run_full_analysis(args):
    """运行完整分析流程"""
    from toolkit.io.dukou_data import load_dukou_data
    from systems.sr.scenarios import ScenarioManager
    from systems.sr.two_step_mc import TwoStepMonteCarlo, calculate_optimal_tolerance
    from toolkit.visualization.sr_plots import (
        SrVisualization, plot_dukou_data
    )

    print("="*70)
    print("中二叠世Sr同位素三步情景分析")
    print("Middle Permian Seawater ⁸⁷Sr/⁸⁶Sr: Three Consecutive Declines")
    print("="*70)

    # 创建输出目录
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    # Step 1: 加载数据
    print("\n[Step 1] 加载Dukou剖面数据...")
    dukou_data = load_dukou_data(
        mn_sr_threshold=args.mn_sr_threshold,
        lowess_frac=args.lowess_frac
    )

    # Step 2: 绘制原始数据
    if args.visualize or args.plot_data:
        print("\n[Step 2] 绘制Dukou Sr演化曲线...")
        fig = plot_dukou_data(
            dukou_data,
            show_target_curve=True,
            save_path=output_dir / '01_dukou_sr_evolution.png',
            show=False
        )
        print(f"   保存: {output_dir / '01_dukou_sr_evolution.png'}")

    # Step 3: 运行三步情景模拟
    manager = ScenarioManager()
    scenarios_to_run = []

    if args.all or args.scenario_a:
        scenarios_to_run.append('A')
    if args.all or args.scenario_b:
        scenarios_to_run.append('B')
    if args.all or args.scenario_c:
        scenarios_to_run.append('C')

    if not scenarios_to_run:
        scenarios_to_run = ['A', 'B', 'C']  # 默认运行所有

    scenario_results = {}

    for scenario_name in scenarios_to_run:
        scenario = manager.get_scenario(scenario_name)

        print("\n" + "="*70)
        print(f"Scenario {scenario_name}: {scenario.name}")
        print("="*70)
        print(f"描述: {scenario.description}")
        print(f"时间: {scenario.age_range[1]:.1f} - {scenario.age_range[0]:.1f} Ma")
        print(f"目标Sr: {scenario.target_sr_range[0]:.5f} - {scenario.target_sr_range[1]:.5f}")

        # 提取目标数据（使用Stage标记而非年龄范围）
        try:
            interval_data = dukou_data.get_stage_data(scenario_name)
            print(f"数据点: {len(interval_data.sample_ids)}个 (Stage {scenario_name})")
            print(f"数据范围: {interval_data.sr_ratios.min():.5f} - {interval_data.sr_ratios.max():.5f}")
            print(f"年龄范围: {interval_data.ages.min():.1f} - {interval_data.ages.max():.1f} Ma")
        except ValueError as e:
            print(f"警告: {e}，尝试使用年龄范围...")
            interval_data = dukou_data.get_interval_data(
                scenario.age_range[0], scenario.age_range[1]
            )
            if len(interval_data.sample_ids) == 0:
                print(f"警告: 该区间没有数据点，跳过")
                continue
            print(f"数据点: {len(interval_data.sample_ids)}个")
            print(f"数据范围: {interval_data.sr_ratios.min():.5f} - {interval_data.sr_ratios.max():.5f}")

        # 显示文献约束（仅作为参考，不作为目标）
        if hasattr(scenario, 'literature_sr_range'):
            lit_range = scenario.literature_sr_range
            print(f"\n文献约束（参考）: {lit_range[0]:.5f} - {lit_range[1]:.5f}")
            print(f"  期望起始Sr: {scenario.literature_sr_start:.5f}")
            print(f"  期望结束Sr: {scenario.literature_sr_end:.5f}")

        # 设置观测数据的起始和结束Sr值及年龄范围（用于模型计算）
        sr_start_observed = interval_data.sr_ratios[np.argmax(interval_data.ages)]  # 最老样本
        sr_end_observed = interval_data.sr_ratios[np.argmin(interval_data.ages)]   # 最新样本
        age_start_observed = interval_data.ages.max()  # 最老年龄
        age_end_observed = interval_data.ages.min()    # 最新年龄
        scenario.set_observed_range(sr_start_observed, sr_end_observed, age_start_observed, age_end_observed)

        # 设置观测曲线（用于模型插值）
        sort_idx = np.argsort(interval_data.ages)  # 按年龄排序（递增）
        scenario.set_observed_curve(
            interval_data.ages[sort_idx],
            interval_data.sr_ratios[sort_idx]
        )

        print(f"\n观测数据锚定:")
        print(f"  起始Sr (最老 {age_start_observed:.1f} Ma): {sr_start_observed:.5f}")
        print(f"  结束Sr (最新 {age_end_observed:.1f} Ma): {sr_end_observed:.5f}")

        # 使用实际观测数据作为MC目标（从老到新排序）
        # 排序：年龄从大到小（从老到新）
        sort_idx = np.argsort(interval_data.ages)[::-1]
        target_ages_raw = interval_data.ages[sort_idx]
        target_sr_raw = interval_data.sr_ratios[sort_idx]

        print(f"\n观测数据目标（从老到新）:")
        print(f"  起始: {target_ages_raw[0]:.1f} Ma, Sr={target_sr_raw[0]:.5f}")
        print(f"  结束: {target_ages_raw[-1]:.1f} Ma, Sr={target_sr_raw[-1]:.5f}")

        # 增加时间分辨率：将观测数据插值到密集网格
        n_time_points = max(50, len(interval_data.ages) * 3)
        print(f"\n时间分辨率优化:")
        print(f"  原始数据点: {len(interval_data.ages)}")
        print(f"  插值后点数: {n_time_points}")

        # 创建密集时间网格（从老到新）
        target_ages_dense = np.linspace(
            target_ages_raw.max(),  # 老
            target_ages_raw.min(),  # 新
            n_time_points
        )

        # 将观测数据插值到密集网格（保持从老到新顺序）
        # 注意：interp需要x递增，所以需要反转
        target_sr_dense = np.interp(
            target_ages_dense,
            target_ages_raw[::-1],  # 递增
            target_sr_raw[::-1]     # 对应值
        )

        # 运行两步蒙特卡洛
        print(f"\n运行两步蒙特卡洛模拟...")
        print(f"  Step 1: {args.n_runs_step1} runs")
        print(f"  Step 2: {args.n_runs_step2} runs")
        print(f"  容差: {args.tolerance}")

        # 使用已设置观测曲线的scenario对象直接运行MC
        tolerance_mc = args.tolerance if args.tolerance else calculate_optimal_tolerance(target_sr_dense)
        mc = TwoStepMonteCarlo(
            scenario=scenario,
            target_ages=target_ages_dense,
            target_sr=target_sr_dense,
            tolerance=tolerance_mc,
            adaptive_tolerance=True
        )
        result = mc.run(
            n_runs_step1=args.n_runs_step1,
            n_runs_step2=args.n_runs_step2,
            tolerance=tolerance_mc,
            verbose=True
        )

        print("\n" + result.summary())

        # 保存结果
        if result.step2 and result.step2.successful_params:
            # 保存参数
            import pandas as pd
            params_df = pd.DataFrame(result.step2.successful_params)
            params_file = output_dir / f'scenario_{scenario_name}_params.csv'
            params_df.to_csv(params_file, index=False)
            print(f"\n参数保存: {params_file}")

            # 存储结果用于可视化（包含原始实测点和插值数据）
            scenario_results[scenario_name] = {
                'ages': target_ages_dense,  # 插值后的密集年龄
                'mean_sr': result.final_mean_sr,
                'ci_lower': result.final_ci_lower,
                'ci_upper': result.final_ci_upper,
                'observed_sr': target_sr_dense,  # 插值后的观测曲线（用于MC目标）
                'observed_ages_raw': target_ages_raw,  # 原始实测年龄（用于绘图）
                'observed_sr_raw': target_sr_raw,      # 原始实测Sr值（用于绘图）
            }

            # 绘制单个情景结果
            if args.visualize:
                viz = SrVisualization()
                fig = viz.plot_monte_carlo_results(
                    ages=target_ages_dense,
                    all_curves=result.step2.successful_curves,
                    target_sr=target_sr_dense,
                    title=f"Scenario {scenario_name}: {scenario.name}",
                    save_path=output_dir / f'scenario_{scenario_name}_mc_results.png',
                    show=False
                )
                print(f"可视化保存: {output_dir / f'scenario_{scenario_name}_mc_results.png'}")

                # 绘制参数密度图
                if len(result.step2.successful_params) > 10:
                    fig = viz.plot_parameter_density(
                        result.step2.successful_params,
                        save_path=output_dir / f'scenario_{scenario_name}_density.png',
                        show=False
                    )
                    print(f"密度图保存: {output_dir / f'scenario_{scenario_name}_density.png'}")

    # Step 4: 综合对比图
    if args.visualize and len(scenario_results) > 1:
        print("\n[Step 4] 生成综合对比图...")
        viz = SrVisualization()
        fig = viz.plot_scenario_comparison(
            scenario_results,
            save_path=output_dir / '04_all_scenarios_comparison.png',
            show=False
        )
        print(f"保存: {output_dir / '04_all_scenarios_comparison.png'}")

    # 生成报告
    print("\n" + "="*70)
    print("分析完成!")
    print("="*70)
    print(f"\n输出文件目录: {output_dir}")
    print("\n生成的文件:")
    for f in sorted(output_dir.glob('*.png')):
        print(f"  - {f.name}")
    for f in sorted(output_dir.glob('*.csv')):
        print(f"  - {f.name}")

    return scenario_results


def main():
    parser = argparse.ArgumentParser(
        description='Middle Permian Sr Isotope Three Declines Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 运行完整分析
  python run_sr_three_declines.py --all --visualize

  # 只运行Scenario A
  python run_sr_three_declines.py --scenario A --n-runs-step1 1000

  # 快速测试（减少运行次数）
  python run_sr_three_declines.py --all --n-runs-step1 100 --n-runs-step2 50

  # 只绘制数据，不运行模拟
  python run_sr_three_declines.py --plot-data
        """
    )

    # 情景选择
    parser.add_argument('--all', action='store_true',
                       help='Run all three scenarios (A, B, C)')
    parser.add_argument('--scenario-a', action='store_true',
                       help='Run Scenario A: Roadian Glaciation')
    parser.add_argument('--scenario-b', action='store_true',
                       help='Run Scenario B: Wordian Rifting')
    parser.add_argument('--scenario-c', action='store_true',
                       help='Run Scenario C: Capitanian LIP')

    # 蒙特卡洛参数
    parser.add_argument('--n-runs-step1', type=int, default=5000,
                       help='Step 1 Monte Carlo runs (default: 5000)')
    parser.add_argument('--n-runs-step2', type=int, default=2000,
                       help='Step 2 Monte Carlo runs (default: 2000)')
    parser.add_argument('--tolerance', type=float, default=0.0002,
                       help='Matching tolerance (default: 0.0002)')

    # 数据参数
    parser.add_argument('--mn-sr-threshold', type=float, default=1.0,
                       help='Mn/Sr screening threshold (default: 1.0)')
    parser.add_argument('--lowess-frac', type=float, default=0.3,
                       help='LOWESS smoothing parameter (default: 0.3)')

    # 输出选项
    parser.add_argument('--output-dir', type=str, default='results/sr_analysis',
                       help='Output directory (default: results/sr_analysis)')
    parser.add_argument('--visualize', action='store_true',
                       help='Generate visualization plots')
    parser.add_argument('--plot-data', action='store_true',
                       help='Only plot input data, skip simulation')

    args = parser.parse_args()

    # 如果没有指定情景，默认运行全部
    if not (args.all or args.scenario_a or args.scenario_b or args.scenario_c or args.plot_data):
        args.all = True

    # 运行分析
    try:
        results = run_full_analysis(args)
        return 0
    except Exception as e:
        print(f"\n错误: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
