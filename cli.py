#!/usr/bin/env python3
"""
同位素质量平衡模型 - 命令行接口
"""

import argparse
import sys
from pathlib import Path
import numpy as np

# 添加项目根目录到路径
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))


def main():
    parser = argparse.ArgumentParser(
        description='Isotope Mass Balance Modeling Framework',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Mg同位素分析 - 碳酸盐体系 (海相碳酸盐岩)
  python cli.py mg --component-type carbonate --file data/Nie_Section_A.xlsx

  # Mg同位素分析 - 碎屑岩体系 (陆源碎屑沉积物)
  python cli.py mg --component-type siliciclastic --delta-sample -0.10

  # Mg同位素批量处理 - 碎屑岩体系
  python cli.py mg --component-type siliciclastic --file data/clay_samples.xlsx

  # C同位素DOC模型
  python cli.py c --scenario dice

  # N同位素正向计算
  python cli.py n --f-assimilator 0.5

  # N同位素反向计算
  python cli.py n --delta15N 3.0 --inverse

  # N同位素关系曲线
  python cli.py n --curve --output n_curve.csv

  # N同位素批量处理Excel文件（含不确定度）
  python cli.py n --file data/nitrogen_data.xlsx --column delta15N --delta-std 0.3 --output results/n_results.xlsx

  # N同位素批量处理（自定义分馏系数不确定度）
  python cli.py n --file data/nitrogen_data.csv --delta-std 0.5 --epsilon-fix-std 0.5 --epsilon-wcd-std 2.0 --output results/n_results.csv

  # N同位素批量处理CSV文件（早三叠世情景）
  python cli.py n --file data/nitrogen_data.csv --scenario early_triassic --output results/n_triassic.csv

  # U同位素单点计算
  python cli.py u --delta-carb -0.65 --steady-state

  # U同位素批量处理Excel文件
  python cli.py u --file data/uranium_data.xlsx --output results/u_results.xlsx

  # U同位素带不确定度分析
  python cli.py u --delta-carb -0.65 --steady-state --uncertainty mc --n-samples 50000

  # U同位素敏感性分析
  python cli.py u --delta-carb -0.65 --sensitivity-analysis

  # U同位素非稳态模拟
  python cli.py u --transient --event-duration 1.0 --peak-f-anox 0.8

  # Sr同位素基本计算
  python cli.py sr --calculate

  # Sr同位素自定义参数计算
  python cli.py sr --calculate --F-highT 20e9 --R-river 0.705

  # Sr同位素蒙特卡洛模拟（二叠纪）
  python cli.py sr --monte-carlo --scenario permian --n-runs 5000

  # Sr同位素情景分析（Wang et al., 2021 Scenario 7）
  python cli.py sr --scenario-analysis scenario7 --n-runs 3000

  # Sr同位素敏感性分析 - 河流同位素
  python cli.py sr --sensitivity R_river --sens-range 0.705 0.725

  # Sr同位素敏感性分析 - 高温热液通量
  python cli.py sr --sensitivity F_highT --sens-range 2e9 35e9

  # 列出支持的体系
  python cli.py list
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # ===== Mg同位素命令 =====
    mg_parser = subparsers.add_parser('mg', help='Mg isotope system (weathering flux model)')

    # 关键参数：区分碳酸盐 vs 碎屑岩体系
    mg_parser.add_argument('--component-type', type=str, default='carbonate',
                          choices=['carbonate', 'siliciclastic', 'silicate', 'detrital'],
                          help='Sample component type (default: carbonate). '
                               'carbonate=海相碳酸盐岩 (Kasemann et al., 2014); '
                               'siliciclastic=陆源碎屑沉积物 (Hu et al., 2023)')
    mg_parser.add_argument('--basin', type=str, default='changjiang',
                          choices=['changjiang', 'global', 'custom'],
                          help='Basin type for siliciclastic system (default: changjiang)')

    mg_parser.add_argument('--file', type=str,
                          help='Input Excel/CSV file with Mg isotope data for batch processing')
    mg_parser.add_argument('--column', type=str, default='delta_26_Mg_iso',
                          choices=['delta_25_Mg_iso', 'delta_26_Mg_iso'],
                          help='Isotope column name')
    mg_parser.add_argument('--sediment-rate', type=float, default=3,
                          help='Reserved for a future age model; row order is used when depth is absent')
    mg_parser.add_argument('--delta-protolith', type=float, default=None,
                          help='Fix silicate protolith delta26Mg; disables its Monte Carlo range')
    mg_parser.add_argument('--delta-protolith-range', type=float, nargs=2,
                          metavar=('LOW', 'HIGH'), default=None,
                          help='Uniform protolith delta26Mg prior '
                               '(default: -0.45 -0.25)')
    mg_parser.add_argument('--delta-fluid-protolith', type=float, default=None,
                          help='Fluid-protolith Mg fractionation in permil (default: -0.25)')
    mg_parser.add_argument('--delta-river-mg', type=float, default=None,
                          help='Dissolved river delta26Mg for flux partitioning (default: -1.14)')
    mg_parser.add_argument('--river-mg-flux', type=float, default=None,
                          help='Total dissolved river Mg flux in mol/yr (default: 3e11)')
    mg_parser.add_argument('--river-mg-flux-rel-std', type=float, default=0.0,
                          help='Relative 1-sigma uncertainty of total river Mg flux (default: 0)')
    mg_parser.add_argument('--mg-monte-carlo', type=int, default=2000,
                          help='Monte Carlo draws per siliciclastic sample (default: 2000)')
    mg_parser.add_argument('--random-seed', type=int, default=42,
                          help='Random seed for siliciclastic uncertainty analysis')
    mg_parser.add_argument('--baseline-samples', type=int, default=5,
                          help=argparse.SUPPRESS)
    mg_parser.add_argument('--baseline-position', choices=['shallow', 'deep'], default='deep',
                          help=argparse.SUPPRESS)
    mg_parser.add_argument('--plot', action='store_true',
                          help='Create a portrait baseline-free Mg-release plot')
    mg_parser.add_argument('--plot-output', type=str, default=None,
                          help='Path for the portrait Mg-release plot (PNG/PDF/SVG)')
    mg_parser.add_argument('--flux-plot-output', type=str, default=None,
                          help='Path for the portrait Changjiang-scaled conditional flux plot')
    mg_parser.add_argument('--time-start', type=float, default=None,
                          help='Time assigned to the first input row (requires --time-end)')
    mg_parser.add_argument('--time-end', type=float, default=None,
                          help='Time assigned to the last input row (requires --time-start)')
    mg_parser.add_argument('--time-unit', type=str, default='Myr',
                          help='Label for user-supplied time anchors (default: Myr)')

    # 基础分析
    mg_parser.add_argument('--weathering-ratio', action='store_true',
                          help='Calculate weathering end-member ratios from sample data')
    mg_parser.add_argument('--delta-sample', type=float,
                          help='Sample δ²⁶Mg value for weathering ratio calculation (‰)')
    mg_parser.add_argument('--delta-seawater', type=float, default=-0.83,
                          help='Seawater δ²⁶Mg value (‰, default: -0.83)')

    # 风化模拟
    mg_parser.add_argument('--weathering-simulation', action='store_true',
                          help='Run weathering regime transition simulation')
    mg_parser.add_argument('--f-initial', type=float, default=0.2,
                          help='Initial silicate weathering fraction (0-1, default: 0.2)')
    mg_parser.add_argument('--f-final', type=float, default=0.8,
                          help='Final silicate weathering fraction (0-1, default: 0.8)')
    mg_parser.add_argument('--transition-mode', type=str, default='linear',
                          choices=['linear', 'exponential'],
                          help='Transition mode (default: linear)')
    mg_parser.add_argument('--transition-start', type=float, default=0.0,
                          help='Transition start time (Myr, default: 0)')
    mg_parser.add_argument('--transition-end', type=float, default=2.0,
                          help='Transition end time (Myr, default: 2)')
    mg_parser.add_argument('--duration', type=float, default=5.0,
                          help='Total simulation duration (Myr, default: 5)')
    mg_parser.add_argument('--flux-multiplier', type=float, default=1.0,
                          help='Weathering flux multiplier relative to modern (default: 1)')

    # Cryogenian论文情景
    mg_parser.add_argument('--cryogenian-scenario', action='store_true',
                          help='Run Cryogenian post-glacial scenario (Kasemann et al., 2014)')
    mg_parser.add_argument('--cryogenian-duration', type=float, default=3.0,
                          help='Cryogenian scenario duration (Myr, default: 3)')

    # 反演计算
    mg_parser.add_argument('--inverse', action='store_true',
                          help='Inverse calculation: infer weathering flux from carbonate data')

    # 通用输出
    mg_parser.add_argument('--output', type=str, help='Output file path (CSV format)')

    mg_parser.add_argument('--n-points', type=int, default=500,
                          help='Number of time points for simulation (default: 500)')

    # 数据解释选项
    mg_parser.add_argument('--data-type', type=str, default='carbonate',
                          choices=['carbonate', 'seawater', 'river'],
                          help='Input data type interpretation (default: carbonate)')
    mg_parser.add_argument('--seawater-correction', type=float, default=None,
                          help='Apply seawater correction factor (‰, default: -2.7 for carbonate)')
    mg_parser.add_argument('--use-raw-values', action='store_true',
                          help='Use raw values directly without correction (for comparison with end-members)')

    # 自定义端元值（用于匹配特定数据集）
    mg_parser.add_argument('--delta-silicate', type=float, default=None,
                          help='Carbonate-model silicate end-member; legacy alias for '
                               '--delta-protolith in siliciclastic mode')
    mg_parser.add_argument('--delta-carbonate', type=float, default=None,
                          help='Custom carbonate end-member (carbonate default: -2.5; '
                               'siliciclastic default: -2.0)')
    mg_parser.add_argument('--apply-offset', type=float, default=None,
                          help='Apply offset correction to all data values (‰)')

    # ===== C同位素命令 =====
    c_parser = subparsers.add_parser('c', help='C isotope system')
    c_parser.add_argument('--scenario', type=str, default='dice',
                         choices=['dice', 'modern'],
                         help='Model scenario')
    c_parser.add_argument('--F-odoc', type=float,
                         help='DOC remineralization flux (mol/Ma)')
    c_parser.add_argument('--target-excursion', type=float,
                         help='Target carbon isotope excursion (negative, in permil)')

    c_parser.add_argument('--output', type=str, help='Output file path')

    # ===== N同位素命令 =====
    n_parser = subparsers.add_parser('n', help='Nitrogen isotope system')
    n_parser.add_argument('--scenario', type=str, default='modern',
                         help='Model scenario: modern, early_triassic, neoproterozoic (default: modern)')
    n_parser.add_argument('--f-assimilator', type=float,
                         help='Nitrate assimilator fraction (0-1) for forward model')
    n_parser.add_argument('--delta15N', type=float,
                         help='Sediment delta15N value for inverse model')
    n_parser.add_argument('--inverse', action='store_true',
                         help='Run inverse model (delta15N -> f_assimilator)')
    n_parser.add_argument('--curve', action='store_true',
                         help='Calculate f_assimilator vs delta15N curve')
    n_parser.add_argument('--f-range', type=float, nargs=2, default=[0.0, 1.0],
                         metavar=('MIN', 'MAX'),
                         help='Range for f_assimilator (default: 0.0 1.0)')
    n_parser.add_argument('--n-points', type=int, default=100,
                         help='Number of points for curve calculation')
    n_parser.add_argument('--monte-carlo', action='store_true',
                         help='Run Monte Carlo uncertainty analysis')
    n_parser.add_argument('--n-samples', type=int, default=10000,
                         help='Number of Monte Carlo samples')
    n_parser.add_argument('--output', type=str,
                         help='Output file path (CSV format)')
    n_parser.add_argument('--file', type=str,
                         help='Input Excel/CSV file with nitrogen isotope data for batch processing')
    n_parser.add_argument('--column', type=str, default='delta15N',
                         help='Column name for delta15N values (default: delta15N)')
    n_parser.add_argument('--delta-std', type=float, default=0.3,
                         help='Measurement uncertainty for delta15N (1σ, ‰, default: 0.3)')
    n_parser.add_argument('--epsilon-fix-std', type=float, default=0.5,
                         help='Uncertainty of epsilon_fix (1σ, ‰, default: 0.5)')
    n_parser.add_argument('--epsilon-wcd-std', type=float, default=2.0,
                         help='Uncertainty of epsilon_wcd (1σ, ‰, default: 2.0)')
    n_parser.add_argument('--n-monte-carlo', type=int, default=2000,
                         help='Number of Monte Carlo samples for uncertainty analysis (default: 2000)')

    # ===== Sr同位素命令 =====
    sr_parser = subparsers.add_parser('sr', help='Strontium isotope system (ocean box model)')
    sr_parser.add_argument('--scenario', type=str, default='modern',
                          choices=['modern', 'permian'],
                          help='Model scenario (default: modern)')

    # 基本计算
    sr_parser.add_argument('--calculate', action='store_true',
                          help='Calculate seawater Sr ratio from fluxes')
    sr_parser.add_argument('--F-river', type=float, default=47.6e9,
                          help='Riverine flux in mol/yr (default: 47.6e9)')
    sr_parser.add_argument('--R-river', type=float, default=0.71107,
                          help='Riverine 87Sr/86Sr (default: 0.71107)')
    sr_parser.add_argument('--F-highT', type=float, default=8.04e9,
                          help='High-T hydrothermal flux in mol/yr (default: 8.04e9)')
    sr_parser.add_argument('--R-highT', type=float, default=0.7037,
                          help='High-T hydrothermal 87Sr/86Sr (default: 0.7037)')
    sr_parser.add_argument('--F-lowT', type=float, default=10e9,
                          help='Low-T hydrothermal flux in mol/yr (default: 10e9)')
    sr_parser.add_argument('--R-lowT', type=float, default=0.7084,
                          help='Low-T hydrothermal 87Sr/86Sr (default: 0.7084)')
    sr_parser.add_argument('--F-dia', type=float, default=3.4e9,
                          help='Diagenetic flux in mol/yr (default: 3.4e9)')
    sr_parser.add_argument('--R-dia', type=float, default=0.7084,
                          help='Diagenetic 87Sr/86Sr (default: 0.7084)')

    # 蒙特卡洛模拟
    sr_parser.add_argument('--monte-carlo', action='store_true',
                          help='Run Monte Carlo stochastic simulation')
    sr_parser.add_argument('--n-runs', type=int, default=5000,
                          help='Number of Monte Carlo runs (default: 5000)')
    sr_parser.add_argument('--time-span', type=float, nargs=2, default=[299, 252],
                          metavar=('START', 'END'),
                          help='Time span in Ma (default: 299 252 for Permian)')
    sr_parser.add_argument('--n-time-points', type=int, default=50,
                          help='Number of time points (default: 50)')

    # 情景分析
    sr_parser.add_argument('--scenario-analysis', type=str,
                          choices=['scenario1', 'scenario2', 'scenario3',
                                  'scenario4', 'scenario5', 'scenario6', 'scenario7'],
                          help='Run scenario analysis (Wang et al., 2021)')

    # 敏感性分析
    sr_parser.add_argument('--sensitivity', type=str,
                          choices=['F_river', 'R_river', 'F_highT', 'R_highT',
                                  'F_lowT', 'R_lowT', 'F_dia', 'R_dia'],
                          help='Parameter for sensitivity analysis')
    sr_parser.add_argument('--sens-range', type=float, nargs=2, default=None,
                          metavar=('MIN', 'MAX'),
                          help='Sensitivity range (auto if not specified)')
    sr_parser.add_argument('--n-steps', type=int, default=50,
                          help='Number of steps for sensitivity (default: 50)')

    # Dukou剖面分析
    sr_parser.add_argument('--dukou-analysis', action='store_true',
                          help='Analyze Dukou section Sr isotope data')
    sr_parser.add_argument('--mn-sr-threshold', type=float, default=1.0,
                          help='Mn/Sr screening threshold (default: 1.0)')
    sr_parser.add_argument('--lowess-frac', type=float, default=0.3,
                          help='LOWESS smoothing parameter (default: 0.3)')

    # 三步情景分析
    sr_parser.add_argument('--scenario-a', action='store_true',
                          help='Run Scenario A: Roadian glaciation (Decline 1)')
    sr_parser.add_argument('--scenario-b', action='store_true',
                          help='Run Scenario B: Wordian rifting (Decline 2)')
    sr_parser.add_argument('--scenario-c', action='store_true',
                          help='Run Scenario C: Capitanian LIP (Decline 3)')
    sr_parser.add_argument('--two-step-mc', action='store_true',
                          help='Run two-step Monte Carlo (Step 1: 5000 + Step 2: 2000)')
    sr_parser.add_argument('--n-runs-step1', type=int, default=5000,
                          help='Step 1 MC runs (default: 5000)')
    sr_parser.add_argument('--n-runs-step2', type=int, default=2000,
                          help='Step 2 MC runs (default: 2000)')

    # 输出
    sr_parser.add_argument('--output', type=str, help='Output CSV file path')
    sr_parser.add_argument('--verbose', '-v', action='store_true',
                          help='Verbose output')

    # ===== U同位素命令 =====
    u_parser = subparsers.add_parser('u', help='Uranium isotope system')
    u_parser.add_argument('--file', type=str,
                         help='Input Excel/CSV file with isotope data')
    u_parser.add_argument('--scenario', type=str, default='modern',
                         choices=['modern', 'oceanic_anoxic_event', 'end_permain',
                                 'frasnian_famennian'],
                         help='Model scenario')
    u_parser.add_argument('--delta-carb', type=float,
                         help='Measured carbonate δ²³⁸U value (single point)')
    u_parser.add_argument('--steady-state', action='store_true',
                         help='Run steady-state model')
    u_parser.add_argument('--transient', action='store_true',
                         help='Run transient (non-steady-state) model')
    u_parser.add_argument('--event-duration', type=float, default=1.0,
                         help='Anoxic event duration (Myr)')
    u_parser.add_argument('--peak-f-anox', type=float, default=0.8,
                         help='Peak anoxic sink fraction')
    u_parser.add_argument('--background-f-anox', type=float, default=0.2,
                         help='Background anoxic sink fraction')
    u_parser.add_argument('--delta-diag', type=float, default=0.4,
                         help='Diagenetic correction factor (‰)')
    u_parser.add_argument('--no-diagenetic-correction', action='store_true',
                         help='Disable diagenetic correction')
    u_parser.add_argument('--uncertainty', type=str,
                         choices=['mc', 'monte-carlo', 'bootstrap'],
                         help='Uncertainty analysis method')
    u_parser.add_argument('--measurement-std', type=float, default=0.05,
                         help='Measurement uncertainty (1σ, ‰)')
    u_parser.add_argument('--n-samples', type=int, default=10000,
                         help='Number of Monte Carlo or bootstrap samples')
    u_parser.add_argument('--sensitivity-analysis', action='store_true',
                         help='Run sensitivity analysis')
    u_parser.add_argument('--confidence-level', type=float, default=0.95,
                         help='Confidence level for uncertainty intervals')
    u_parser.add_argument('--no-uncertainty', action='store_true',
                         help='Disable uncertainty calculation for batch processing')
    u_parser.add_argument('--output', type=str,
                         help='Output file path (for batch processing)')

    # ===== 体系列表命令 =====
    subparsers.add_parser('list', help='List available isotope systems')

    # ===== 信息命令 =====
    info_parser = subparsers.add_parser('info', help='Show system information')
    info_parser.add_argument('element', type=str,
                            choices=['mg', 'c', 'n', 'u', 's', 'sr', 'nd'],
                            help='Element symbol')

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    # 执行命令
    if args.command == 'list':
        list_systems()
    elif args.command == 'info':
        show_info(args.element)
    elif args.command == 'mg':
        run_mg_analysis(args)
    elif args.command == 'c':
        run_c_analysis(args)
    elif args.command == 'u':
        run_u_analysis(args)
    elif args.command == 'n':
        run_n_analysis(args)
    elif args.command == 'sr':
        run_sr_analysis(args)


def list_systems():
    """列出可用的同位素体系"""
    print("\n=== Available Isotope Systems ===\n")

    systems = [
        ('mg', 'Magnesium', 'Mg同位素风化体系'),
        ('c', 'Carbon', '碳循环，DOC氧化与碳同位素负漂'),
        ('n', 'Nitrogen', '氮循环，硝酸盐可利用性示踪'),
        ('u', 'Uranium', '海洋铀循环，氧化还原条件示踪'),
        ('s', 'Sulfur', '硫循环，硫酸盐还原（计划中）'),
        ('sr', 'Strontium', 'Sr同位素海洋箱模型（Wang et al., 2021）'),
        ('nd', 'Neodymium', 'Nd同位素，洋流循环（计划中）'),
    ]

    for element, name, description in systems:
        status = "✓" if element in ['mg', 'c', 'n', 'u', 'sr'] else "○"
        print(f"  {status} {element.upper():2} - {name:12} : {description}")

    print("\n✓ = Implemented, ○ = Planned")

    # Mg 子体系详细说明
    print("\n--- Mg Isotope Sub-systems ---")
    print("\n  1. Carbonate System (--component-type carbonate)")
    print("     Reference: Kasemann et al. (2014)")
    print("     Input: 海相碳酸盐岩 δ²⁶Mg")
    print("     Model: 海水沉淀分馏")
    print("     Output: 风化比例、海水演化")

    print("\n  2. Siliciclastic System (--component-type siliciclastic)")
    print("     Reference: Hu et al. (2023)")
    print("     Input: 陆源碎屑沉积物 δ²⁶Mg (黏土矿物)")
    print("     Model: 风化残余分馏 (Rayleigh)")
    print("     Output: 硅酸盐风化通量、SWI指数")

    print()


def show_info(element: str):
    """显示体系详细信息"""
    print(f"\n=== {element.upper()} Isotope System ===\n")

    if element == 'mg':
        print("Mg同位素体系包含两种子体系，使用 --component-type 参数选择：\n")

        # 碳酸盐体系
        print("[1] Carbonate System (--component-type carbonate)")
        print("-" * 50)
        from systems.mg import get_mg_parameters
        params = get_mg_parameters()

        print(f"Reference: Kasemann et al. (2014) EPSL")
        print(f"Reference standard: {params.reference_standard}")
        print(f"\nEnd-members (δ²⁶Mg, ‰):")
        for name, data in params.end_members.items():
            print(f"  {name:12}: {data['delta26']:+.2f} ± {data.get('uncertainty', 0):.2f}")

        print(f"\nReservoir mass: {params.reservoir_mass:.2e} mol")
        print(f"Input fluxes (mol/Ma):")
        for name, value in params.input_fluxes.items():
            print(f"  {name}: {value:.2e}")

        print(f"\nFractionation factors:")
        for name, value in params.fractionation_factors.items():
            print(f"  {name}: {value:+.2f}‰")

        # 碎屑岩体系
        print("\n[2] Siliciclastic System (--component-type siliciclastic)")
        print("-" * 50)
        from systems.mg.silicate import SilicateWeatheringParams
        sil_params = SilicateWeatheringParams()

        print(f"Reference: Hu et al. (2023) Global and Planetary Change")
        print(f"Reference standard: DSM3")
        print(f"\nEnd-members (δ²⁶Mg, ‰):")
        print(f"  UCC:           {sil_params.d26Mg_UCC:+.2f}‰")
        print(f"  Carbonate:     {sil_params.d26Mg_carbonate:+.2f}‰")
        print(f"  River water:   {sil_params.d26Mg_river_water:+.2f}‰")

        print(f"\nFractionation factors:")
        print(f"  Δ(fluid-protolith): {sil_params.Delta_fluid_protolith:+.2f}‰")
        print(f"  Literature range:   {sil_params.Delta_fluid_protolith_range[0]:+.2f} "
              f"to {sil_params.Delta_fluid_protolith_range[1]:+.2f}‰")

        print(f"\nFluxes:")
        print(f"  River total: {sil_params.F_river_total/1e10:.1f} × 10¹⁰ mol/yr")

    elif element == 'c':
        from systems.c import CIsotopeSystem, get_c_parameters
        params = get_c_parameters('dice')

        print(f"Element: {params.name}")
        print(f"Reference: {params.reference_standard}")
        print(f"\nEnd-members (δ¹³C, ‰):")
        for name, data in params.end_members.items():
            print(f"  {name:12}: {data['delta13']:+.1f}")

        print(f"\nDIC reservoir: {params.reservoir_mass:.2e} mol")
        print(f"Organic burial fraction: {params.fractionation_factors.get('organic_burial_fraction', 0.14)}")

    elif element == 'u':
        from systems.u import UIsotopeSystem, get_u_parameters
        params = get_u_parameters('modern')

        print(f"Element: {params.name}")
        print(f"Reference: {params.reference_standard}")
        print(f"\nEnd-members (δ²³⁸U, ‰):")
        for name, data in params.end_members.items():
            if isinstance(data, dict) and 'delta238' in data:
                print(f"  {name:15}: {data['delta238']:+.2f}")

        print(f"\nReservoir mass: {params.reservoir_mass:.2e} mol")
        print(f"Residence time: ~0.5 Myr")
        print(f"\nFractionation factors:")
        print(f"  Δ(ox):   {params.fractionation_factors.get('delta_sw_ox', 0):.2f}‰")
        print(f"  Δ(anox): +{params.fractionation_factors.get('delta_sw_anox', 0):.2f}‰")
        print(f"  Δ(diag): +{params.fractionation_factors.get('delta_diag', 0):.2f}‰")

    elif element == 'n':
        from systems.n import NIsotopeSystem, get_n_parameters
        params = get_n_parameters('modern')

        print(f"Element: {params.name}")
        print(f"Reference: {params.reference_standard}")
        print(f"\nModel: Two-box steady-state nitrogen cycle model")
        print(f"  Reference: Kang et al. (2023) & Ma et al. (2025)")

        print(f"\nReservoir mass: {params.reservoir_mass:.2e} mol N")
        print(f"\nInput fluxes (Tg N/a):")
        for name, value in params.input_fluxes.items():
            print(f"  {name}: {value:.1f}")

        print(f"\nOutput fluxes (Tg N/a):")
        for name, value in params.output_fluxes.items():
            print(f"  {name}: {value:.1f}")

        print(f"\nFractionation factors:")
        ff = params.fractionation_factors
        print(f"  ε_fix: {ff.get('epsilon_fixation', -0.5):+.1f}‰ (range: {ff.get('epsilon_fixation_min', -2.0):+.1f} to {ff.get('epsilon_fixation_max', 1.0):+.1f})")
        print(f"  ε_wcd: {ff.get('epsilon_wcd', -26.0):+.1f}‰ (range: {ff.get('epsilon_wcd_min', -30.0):+.1f} to {ff.get('epsilon_wcd_max', -22.0):+.1f})")
        print(f"  ε_sd:  {ff.get('epsilon_sd', 0.0):+.1f}‰")

        print(f"\nKey relationship:")
        print(f"  δ¹⁵N_sed = (1-f) × δ¹⁵N_NH4 + f × δ¹⁵N_NO3")
        print(f"  where f = nitrate assimilator fraction")
        print(f"  Peak δ¹⁵N occurs at f ≈ 0.48")

    elif element == 'sr':
        from systems.sr import SrIsotopeSystem, get_sr_parameters
        params = get_sr_parameters()

        print(f"Element: {params.name}")
        print(f"Reference: Wang et al. (2021) Earth-Science Reviews")
        print(f"Reference standard: {params.reference_standard}")

        print(f"\nModel: Four-component oceanic box model")
        print(f"  R_sw = (F_riv×R_riv + F_highT×R_highT + F_lowT×R_lowT + F_dia×R_dia) / F_total")

        print(f"\nEnd-members (⁸⁷Sr/⁸⁶Sr):")
        for name, data in params.end_members.items():
            if isinstance(data, dict) and 'ratio' in data:
                print(f"  {name:20}: {data['ratio']:.5f}")

        print(f"\nReservoir mass: {params.reservoir_mass:.2e} mol")
        print(f"Residence time: ~2.6 Myr")

        print(f"\nModern fluxes (10⁹ mol/yr):")
        for name, value in params.input_fluxes.items():
            print(f"  {name:20}: {value/1e9:.2f}")

        print(f"\nKey findings from Wang et al. (2021):")
        print(f"  - Hydrothermal activity is the main driver of Permian Sr isotope trends")
        print(f"  - Low-temperature hydrothermal flux increases with global warming")
        print(f"  - Continental weathering changes alone cannot explain the record")

    print()


def run_mg_analysis(args):
    """
    运行Mg同位素风化分析

    根据 --component-type 参数自动选择体系：
    - carbonate: 碳酸盐体系 (Kasemann et al., 2014)
    - siliciclastic: 碎屑岩体系 (Hu et al., 2023)
    """
    # 根据体系类型分派
    component_type = args.component_type.lower()

    if component_type in ('siliciclastic', 'silicate', 'detrital'):
        run_mg_siliciclastic_analysis(args)
    else:
        run_mg_carbonate_analysis(args)


def run_mg_carbonate_analysis(args):
    """运行碳酸盐体系Mg同位素分析 (基于Kasemann et al., 2014)"""
    print("\n" + "="*70)
    print("Mg Isotope Weathering Flux Model")
    print("Component Type: CARBONATE (海相碳酸盐岩)")
    print("Based on: Kasemann et al. (2014) EPSL")
    print("="*70 + "\n")

    from systems.mg import create_mg_system

    system = create_mg_system('carbonate', scenario='modern')

    # 应用自定义端元值
    if args.delta_silicate is not None:
        system._delta_sil = args.delta_silicate
        print(f"[Custom End-Members]")
        print(f"  Silicate δ²⁶Mg: {args.delta_silicate:+.2f}‰ (user-defined)")
    if args.delta_carbonate is not None:
        system._delta_carb = args.delta_carbonate
        print(f"  Carbonate δ²⁶Mg: {args.delta_carbonate:+.2f}‰ (user-defined)")
    if args.delta_silicate is not None or args.delta_carbonate is not None:
        print()

    # ===== 风化比例计算 =====
    if args.weathering_ratio or args.delta_sample is not None:
        print("[Weathering End-Member Calculation]")
        print("-" * 50)

        delta_sample = args.delta_sample if args.delta_sample is not None else -2.0

        print(f"  Input:")
        print(f"    Sample δ²⁶Mg: {delta_sample:.2f}‰")
        print(f"    Seawater δ²⁶Mg: {args.delta_seawater:.2f}‰")
        print(f"    Silicate end-member: {system._delta_sil:.2f}‰")
        print(f"    Carbonate end-member: {system._delta_carb:.2f}‰")

        ratios = system.calculate_weathering_ratio(
            delta_sample=delta_sample,
            delta_seawater=args.delta_seawater
        )

        print(f"\n  Results:")
        print(f"    Silicate weathering fraction (f_sil): {ratios['f_silicate']:.2%}")
        print(f"    Carbonate weathering fraction (f_carb): {ratios['f_carbonate']:.2%}")
        print(f"    Inferred river δ²⁶Mg: {ratios['delta_river']:+.2f}‰")
        print()

    # ===== 风化转变模拟 =====
    if args.weathering_simulation:
        print("[Weathering Regime Transition Simulation]")
        print("-" * 50)

        print(f"  Configuration:")
        print(f"    Initial f_silicate: {args.f_initial:.2f}")
        print(f"    Final f_silicate: {args.f_final:.2f}")
        print(f"    Transition mode: {args.transition_mode}")
        print(f"    Transition period: {args.transition_start}-{args.transition_end} Myr")
        print(f"    Total duration: {args.duration} Myr")
        print(f"    Flux multiplier: {args.flux_multiplier}× modern")
        print()

        transition = {
            't_start': args.transition_start,
            't_end': args.transition_end,
            'f_initial': args.f_initial,
            'f_final': args.f_final,
            'mode': args.transition_mode
        }

        # 估算初始海水δ²⁶Mg
        from systems.mg.carbonate import WeatheringFluxConfig
        initial_config = WeatheringFluxConfig(
            f_silicate=args.f_initial,
            F_riv_multiplier=args.flux_multiplier
        )
        initial_delta = system.weathering_model.steady_state_seawater(initial_config)

        print(f"  Running simulation (initial δ²⁶Mg_sw = {initial_delta:+.2f}‰)...")

        result = system.simulate_weathering_transition(
            time_span=(0, args.duration),
            transition=transition,
            initial_delta=initial_delta,
            n_points=args.n_points
        )

        if result.success:
            print(f"\n  Results at key time points:")
            print(f"  {'Time (Myr)':<12} {'f_sil':<10} {'f_carb':<10} {'δ²⁶Mg_sw':<12}")
            print("  " + "-" * 44)

            times_ma = result.time / 1e6
            delta_sw = result.data['delta_sw']
            f_sil = result.data['f_silicate']

            key_times = [0, args.transition_start,
                        (args.transition_start + args.transition_end)/2,
                        args.transition_end, args.duration]

            for t in sorted(set(key_times)):
                if t <= times_ma[-1]:
                    idx = np.argmin(np.abs(times_ma - t))
                    print(f"  {t:<12.2f} {f_sil[idx]:<10.2f} {1-f_sil[idx]:<10.2f} {delta_sw[idx]:<+12.2f}")

            # 保存结果
            if args.output:
                save_mg_simulation_results(result, args.output)
                print(f"\n  Results saved to: {args.output}")
        else:
            print(f"  Error: {result.message}")
        print()

    # ===== Cryogenian情景 =====
    if args.cryogenian_scenario:
        print("[Cryogenian Post-Glacial Scenario]")
        print("-" * 50)
        print(f"  Based on: Kasemann et al. (2014)")
        print(f"  Duration: {args.cryogenian_duration} Myr")
        print()
        print(f"  Scenario parameters:")
        print(f"    Phase 1 (0-0.5 Myr): Mixed weathering, 9× modern flux")
        print(f"    Phase 2 (0.5-{args.cryogenian_duration} Myr): Silicate-dominated, 6× modern flux")
        print()

        result = system.simulate_cryogenian_scenario(
            duration_ma=args.cryogenian_duration,
            n_points=args.n_points
        )

        if result.success:
            print(f"  Simulation results:")
            print(f"  {'Time (Myr)':<12} {'Flux (×)':<12} {'f_sil':<10} {'δ²⁶Mg_sw':<12}")
            print("  " + "-" * 46)

            times_ma = result.time / 1e6
            delta_sw = result.data['delta_sw']
            f_sil = result.data['f_silicate']
            flux_mult = result.data['flux_multiplier']

            for t in [0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]:
                if t <= times_ma[-1]:
                    idx = np.argmin(np.abs(times_ma - t))
                    print(f"  {t:<12.2f} {flux_mult[idx]:<12.1f} {f_sil[idx]:<10.2f} {delta_sw[idx]:<+12.2f}")

            if args.output:
                save_mg_simulation_results(result, args.output)
                print(f"\n  Results saved to: {args.output}")
        else:
            print(f"  Error: {result.message}")
        print()

    # ===== 批量处理模式 =====
    if args.file:
        print("[Batch Processing Mode]")
        print("-" * 50)
        print(f"  Input file: {args.file}")

        from toolkit.io import BatchProcessor
        processor = BatchProcessor(element='mg')

        try:
            results_df = processor.process_file(
                args.file,
                output_path=args.output,
                show_progress=True
            )

            print(f"\n  Summary:")
            print(f"    Total samples: {len(results_df)}")
            success_count = results_df['processing_success'].sum()
            print(f"    Successful: {success_count}")

            if 'f_carbonate' in results_df.columns:
                print(f"\n  Weathering fraction statistics:")
                print(f"    Mean carbonate: {results_df['f_carbonate'].mean():.2%}")
                print(f"    Mean silicate: {results_df['f_silicate'].mean():.2%}")

        except Exception as e:
            print(f"  Error: {e}")
            import traceback
            traceback.print_exc()
        print()

    # ===== 反演计算 =====
    if args.inverse:
        print("[Inverse Calculation - Weathering Flux Inference]")
        print("-" * 50)

        if not args.file:
            print("  Error: --inverse requires --file to specify input data")
            print()
        else:
            run_mg_inverse(system, args)

    # ===== 剖面数据处理 =====
    if args.file and not args.inverse and not hasattr(args, 'processed_batch'):
        run_mg_section_analysis(system, args)

    print("="*70)
    print("Analysis complete")
    print("="*70 + "\n")


def run_mg_inverse(system, args):
    """运行Mg同位素反演计算"""
    import pandas as pd

    try:
        df = pd.read_excel(args.file)
        print(f"  Loaded {len(df)} samples from {args.file}")

        # 获取δ²⁶Mg数据
        if 'delta_26_Mg_iso' in df.columns:
            delta_data = df['delta_26_Mg_iso'].values
        else:
            print("  Error: Column 'delta_26_Mg_iso' not found")
            return

        # 获取年龄数据（如果有）
        age_data = None
        for col in ['age', 'age_ma', 'Age', 'Age_Ma']:
            if col in df.columns:
                age_data = df[col].values
                break

        if age_data is None:
            # 如果没有年龄，假设等间距分布
            print("  Note: No age column found, assuming uniform spacing")
            age_data = np.linspace(0, args.duration, len(df))

        # 运行反演
        result = system.inverse_weathering_flux(
            age_data=age_data,
            delta_carb_data=delta_data,
            assume_steady_state=True
        )

        if result.success:
            print(f"\n  Inversion Results:")
            print(f"  {'Age (Myr)':<12} {'δ²⁶Mg':<10} {'f_silicate':<12} {'f_carbonate':<12}")
            print("  " + "-" * 46)

            ages = result.get('age_ma')
            f_sil = result.get('f_silicate')
            f_carb = result.get('f_carbonate')

            for i in range(min(10, len(ages))):  # 显示前10个
                print(f"  {ages[i]:<12.2f} {delta_data[i]:<+10.2f} {f_sil[i]:<12.2f} {f_carb[i]:<12.2f}")

            if len(ages) > 10:
                print(f"  ... and {len(ages)-10} more samples")

            # 统计信息
            print(f"\n  Statistics:")
            print(f"    Mean f_silicate: {np.mean(f_sil):.2%}")
            print(f"    Range: {np.min(f_sil):.2f} - {np.max(f_sil):.2f}")

            # 保存结果
            if args.output:
                output_df = pd.DataFrame({
                    'age_myr': ages,
                    'delta_26Mg': delta_data,
                    'f_silicate': f_sil,
                    'f_carbonate': f_carb,
                    'delta_river_inferred': result.get('delta_river_inferred')
                })
                output_df.to_csv(args.output, index=False)
                print(f"\n  Results saved to: {args.output}")

    except Exception as e:
        print(f"  Error: {e}")
        import traceback
        traceback.print_exc()
    print()


def run_mg_section_analysis(system, args):
    """分析Mg同位素剖面数据（如Nie Section）"""
    import pandas as pd

    print("[Section Data Analysis]")
    print("-" * 50)
    print(f"  Reference Standard: D3MS = DSM3 (Dead Sea Magnesium)")

    try:
        df = pd.read_excel(args.file)

        # 移除完全为空的列
        df = df.dropna(axis=1, how='all')

        print(f"  Loaded {len(df)} samples from {args.file}")
        print(f"  Columns: {list(df.columns)}")

        # 获取δ²⁶Mg数据
        if 'delta_26_Mg_iso' not in df.columns:
            print(f"  Error: Required column 'delta_26_Mg_iso' not found")
            return

        delta_raw = df['delta_26_Mg_iso'].values

        # 应用偏移校正（如果有）
        if args.apply_offset is not None:
            delta_raw = delta_raw + args.apply_offset
            print(f"  Applied offset correction: {args.apply_offset:+.2f}‰")

        # 获取不确定度（如果有）
        delta_err = None
        if 'delta_26_Mg_iso_2sd' in df.columns:
            delta_err = df['delta_26_Mg_iso_2sd'].values / 2  # 转换为1σ
            print(f"  Uncertainty data found (2σ)")

        # 根据数据类型进行校正
        correction = 0
        if args.use_raw_values:
            delta_data = delta_raw
            print(f"  Using raw values (no correction)")
        else:
            # 根据数据类型应用校正
            if args.data_type == 'carbonate':
                # 碳酸盐沉积物 → 海水: δ_sw = δ_carb - Δ_carb
                correction = args.seawater_correction if args.seawater_correction is not None else -2.7
                delta_data = delta_raw - correction  # 反推海水值
                print(f"  Data type: Carbonate (applying correction: {correction:+.2f}‰)")
                print(f"    δ_sw = δ_carb - ({correction:+.2f})")
            elif args.data_type == 'seawater':
                delta_data = delta_raw
                print(f"  Data type: Seawater (no correction)")
            elif args.data_type == 'river':
                delta_data = delta_raw
                print(f"  Data type: River water (no correction)")

        print(f"\n  Raw data range: {np.min(delta_raw):+.2f} to {np.max(delta_raw):+.2f}‰")
        if not args.use_raw_values and args.data_type == 'carbonate':
            print(f"  Corrected range: {np.min(delta_data):+.2f} to {np.max(delta_data):+.2f}‰")

        # 计算每个样品的风化比例
        print(f"\n  Calculating weathering ratios...")
        print(f"    Using end-members: silicate={system._delta_sil:+.2f}‰, carbonate={system._delta_carb:+.2f}‰")

        f_sil_list = []
        f_carb_list = []
        delta_riv_list = []

        # 根据数据类型选择计算方法
        if args.data_type == 'river' or args.use_raw_values:
            # 直接计算风化比例（数据已经是河流或可直接对比端元）
            from systems.mg.parameters import solve_f_silicate
            for delta in delta_data:
                f_sil = solve_f_silicate(delta, system._delta_sil, system._delta_carb)
                delta_riv = delta  # 数据本身就是河流值
                f_sil_list.append(f_sil)
                f_carb_list.append(1 - f_sil)
                delta_riv_list.append(delta_riv)
        else:
            # 使用标准方法（数据→风化比例）
            for delta in delta_data:
                ratios = system.calculate_weathering_ratio(
                    delta_sample=delta,
                    delta_seawater=args.delta_seawater
                )
                f_sil_list.append(ratios['f_silicate'])
                f_carb_list.append(ratios['f_carbonate'])
                delta_riv_list.append(ratios['delta_river'])

        f_sil_array = np.array(f_sil_list)
        f_carb_array = np.array(f_carb_list)

        print(f"\n  Weathering Ratio Results:")
        print(f"  {'Sample':<8} {'δ²⁶Mg':<10} {'f_silicate':<12} {'f_carbonate':<12}")
        print("  " + "-" * 42)

        for i in range(min(5, len(delta_data))):
            print(f"  {i+1:<8} {delta_data[i]:<+10.2f} {f_sil_array[i]:<12.2f} {f_carb_array[i]:<12.2f}")

        if len(delta_data) > 5:
            print(f"  ... and {len(delta_data)-5} more samples")

        print(f"\n  Summary Statistics:")
        print(f"    Mean δ²⁶Mg: {np.mean(delta_data):+.2f} ± {np.std(delta_data):.2f}‰")
        print(f"    Range: {np.min(delta_data):+.2f} to {np.max(delta_data):+.2f}‰")
        print(f"    Mean f_silicate: {np.mean(f_sil_array):.2%} (range: {np.min(f_sil_array):.0%} - {np.max(f_sil_array):.0%})")
        print(f"    Mean f_carbonate: {np.mean(f_carb_array):.2%}")

        # 识别风化转变趋势
        if len(delta_data) >= 5:
            # 简单线性趋势分析
            x = np.arange(len(delta_data))
            slope, intercept = np.polyfit(x, f_sil_array, 1)
            print(f"\n  Trend Analysis:")
            print(f"    f_silicate trend slope: {slope:+.4f} per sample")
            if slope > 0.05:
                trend_desc = "increasing (carbonate → silicate transition)"
            elif slope < -0.05:
                trend_desc = "decreasing (silicate → carbonate transition)"
            else:
                trend_desc = "relatively stable"
            print(f"    Interpretation: {trend_desc}")

        # 保存结果
        if args.output:
            output_df = pd.DataFrame({
                'sample_index': np.arange(1, len(df) + 1),
                'delta_26Mg_raw': delta_raw,
                'delta_26Mg_corrected': delta_data if not args.use_raw_values else delta_raw,
                'delta_26Mg_1sigma': delta_err if delta_err is not None else np.nan,
                'f_silicate': f_sil_array,
                'f_carbonate': f_carb_array,
                'delta_river_inferred': delta_riv_list
            })
            output_df.to_csv(args.output, index=False)
            print(f"\n  Results saved to: {args.output}")

        # 如果用户要求，基于数据趋势进行模拟
        if args.weathering_simulation:
            print(f"\n  Running simulation based on section data trend...")
            # 使用数据中的最小和最大f_silicate作为模拟边界
            f_min = max(0.1, np.min(f_sil_array))
            f_max = min(0.9, np.max(f_sil_array))

            print(f"    Simulating transition: f_sil {f_min:.2f} → {f_max:.2f}")

            transition = {
                't_start': 0,
                't_end': args.duration * 0.6,
                'f_initial': f_min,
                'f_final': f_max,
                'mode': 'linear'
            }

            from systems.mg.model import WeatheringFluxConfig
            initial_config = WeatheringFluxConfig(f_silicate=f_min)
            initial_delta = system.weathering_model.steady_state_seawater(initial_config)

            result = system.simulate_weathering_transition(
                time_span=(0, args.duration),
                transition=transition,
                initial_delta=initial_delta,
                n_points=args.n_points
            )

            if result.success:
                print(f"    Simulation completed successfully")
                print(f"    Initial δ²⁶Mg_sw: {initial_delta:+.2f}‰")
                print(f"    Final δ²⁶Mg_sw: {result.data['delta_sw'][-1]:+.2f}‰")

                # 保存模拟结果
                if args.output:
                    sim_output = args.output.replace('.csv', '_simulation.csv')
                    save_mg_simulation_results(result, sim_output)
                    print(f"    Simulation saved to: {sim_output}")

        args.processed_batch = True  # 标记已处理

    except Exception as e:
        print(f"  Error: {e}")
        import traceback
        traceback.print_exc()
    print()


def save_mg_simulation_results(result, output_path):
    """保存Mg模拟结果到CSV文件"""
    import pandas as pd

    data = {
        'time_myr': result.time / 1e6,
        'M_sw_mol': result.data['M_sw'],
        'delta_sw_permil': result.data['delta_sw'],
        'f_silicate': result.data['f_silicate'],
        'f_carbonate': result.data['f_carbonate'],
        'flux_multiplier': result.data['flux_multiplier'],
        'delta_river_permil': result.data['delta_river']
    }

    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)


def run_c_analysis(args):
    """运行C同位素分析"""
    print(f"\n=== C Isotope Analysis ({args.scenario.upper()}) ===\n")

    from systems.c import CIsotopeSystem

    system = CIsotopeSystem(scenario=args.scenario)

    if args.F_odoc is not None:
        print(f"DOC remineralization flux: {args.F_odoc:.2e} mol/Ma")

        result = system.solve_steady_state(F_odoc=args.F_odoc)

        if result.success:
            print(f"\nSteady-state results:")
            print(f"  δ¹³C_carb: {result.get('delta13C_carb'):.2f}‰")
            print(f"  δ¹³C_org:  {result.get('delta13C_org'):.2f}‰")

            initial = system.params.end_members['seawater_dic']['delta13']
            excursion = result.get('delta13C_carb') - initial
            print(f"  Carbon isotope excursion: {excursion:.2f}‰")

    elif args.target_excursion is not None:
        print(f"Finding DOC flux for {args.target_excursion}‰ excursion...")

        result = system.find_doc_for_excursion(args.target_excursion)

        print(f"\nResult:")
        print(f"  Required F_odoc: {result['F_odoc']:.2e} mol/Ma")
        print(f"  Achieved excursion: {result['excursion']:.2f}‰")

    else:
        # 运行完整模型
        print("Running complete DOC excursion model...")

        result = system.doc_excursion_model(F_odoc_range=(0, 10e18), n_points=100)

        if result.success:
            F_odoc = result.get('F_odoc')
            delta_delta13C = result.get('delta_delta13C')

            # 找到关键值
            idx_2 = np.argmin(np.abs(delta_delta13C + 2.0))
            idx_4 = np.argmin(np.abs(delta_delta13C + 4.0))

            print(f"\nKey results:")
            print(f"  ~2‰ excursion: F_odoc ≈ {F_odoc[idx_2]/1e18:.2f}×10¹⁸ mol/Ma")
            print(f"  ~4‰ excursion: F_odoc ≈ {F_odoc[idx_4]/1e18:.2f}×10¹⁸ mol/Ma")



    print()


def run_u_analysis(args):
    """运行U同位素分析"""
    print(f"\n=== U Isotope Analysis ({args.scenario}) ===\n")

    # 批量处理模式
    if args.file:
        from toolkit.io import BatchProcessor

        print(f"Batch processing mode")
        print(f"Input file: {args.file}")

        processor = BatchProcessor(
            element='u',
            scenario=args.scenario,
            apply_diagenetic_correction=not args.no_diagenetic_correction,
            delta_diag=args.delta_diag,
            include_uncertainty=not args.no_uncertainty,
            n_monte_carlo=args.n_samples if args.n_samples < 5000 else 1000
        )

        try:
            results_df = processor.process_file(
                args.file,
                output_path=args.output,
                show_progress=True
            )

            # 显示摘要
            print("\n" + "=" * 50)
            print("Processing Summary:")
            print("=" * 50)
            print(f"Total samples: {len(results_df)}")
            success_count = results_df['processing_success'].sum()
            print(f"Successful: {success_count}")
            print(f"Failed: {len(results_df) - success_count}")

            if 'f_anox' in results_df.columns:
                print(f"\nf_anox statistics:")
                print(f"  Mean: {results_df['f_anox'].mean():.1%}")
                print(f"  Range: [{results_df['f_anox'].min():.1%}, {results_df['f_anox'].max():.1%}]")

            return
        except Exception as e:
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()
            return

    # 单点计算模式
    from systems.u import UIsotopeSystem

    system = UIsotopeSystem(scenario=args.scenario)
    info = system.get_model_info()

    print(f"Model type: {info['model_type']}")
    print(f"Reservoir mass: {info['reservoir_mass']:.2e} mol")
    print(f"Residence time: {info['residence_time']:.2f} Myr")

    if args.steady_state and args.delta_carb is not None:
        # 稳态计算
        print(f"\n--- Steady-State Calculation ---")

        apply_diag = not args.no_diagenetic_correction
        result = system.calculate_f_anox_steady_state(
            delta238_carb=args.delta_carb,
            apply_diagenetic_correction=apply_diag,
            delta_diag=args.delta_diag
        )

        # === 输入参数 ===
        print(f"\n[Input Parameters]")
        print(f"  Measured δ²³⁸U_carb:     {args.delta_carb:+.2f}‰")
        print(f"  Diagenetic correction:   {'Yes' if apply_diag else 'No'}")
        if apply_diag:
            print(f"    Δ_diag:                {result['delta_diag']:+.2f}‰")
            print(f"    Corrected δ²³⁸U_carb:  {result['delta238_carb_corrected']:+.2f}‰")

        # === 同位素结果 ===
        print(f"\n[Isotope Results]")
        print(f"  Seawater δ²³⁸U:          {result['delta238_seawater']:+.2f}‰")
        print(f"  Oxic sink δ²³⁸U:         {result['delta238_oxic_sink']:+.2f}‰")
        print(f"  Anoxic sink δ²³⁸U:       {result['delta238_anoxic_sink']:+.2f}‰")

        # === 汇比例 ===
        print(f"\n[Sink Fractions]")
        print(f"  f_anox (anoxic):         {result['f_anox']:+.1%}")
        print(f"  f_oxic (oxic):           {result['f_oxic']:+.1%}")

        # === 缺氧面积估算 ===
        anoxic_area = system.estimate_anoxic_area(result['f_anox'])
        print(f"\n[Anoxic Seafloor Area]")
        print(f"  Estimated:               ~{anoxic_area:.1f}%")
        print(f"  (Based on Tissot & Dauphas 2015; Song et al. 2017)")

        # 不确定度分析
        if args.uncertainty or args.sensitivity_analysis:
            from systems.u import UncertaintyAnalyzer
            analyzer = UncertaintyAnalyzer(system)

            if args.sensitivity_analysis:
                print("\n[Sensitivity Analysis]")
                sens_result = analyzer.sensitivity_analysis(args.delta_carb)

                print(f"  Baseline f_anox:         {sens_result['baseline_f_anox']:+.1%}")
                print("  Parameter sensitivities (ranked):")
                for item in sens_result['tornado_data'][:5]:
                    print(f"    {item['parameter']:15}: range [{item['min_effect']:+.1%}, {item['max_effect']:+.1%}]")

            if args.uncertainty in ['mc', 'monte-carlo']:
                print(f"\n[Monte Carlo Uncertainty Analysis]")
                print(f"  Samples:                 {args.n_samples}")
                print(f"  Confidence level:        {args.confidence_level*100:.0f}%")
                print("\n  Uncertain parameters:")
                print(f"    Δ_sw-anox:             0.77 ± 0.04 ‰")
                print(f"    Δ_diag:                0.40 (range 0.30-0.50) ‰")
                print(f"    δ_river:               -0.29 ± 0.16 ‰")
                print(f"    Measurement:           ±{args.measurement_std} ‰")

                mc_result = analyzer.monte_carlo_steady_state(
                    delta238_carb=args.delta_carb,
                    n_samples=args.n_samples,
                    confidence_level=args.confidence_level,
                    apply_diagenetic_correction=apply_diag
                )

                ci_lower, ci_upper = mc_result['f_anox_ci']
                print(f"\n  [Uncertainty Results]")
                print(f"    f_anox mean:           {mc_result['f_anox_mean']:+.1%}")
                print(f"    f_anox median:         {mc_result['f_anox_median']:+.1%}")
                print(f"    Std dev:               {mc_result['f_anox_std']:+.1%}")
                print(f"    {args.confidence_level*100:.0f}% CI:                [{ci_lower:+.1%}, {ci_upper:+.1%}]")

                # 添加缺氧面积的不确定度估算
                area_mean = system.estimate_anoxic_area(mc_result['f_anox_mean'])
                area_lower = system.estimate_anoxic_area(ci_lower)
                area_upper = system.estimate_anoxic_area(ci_upper)
                print(f"\n    Anoxic area:           ~{area_mean:.1f}%")
                print(f"    Area {args.confidence_level*100:.0f}% CI:          [{area_lower:.1f}%, {area_upper:.1f}%]")

                if mc_result['convergence']['converged']:
                    print(f"\n  ✓ MCMC converged (R̂ = {mc_result['convergence']['r_hat']:.3f})")
                else:
                    print(f"\n  ⚠ Consider increasing n_samples (R̂ = {mc_result['convergence']['r_hat']:.3f})")

    elif args.transient:
        # 非稳态模拟
        print(f"\n--- Transient Model Simulation ---")
        print(f"Event duration: {args.event_duration} Myr")
        print(f"Peak f_anox: {args.peak_f_anox}")
        print(f"Background f_anox: {args.background_f_anox}")

        result = system.simulate_anoxic_event(
            event_duration=args.event_duration,
            peak_f_anox=args.peak_f_anox,
            background_f_anox=args.background_f_anox
        )

        if result.success:
            data = result.data
            delta_sw = data['delta_seawater']
            f_anox = data['f_anox']
            time = data['time_myr']

            # 找到最小值和事件结束后的恢复情况
            min_idx = np.argmin(delta_sw)
            end_idx = -1

            print(f"\nResults:")
            print(f"  Background δ²³⁸U_sw: {delta_sw[0]:+.2f}‰")
            print(f"  Minimum δ²³⁸U_sw:    {delta_sw[min_idx]:+.2f}‰ (at t={time[min_idx]:.2f} Myr)")
            print(f"  Final δ²³⁸U_sw:      {delta_sw[end_idx]:+.2f}‰")
            print(f"  Excursion amplitude: {delta_sw[min_idx] - delta_sw[0]:.2f}‰")

            # 检查是否达到平衡
            equil = system.solve_equilibration_time(
                target_f_anox=args.peak_f_anox,
                initial_f_anox=args.background_f_anox
            )
            print(f"\n  Equilibration time: ~{equil['equilibration_time']:.2f} Myr")
            print(f"  Theoretical minimum: {equil['final_delta']:+.2f}‰")
        else:
            print(f"\nError: {result.message}")

    else:
        # 显示示例计算
        print("\n--- Example Calculations ---")

        # 示例1: 现代条件
        print("\n1. Modern oxic ocean:")
        result1 = system.calculate_f_anox_steady_state(-0.45)
        print(f"   δ²³⁸U_carb = -0.45‰ → f_anox = {result1['f_anox']:.1%}")

        # 示例2: 缺氧条件
        print("\n2. Anoxic event:")
        result2 = system.calculate_f_anox_steady_state(-0.85)
        print(f"   δ²³⁸U_carb = -0.85‰ → f_anox = {result2['f_anox']:.1%}")

        # 示例3: 反向计算
        print("\n3. Reverse calculation (f_anox → δ²³⁸U):")
        result3 = system.calculate_seawater_delta_steady_state(f_anox=0.5)
        print(f"   f_anox = 50% → δ²³⁸U_seawater = {result3['delta238_seawater']:+.2f}‰")

    print()


def run_n_analysis(args):
    """
    运行N同位素分析

    基于Kang et al. (2023)和Ma et al. (2025)的双箱稳态模型
    """
    print(f"\n=== N Isotope Analysis ({args.scenario}) ===\n")

    from systems.n import NIsotopeSystem

    system = NIsotopeSystem(scenario=args.scenario)

    # 显示模型信息
    info = system.get_model_info()
    print(f"Model: {info['name']} Isotope System")
    print(f"Scenario: {args.scenario}")
    print(f"Reference: Atmospheric N₂ (δ¹⁵N = 0‰)")
    print()

    # ===== 批量处理模式 =====
    if args.file:
        print("[Batch Processing Mode]")
        print("-" * 50)
        print(f"  Input file: {args.file}")
        print(f"  Scenario: {args.scenario}")
        print(f"  Delta15N column: {args.column}")
        print(f"  Measurement uncertainty (1σ): {args.delta_std}‰")
        print(f"  Epsilon_fix uncertainty: ±{args.epsilon_fix_std}‰")
        print(f"  Epsilon_wcd uncertainty: ±{args.epsilon_wcd_std}‰")
        print(f"  Monte Carlo samples: {args.n_monte_carlo}")
        print()

        from toolkit.io import BatchProcessor

        try:
            # 创建批量处理器（启用不确定度计算）
            processor = BatchProcessor(
                element='n',
                scenario=args.scenario,
                include_uncertainty=True,
                n_monte_carlo=args.n_monte_carlo,
                epsilon_fix_std=args.epsilon_fix_std,
                epsilon_wcd_std=args.epsilon_wcd_std
            )

            # 处理文件
            results_df = processor.process_file(
                args.file,
                output_path=args.output,
                show_progress=True
            )

            # 显示摘要
            print("\n" + "=" * 50)
            print("Processing Summary:")
            print("=" * 50)
            print(f"Total samples: {len(results_df)}")

            if 'f_assimilator' in results_df.columns:
                success_count = (~results_df['f_assimilator'].isna()).sum()
                print(f"Successful: {success_count}")
                print(f"Failed: {len(results_df) - success_count}")
                print(f"\nf_assimilator statistics:")
                print(f"  Mean: {results_df['f_assimilator'].mean():.3f}")
                print(f"  Range: [{results_df['f_assimilator'].min():.3f}, {results_df['f_assimilator'].max():.3f}]")

                # 显示不确定度统计（如果有）
                if 'f_assimilator_std' in results_df.columns:
                    mean_std = results_df['f_assimilator_std'].mean()
                    print(f"\nUncertainty (Monte Carlo, n={args.n_monte_carlo}):")
                    print(f"  Mean 1σ uncertainty: ±{mean_std:.3f}")
                    print(f"  Typical 95% CI width: ±{mean_std * 2:.3f}")

            if args.output:
                print(f"\nResults saved to: {args.output}")
                print("\nOutput columns:")
                print("  - f_assimilator: 反演结果")
                print("  - f_assimilator_std: 标准差 (1σ)")
                print("  - f_assimilator_ci68: 68% 置信区间")
                print("  - f_assimilator_ci95: 95% 置信区间")

        except Exception as e:
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()

        print()
        return

    # ===== 关系曲线计算 =====
    if args.curve:
        print("[f_assimilator vs δ¹⁵N Curve Calculation]")
        print("-" * 50)

        f_min, f_max = args.f_range
        print(f"  Range: f_assimilator = [{f_min:.2f}, {f_max:.2f}]")
        print(f"  Points: {args.n_points}")

        if args.monte_carlo:
            print(f"  Monte Carlo samples: {args.n_samples}")
        print()

        curve = system.calculate_f_assimilator_curve(
            f_range=(f_min, f_max),
            n_points=args.n_points,
            n_monte_carlo=args.n_samples if args.monte_carlo else 1
        )

        print("  Results (first 5 points):")
        print(f"  {'f_assimilator':<15} {'δ¹⁵N_sed':<12}")
        print("  " + "-" * 27)
        for i in range(min(5, len(curve['f_assimilator']))):
            f = curve['f_assimilator'][i]
            delta = curve['delta15N_sed_mean'][i]
            print(f"  {f:<15.3f} {delta:<+12.2f}")

        if len(curve['f_assimilator']) > 5:
            print(f"  ... and {len(curve['f_assimilator'])-5} more points")

        # 找到峰值点
        max_idx = np.argmax(curve['delta15N_sed_mean'])
        print(f"\n  Peak δ¹⁵N: {curve['delta15N_sed_mean'][max_idx]:+.2f}‰")
        print(f"  At f_assimilator: {curve['f_assimilator'][max_idx]:.3f}")

        # 保存结果
        if args.output:
            import pandas as pd
            df = pd.DataFrame({
                'f_assimilator': curve['f_assimilator'],
                'delta15N_sed_mean': curve['delta15N_sed_mean']
            })
            if args.monte_carlo:
                df['delta15N_sed_ci68_lower'] = curve['delta15N_sed_ci68_lower']
                df['delta15N_sed_ci68_upper'] = curve['delta15N_sed_ci68_upper']
                df['delta15N_sed_ci95_lower'] = curve['delta15N_sed_ci95_lower']
                df['delta15N_sed_ci95_upper'] = curve['delta15N_sed_ci95_upper']
            df.to_csv(args.output, index=False)
            print(f"\n  Results saved to: {args.output}")

        print()
        return

    # ===== 反向模型：从 delta15N 反演 f_assimilator =====
    if args.inverse and args.delta15N is not None:
        print("[Inverse Model: δ¹⁵N → f_assimilator]")
        print("-" * 50)

        delta15N = args.delta15N
        f_min, f_max = args.f_range

        print(f"  Input δ¹⁵N: {delta15N:+.2f}‰")
        print(f"  Search range: f_assimilator = [{f_min:.2f}, {f_max:.2f}]")
        print()

        result = system.inverse_model(
            delta15N_sed=delta15N,
            f_range=(f_min, f_max)
        )

        print(f"  Results:")
        print(f"    f_assimilator: {result['f_assimilator']:.4f}")
        print(f"    Calculated δ¹⁵N: {result['delta15N_sed_calculated']:+.2f}‰")
        print(f"    Residual: {result['residual']:+.4f}‰")

        # 解释
        f = result['f_assimilator']
        if f < 0.1:
            interp = "Nitrate-depleted (anoxic)"
        elif f < 0.3:
            interp = "Nitrate-limited"
        elif f < 0.5:
            interp = "Moderate nitrate availability"
        else:
            interp = "Nitrate-sufficient (oxic)"
        print(f"    Interpretation: {interp}")

        # 蒙特卡洛不确定性
        if args.monte_carlo:
            print(f"\n  Monte Carlo Uncertainty Analysis ({args.n_samples} samples)...")
            mc_result = system.monte_carlo_simulation(
                f_assimilator=f,
                n_samples=args.n_samples
            )
            print(f"    δ¹⁵N_sed: {mc_result['delta15N_sed_mean']:.2f} ± {mc_result['delta15N_sed_std']:.2f}‰")
            print(f"    95% CI: [{mc_result['delta15N_sed_ci95'][0]:.2f}, {mc_result['delta15N_sed_ci95'][1]:.2f}]‰")

        print()
        return

    # ===== 正向模型：从 f_assimilator 计算 delta15N =====
    if args.f_assimilator is not None:
        print("[Forward Model: f_assimilator → δ¹⁵N]")
        print("-" * 50)

        f = args.f_assimilator
        print(f"  Input f_assimilator: {f:.4f}")
        print()

        # 计算储库同位素
        reservoirs = system.calculate_reservoir_isotopes(f)
        delta15N = system.forward_model(f)

        print(f"  Results:")
        print(f"    δ¹⁵N_ammonium: {reservoirs['delta15N_ammonium']:+.2f}‰")
        print(f"    δ¹⁵N_nitrate: {reservoirs['delta15N_nitrate']:+.2f}‰")
        print(f"    δ¹⁵N_sediment: {delta15N:+.2f}‰")

        # 蒙特卡洛不确定性
        if args.monte_carlo:
            print(f"\n  Monte Carlo Uncertainty Analysis ({args.n_samples} samples)...")
            mc_result = system.monte_carlo_simulation(
                f_assimilator=f,
                n_samples=args.n_samples
            )
            print(f"    Mean δ¹⁵N_sed: {mc_result['delta15N_sed_mean']:.2f}‰")
            print(f"    Std dev: {mc_result['delta15N_sed_std']:.2f}‰")
            print(f"    95% CI: [{mc_result['delta15N_sed_ci95'][0]:.2f}, {mc_result['delta15N_sed_ci95'][1]:.2f}]‰")

        print()
        return

    # ===== 默认：显示示例 =====
    print("[Example Calculations]")
    print("-" * 50)

    f_values = [0.0, 0.11, 0.25, 0.48, 0.7, 1.0]

    print(f"  {'f_assimilator':<15} {'δ¹⁵N_sed':<12} {'Nitrate Status'}")
    print("  " + "-" * 50)

    for f in f_values:
        delta15N = system.forward_model(f)
        if f < 0.1:
            status = "Depleted (Anoxic)"
        elif f < 0.3:
            status = "Limited"
        elif f < 0.5:
            status = "Moderate"
        else:
            status = "Sufficient (Oxic)"
        print(f"  {f:<15.2f} {delta15N:<+12.2f} {status}")

    print(f"\n  Note: Peak δ¹⁵N occurs at f ≈ 0.48")
    print(f"        Lower f = more water-column denitrification")
    print(f"        Higher f = more nitrate assimilation")
    print()


def run_mg_siliciclastic_analysis(args):
    """Run the Hu et al. (2023) siliciclastic Mg-isotope workflow."""
    print("\n" + "="*70)
    print("Mg Isotope Weathering Flux Model")
    print("Component Type: SILICICLASTIC (陆源碎屑沉积物)")
    print("Based on: Hu et al. (2023) Global and Planetary Change")
    print("="*70 + "\n")

    from systems.mg import create_mg_system
    import pandas as pd

    system = create_mg_system('siliciclastic', basin=args.basin)
    params = system.params

    # Explicit siliciclastic options take precedence. The old
    # --delta-silicate flag is retained as a protolith alias.
    protolith_override = args.delta_protolith
    if protolith_override is None and args.delta_silicate is not None:
        protolith_override = args.delta_silicate
    if protolith_override is not None and args.delta_protolith_range is not None:
        raise ValueError(
            "Use either --delta-protolith or --delta-protolith-range, not both"
        )
    if protolith_override is not None:
        params.d26Mg_UCC = protolith_override
        params.d26Mg_UCC_range = None
    elif args.delta_protolith_range is not None:
        params.d26Mg_UCC_range = tuple(args.delta_protolith_range)
    if args.delta_fluid_protolith is not None:
        params.Delta_fluid_protolith = args.delta_fluid_protolith
    if args.delta_river_mg is not None:
        params.d26Mg_river_water = args.delta_river_mg
    if args.delta_carbonate is not None:
        params.d26Mg_carbonate = args.delta_carbonate
    if args.river_mg_flux is not None:
        params.F_river_total = args.river_mg_flux
    params.uncertainty_river_flux_fraction = args.river_mg_flux_rel_std
    system._validate_parameters()

    print(f"[System Configuration]")
    print(f"  Basin type: {args.basin}")
    print(f"  Protolith δ²⁶Mg reference: {params.d26Mg_UCC:+.2f}‰")
    if params.d26Mg_UCC_range is None:
        print("  Protolith Monte Carlo prior: fixed at reference value")
    else:
        print(
            f"  Protolith Monte Carlo prior: Uniform("
            f"{params.d26Mg_UCC_range[0]:+.2f}, "
            f"{params.d26Mg_UCC_range[1]:+.2f})‰; shared across profile"
        )
    print(f"  Δ(fluid-protolith): {params.Delta_fluid_protolith:+.2f}‰")
    print(f"  River water δ²⁶Mg: {params.d26Mg_river_water:+.2f}‰")
    print(f"  Carbonate flux end-member: {params.d26Mg_carbonate:+.2f}‰")
    print(f"  Total river Mg flux prior: {params.F_river_total:.3e} mol/yr")
    print("  Absolute fluxes are conditional on the total river Mg flux prior.")
    print()

    def _format_percent(value):
        return "n/a" if not np.isfinite(value) else f"{value:.1%}"

    def _format_flux(value):
        return "n/a" if not np.isfinite(value) else f"{value:.3e} mol/yr"

    # ===== 单点计算模式 =====
    if args.delta_sample is not None:
        print("[Single Sample Analysis]")
        print("-" * 50)
        d26Mg_clay = args.delta_sample
        print(f"  Input δ²⁶Mg (clay): {d26Mg_clay:.2f}‰")
        result = system.calculate_weathering_flux(d26Mg_clay)

        if result.success:
            print(f"  Results:")
            print(f"    Model status: {result.model_status}")
            print(f"    Weathering stage: {result.weathering_stage}")
            print(f"    Retained Mg (f_Mg): {_format_percent(result.f_Mg)}")
            print(f"    Mg loss fraction: {_format_percent(result.Mg_loss_fraction)}")
            print(f"    Cumulative silicate flux δ²⁶Mg: {result.d26Mg_silicate:+.3f}‰")
            print(f"    Silicate fraction of dissolved Mg: {_format_percent(result.silicate_flux_fraction)}")
            print(f"    Conditional silicate Mg flux: {_format_flux(result.F_silicate)}")
            if result.message:
                print(f"    Caution: {result.message}")

            mc = system.monte_carlo_analysis(
                d26Mg_clay,
                n_iterations=args.mg_monte_carlo,
                random_seed=args.random_seed,
            )
            ci = mc['F_silicate']['ci_95']
            print(f"    Monte Carlo median F_silicate: {_format_flux(mc['F_silicate']['median'])}")
            if np.all(np.isfinite(ci)):
                print(f"    95% interval: [{ci[0]:.3e}, {ci[1]:.3e}] mol/yr")
            print(f"    Valid Monte Carlo fraction: {mc['valid_flux_fraction']:.1%}")
        else:
            print(f"  Error: {result.message}")
        print()

    # ===== 批量处理模式 =====
    if args.file:
        print("[Batch Processing Mode]")
        print("-" * 50)
        print(f"  Input file: {args.file}")

        try:
            if args.file.lower().endswith('.csv'):
                df = pd.read_csv(args.file)
            else:
                df = pd.read_excel(args.file)

            # Remove spreadsheet padding while preserving all analytical data.
            df = df.dropna(axis=1, how='all')
            unnamed = [c for c in df.columns if str(c).startswith('Unnamed:')]
            df = df.drop(columns=unnamed, errors='ignore')
            n_samples = len(df)
            if n_samples == 0:
                raise ValueError("Input file contains no samples")

            iso_col = args.column if args.column in df.columns else None
            if iso_col is None:
                candidates = ['delta_26_Mg_iso', 'd26Mg', 'δ26Mg', 'delta_26Mg']
                iso_col = next((c for c in candidates if c in df.columns), None)
            if iso_col is None:
                raise ValueError("No supported delta26Mg column found")

            df.insert(0, 'sample_index', np.arange(1, n_samples + 1))
            relative_depth = (
                np.zeros(1)
                if n_samples == 1
                else np.linspace(0.0, 1.0, n_samples)
            )
            df.insert(1, 'relative_depth', relative_depth)

            if iso_col == 'delta_25_Mg_iso':
                d26_values = pd.to_numeric(df[iso_col], errors='coerce') / 0.521
                input_note = "delta26Mg inferred from delta25Mg / 0.521"
            else:
                d26_values = pd.to_numeric(df[iso_col], errors='coerce')
                input_note = f"delta26Mg read from {iso_col}"
            df['delta26Mg_model_input'] = d26_values

            std_col = next(
                (c for c in ['delta_26_Mg_iso_2sd', 'd26Mg_2sd', 'delta26Mg_2sd'] if c in df.columns),
                None,
            )
            if std_col is not None:
                d26_std = pd.to_numeric(df[std_col], errors='coerce') / 2.0
            else:
                d26_std = pd.Series(params.uncertainty_clay, index=df.index)

            # Mass-dependent isotope quality control. The residual and z score
            # are reported rather than silently removing samples.
            if {'delta_25_Mg_iso', 'delta_26_Mg_iso'}.issubset(df.columns):
                residual = (
                    pd.to_numeric(df['delta_25_Mg_iso'], errors='coerce')
                    - 0.521 * pd.to_numeric(df['delta_26_Mg_iso'], errors='coerce')
                )
                df['d25_d26_residual'] = residual
                if {'delta_25_Mg_iso_2sd', 'delta_26_Mg_iso_2sd'}.issubset(df.columns):
                    sigma25 = pd.to_numeric(df['delta_25_Mg_iso_2sd'], errors='coerce') / 2.0
                    sigma26 = pd.to_numeric(df['delta_26_Mg_iso_2sd'], errors='coerce') / 2.0
                    residual_sigma = np.sqrt(sigma25**2 + (0.521 * sigma26)**2)
                    df['d25_d26_zscore'] = residual / residual_sigma
                    df['mass_dependent_qc'] = df['d25_d26_zscore'].abs() <= 2.0

            rows = []
            for idx, row in df.iterrows():
                d26Mg = row['delta26Mg_model_input']
                if not np.isfinite(d26Mg):
                    rows.append({'calculation_success': False, 'message': 'Missing isotope data'})
                    continue

                result = system.calculate_weathering_flux(float(d26Mg))
                sample_std = d26_std.iloc[idx]
                if not np.isfinite(sample_std):
                    sample_std = params.uncertainty_clay
                mc = system.monte_carlo_analysis(
                    float(d26Mg),
                    d26Mg_clay_std=float(sample_std),
                    n_iterations=args.mg_monte_carlo,
                    random_seed=args.random_seed + int(idx),
                )
                f_ci = mc['f_Mg']['ci_95']
                flux_ci = mc['F_silicate']['ci_95']
                rows.append({
                    'calculation_success': result.success,
                    'model_status': result.model_status,
                    'within_hu_calibration': result.within_hu_calibration,
                    'weathering_stage': result.weathering_stage,
                    'delta26Mg_model_1sigma': sample_std,
                    'delta26Mg_protolith_reference': params.d26Mg_UCC,
                    'delta26Mg_protolith_prior_distribution': (
                        'fixed' if params.d26Mg_UCC_range is None else 'uniform'
                    ),
                    'delta26Mg_protolith_prior_low': (
                        params.d26Mg_UCC
                        if params.d26Mg_UCC_range is None
                        else params.d26Mg_UCC_range[0]
                    ),
                    'delta26Mg_protolith_prior_high': (
                        params.d26Mg_UCC
                        if params.d26Mg_UCC_range is None
                        else params.d26Mg_UCC_range[1]
                    ),
                    'delta26Mg_protolith_prior_draw_median': mc['d26Mg_protolith']['median'],
                    'Delta_fluid_protolith': params.Delta_fluid_protolith,
                    'delta26Mg_river_prior': params.d26Mg_river_water,
                    'delta26Mg_carbonate_prior': params.d26Mg_carbonate,
                    'f_Mg': result.f_Mg,
                    'Mg_loss_fraction': result.Mg_loss_fraction,
                    'd26Mg_silicate_flux': result.d26Mg_silicate,
                    'silicate_flux_fraction': result.silicate_flux_fraction,
                    'F_silicate_mol_yr': result.F_silicate,
                    'F_carbonate_mol_yr': result.F_carbonate,
                    'F_river_total_prior_mol_yr': params.F_river_total,
                    'f_Mg_mc_median': mc['f_Mg']['median'],
                    'f_Mg_mc_ci95_low': f_ci[0],
                    'f_Mg_mc_ci95_high': f_ci[1],
                    'F_silicate_mc_median_mol_yr': mc['F_silicate']['median'],
                    'F_silicate_mc_ci95_low_mol_yr': flux_ci[0],
                    'F_silicate_mc_ci95_high_mol_yr': flux_ci[1],
                    'mc_valid_flux_fraction': mc['valid_flux_fraction'],
                    'mc_iterations': args.mg_monte_carlo,
                    'absolute_flux_is_conditional': True,
                    'message': result.message,
                })

            results_df = pd.DataFrame(rows, index=df.index)

            if args.baseline_samples < 1:
                raise ValueError("--baseline-samples must be at least 1")
            baseline_count = min(args.baseline_samples, n_samples)
            profile_mc = system.monte_carlo_profile(
                d26_values.to_numpy(dtype=float),
                d26Mg_clay_std=d26_std.to_numpy(dtype=float),
                baseline_count=baseline_count,
                baseline_position='start' if args.baseline_position == 'shallow' else 'end',
                n_iterations=args.mg_monte_carlo,
                random_seed=args.random_seed,
            )
            results_df['Mg_release_fraction_mc_median'] = profile_mc['Mg_release_fraction_median']
            results_df['Mg_release_fraction_mc_ci95_low'] = profile_mc['Mg_release_fraction_ci95_low']
            results_df['Mg_release_fraction_mc_ci95_high'] = profile_mc['Mg_release_fraction_ci95_high']
            results_df['Mg_release_fraction_mc_valid_fraction'] = profile_mc['valid_Mg_release_fraction']
            results_df['conditional_F_silicate_profile_mc_median_mol_yr'] = profile_mc['conditional_F_silicate_median']
            results_df['conditional_F_silicate_profile_mc_ci95_low_mol_yr'] = profile_mc['conditional_F_silicate_ci95_low']
            results_df['conditional_F_silicate_profile_mc_ci95_high_mol_yr'] = profile_mc['conditional_F_silicate_ci95_high']
            results_df['conditional_F_silicate_profile_mc_valid_fraction'] = profile_mc['valid_conditional_F_silicate']
            results_df['weathering_flux_assumes_constant_rock_supply_and_Mg'] = True

            change = profile_mc['change_point']
            split_ci = change['after_sample_ci95']
            shallow_ci = change['shallow_state']['ci_95']
            deep_ci = change['deep_state']['ci_95']
            ratio_ci = change['deep_to_shallow_ratio']['ci_95']
            improvement_ci = change['cost_improvement']['ci_95']
            profile_summary = {
                'profile_protolith_prior_draw_median': profile_mc['d26Mg_protolith']['median'],
                'profile_protolith_prior_draw_ci95_low': profile_mc['d26Mg_protolith']['ci_95'][0],
                'profile_protolith_prior_draw_ci95_high': profile_mc['d26Mg_protolith']['ci_95'][1],
                'protolith_is_shared_across_profile': profile_mc['protolith_is_shared_across_profile'],
                'automatic_change_point_after_sample': change['after_sample_best'],
                'automatic_change_point_posterior_mode': change['after_sample_mode'],
                'automatic_change_point_posterior_median': change['after_sample_median'],
                'automatic_change_point_ci95_low': split_ci[0],
                'automatic_change_point_ci95_high': split_ci[1],
                'automatic_change_point_valid_fraction': change['n_valid'] / args.mg_monte_carlo,
                'shallow_Mg_release_fraction_mc_median': change['shallow_state']['median'],
                'shallow_Mg_release_fraction_mc_ci95_low': shallow_ci[0],
                'shallow_Mg_release_fraction_mc_ci95_high': shallow_ci[1],
                'deep_Mg_release_fraction_mc_median': change['deep_state']['median'],
                'deep_Mg_release_fraction_mc_ci95_low': deep_ci[0],
                'deep_Mg_release_fraction_mc_ci95_high': deep_ci[1],
                'deep_to_shallow_weathering_flux_ratio_median': change['deep_to_shallow_ratio']['median'],
                'deep_to_shallow_weathering_flux_ratio_ci95_low': ratio_ci[0],
                'deep_to_shallow_weathering_flux_ratio_ci95_high': ratio_ci[1],
                'probability_deep_weathering_flux_greater': change['probability_deep_greater'],
                'change_point_l1_cost_improvement': change['central_cost_improvement'],
                'change_point_l1_cost_improvement_mc_median': change['cost_improvement']['median'],
                'change_point_l1_cost_improvement_mc_ci95_low': improvement_ci[0],
                'change_point_l1_cost_improvement_mc_ci95_high': improvement_ci[1],
            }
            for column, value in profile_summary.items():
                results_df[column] = value

            if (args.time_start is None) != (args.time_end is None):
                raise ValueError("--time-start and --time-end must be supplied together")
            if args.time_start is not None:
                df['model_time'] = np.linspace(args.time_start, args.time_end, n_samples)
                plot_y_column = 'model_time'
                plot_y_label = f"Time ({args.time_unit})"
                time_axis_note = (
                    f"linear time anchors {args.time_start:g} to {args.time_end:g} "
                    f"{args.time_unit}"
                )
            else:
                df['relative_stratigraphic_time'] = relative_depth
                plot_y_column = 'relative_stratigraphic_time'
                plot_y_label = 'Relative stratigraphic position'
                time_axis_note = "relative row-order axis (no absolute ages)"

            output_df = pd.concat([df, results_df], axis=1)

            valid_flux = output_df[
                'conditional_F_silicate_profile_mc_median_mol_yr'
            ].dropna()
            status_counts = output_df['model_status'].value_counts(dropna=False)
            print(f"  Loaded {n_samples} samples; row order treated as increasing depth")
            print(f"  Model input: {input_note}")
            print("  Deterministic status counts at the reference protolith:")
            for status, count in status_counts.items():
                print(f"    {status}: {count}")
            if len(valid_flux):
                print(f"  Conditional F_silicate median range: {valid_flux.min():.3e} - {valid_flux.max():.3e} mol/yr")
                x = output_df.loc[valid_flux.index, 'relative_depth'].to_numpy()
                if len(valid_flux) > 1:
                    slope = np.polyfit(x, valid_flux.to_numpy(), 1)[0]
                    print(f"  Flux trend over normalized depth: {slope:+.3e} mol/yr per relative-depth unit")

            release = output_df['Mg_release_fraction_mc_median']
            reliable_release = output_df[
                np.isfinite(release)
                & (output_df['Mg_release_fraction_mc_valid_fraction'] >= 0.5)
            ]
            print(f"  Plot time axis: {time_axis_note}")
            print("  Flux proxy: Mg release fraction (1 - f_Mg); no manual baseline")
            low_valid_count = int(
                (output_df['Mg_release_fraction_mc_valid_fraction'] < 0.5).sum()
            )
            print(f"  Samples with <50% physical Mg-release draws: {low_valid_count}")
            if len(reliable_release):
                print(
                    f"  Reliable Mg-release range: {reliable_release['Mg_release_fraction_mc_median'].min():.3f}"
                    f"-{reliable_release['Mg_release_fraction_mc_median'].max():.3f}"
                )
            if change['after_sample_best'] is not None:
                split = int(change['after_sample_best'])
                print(
                    f"  Automatic transition: between samples {split} and {split + 1} "
                    f"(per-draw split 95% interval after samples "
                    f"{split_ci[0]:.0f}-{split_ci[1]:.0f})"
                )
                print(
                    f"  Shallow/deep Mg-release states: "
                    f"{change['shallow_state']['median']:.3f} / "
                    f"{change['deep_state']['median']:.3f}"
                )
                print(
                    f"  Deep/shallow weathering-flux ratio: "
                    f"{change['deep_to_shallow_ratio']['median']:.2f}x "
                    f"(95% CI {ratio_ci[0]:.2f}-{ratio_ci[1]:.2f}x)"
                )
                print("  Direction labels describe row order, not absolute geological time.")

            if args.output:
                output_path = Path(args.output)
                output_path.parent.mkdir(parents=True, exist_ok=True)
                if output_path.suffix.lower() in ('.xlsx', '.xls'):
                    output_df.to_excel(output_path, index=False)
                else:
                    output_df.to_csv(output_path, index=False)
                print(f"\n  Results saved to: {args.output}")

            if args.plot or args.plot_output or args.flux_plot_output:
                if args.plot_output:
                    plot_path = Path(args.plot_output)
                elif args.output:
                    plot_path = Path(args.output).with_suffix('.png')
                else:
                    plot_path = Path('results') / f"{Path(args.file).stem}_Mg_release.png"
                plot_path.parent.mkdir(parents=True, exist_ok=True)
                plot_mg_flux_profile(
                    output_df,
                    y_column=plot_y_column,
                    y_label=plot_y_label,
                    output_path=plot_path,
                    title=f"{Path(args.file).stem.replace('_', ' ')}: Mg weathering proxy",
                )
                print(f"  Portrait Mg-release plot saved to: {plot_path}")

                if args.flux_plot_output:
                    flux_plot_path = Path(args.flux_plot_output)
                elif args.output:
                    flux_plot_path = Path(args.output).with_name(
                        f"{Path(args.output).stem}_Changjiang_flux.png"
                    )
                else:
                    flux_plot_path = Path('results') / (
                        f"{Path(args.file).stem}_Changjiang_flux.png"
                    )
                flux_plot_path.parent.mkdir(parents=True, exist_ok=True)
                plot_mg_conditional_flux_profile(
                    output_df,
                    y_column=plot_y_column,
                    y_label=plot_y_label,
                    output_path=flux_plot_path,
                    title=(
                        f"{Path(args.file).stem.replace('_', ' ')}: "
                        "conditional silicate Mg flux"
                    ),
                )
                print(
                    f"  Portrait Changjiang-scaled flux plot saved to: "
                    f"{flux_plot_path}"
                )
            print()

        except Exception as e:
            print(f"  Error: {e}")
            import traceback
            traceback.print_exc()

    # ===== 如果没有指定操作，显示示例 =====
    if args.delta_sample is None and not args.file:
        print("[Example Calculation - Hu et al., 2023]")
        print("-" * 50)
        print("  Using typical Changjiang River values:")
        example_values = [-0.15, -0.10, -0.05, 0.00]
        print(f"  {'δ²⁶Mg_clay':<12} {'Stage':<15} {'f_Mg':<10} {'F_sili (mol/yr)':<18}")
        print("  " + "-" * 59)
        for d26Mg in example_values:
            result = system.calculate_weathering_flux(d26Mg)
            if result.success:
                print(f"  {d26Mg:<+12.2f} {result.weathering_stage:<15} "
                      f"{result.f_Mg:<10.3f} {result.F_silicate:<18.3e}")
        print()
        print("  Note: Higher δ²⁶Mg_clay indicates stronger weathering")
        print("        (more ²⁴Mg leached, ²⁶Mg enriched in residue)")
        print("        Absolute flux scales with the total river Mg flux prior.")
        print()

    print("="*70)
    print("Analysis complete")
    print("="*70 + "\n")


def plot_mg_flux_profile(df, y_column, y_label, output_path, title):
    """Plot a narrow portrait profile of baseline-free Mg release."""
    import pandas as pd

    protolith_note = _format_protolith_note(df)
    ratio = pd.to_numeric(
        df.get('deep_to_shallow_weathering_flux_ratio_median'), errors='coerce'
    ).dropna()
    ratio_low = pd.to_numeric(
        df.get('deep_to_shallow_weathering_flux_ratio_ci95_low'), errors='coerce'
    ).dropna()
    ratio_high = pd.to_numeric(
        df.get('deep_to_shallow_weathering_flux_ratio_ci95_high'), errors='coerce'
    ).dropna()
    note = None
    if len(ratio) and len(ratio_low) and len(ratio_high):
        note = (
            f"{protolith_note}\n"
            f"Deep / shallow = {ratio.iloc[0]:.2f}×\n"
            f"95% CI: {ratio_low.iloc[0]:.2f}-{ratio_high.iloc[0]:.2f}×"
        )
    _plot_geochemical_profile(
        df=df,
        y_column=y_column,
        y_label=y_label,
        median_column='Mg_release_fraction_mc_median',
        low_column='Mg_release_fraction_mc_ci95_low',
        high_column='Mg_release_fraction_mc_ci95_high',
        valid_column='Mg_release_fraction_mc_valid_fraction',
        x_label=r'Mg release fraction, $1-f_{\mathrm{Mg}}$',
        output_path=output_path,
        title=title,
        x_scale=1.0,
        x_limits=(-0.04, 1.02),
        note=note,
    )


def plot_mg_conditional_flux_profile(df, y_column, y_label, output_path, title):
    """Plot dissolved silicate Mg flux conditional on the Changjiang prior."""
    river_flux = float(
        df['F_river_total_prior_mol_yr'].dropna().iloc[0]
    )
    note = (
        f"{_format_protolith_note(df)}\n"
        f"River Mg prior: {river_flux / 1e11:.2f} × $10^{{11}}$ mol yr$^{{-1}}$\n"
        "Absolute scale is conditional"
    )
    _plot_geochemical_profile(
        df=df,
        y_column=y_column,
        y_label=y_label,
        median_column='conditional_F_silicate_profile_mc_median_mol_yr',
        low_column='conditional_F_silicate_profile_mc_ci95_low_mol_yr',
        high_column='conditional_F_silicate_profile_mc_ci95_high_mol_yr',
        valid_column='conditional_F_silicate_profile_mc_valid_fraction',
        x_label=r'Conditional silicate Mg flux ($10^{11}$ mol yr$^{-1}$)',
        output_path=output_path,
        title=title,
        x_scale=1e11,
        note=note,
    )


def _format_protolith_note(df):
    distribution = str(
        df['delta26Mg_protolith_prior_distribution'].dropna().iloc[0]
    )
    low = float(df['delta26Mg_protolith_prior_low'].dropna().iloc[0])
    high = float(df['delta26Mg_protolith_prior_high'].dropna().iloc[0])
    if distribution == 'fixed':
        return rf'Protolith $\delta^{{26}}$Mg = {low:+.2f}‰'
    return rf'Protolith $\delta^{{26}}$Mg ~ U({low:+.2f}, {high:+.2f})‰'


def _plot_geochemical_profile(
    df,
    y_column,
    y_label,
    median_column,
    low_column,
    high_column,
    valid_column,
    x_label,
    output_path,
    title,
    x_scale=1.0,
    x_limits=None,
    note=None,
):
    """Render a narrow publication-style geochemical profile."""
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.ticker import AutoMinorLocator

    required = [
        y_column,
        median_column,
        low_column,
        high_column,
        valid_column,
    ]
    plot_data = df[required].apply(pd.to_numeric, errors='coerce')
    plot_data = plot_data[np.isfinite(plot_data[y_column])]
    if plot_data.empty:
        raise ValueError("No profile positions are available for plotting")

    y = plot_data[y_column].to_numpy()
    median = plot_data[median_column].to_numpy() / x_scale
    low = plot_data[low_column].to_numpy() / x_scale
    high = plot_data[high_column].to_numpy() / x_scale
    valid_fraction = plot_data[valid_column].to_numpy()
    reliable = (valid_fraction >= 0.5) & np.isfinite(median)
    if not np.any(np.isfinite(median)):
        raise ValueError("No finite profile estimates are available for plotting")

    if x_limits is None:
        finite_limits = np.concatenate(
            [low[np.isfinite(low)], high[np.isfinite(high)]]
        )
        x_min = float(np.min(finite_limits))
        x_max = float(np.max(finite_limits))
        padding = max(0.06 * (x_max - x_min), 0.02 * abs(x_max), 0.01)
        x_limits = (max(0.0, x_min - padding), x_max + padding)

    fig, ax = plt.subplots(figsize=(4.1, 8.2), dpi=300)
    ax.fill_betweenx(
        y,
        np.where(reliable, low, np.nan),
        np.where(reliable, high, np.nan),
        color='#abc8d2',
        alpha=0.48,
        linewidth=0,
    )
    ax.plot(
        np.where(reliable, median, np.nan),
        y,
        color='#135f76',
        linewidth=1.25,
        zorder=2,
    )
    ax.scatter(
        median[reliable],
        y[reliable],
        color='#135f76',
        edgecolor='white',
        linewidth=0.4,
        s=19,
        zorder=3,
    )
    if np.any(~reliable):
        boundary_x = x_limits[0] + 0.025 * (x_limits[1] - x_limits[0])
        invalid_x = np.where(
            np.isfinite(median[~reliable]),
            median[~reliable],
            boundary_x,
        )
        ax.scatter(
            invalid_x,
            y[~reliable],
            facecolor='none',
            edgecolor='#555555',
            linewidth=0.9,
            s=25,
            marker='s',
            label='<50% valid draws',
            zorder=3,
        )

    split = pd.to_numeric(
        df.get('automatic_change_point_after_sample'), errors='coerce'
    ).dropna()
    if len(split):
        split_index = int(split.iloc[0])
        if 0 < split_index < len(y):
            split_y = 0.5 * (y[split_index - 1] + y[split_index])
            ax.axhline(
                split_y,
                color='#9b303b',
                linestyle=(0, (4, 2)),
                linewidth=1.0,
                label=f'Split {split_index}|{split_index + 1}',
            )

    if note:
        ax.text(
            0.04,
            0.025,
            note,
            transform=ax.transAxes,
            fontsize=7.3,
            va='bottom',
            ha='left',
            color='#222222',
        )

    ax.set_xlabel(x_label, fontsize=9)
    ax.set_ylabel(y_label, fontsize=9)
    ax.set_title(title.replace(': ', '\n'), fontsize=9.5, pad=8)
    ax.set_ylim(y[-1], y[0])
    ax.set_xlim(*x_limits)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(
        axis='both',
        which='major',
        direction='in',
        top=True,
        right=True,
        length=4,
        width=0.8,
        labelsize=8,
    )
    ax.tick_params(
        axis='both',
        which='minor',
        direction='in',
        top=True,
        right=True,
        length=2,
        width=0.6,
    )
    for spine in ax.spines.values():
        spine.set_color('#222222')
        spine.set_linewidth(0.8)
    ax.grid(False)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(
            handles,
            labels,
            frameon=False,
            loc='best',
            fontsize=7.2,
            handlelength=2.2,
        )
    fig.tight_layout(pad=0.8)
    fig.savefig(output_path, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def run_sr_analysis(args):
    """
    运行Sr同位素分析
    基于Wang et al. (2021)随机海洋箱模型
    """
    from systems.sr import SrIsotopeSystem, SrFluxConfig
    import pandas as pd
    import numpy as np

    print("\n" + "="*70)
    print("Sr Isotope Ocean Box Model")
    print("Based on: Wang et al. (2021) Earth-Science Reviews")
    print("="*70 + "\n")

    # ===== Dukou剖面数据分析 =====
    if args.dukou_analysis:
        print("[Dukou Section Sr Isotope Analysis]")
        print("-" * 50)

        from toolkit.io.dukou_data import load_dukou_data
        from systems.sr.age_model import INTERVALS

        # 加载数据
        data = load_dukou_data(
            mn_sr_threshold=args.mn_sr_threshold,
            lowess_frac=args.lowess_frac
        )

        print(f"\n数据摘要:")
        print(f"  样本数: {len(data.sample_ids)}")
        print(f"  年龄范围: {data.ages.min():.1f} - {data.ages.max():.1f} Ma")
        print(f"  Sr范围: {data.sr_ratios.min():.5f} - {data.sr_ratios.max():.5f}")

        # 显示三个下降期数据
        print(f"\n三个下降期数据:")
        for interval_key in ['interval_2', 'interval_3', 'interval_4']:
            interval_info = INTERVALS[interval_key]
            age_min, age_max = interval_info['age_range']
            stats = data.get_interval_statistics(age_min, age_max)
            print(f"\n  {interval_info['name']}:")
            print(f"    年龄: {age_max:.1f} - {age_min:.1f} Ma")
            print(f"    样本: {stats['n_samples']}个")
            print(f"    Sr: {stats['sr_range'][0]:.5f} - {stats['sr_range'][1]:.5f}")

        # 保存结果
        if args.output:
            df = pd.DataFrame({
                'Sample_ID': data.sample_ids,
                'Age_Ma': data.ages,
                'Position_m': data.positions,
                'Sr_87_86': data.sr_ratios,
                'Mn_Sr': data.mn_sr,
                'Uncertainty_2s': data.uncertainties
            })

            if data.target_sr is not None:
                # 添加目标曲线
                target_df = pd.DataFrame({
                    'Target_Age_Ma': data.target_ages,
                    'Target_Sr': data.target_sr,
                    'CI_Lower': data.target_ci_lower,
                    'CI_Upper': data.target_ci_upper
                })

                # 保存到Excel的两个sheet
                with pd.ExcelWriter(args.output) as writer:
                    df.to_excel(writer, sheet_name='Data', index=False)
                    target_df.to_excel(writer, sheet_name='Target_Curve', index=False)
            else:
                df.to_csv(args.output, index=False)

            print(f"\n结果保存到: {args.output}")

        print()
        return

    # ===== 三步情景分析 =====
    if args.scenario_a or args.scenario_b or args.scenario_c:
        from systems.sr.scenarios import ScenarioManager
        from systems.sr.two_step_mc import run_two_step_mc_for_scenario
        from toolkit.io.dukou_data import load_dukou_data

        # 确定要运行的情景
        scenarios_to_run = []
        if args.scenario_a:
            scenarios_to_run.append('A')
        if args.scenario_b:
            scenarios_to_run.append('B')
        if args.scenario_c:
            scenarios_to_run.append('C')

        # 加载Dukou数据作为目标
        print("加载Dukou数据作为目标曲线...")
        dukou_data = load_dukou_data(
            mn_sr_threshold=args.mn_sr_threshold,
            lowess_frac=args.lowess_frac
        )

        manager = ScenarioManager()

        for scenario_name in scenarios_to_run:
            scenario = manager.get_scenario(scenario_name)

            print("\n" + "="*70)
            print(f"Scenario {scenario_name}: {scenario.name}")
            print("="*70)
            print(f"描述: {scenario.description}")
            print(f"年龄范围: {scenario.age_range[1]:.1f} - {scenario.age_range[0]:.1f} Ma")
            print(f"目标Sr: {scenario.target_sr_range[0]:.5f} - {scenario.target_sr_range[1]:.5f}")

            # 提取目标数据
            interval_data = dukou_data.get_interval_data(
                scenario.age_range[0], scenario.age_range[1]
            )

            if len(interval_data.sample_ids) == 0:
                print(f"警告: 该区间没有数据点，跳过")
                continue

            # 运行两步蒙特卡洛
            if args.two_step_mc:
                print(f"\n运行两步蒙特卡洛模拟...")
                print(f"  Step 1: {args.n_runs_step1} runs")
                print(f"  Step 2: {args.n_runs_step2} runs")

                result = run_two_step_mc_for_scenario(
                    scenario_name=scenario_name,
                    target_ages=interval_data.ages,
                    target_sr=interval_data.sr_ratios,
                    n_runs_step1=args.n_runs_step1,
                    n_runs_step2=args.n_runs_step2,
                    verbose=args.verbose
                )

                print("\n" + result.summary())

                # 保存结果
                if args.output and result.step2:
                    output_file = args.output.replace('.csv', f'_scenario_{scenario_name}.csv')

                    # 创建结果DataFrame
                    results_df = pd.DataFrame(result.step2.successful_params)
                    results_df.to_csv(output_file, index=False)
                    print(f"\n参数结果保存到: {output_file}")
            else:
                # 简单运行情景
                print(f"\n运行简单蒙特卡洛模拟 (100 runs)...")
                results = manager.run_all_scenarios(n_runs=100, ages_per_scenario=10)

                if scenario_name in results:
                    result = results[scenario_name]
                    if result.success:
                        print(f"成功: {len(result.successful_params)} 个参数组合")

        print()
        return

    # 创建Sr同位素体系
    sr = SrIsotopeSystem(scenario=args.scenario)

    # ===== 基本计算 =====
    if args.calculate or not (args.monte_carlo or args.scenario_analysis or args.sensitivity):
        print("[Seawater Sr Isotope Calculation]")
        print("-" * 50)

        # 创建通量配置
        config = SrFluxConfig(
            F_river=args.F_river,
            R_river=args.R_river,
            F_hydrothermal_highT=args.F_highT,
            R_hydrothermal_highT=args.R_highT,
            F_hydrothermal_lowT=args.F_lowT,
            R_hydrothermal_lowT=args.R_lowT,
            F_diagenetic=args.F_dia,
            R_diagenetic=args.R_dia,
        )

        # 计算海水Sr同位素
        ratio = sr.calculate_seawater_ratio(config)

        print("Input parameters:")
        print(f"  Riverine:     F = {args.F_river/1e9:.2f}×10⁹ mol/yr, R = {args.R_river:.5f}")
        print(f"  High-T hydro: F = {args.F_highT/1e9:.2f}×10⁹ mol/yr, R = {args.R_highT:.5f}")
        print(f"  Low-T hydro:  F = {args.F_lowT/1e9:.2f}×10⁹ mol/yr, R = {args.R_lowT:.5f}")
        print(f"  Diagenetic:   F = {args.F_dia/1e9:.2f}×10⁹ mol/yr, R = {args.R_dia:.5f}")

        print(f"\nResult:")
        print(f"  Seawater ⁸⁷Sr/⁸⁶Sr = {ratio:.5f}")

        # 与现代海水比较
        modern_observed = 0.70917
        diff = (ratio - modern_observed) * 10000  # ppm
        print(f"  Modern observed:    {modern_observed:.5f}")
        print(f"  Difference:         {diff:+.2f} ppm")
        print()

    # ===== 敏感性分析 =====
    if args.sensitivity:
        print(f"[Sensitivity Analysis: {args.sensitivity}]")
        print("-" * 50)

        # 确定参数范围
        if args.sens_range:
            param_range = tuple(args.sens_range)
        else:
            # 自动范围
            default_ranges = {
                'F_river': (10e9, 190e9),
                'R_river': (0.705, 0.725),
                'F_highT': (2e9, 35e9),
                'R_highT': (0.703, 0.704),
                'F_lowT': (2.5e9, 40e9),
                'R_lowT': (0.7025, 0.7084),
                'F_dia': (1e9, 10e9),
                'R_dia': (0.707, 0.710),
            }
            param_range = default_ranges.get(args.sensitivity, (0, 1))

        result = sr.sensitivity_analysis(
            param_name=args.sensitivity,
            param_range=param_range,
            n_steps=args.n_steps
        )

        if result.success:
            param_vals = result.data['param_values']
            sr_vals = result.data['sr_ratios']

            print(f"Parameter range: [{param_vals.min():.4e}, {param_vals.max():.4e}]")
            print(f"Sr ratio range:  [{sr_vals.min():.5f}, {sr_vals.max():.5f}]")

            # 显示关键点
            print(f"\nKey points:")
            for i in [0, len(param_vals)//2, -1]:
                print(f"  {args.sensitivity}={param_vals[i]:.4e}: R_sw={sr_vals[i]:.5f}")

            # 保存结果
            if args.output:
                df = pd.DataFrame({
                    args.sensitivity: param_vals,
                    'Sr_ratio': sr_vals,
                    'sensitivity': result.data['sensitivity']
                })
                df.to_csv(args.output, index=False)
                print(f"\nResults saved to: {args.output}")

        print()

    # ===== 情景分析 =====
    if args.scenario_analysis:
        print(f"[Scenario Analysis: {args.scenario_analysis}]")
        print("-" * 50)
        print(f"Running {args.n_runs} Monte Carlo simulations...")
        print("Note: No observation-based filtering in scenario analysis")

        # 创建自定义的随机模型，不过滤
        from systems.sr.model import StochasticSrModel
        import numpy as np

        time_points = np.linspace(args.time_span[0], args.time_span[1], args.n_time_points)
        stochastic = StochasticSrModel(
            time_points=time_points,
            observed_data=None  # 不过滤
        )

        # 获取情景参数范围
        from systems.sr.parameters import STOCHASTIC_RANGES
        param_ranges = STOCHASTIC_RANGES.copy()

        # 根据情景移除固定参数
        fixed_map = {
            'scenario1': ['F_diagenetic', 'R_diagenetic'],
            'scenario2': ['R_river', 'F_diagenetic', 'R_diagenetic'],
            'scenario3': ['F_river', 'F_diagenetic', 'R_diagenetic'],
            'scenario4': ['F_hydrothermal_highT', 'F_diagenetic', 'R_diagenetic'],
            'scenario5': ['F_hydrothermal_lowT', 'F_diagenetic', 'R_diagenetic'],
            'scenario6': ['F_hydrothermal_highT', 'F_hydrothermal_lowT', 'F_diagenetic', 'R_diagenetic'],
            'scenario7': ['F_river', 'R_river', 'F_diagenetic', 'R_diagenetic'],
        }

        for param in fixed_map.get(args.scenario_analysis, []):
            if param in param_ranges:
                del param_ranges[param]

        result = stochastic.run_monte_carlo(
            n_runs=args.n_runs,
            param_ranges=param_ranges,
            verbose=args.verbose
        )

        if result.success:
            print(f"\nResults:")
            print(f"  Successful runs: {result.data['n_successful']}")
            print(f"  Success rate: {result.data['success_rate']:.2%}")

            # 显示参数统计
            print(f"\nParameter statistics:")
            for param_name, stats in result.statistics.items():
                print(f"  {param_name}:")
                print(f"    Mean:   {stats['mean']:.3e}")
                print(f"    Median: {stats['median']:.3e}")
                print(f"    Range:  [{stats['min']:.3e}, {stats['max']:.3e}]")

            # 保存结果
            if args.output:
                output_data = {
                    'time_Ma': result.time,
                    'mean_Sr_ratio': result.data['mean_ratio'],
                    'std_Sr_ratio': result.data['std_ratio'],
                    'p2.5': result.data['percentile_2.5'],
                    'p97.5': result.data['percentile_97.5'],
                }

                # 添加参数统计
                for param_name in result.statistics.keys():
                    param_vals = result.get_parameter_array(param_name)
                    output_data[f'{param_name}_values'] = param_vals[:len(result.time)]

                df = pd.DataFrame(output_data)
                df.to_csv(args.output, index=False)
                print(f"\nResults saved to: {args.output}")
        else:
            print(f"  Error: {result.message}")

        print()

    # ===== 蒙特卡洛模拟 =====
    if args.monte_carlo and not args.scenario_analysis:
        print("[Monte Carlo Stochastic Simulation]")
        print("-" * 50)
        print(f"Scenario: {args.scenario}")
        print(f"Time span: {args.time_span[0]} - {args.time_span[1]} Ma")
        print(f"Number of runs: {args.n_runs}")
        print("Note: Running without observation-based filtering")
        print()

        # 使用不过滤的随机模型
        from systems.sr.model import StochasticSrModel
        import numpy as np
        from systems.sr.parameters import STOCHASTIC_RANGES

        time_points = np.linspace(args.time_span[0], args.time_span[1], args.n_time_points)
        stochastic = StochasticSrModel(
            time_points=time_points,
            observed_data=None  # 不过滤
        )

        result = stochastic.run_monte_carlo(
            n_runs=args.n_runs,
            param_ranges=STOCHASTIC_RANGES,
            verbose=args.verbose
        )

        if result.success:
            print(f"\nSimulation Results:")
            print(f"  Success rate: {result.data['success_rate']:.2%}")
            print(f"  Successful runs: {result.data['n_successful']}/{args.n_runs}")

            print(f"\nSeawater Sr isotope evolution:")
            ages = result.time
            mean_ratio = result.data['mean_ratio']
            std_ratio = result.data['std_ratio']

            # 显示关键时间点
            key_ages = [ages[0], ages[len(ages)//2], ages[-1]]
            print(f"  {'Age (Ma)':<12} {'Mean':<10} {'Std':<10}")
            print(f"  {'-'*32}")
            for age in key_ages:
                idx = np.argmin(np.abs(ages - age))
                print(f"  {age:<12.1f} {mean_ratio[idx]:<10.5f} {std_ratio[idx]:<10.5f}")

            print(f"\nParameter statistics:")
            for param_name, stats in list(result.statistics.items())[:3]:
                print(f"  {param_name}: {stats['mean']:.3e} ± {stats['std']:.3e}")

            # 保存结果
            if args.output:
                df = pd.DataFrame({
                    'time_Ma': ages,
                    'mean_Sr_ratio': mean_ratio,
                    'std_Sr_ratio': std_ratio,
                    'p2.5': result.data['percentile_2.5'],
                    'p97.5': result.data['percentile_97.5'],
                })
                df.to_csv(args.output, index=False)
                print(f"\nResults saved to: {args.output}")
        else:
            print(f"  Simulation failed: {result.message}")

        print()

    print("="*70)
    print("Analysis complete")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
