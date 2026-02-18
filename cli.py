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
  # Mg同位素分析
  python cli.py mg --file data/Nie_Section_A.xlsx --column delta_26_Mg_iso
  
  # C同位素DOC模型
  python cli.py c --scenario dice --plot
  
  # 列出支持的体系
  python cli.py list
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # ===== Mg同位素命令 =====
    mg_parser = subparsers.add_parser('mg', help='Mg isotope system')
    mg_parser.add_argument('--file', type=str, help='Input data file (xlsx/csv)')
    mg_parser.add_argument('--column', type=str, default='delta_26_Mg_iso',
                          choices=['delta_25_Mg_iso', 'delta_26_Mg_iso'],
                          help='Isotope column name')
    mg_parser.add_argument('--sediment-rate', type=float, default=3,
                          help='Sedimentation rate (m/Ma)')
    mg_parser.add_argument('--weathering-ratio', action='store_true',
                          help='Calculate weathering end-member ratios')
    mg_parser.add_argument('--evolution', action='store_true',
                          help='Run seawater evolution model')
    mg_parser.add_argument('--output', type=str, help='Output file path')
    
    # ===== C同位素命令 =====
    c_parser = subparsers.add_parser('c', help='C isotope system')
    c_parser.add_argument('--scenario', type=str, default='dice',
                         choices=['dice', 'modern'],
                         help='Model scenario')
    c_parser.add_argument('--F-odoc', type=float,
                         help='DOC remineralization flux (mol/Ma)')
    c_parser.add_argument('--target-excursion', type=float,
                         help='Target carbon isotope excursion (negative, in permil)')
    c_parser.add_argument('--plot', action='store_true',
                         help='Generate plots')
    c_parser.add_argument('--output', type=str, help='Output file path')
    
    # ===== 体系列表命令 =====
    subparsers.add_parser('list', help='List available isotope systems')
    
    # ===== 信息命令 =====
    info_parser = subparsers.add_parser('info', help='Show system information')
    info_parser.add_argument('element', type=str,
                            choices=['mg', 'c', 's', 'sr', 'nd'],
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


def list_systems():
    """列出可用的同位素体系"""
    print("\n=== Available Isotope Systems ===\n")
    
    systems = [
        ('mg', 'Magnesium', '风化-沉积体系，估算碳酸盐/硅酸盐风化比例'),
        ('c', 'Carbon', '碳循环，DOC氧化与碳同位素负漂'),
        ('s', 'Sulfur', '硫循环，硫酸盐还原（计划中）'),
        ('sr', 'Strontium', 'Sr同位素，风化示踪（计划中）'),
        ('nd', 'Neodymium', 'Nd同位素，洋流循环（计划中）'),
    ]
    
    for element, name, description in systems:
        status = "✓" if element in ['mg', 'c'] else "○"
        print(f"  {status} {element.upper():2} - {name:12} : {description}")
    
    print("\n✓ = Implemented, ○ = Planned")
    print()


def show_info(element: str):
    """显示体系详细信息"""
    print(f"\n=== {element.upper()} Isotope System ===\n")
    
    if element == 'mg':
        from systems.mg import MgIsotopeSystem, get_mg_parameters
        params = get_mg_parameters()
        
        print(f"Element: {params.name}")
        print(f"Reference: {params.reference_standard}")
        print(f"\nEnd-members (δ²⁶Mg, ‰):")
        for name, data in params.end_members.items():
            print(f"  {name:12}: {data['delta26']:+.2f} ± {data.get('uncertainty', 0):.2f}")
        
        print(f"\nReservoir mass: {params.reservoir_mass:.2e} mol")
        print(f"Input fluxes (mol/Ma):")
        for name, value in params.input_fluxes.items():
            print(f"  {name}: {value:.2e}")
    
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
    
    print()


def run_mg_analysis(args):
    """运行Mg同位素分析"""
    print("\n=== Mg Isotope Analysis ===\n")
    
    from systems.mg import MgIsotopeSystem
    
    system = MgIsotopeSystem()
    
    if args.file:
        # 数据加载和分析
        print(f"Loading data from: {args.file}")
        print(f"Using column: {args.column}")
        print(f"Sedimentation rate: {args.sediment_rate} m/Ma")
        
        # TODO: 实现数据加载和完整分析流程
        print("\n(Data loading not yet fully implemented in new architecture)")
    
    if args.weathering_ratio:
        print("\nCalculating weathering end-member ratios...")
        # 示例计算
        delta_sample = -2.0  # 示例值
        ratios = system.calculate_weathering_ratio(delta_sample, delta_seawater=-0.83)
        print(f"  Sample δ²⁶Mg: {delta_sample}‰")
        print(f"  Carbonate fraction: {ratios['f_carbonate']:.2%}")
        print(f"  Silicate fraction: {ratios['f_silicate']:.2%}")
    
    if args.evolution:
        print("\nRunning seawater evolution model...")
        result = system.seawater_evolution(
            time_span=(0, 100),  # 100 Ma
            initial_delta=-0.5,
            flux_scenario='modern'
        )
        
        if result.success:
            print(f"  Time span: 0-100 Ma")
            print(f"  Final δ²⁶Mg_seawater: {result.values[-1, 1]:.3f}‰")
        else:
            print(f"  Error: {result.message}")
    
    print()


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
            
            if args.plot:
                print("\nGenerating plots...")
                # TODO: 实现绘图
                print("  (Plotting not yet fully implemented)")
    
    print()


if __name__ == '__main__':
    main()
