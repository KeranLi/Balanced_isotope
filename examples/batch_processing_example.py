#!/usr/bin/env python3
"""
批量处理功能示例

演示如何使用 toolkit.io 模块进行批量数据处理
"""

import numpy as np
import pandas as pd
import sys
from pathlib import Path

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from toolkit.io import DataHandler, BatchProcessor


def example_1_basic_batch_processing():
    """
    示例1: 基本批量处理
    
    使用BatchProcessor处理整个Excel文件
    """
    print("=" * 70)
    print("示例1: 基本批量处理")
    print("=" * 70)
    
    # 创建处理器
    processor = BatchProcessor(
        element='u',
        scenario='modern',
        apply_diagenetic_correction=True,
        delta_diag=0.4,
        include_uncertainty=True,
        n_monte_carlo=1000
    )
    
    # 处理文件
    input_file = 'data/uranium_example.xlsx'
    output_file = 'results/batch_example_results.xlsx'
    
    print(f"\n输入文件: {input_file}")
    print(f"输出文件: {output_file}")
    print("\n处理中...")
    
    results_df = processor.process_file(
        input_file,
        output_path=output_file,
        show_progress=True
    )
    
    # 显示摘要
    print("\n" + "-" * 50)
    print("处理结果摘要:")
    print("-" * 50)
    print(f"总样品数: {len(results_df)}")
    print(f"成功: {results_df['processing_success'].sum()}")
    print(f"失败: {(~results_df['processing_success']).sum()}")
    
    print(f"\nf_anox 统计:")
    print(f"  均值: {results_df['f_anox'].mean():.1%}")
    print(f"  标准差: {results_df['f_anox'].std():.1%}")
    print(f"  最小值: {results_df['f_anox'].min():.1%}")
    print(f"  最大值: {results_df['f_anox'].max():.1%}")
    
    print(f"\n缺氧面积统计:")
    print(f"  均值: {results_df['anoxic_area_percent'].mean():.1f}%")
    print(f"  范围: [{results_df['anoxic_area_percent'].min():.1f}%, "
          f"{results_df['anoxic_area_percent'].max():.1f}%]")


def example_2_data_reading():
    """
    示例2: 数据读取和验证
    
    使用DataHandler读取和验证数据
    """
    print("\n" + "=" * 70)
    print("示例2: 数据读取和验证")
    print("=" * 70)
    
    # 创建数据处理器
    handler = DataHandler('data/uranium_example.xlsx')
    
    # 读取数据
    print("\n读取数据...")
    data = handler.read()
    
    # 获取数据信息
    info = handler.get_info()
    print("\n数据信息:")
    print(f"  行数: {info['n_rows']}")
    print(f"  列数: {info['n_columns']}")
    print(f"  检测到的同位素: {info['detected_isotopes']}")
    print(f"  检测到的列: {list(info['detected_columns'].keys())}")
    
    # 验证数据
    print("\n验证数据...")
    is_valid, errors = handler.validate_data('u')
    if is_valid:
        print("  ✓ 数据验证通过")
    else:
        print(f"  ✗ 验证错误: {errors}")
    
    # 获取U同位素数据
    u_data = handler.get_isotope_data('u')
    print(f"\nU同位素数据:")
    print(f"  δ²³⁸U值数量: {len(u_data['delta'])}")
    if 'delta_std' in u_data:
        print(f"  不确定度数量: {len(u_data['delta_std'])}")
    if 'sample_id' in u_data:
        print(f"  样品ID: {u_data['sample_id'][:3]}...")


def example_3_custom_processing():
    """
    示例3: 自定义处理逻辑
    
    遍历数据并自定义处理每个样品
    """
    print("\n" + "=" * 70)
    print("示例3: 自定义处理逻辑")
    print("=" * 70)
    
    from systems.u import UIsotopeSystem
    
    # 读取数据
    handler = DataHandler('data/uranium_example.xlsx')
    data = handler.read()
    u_data = handler.get_isotope_data('u')
    
    # 创建U同位素体系
    u_system = UIsotopeSystem(scenario='modern')
    
    # 自定义处理每个样品
    results = []
    print("\n处理样品...")
    
    for i in range(len(u_data['delta'])):
        delta238 = u_data['delta'][i]
        sample_id = u_data.get('sample_id', [f'Sample_{i}'])[i] if 'sample_id' in u_data else f'Sample_{i}'
        elevation = u_data.get('elevation_m', [np.nan])[i] if 'elevation_m' in u_data else np.nan
        
        # 稳态计算
        result = u_system.calculate_f_anox_steady_state(delta238_carb=delta238)
        
        # 保存结果
        results.append({
            'sample_id': sample_id,
            'elevation': elevation,
            'delta238_u': delta238,
            'f_anox': result['f_anox'],
            'anoxic_area': u_system.estimate_anoxic_area(result['f_anox'])
        })
    
    # 转换为DataFrame
    results_df = pd.DataFrame(results)
    
    print(f"\n处理完成，共 {len(results_df)} 个样品")
    print("\n结果预览:")
    print(results_df.head())
    
    # 保存结果
    output_file = 'results/custom_processing_results.xlsx'
    handler.write(results_df, output_file)
    print(f"\n结果已保存到: {output_file}")


def example_4_multiple_files():
    """
    示例4: 批量处理多个文件
    
    使用shell命令或Python批量处理多个文件
    """
    print("\n" + "=" * 70)
    print("示例4: 批量处理多个文件")
    print("=" * 70)
    
    import glob
    
    # 查找所有数据文件
    data_files = glob.glob('data/*.xlsx')
    print(f"\n找到 {len(data_files)} 个数据文件:")
    for f in data_files:
        print(f"  - {Path(f).name}")
    
    # 创建处理器
    processor = BatchProcessor(element='mg')
    
    print("\n处理Mg同位素数据文件...")
    for file_path in data_files:
        if 'Nie_Section' in file_path:
            print(f"\n处理: {file_path}")
            try:
                results_df = processor.process_file(file_path, show_progress=False)
                
                # 保存结果
                output_name = f"results/{Path(file_path).stem}_results.xlsx"
                processor.save_results(results_df, output_name)
                print(f"  ✓ 结果保存到: {output_name}")
                print(f"  样品数: {len(results_df)}")
                if 'f_carbonate' in results_df.columns:
                    print(f"  平均碳酸盐比例: {results_df['f_carbonate'].mean():.1%}")
            except Exception as e:
                print(f"  ✗ 错误: {e}")


def example_5_summary_statistics():
    """
    示例5: 生成处理摘要统计
    
    批量处理后生成摘要统计信息
    """
    print("\n" + "=" * 70)
    print("示例5: 生成处理摘要统计")
    print("=" * 70)
    
    # 先处理数据
    processor = BatchProcessor(element='u')
    results = processor.process_file('data/uranium_example.xlsx', show_progress=False)
    
    # 转换为ProcessingResult列表格式
    processing_results = []
    for _, row in results.iterrows():
        from toolkit.io.batch_processor import ProcessingResult
        
        input_data = {col: row[col] for col in ['sample_id', 'elevation_m', 'delta_238_u'] 
                     if col in row}
        output_data = {col: row[col] for col in ['f_anox', 'anoxic_area_percent'] 
                      if col in row}
        
        processing_results.append(ProcessingResult(
            success=row['processing_success'],
            index=len(processing_results),
            input_data=input_data,
            output_data=output_data,
            error_message=row.get('processing_error', '')
        ))
    
    # 获取摘要
    summary = processor.get_summary(processing_results)
    
    print("\n处理摘要:")
    print("-" * 50)
    for key, value in summary.items():
        if isinstance(value, float):
            if 'rate' in key or 'f_' in key:
                print(f"  {key}: {value:.1%}")
            else:
                print(f"  {key}: {value:.3f}")
        else:
            print(f"  {key}: {value}")


def create_sample_data():
    """
    创建示例数据文件
    """
    print("\n" + "=" * 70)
    print("创建示例数据文件")
    print("=" * 70)
    
    # 创建模拟数据
    np.random.seed(42)
    n = 20
    
    data = {
        'sample_id': [f'Sample_{i+1:02d}' for i in range(n)],
        'depth_m': np.linspace(0, 100, n),
        'delta_238_u': np.random.normal(-0.45, 0.1, n),
        'delta_238_u_std': np.random.uniform(0.03, 0.08, n)
    }
    
    df = pd.DataFrame(data)
    
    # 保存
    output_file = 'data/my_uranium_data.xlsx'
    df.to_excel(output_file, index=False)
    print(f"\n示例数据已创建: {output_file}")
    print(f"样品数: {len(df)}")
    print("\n数据预览:")
    print(df.head())


def main():
    """运行所有示例"""
    print("\n" + "=" * 70)
    print("批量处理功能完整示例")
    print("=" * 70)
    
    # 创建必要的目录
    Path('results').mkdir(exist_ok=True)
    
    # 运行示例
    example_1_basic_batch_processing()
    example_2_data_reading()
    example_3_custom_processing()
    example_4_multiple_files()
    example_5_summary_statistics()
    create_sample_data()
    
    print("\n" + "=" * 70)
    print("所有示例运行完成!")
    print("=" * 70)
    print("\n更多用法请参考:")
    print("  - BATCH_PROCESSING.md")
    print("  - python cli.py u --help")
    print()


if __name__ == '__main__':
    main()
