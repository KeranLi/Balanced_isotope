from utils.data_processor import DataProcessor
from utils.plotter import Plotter
from calculation.mass_balance_model import MassBalanceModel
from calculation.monte_carlo import MonteCarloSimulation
from calculation.isotope import MgIsotopeModel

import numpy as np
import argparse
import matplotlib.pyplot as plt

def main():
    # 设置命令行参数解析器
    parser = argparse.ArgumentParser(description="Run the Mg isotope weathering model")
    
    # 添加命令行参数：文件路径、沉积速率、同位素列名和蒙特卡洛模拟次数
    parser.add_argument('file_path', type=str, help="Path to the data file (xlsx or csv)")
    parser.add_argument('rm_col', type=str, choices=['delta_25_Mg_iso', 'delta_26_Mg_iso'], help="Column name for the isotope data")
    parser.add_argument('--sedimentation_rate', type=float, default=3, help="Sedimentation rate in meters per million years (default: 3)")
    parser.add_argument('--n_simulations', type=int, default=100, help="Number of Monte Carlo simulations (default: 100)")

    # 解析命令行参数
    args = parser.parse_args()
    
    # 打印参数
    print(f"Using file: {args.file_path}")
    print(f"Using isotope column: {args.rm_col}")
    print(f"Using sedimentation rate: {args.sedimentation_rate}")
    print(f"Using {args.n_simulations} Monte Carlo simulations")

    # 1. 数据加载
    data_processor = DataProcessor(args.file_path, data_format='xlsx', rm_col=args.rm_col, sediment_rate=args.sedimentation_rate)
    
    # 2. 数据处理
    # 获取处理后的数据（如果没有深度数据，自动生成）
    age_data, depth_data, rm_data = data_processor.process_data()

    # 插值操作：设定插值年龄范围
    age_interpolation_range = np.arange(0, 2000, 0.01)  # 插值范围
    RMg_interpolated = data_processor.interpolate_data(age_interpolation_range)

    # 3. 质量平衡模型计算
    mass_balance_model = MassBalanceModel(RMg_interpolated, age_interpolation_range)
    swpre = mass_balance_model.calculate_swpre()  # 计算swpre

    # 4. 海水Mg同位素计算
    model = MgIsotopeModel(RMg_interpolated, age_interpolation_range, sediment_rate=args.sedimentation_rate, 
                           mass_balance_model=mass_balance_model, monte_carlo_simulation=None, plotter=None)

    # 假设碳酸盐和硅酸盐风化输入通量
    Fcarb_in = 0.5  # 示例碳酸盐风化输入通量
    Fsilicate_in = 0.3  # 示例硅酸盐风化输入通量
    carb_isotope = -0.2  # 示例碳酸盐Mg同位素值
    silicate_isotope = -0.5  # 示例硅酸盐Mg同位素值
    seawater_isotope = model.calculate_seawater_isotope(Fcarb_in, Fsilicate_in, carb_isotope, silicate_isotope)

    # 5. 蒙特卡洛模拟
    monte_carlo = MonteCarloSimulation(swpre, RMg_interpolated)
    simulation_results = monte_carlo.run_simulation(n_simulations=args.n_simulations)

    # 6. 生成绘图数据
    X = age_interpolation_range
    Y = swpre

    # 生成网格数据
    Xedge = np.linspace(min(X), max(X), 300)
    Yedge = np.linspace(min(Y), max(Y), 300)
    N, binX, binY = np.histogram2d(X, Y, bins=[Xedge, Yedge])

    Xedgemesh, Yedgemesh = np.meshgrid(Xedge[:-1], Yedge[:-1])

    # 计算排序后的统计数据
    Agesort2 = np.sort(X)
    mmRsort = [np.mean(Y[X == age]) for age in Agesort2]
    eeRsort = [np.std(Y[X == age]) for age in Agesort2]

    # 7. 使用simulation_results更新绘图数据
    # 假设simulation_results返回的是多个模拟结果列表，例如fcarbin, fcarbout等
    fcarbin_all = np.concatenate([result[0] for result in simulation_results])
    fcarbout_all = np.concatenate([result[1] for result in simulation_results])

    # 8. 绘图：展示深度和同位素的关系图
    plotter = Plotter()
    plotter.plot_depth_vs_isotope(depth_data, rm_data)  # 绘制深度和同位素的关系图

    # 绘制海水Mg同位素随时间变化的图
    plt.figure(figsize=(10, 6))
    plt.plot(age_interpolation_range, seawater_isotope, label='Seawater Mg Isotope', color='blue')
    plt.xlabel('Age (Million years)')
    plt.ylabel('Mg Isotope')
    plt.title('Seawater Mg Isotope Changes Over Time')
    plt.legend()
    plt.show()

    # 绘制结果
    plotter.plot_results(X, Y, Xedgemesh, Yedgemesh, N, Agesort2, mmRsort, eeRsort, fcarbin_all, fcarbout_all, binX, binY)

if __name__ == "__main__":
    main()
