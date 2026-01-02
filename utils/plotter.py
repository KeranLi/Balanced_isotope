import matplotlib.pyplot as plt
import numpy as np  # 确保导入 numpy

class Plotter:
    def __init__(self):
        """
        初始化绘图类。
        """
        pass

    def plot_depth_vs_isotope(self, depth_data, rm_data):
        """
        绘制深度与同位素的关系图。
        :param depth_data: 深度数据
        :param rm_data: 同位素数据（如 delta_25_Mg_iso 或 delta_26_Mg_iso）
        """
        plt.figure(figsize=(8, 6))
        plt.scatter(depth_data, rm_data, color='green', alpha=0.5)
        plt.xlabel('Depth (meters)')
        plt.ylabel('Isotope Data')
        plt.title('Depth vs Isotope Data')
        plt.show()

    def plot_results(self, X, Y, Xedgemesh, Yedgemesh, N, Agesort2, mmRsort, eeRsort, fcarbin_all, fcarbout_all, binX, binY):
        """
        绘制结果图表，包括模拟数据。
        :param X: X轴数据
        :param Y: Y轴数据
        :param Xedgemesh: X轴边缘网格数据
        :param Yedgemesh: Y轴边缘网格数据
        :param N: 数据矩阵
        :param Agesort2: 排序后的年龄数据
        :param mmRsort: 平均值
        :param eeRsort: 标准差
        :param fcarbin_all: 模拟的 fcarbin 数据
        :param fcarbout_all: 模拟的 fcarbout 数据
        :param binX: X轴的 bin 边缘
        :param binY: Y轴的 bin 边缘
        """
        # 转换 mmRsort 和 eeRsort 为 numpy 数组，以便进行数学运算
        mmRsort = np.array(mmRsort)
        eeRsort = np.array(eeRsort)

        # 绘制热图
        plt.pcolor(Xedgemesh, Yedgemesh, N)
        plt.colorbar(label="Frequency")

        # 绘制模拟的 fcarbin 和 fcarbout 数据
        plt.scatter(fcarbin_all, fcarbout_all, color='red', alpha=0.3, label='Simulation Results')

        # 绘制平均值和标准差曲线
        plt.plot(Agesort2, mmRsort, label='Mean', color='black')
        plt.plot(Agesort2, mmRsort - eeRsort, label='Lower Bound', color='blue')
        plt.plot(Agesort2, mmRsort + eeRsort, label='Upper Bound', color='blue')

        # 使用 binX 和 binY 来调整热图的范围
        plt.xlim(binX[0], binX[-1])  # 设置 X 轴范围
        plt.ylim(binY[0], binY[-1])  # 设置 Y 轴范围

        # 设置标签和标题
        plt.xlabel('Age (Million years)')
        plt.ylabel('swpre')
        plt.title('Mass Balance Model with Monte Carlo Simulation')

        plt.legend()
        plt.show()
