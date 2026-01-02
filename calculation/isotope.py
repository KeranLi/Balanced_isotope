import numpy as np
import matplotlib.pyplot as plt

class MgIsotopeModel:
    def __init__(self, RMg_interpolated, age_interpolation_range, sediment_rate, mass_balance_model, monte_carlo_simulation, plotter):
        """
        初始化模型，设置初始值。
        :param RMg_interpolated: 插值后的Mg同位素数据
        :param age_interpolation_range: 年龄范围
        :param sediment_rate: 沉积速率
        :param mass_balance_model: 风化质量平衡模型实例
        :param monte_carlo_simulation: 传入蒙特卡洛模拟实例
        :param plotter: 传入绘图实例
        """
        self.RMg_interpolated = RMg_interpolated
        self.age_interpolation_range = age_interpolation_range
        self.sediment_rate = sediment_rate
        self.mass_balance_model = mass_balance_model  # 传递质量平衡模型实例
        self.moncarlo_simulation = monte_carlo_simulation
        self.plotter = plotter
        self.seawater_isotope = None  # 海水中的Mg同位素变化

    def calculate_seawater_isotope(self, Fcarb_in, Fsilicate_in, carb_isotope, silicate_isotope, dt=1e-10):
        """
        根据风化输入和质量平衡方程计算海水中的Mg同位素（考虑微分）。
        :param Fcarb_in: 碳酸盐风化输入通量
        :param Fsilicate_in: 硅酸盐风化输入通量
        :param carb_isotope: 碳酸盐的Mg同位素组成
        :param silicate_isotope: 硅酸盐的Mg同位素组成
        :param dt: 时间步长
        :return: 海水中的Mg同位素
        """
        # 获取风化过程的Mg同位素（swpre）
        swpre = self.mass_balance_model.calculate_swpre()

        # 初始化海水同位素数组
        seawater_isotope = np.zeros_like(swpre)
        seawater_isotope[0] = silicate_isotope  # 初始海水同位素（可以设置为初始值，例如硅酸盐的同位素值）

        # 使用微分方程求解海水Mg同位素随时间变化
        for i in range(1, len(swpre)):
            # 质量平衡公式：风化输入通量和风化同位素的组合
            delta_seawater_isotope = (Fsilicate_in * silicate_isotope + Fcarb_in * carb_isotope) / (Fsilicate_in + Fcarb_in)
            
            # 微分方程的数值解法（Euler方法）
            seawater_isotope[i] = seawater_isotope[i-1] + delta_seawater_isotope * dt

        return seawater_isotope



