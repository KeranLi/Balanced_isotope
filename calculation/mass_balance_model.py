import numpy as np
from scipy.interpolate import interp1d

class MassBalanceModel:
    def __init__(self, RMg_interpolated, age_interpolation_range):
        """
        初始化质量平衡模型类。
        :param RMg_interpolated: 插值后的Mg同位素数据
        :param age_interpolation_range: 年龄范围
        """
        self.RMg_interpolated = RMg_interpolated
        self.age_interpolation_range = age_interpolation_range
        self.swpre = None  # 确保swpre初始化为None
    
    def calculate_swpre(self):
        """
        计算swpre值，代表风化质量平衡模型的关键计算步骤。
        :return: swpre计算结果
        """
        # 如果self.swpre尚未计算，则进行计算
        if self.swpre is None:
            self.swpre = np.zeros(len(self.RMg_interpolated))
            for i in range(len(self.RMg_interpolated)):
                decay = np.arange(i, len(self.RMg_interpolated))
                Rdecay = np.exp(-0.001 * (decay - i) / 100)
                self.swpre[i] = np.sum(self.RMg_interpolated[i:] * Rdecay) / np.sum(Rdecay)
        return self.swpre
    
    def update_swpre_interpolation(self, new_age_range):
        """
        更新swpre的插值。
        :param new_age_range: 新的年龄范围
        :return: 插值后的swpre
        """
        # 如果swpre未计算，先计算它
        if self.swpre is None:
            self.calculate_swpre()

        # 使用 fill_value="extrapolate" 进行外推处理
        return interp1d(self.age_interpolation_range, self.swpre, kind='linear', fill_value="extrapolate")(new_age_range)
