import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

class DataProcessor:
    def __init__(self, file_path, data_format='xlsx', depth_col=None, age_col=None, rm_col='RMgPh', sediment_rate=3):
        """
        初始化数据处理类，处理没有深度数据的情况。
        :param file_path: 数据文件路径
        :param data_format: 数据格式，默认xlsx，支持xlsx, csv等
        :param depth_col: 深度列名（如果有的话）
        :param age_col: 年龄列名（如果有的话）
        :param rm_col: RMgPh列名
        :param sediment_rate: 沉积速率（米/百万年）
        """
        self.file_path = file_path
        self.data_format = data_format
        self.depth_col = depth_col
        self.age_col = age_col
        self.rm_col = rm_col
        self.sediment_rate = sediment_rate
        
        self.data = self.load_data()
    
    def load_data(self):
        """
        根据文件格式加载数据。
        :return: 加载的数据
        """
        if self.data_format == 'xlsx':
            data = pd.read_excel(self.file_path)
        elif self.data_format == 'csv':
            data = pd.read_csv(self.file_path)
        else:
            raise ValueError("Unsupported file format. Use 'xlsx' or 'csv'.")
        return data
    
    def process_data(self):
        """
        处理数据，生成虚拟的深度和年龄数据（如果没有的话）。
        :return: 处理后的深度、年龄数据，以及RMgPh数据
        """
        # 如果没有深度数据，我们根据行号生成深度数据
        if self.depth_col is None:
            self.data['Depth'] = np.arange(0, len(self.data) * self.sediment_rate, self.sediment_rate)
        
        # 如果没有年龄数据，我们根据深度推算年龄
        if self.age_col is None:
            self.data['Age'] = self.data['Depth'] / self.sediment_rate
        
        age_data = self.data['Age'].values
        depth_data = self.data['Depth'].values
        rm_data = self.data[self.rm_col].values
        
        return age_data, depth_data, rm_data
    
    def interpolate_data(self, age_interpolation_range, kind='linear'):
        """
        对原始数据进行插值处理。
        :param age_interpolation_range: 插值的年龄范围
        :param kind: 插值方式（默认线性插值）
        :return: 插值后的数据
        """
        age_data, _, rm_data = self.process_data()
        interp_func = interp1d(age_data, rm_data, kind=kind, fill_value="extrapolate")
        return interp_func(age_interpolation_range)
