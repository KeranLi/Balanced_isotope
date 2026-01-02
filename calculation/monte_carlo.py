import numpy as np
from scipy.stats import norm

class MonteCarloSimulation:
    def __init__(self, swpreinter, RMg_interpolated):
        """
        初始化蒙特卡洛模拟类。
        :param swpreinter: 插值后的swpre值
        :param RMg_interpolated: 插值后的RMg数据
        """
        self.swpreinter = swpreinter
        self.RMg_interpolated = RMg_interpolated

    def run_simulation(self, n_simulations=100, size=500000000):
        """
        执行蒙特卡洛模拟来估算风化平衡。
        :param n_simulations: 模拟次数
        :param size: 每次模拟的大小
        :return: 模拟结果
        """
        
        results = []
        
        # 将size限制为swpreinter的大小，避免超出索引范围
        size = min(size, len(self.RMg_interpolated))

        for _ in range(n_simulations):
            # 生成随机索引，确保所有索引都在有效范围内
            tt_index = np.random.randint(0, len(self.RMg_interpolated), size=size)  # 索引范围限定为 [0, len(self.RMg_interpolated) - 1]

            # 生成随机数据，fcarbin 和 fcarbout 应为一维数组
            fcarbin = np.random.rand(len(tt_index))  # 生成长度为 size 的一维数组
            fcarbout = np.random.rand(len(tt_index))  # 生成长度为 size 的一维数组
            carbin = self.swpreinter[tt_index] + norm.rvs(loc=-1.5, scale=0.2, size=len(fcarbin))
            silin = norm.rvs(loc=-0.4, scale=0.1, size=len(fcarbin))
            frsil = norm.rvs(loc=0.8, scale=0.4, size=len(fcarbin))
            frcarb = norm.rvs(loc=-1.5, scale=0.2, size=len(fcarbin))
            RMgsolv = self.RMg_interpolated[tt_index]
            sw = RMgsolv + norm.rvs(loc=0, scale=0.1, size=len(fcarbin))

            misfit = fcarbin * (carbin - silin) - fcarbout * (frcarb - frsil) - frsil - sw + silin
            Df = np.abs(misfit)
            sindex = np.where(Df < 0.001)  # 获取满足条件的索引

            # 使用sindex来索引fcarbin和fcarbout
            fcarbinsolved = fcarbin[sindex]
            fcarboutsolved = fcarbout[sindex]

            results.append((fcarbinsolved, fcarboutsolved))
        
        return results