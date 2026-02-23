"""
Sr同位素体系模块 - 海洋Sr循环模型

基于Wang等(2021)论文实现：
- 随机海洋箱式模型（Stochastic Oceanic Box Model）
- 四端元混合模型（河流、高温热液、低温热液、成岩作用）
- 用于模拟二叠纪海水⁸⁷Sr/⁸⁶Sr演化

参考：
Wang et al. (2021) Earth-Science Reviews 218: 103679
"Revisiting the Permian seawater Sr/86Sr record: New perspectives from 
brachiopod proxy data and stochastic oceanic box models"

使用方式:
    from systems.sr import SrIsotopeSystem, SrFluxConfig
    
    # 创建Sr同位素体系
    sr = SrIsotopeSystem(scenario='permian')
    
    # 运行随机模拟
    result = sr.monte_carlo_simulation(
        time_span=(299, 252),  # 二叠纪 (Ma)
        n_runs=5000,
        filter_by_data=True
    )

主要功能：
1. 质量平衡计算 - 四端元混合模型
2. 随机蒙特卡洛模拟 - 探索参数空间
3. 情景分析 - 测试不同约束条件
4. 与观测数据对比 - 筛选符合经验记录的参数组合
"""

from systems.sr.model import (
    SrIsotopeSystem,
    SrFluxConfig,
    SrModelResult,
    StochasticSrModel,
)
from systems.sr.parameters import (
    get_sr_parameters,
    get_permian_parameters,
    get_modern_parameters,
)

__all__ = [
    # 主类
    'SrIsotopeSystem',
    'SrFluxConfig',
    'SrModelResult',
    'StochasticSrModel',
    # 参数函数
    'get_sr_parameters',
    'get_permian_parameters',
    'get_modern_parameters',
]
