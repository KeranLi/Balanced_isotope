"""
Mg同位素体系模块 - 统一接口

提供两种分析体系：
1. 碳酸盐体系 (carbonate): 基于海相碳酸盐岩 δ²⁶Mg 反演风化
   - 模块: systems.mg.carbonate
   - 参考: Kasemann et al. (2014)
   - 适用: 海相碳酸盐岩、海水演化模拟

2. 碎屑岩体系 (siliciclastic): 基于陆源碎屑沉积物 δ²⁶Mg 反演风化
   - 模块: systems.mg.silicate
   - 参考: Hu et al. (2023)
   - 适用: 黏土矿物、河流沉积物、大陆边缘泥质岩

使用方式:
    # 方式1: 直接导入特定体系
    from systems.mg.carbonate import MgIsotopeSystem as CarbonateMgSystem
    from systems.mg.silicate import SilicateMgSystem
    
    # 方式2: 使用工厂函数
    from systems.mg import create_mg_system
    system = create_mg_system(component_type='carbonate')  # 或 'siliciclastic'

体系对比:
    | 特征           | 碳酸盐体系              | 碎屑岩体系              |
    |----------------|------------------------|------------------------|
    | 输入数据       | 碳酸盐岩 δ²⁶Mg         | 碎屑沉积物 δ²⁶Mg       |
    | 模型原理       | 海水沉淀分馏            | 风化残余分馏            |
    | 应用对象       | 海相碳酸盐岩            | 陆源碎屑沉积物          |
    | 主要端元       | 硅酸盐风化、碳酸盐风化   | UCC、硅酸盐风化         |
    | 输出结果       | 风化比例、海水演化      | 风化通量、SWI指数       |
"""

# 碳酸盐体系 (Kasemann et al., 2014)
from systems.mg.carbonate import (
    MgIsotopeSystem as CarbonateMgSystem,
    MgWeatheringModel,
    WeatheringFluxConfig,
)
from systems.mg.parameters import (
    get_mg_parameters,
    get_cryogenian_parameters,
    calculate_river_delta26,
    solve_f_silicate,
)

# 碎屑岩体系 (Hu et al., 2023)
from systems.mg.silicate import (
    SilicateMgSystem,
    SilicateWeatheringParams,
    SilicateModelResult,
)

# 为兼容性保留原名
MgIsotopeSystem = CarbonateMgSystem


def create_mg_system(component_type: str = 'carbonate', **kwargs):
    """
    工厂函数：创建 Mg 同位素体系实例
    
    Parameters:
        component_type: 体系类型
            - 'carbonate': 碳酸盐体系 (Kasemann et al., 2014)
            - 'siliciclastic': 碎屑岩体系 (Hu et al., 2023)
        **kwargs: 传递给体系构造函数的参数
        
    Returns:
        对应的体系实例
        
    Raises:
        ValueError: 如果 component_type 不支持
        
    Examples:
        >>> # 碳酸盐体系
        >>> system = create_mg_system('carbonate', scenario='modern')
        >>> 
        >>> # 碎屑岩体系
        >>> system = create_mg_system('siliciclastic', basin='changjiang')
    """
    component_type = component_type.lower().strip()
    
    if component_type in ('carbonate', 'carb'):
        return CarbonateMgSystem(**kwargs)
    
    elif component_type in ('siliciclastic', 'silicate', 'detrital', 'clastic'):
        return SilicateMgSystem(**kwargs)
    
    else:
        raise ValueError(
            f"不支持的体系类型: '{component_type}'\n"
            f"支持的类型: 'carbonate' (碳酸盐), 'siliciclastic' (碎屑岩)"
        )


def list_mg_systems():
    """列出可用的 Mg 同位素体系"""
    systems = [
        {
            'type': 'carbonate',
            'aliases': ['carb'],
            'name': '碳酸盐体系',
            'reference': 'Kasemann et al. (2014)',
            'description': '基于海相碳酸盐岩 δ²⁶Mg 反演风化',
            'applicability': ['海相碳酸盐岩', '海水演化模拟'],
        },
        {
            'type': 'siliciclastic',
            'aliases': ['silicate', 'detrital', 'clastic'],
            'name': '碎屑岩体系',
            'reference': 'Hu et al. (2023)',
            'description': '基于陆源碎屑沉积物 δ²⁶Mg 反演风化',
            'applicability': ['黏土矿物', '河流沉积物', '泥质岩'],
        },
    ]
    return systems


__all__ = [
    # 碳酸盐体系
    'CarbonateMgSystem',
    'MgWeatheringModel',
    'WeatheringFluxConfig',
    # 碎屑岩体系
    'SilicateMgSystem',
    'SilicateWeatheringParams',
    'SilicateModelResult',
    # 参数函数
    'get_mg_parameters',
    'get_cryogenian_parameters',
    'calculate_river_delta26',
    'solve_f_silicate',
    # 工厂函数
    'create_mg_system',
    'list_mg_systems',
    # 兼容性保留
    'MgIsotopeSystem',
]
