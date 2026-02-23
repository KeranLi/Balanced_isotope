# Sr同位素海洋箱模型使用文档

基于 Wang et al. (2021) Earth-Science Reviews 218: 103679

## 1. 模型概述

本模型实现了Wang等人(2021)提出的**随机海洋箱式模型**（Stochastic Oceanic Box Model），用于模拟海水Sr同位素（⁸⁷Sr/⁸⁶Sr）的演化。

### 1.1 核心特点

- **四端元混合模型**：河流、高温热液、低温热液、成岩作用
- **蒙特卡洛随机模拟**：探索参数空间，减少模型偏差
- **情景分析**：测试不同约束条件下的解空间
- **与观测数据对比**：筛选符合经验记录的参数组合

### 1.2 质量平衡方程

```
F_in = F_riv + F_highT + F_lowT + F_dia = F_out

R_ocean = (F_riv×R_riv + F_highT×R_highT + F_lowT×R_lowT + F_dia×R_dia) / 
          (F_riv + F_highT + F_lowT + F_dia)
```

## 2. 快速开始

### 2.1 基本使用

```python
from systems.sr import SrIsotopeSystem

# 创建Sr同位素体系
sr = SrIsotopeSystem(scenario='modern')

# 计算海水Sr同位素比值
ratio = sr.calculate_seawater_ratio()
print(f"现代海水⁸⁷Sr/⁸⁶Sr: {ratio:.5f}")
```

### 2.2 自定义参数

```python
from systems.sr import SrFluxConfig

# 创建自定义通量配置
config = SrFluxConfig(
    F_river=100e9,          # 河流通量 (mol/yr)
    R_river=0.712,          # 河流⁸⁷Sr/⁸⁶Sr
    F_hydrothermal_highT=20e9,
    R_hydrothermal_highT=0.7037,
    F_hydrothermal_lowT=15e9,
    R_hydrothermal_lowT=0.7084,
)

ratio = sr.calculate_seawater_ratio(config)
```

## 3. 蒙特卡洛随机模拟

### 3.1 基本模拟

```python
sr = SrIsotopeSystem(scenario='permian')

result = sr.monte_carlo_simulation(
    time_span=(299, 252),    # 二叠纪 (Ma)
    n_runs=5000,             # 蒙特卡洛运行次数
    n_time_points=50,        # 时间点数量
    filter_by_data=True,     # 根据观测数据过滤
)

if result.success:
    print(f"成功率: {result.data['success_rate']:.2%}")
    print(f"成功运行次数: {result.data['n_successful']}")
```

### 3.2 结果分析

```python
# 获取参数统计信息
for param_name, stats in result.statistics.items():
    print(f"{param_name}:")
    print(f"  平均值: {stats['mean']:.3e}")
    print(f"  中位数: {stats['median']:.3e}")
    print(f"  范围: [{stats['min']:.3e}, {stats['max']:.3e}]")

# 获取Sr同位素演化
mean_ratio = result.data['mean_ratio']      # 平均值
std_ratio = result.data['std_ratio']        # 标准差
p25 = result.data['percentile_2.5']         # 2.5%分位数
p975 = result.data['percentile_97.5']       # 97.5%分位数
```

## 4. 情景分析

Wang et al. (2021) 论文中测试了7种情景：

```python
sr = SrIsotopeSystem(scenario='permian')

# 运行特定情景
result = sr.scenario_analysis(
    scenario_name='scenario1',   # 选择情景
    n_runs=5000,
    time_span=(299, 252),
)

# 支持的情景：
# - 'scenario1': 所有参数可变（除成岩）
# - 'scenario2': 固定河流同位素
# - 'scenario3': 固定河流通量
# - 'scenario4': 固定高温热液通量
# - 'scenario5': 固定低温热液通量
# - 'scenario6': 固定所有热液通量
# - 'scenario7': 固定河流参数
```

## 5. 参数敏感性分析

```python
sr = SrIsotopeSystem(scenario='modern')

# 分析河流同位素的敏感性
result = sr.sensitivity_analysis(
    param_name='R_river',
    param_range=(0.705, 0.725),
    n_steps=50,
)

# 结果
param_values = result.data['param_values']
sr_ratios = result.data['sr_ratios']
sensitivity = result.data['sensitivity']  # d(R_sw)/d(param)
```

## 6. 端元参数

### 6.1 现代参数（默认值）

| 端元 | 通量 (10⁹ mol/yr) | ⁸⁷Sr/⁸⁶Sr | 来源 |
|------|------------------|-----------|------|
| 河流 | 47.6 | 0.71107 | Peucker-Ehrenbrink & Fiske, 2019 |
| 高温热液 | 8.04 | 0.7037 | Li & Elderfield, 2013 |
| 低温热液 | 10 | 0.7084 | Coogan & Dosso, 2015 |
| 成岩作用 | 3.4 | 0.7084 | 估算 |

### 6.2 随机模拟参数范围

基于Wang et al. (2021) Table 1：

| 参数 | 范围 | 不确定性 |
|------|------|----------|
| F_river | 10-190 × 10⁹ mol/yr | ±10% |
| R_river | 0.7025-0.75 | ±0.0005 |
| F_hydrothermal_lowT | 2.5-40 × 10⁹ mol/yr | ±10% |
| R_hydrothermal_lowT | 0.7025-0.7084 | ±0.0005 |
| F_hydrothermal_highT | 2-35 × 10⁹ mol/yr | ±10% |
| R_hydrothermal_highT | 0.703-0.704 | ±0.0001 |

## 7. 二叠纪应用示例

```python
from systems.sr import SrIsotopeSystem
import numpy as np

sr = SrIsotopeSystem(scenario='permian')

# 二叠纪关键时间点
ages_ma = [299, 273, 265, 259, 252]
stages = ["Asselian", "Kungurian", "Wordian", "Capitanian", "Changhsingian"]

# 典型观测值（基于Wang et al., 2021）
observed_ratios = [0.70827, 0.7070, 0.70683, 0.7069, 0.70717]

print("二叠纪Sr同位素演化:")
for age, stage, ratio in zip(ages_ma, stages, observed_ratios):
    print(f"{stage} ({age} Ma): {ratio:.5f}")

# 运行蒙特卡洛模拟探索驱动机制
result = sr.monte_carlo_simulation(
    time_span=(299, 252),
    n_runs=5000,
)
```

## 8. 主要结论（基于Wang et al., 2021）

1. **热液活动是主要驱动因素**：二叠纪Sr同位素的持续下降最可能由增强的热液活动引起

2. **低温热液的贡献**：全球变暖可能导致低温热液系统的Sr同位素组成向更非放射成因方向偏移

3. **大陆风化的作用**：传统解释强调大陆风化，但模型结果表明热液输入的变化更能解释观测记录

4. **Neotethys洋打开的影响**：与Neotethys洋打开相关的构造活动可能增强了热液活动

## 9. 运行示例

```bash
# 运行完整演示
python examples/sr_isotope_demo.py

# 运行主程序（包含所有同位素系统）
python run.py
```

## 10. 参考文献

- Wang et al. (2021). Revisiting the Permian seawater Sr/86Sr record: New perspectives from brachiopod proxy data and stochastic oceanic box models. *Earth-Science Reviews*, 218, 103679.

- Peucker-Ehrenbrink, B., & Fiske, G. J. (2019). A continental perspective of the seawater 87Sr/86Sr record: A review. *Chemical Geology*, 510, 140-165.

- Coogan, L. A., & Dosso, S. E. (2015). Alteration of ocean crust provides a strong temperature dependent feedback on the geological carbon cycle and is a primary driver of the Sr-isotopic composition of seawater. *Earth and Planetary Science Letters*, 415, 38-46.

- Li, G. J., & Elderfield, H. (2013). Evolution of carbon cycle over the past 100 million years. *Geochimica et Cosmochimica Acta*, 103, 11-25.
