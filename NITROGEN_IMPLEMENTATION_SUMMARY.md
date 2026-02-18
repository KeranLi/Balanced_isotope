# 氮同位素体系实现总结

## 已完成的工作

### 1. 文献学习与模拟思路总结

深入学习了 `reference/nitrogen/` 目录中的两篇核心文献：

#### Kang et al. (2023) - 新元古代氮循环
- **核心发现**: 800 Ma 前后沉积物 δ¹⁵N 出现阶跃上升（~1.5‰ → ~4.5‰）
- **模型解释**: 硝酸盐占比从 ~0.11 增加到 ~0.16，增加约50%
- **生态意义**: 硝酸盐可利用性限制了真核生物的生态崛起

#### Ma et al. (2025) - 早三叠世氮循环
- **核心发现**: 早三叠世海洋长期处于硝酸盐匮乏状态
- **阶段演化**: 
  - 第I阶段 (PTB-SSB): f_assimilator ≈ 0.06 (极端缺氧高温)
  - 第II阶段 (Spathian早期): f_assimilator ≈ 0.20 (氧化降温)
  - 第III阶段 (Spathian晚期): f_assimilator ≈ 0.09 (再缺氧)
- **生态影响**: 硝酸盐匮乏 + 铵盐毒性共同延迟生态系统复苏

### 2. 双箱稳态模型实现

#### 核心方程
```
储库1: N_fixer/ammonium (固氮生物-铵)
  输入: F_fix (固氮作用)
  输出: F_remin (再矿化) + F_fixer_burial (埋藏)
  同位素: δ¹⁵N_ammonium = δ¹⁵N_atmosphere + ε_fix

储库2: N_assimilator/nitrate (硝酸盐同化生物-硝酸盐)  
  输入: F_remin (再矿化)
  输出: F_wcd (水柱反硝化) + F_sd (沉积反硝化) + F_assimilator_burial (埋藏)
  同位素: δ¹⁵N_nitrate = δ¹⁵N_ammonium - (F_wcd×ε_wcd + F_sd×ε_sd)/F_remin

沉积物: 
  δ¹⁵N_sed = (1-f_assimilator)×δ¹⁵N_ammonium + f_assimilator×δ¹⁵N_nitrate
  f_assimilator = F_assimilator_burial / F_total_burial
```

#### 关键非线性关系
- f = 0 时: δ¹⁵N_sed ≈ -0.5‰ (纯固氮)
- f ≈ 0.47 时: δ¹⁵N_sed 达到最大值 (~6.2‰)
- f = 1 时: δ¹⁵N_sed ≈ -0.5‰ (纯硝酸盐同化)

### 3. 新增文件

```
systems/
├── n/
│   ├── __init__.py           # 模块初始化
│   ├── parameters.py         # 氮同位素参数定义
│   ├── model.py              # 核心模型实现
│   ├── DESIGN.md             # 设计文档
│   └── README.md             # 使用说明
toolkit/
└── math/
    └── statistics.py         # 统计工具 (Bootstrap, LOWESS, Changepoint)
examples/
└── nitrogen_example.py       # 使用示例
```

### 4. 新增功能

#### 氮同位素体系 (`systems.n`)
- `NIsotopeSystem` 类: 核心模型
  - `forward_model()`: 正向计算 δ¹⁵N_sed
  - `inverse_model()`: 反向反演 f_assimilator
  - `monte_carlo_simulation()`: 不确定性分析
  - `calculate_f_assimilator_curve()`: 关系曲线

#### 统计工具 (`toolkit.math.statistics`)
- `Bootstrap`: Bootstrap 重采样和置信区间
- `LOWESS`: 局部加权回归平滑
- `ChangepointDetector`: 变点检测 (AMOC, PELT)
- `TimeSeriesAnalysis`: 时间序列分析工具

### 5. 使用示例

```python
from systems.n import NIsotopeSystem
from toolkit.math.statistics import Bootstrap, LOWESS

# 创建模型
n_system = NIsotopeSystem(scenario='early_triassic')

# 正向计算
delta15N = n_system.forward_model(f_assimilator=0.2)

# 反向反演
result = n_system.inverse_model(delta15N_sed=3.0)

# 蒙特卡洛不确定性
mc_result = n_system.monte_carlo_simulation(f_assimilator=0.2)

# Bootstrap 分析
bootstrap_result = Bootstrap.confidence_interval(data)

# LOWESS 平滑
lowess_result = LOWESS.fit(ages, delta15N_values)
```

### 6. 验证结果

示例运行输出显示模型正确实现了：
1. 正向模型的非线性关系（f ≈ 0.47 时 δ¹⁵N 达到峰值）
2. 反向模型的合理反演范围
3. 蒙特卡洛不确定性的正确传播
4. 统计工具的正确功能

## 后续可扩展方向

1. **硫同位素体系**: 完善 `systems/s` 的多同位素质量无关分馏
2. **锶同位素体系**: 实现放射成因同位素体系
3. **耦合模型**: 碳-氮-硫耦合循环模型
4. **可视化工具**: 添加专门的绘图模块
5. **数据导入**: 支持从 Excel/CSV 导入分析数据
