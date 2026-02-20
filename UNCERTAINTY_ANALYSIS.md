# 铀同位素模型 - 不确定度分析指南

本文档介绍U同位素海洋循环模型中的不确定度分析方法。

## 背景

海洋铀循环模型涉及多个参数，这些参数都存在不同程度的不确定度：

| 参数 | 最佳估计 | 不确定度 | 来源 |
|------|----------|----------|------|
| δ²³⁸U_river | -0.29‰ | ±0.16‰ | Andersen et al. (2016) |
| δ²³⁸U_seawater | -0.392‰ | ±0.005‰ | Tissot & Dauphas (2015) |
| Δ_sw-anox | +0.77‰ | ±0.04‰ | Stirling et al. (2015) |
| Δ_diag | +0.40‰ | 0.30-0.50‰ | Elrick et al. (2017) |

这些不确定度会通过模型传播到最终的f_anox估计中。

## 可用的不确定度分析方法

### 1. 蒙特卡洛模拟 (Monte Carlo)

通过随机采样参数分布，评估f_anox的统计分布。

```bash
python cli.py u --delta-carb -0.65 --steady-state --uncertainty mc --n-samples 50000
```

**Python API:**
```python
from systems.u import UIsotopeSystem, UncertaintyAnalyzer

u_system = UIsotopeSystem()
analyzer = UncertaintyAnalyzer(u_system)

result = analyzer.monte_carlo_steady_state(
    delta238_carb=-0.65,
    n_samples=50000,
    confidence_level=0.95
)

print(f"f_anox = {result['f_anox_mean']:.1%} ± {result['f_anox_std']:.1%}")
print(f"95% CI: [{result['f_anox_ci'][0]:.1%}, {result['f_anox_ci'][1]:.1%}]")
```

**输出解读:**
- `f_anox_mean`: f_anox的均值
- `f_anox_median`: f_anox的中位数（更稳健）
- `f_anox_std`: 标准差
- `f_anox_ci`: 置信区间（默认95%）
- `convergence/r_hat`: Gelman-Rubin统计量（<1.1表示收敛）

### 2. Bootstrap分析

适用于有多个重复测量样本的情况。

```python
measurements = np.array([-0.62, -0.58, -0.65, -0.60])  # 多个样品

result = analyzer.bootstrap_analysis(
    delta238_measurements=measurements,
    n_bootstrap=10000
)
```

### 3. 敏感性分析

识别对f_anox计算影响最大的参数。

```bash
python cli.py u --delta-carb -0.65 --steady-state --sensitivity-analysis
```

**输出示例:**
```
Parameter sensitivities (ranked by importance):
  delta_river    : range [-20.8%, +20.8%]
  delta_diag     : range [-13.0%, +13.0%]
  delta_sw_anox  : range [-7.5%, +9.2%]
```

**解读:**
- δ_river（河流输入同位素）是影响最大的参数
- 为提高f_anox估计精度，应优先改进该参数的测定

### 4. 自定义参数不确定度

根据具体研究情况调整参数不确定度：

```python
from systems.u import ParameterUncertainty

# 自定义参数分布
analyzer.config.delta_river = ParameterUncertainty(
    mean=-0.29, std=0.10, distribution='normal'
)
analyzer.config.delta_diag = ParameterUncertainty(
    mean=0.40, lower=0.20, upper=0.60, distribution='uniform'
)
```

支持的分布类型:
- `'normal'`: 正态分布，需提供mean和std
- `'uniform'`: 均匀分布，需提供lower和upper

### 5. 仅测量不确定度

当模型参数确定，仅评估测量误差时：

```python
result = analyzer.measurement_uncertainty_only(
    delta238_carb=-0.65,
    measurement_std=0.05
)
```

## 典型结果示例

### 示例1: 现代氧化海洋条件
```
输入: δ²³⁸U_carb = -0.45‰

蒙特卡洛结果 (10,000 samples):
  f_anox = 20.3% ± 15.2%
  95% CI: [0.0%, 47.5%]
```

### 示例2: 缺氧事件条件
```
输入: δ²³⁸U_carb = -0.85‰

蒙特卡洛结果 (10,000 samples):
  f_anox = 84.7% ± 8.3%
  95% CI: [67.2%, 96.1%]
```

## 不确定度来源分解

通过比较不同情况下的不确定度，可以识别主要误差来源：

| 情况 | f_anox 不确定度 | 说明 |
|------|----------------|------|
| 仅测量误差 | ±5% | 理想情况 |
| 仅Δ_sw-anox | ±8% | 分馏系数误差 |
| 仅δ_river | ±20% | 端元值误差 |
| 全参数 | ±22% | 综合效应 |

## 最佳实践建议

1. **报告不确定度**: 总是报告95%置信区间，而不仅是点估计
2. **检查收敛**: 确保蒙特卡洛模拟的R̂ < 1.1
3. **样本数量**: 至少使用10,000次采样以获得稳定结果
4. **敏感性优先**: 根据敏感性分析结果，优先改进关键参数的测定精度
5. **测量精度**: 提高δ²³⁸U测量精度可以显著降低不确定度

## 参考文献

- Andersen, M.B., et al. (2016). Closing in on the marine ²³⁸U/²³⁵U budget. *Chemical Geology*, 420, 11-22.
- Tissot, F.L., & Dauphas, N. (2015). Uranium isotopic compositions of the crust and ocean. *GCA*, 167, 113-143.
- Stirling, C.H., et al. (2015). Isotope fractionation of ²³⁸U and ²³⁵U during biologically-mediated uranium reduction. *GCA*, 163, 200-218.
- Elrick, M., et al. (2017). Global-ocean redox variation based on uranium isotope trends. *Geology*, 45(2), 163-166.
