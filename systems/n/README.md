# 氮同位素体系 (Nitrogen Isotope System)

## 概述

本模块实现了基于 Kang et al. (2023) 和 Ma et al. (2025) 的双箱稳态氮循环模型，用于定量重建地质历史时期海洋硝酸盐可利用性。

## 核心模型

### 双箱模型结构

```
大气 N₂ (δ¹⁵N = 0‰)
    ↓ F_fix (固氮作用, ε_fix = -2‰ ~ +1‰)
    
N_fixer/ammonium 储库 ──F_remin──→ N_assimilator/nitrate 储库
(固氮生物-铵)              (再矿化)      (硝酸盐同化生物-硝酸盐)
    │                                        │
    ↓ F_fixer_burial                         ↓ F_assimilator_burial
    │              沉积物 δ¹⁵N_sed            │
    └─────────────── (混合) ─────────────────┘
                              ↑
                    F_wcd (水柱反硝化, ε_wcd = -30‰ ~ -22‰)
                    F_sd (沉积反硝化, ε_sd = 0‰)
```

### 关键方程

**质量平衡：**
- 固氮生物-铵储库: F_fix = F_remin + F_fixer_burial
- 硝酸盐同化生物-硝酸盐储库: F_remin = F_wcd + F_sd + F_assimilator_burial

**同位素平衡：**
- δ¹⁵N_ammonium = δ¹⁵N_atmosphere + ε_fix
- δ¹⁵N_nitrate = δ¹⁵N_ammonium - (F_wcd×ε_wcd + F_sd×ε_sd) / F_remin

**沉积物同位素：**
- δ¹⁵N_sed = (1-f_assimilator) × δ¹⁵N_ammonium + f_assimilator × δ¹⁵N_nitrate
- f_assimilator = F_assimilator_burial / F_total_burial (硝酸盐占比)

### 非线性关系

δ¹⁵N_sed 与 f_assimilator 存在非线性关系：
- f = 0: δ¹⁵N_sed ≈ -0.4‰ (纯固氮)
- f ≈ 0.48: δ¹⁵N_sed 达到最大值 (~5.7‰)
- f = 1: δ¹⁵N_sed ≈ -0.7‰ (纯硝酸盐同化)

## 快速开始

```python
from systems.n import NIsotopeSystem

# 创建模型实例
n_system = NIsotopeSystem(scenario='modern')

# 正向模型: 从硝酸盐占比计算沉积物同位素
delta15N = n_system.forward_model(f_assimilator=0.3)
print(f"δ¹⁵N_sed = {delta15N:.2f}‰")

# 反向模型: 从沉积物同位素反演硝酸盐占比
result = n_system.inverse_model(delta15N_sed=3.0)
print(f"f_assimilator = {result['f_assimilator']:.3f}")

# 蒙特卡洛不确定性分析
mc_result = n_system.monte_carlo_simulation(f_assimilator=0.2)
print(f"δ¹⁵N_sed = {mc_result['delta15N_sed_mean']:.2f} ± {mc_result['delta15N_sed_std']:.2f}‰")
```

## 主要功能

### 1. 正向模型 (`forward_model`)
从硝酸盐占比计算沉积物氮同位素组成。

### 2. 反向模型 (`inverse_model`)
从沉积物氮同位素反演硝酸盐占比。注意：由于非线性关系，同一δ¹⁵N值可能对应两个f值，在缺氧环境下通常选择较小的值。

### 3. 蒙特卡洛模拟 (`monte_carlo_simulation`)
评估分馏系数不确定性对结果的影响。

### 4. 关系曲线计算 (`calculate_f_assimilator_curve`)
计算 f_assimilator 与 δ¹⁵N_sed 的完整关系曲线，用于绘图。

### 5. 时间序列分析 (配合 `toolkit.math.statistics`)
- Bootstrap 置信区间估计
- LOWESS 平滑
- 变点检测

## 应用场景

### 1. 新元古代真核生物崛起 (Kang et al. 2023)

```python
n_system = NIsotopeSystem(scenario='neoproterozoic')

# 800 Ma 前 (缺氧环境)
result_before = n_system.inverse_model(delta15N_sed=1.5, f_range=(0, 0.2))
print(f"800 Ma 前 f_assimilator ≈ {result_before['f_assimilator']:.2f}")

# 800 Ma 后 (氧化环境)
result_after = n_system.inverse_model(delta15N_sed=4.5, f_range=(0.1, 0.5))
print(f"800 Ma 后 f_assimilator ≈ {result_after['f_assimilator']:.2f}")
```

### 2. 早三叠世生态系统复苏 (Ma et al. 2025)

```python
n_system = NIsotopeSystem(scenario='early_triassic')

# 分阶段分析
stages = {
    'I (Griesbachian-Smithian)': (0, 0.1),
    'II (Spathian early)': (0.15, 0.25),
    'III (Spathian late)': (0.05, 0.15)
}

for stage_name, f_range in stages.items():
    result = n_system.inverse_model(delta15N_sed=2.5, f_range=f_range)
    print(f"{stage_name}: f = {result['f_assimilator']:.3f}")
```

## 参数说明

### 通量参数 (Tg N/a)
| 参数 | 现代值 | 说明 |
|------|--------|------|
| F_fix | 205 | 固氮通量 |
| F_total_burial | 25 | 总埋藏通量 |
| F_wcd | 140 | 水柱反硝化通量 |
| F_sd | 40 | 沉积反硝化通量 |

### 分馏系数 (‰)
| 过程 | ε 范围 | 说明 |
|------|--------|------|
| 固氮 | -2 ~ +1 | 生物固氮 |
| 水柱反硝化 | -30 ~ -22 | 缺氧水体中 |
| 沉积反硝化 | 0 | 沉积物中 |

## 参考文献

1. **Kang et al. (2023)** Nitrate limitation in early Neoproterozoic oceans 
   delayed the ecological rise of eukaryotes. *Science Advances*, 9, eade9647.

2. **Ma et al. (2025)** Prolonged nitrate depletion delayed marine ecosystem 
   recovery after the end-Permian mass extinction. *Science China Earth Sciences*, 
   68(9), 3035-3049.

3. **Stüeken et al. (2016)** Nitrogen isotope evidence for anoxic deep ocean 
   during the Late Archean. *Earth and Planetary Science Letters*, 447, 130-138.

4. **Kipp et al. (2018)** Pervasive aerobic nitrogen cycling in the surface 
   ocean across the Paleoproterozoic Era. *Earth and Planetary Science Letters*, 
   500, 117-126.
