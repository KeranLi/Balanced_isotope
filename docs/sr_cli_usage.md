# Sr同位素CLI使用指南

## 快速开始

### 1. 基本计算

计算现代海水Sr同位素比值：

```bash
python cli.py sr --calculate
```

输出：
```
Seawater ⁸⁷Sr/⁸⁶Sr = 0.70969
Modern observed:    0.70917
Difference:         +5.24 ppm
```

### 2. 自定义参数计算

```bash
python cli.py sr --calculate --F-highT 20e9 --R-river 0.705
```

参数说明：
- `--F-river`: 河流通量 (mol/yr, 默认: 47.6e9)
- `--R-river`: 河流⁸⁷Sr/⁸⁶Sr (默认: 0.71107)
- `--F-highT`: 高温热液通量 (mol/yr, 默认: 8.04e9)
- `--R-highT`: 高温热液⁸⁷Sr/⁸⁶Sr (默认: 0.7037)
- `--F-lowT`: 低温热液通量 (mol/yr, 默认: 10e9)
- `--R-lowT`: 低温热液⁸⁷Sr/⁸⁶Sr (默认: 0.7084)
- `--F-dia`: 成岩通量 (mol/yr, 默认: 3.4e9)
- `--R-dia`: 成岩⁸⁷Sr/⁸⁶Sr (默认: 0.7084)

## 蒙特卡洛随机模拟

运行5000次随机模拟（二叠纪）：

```bash
python cli.py sr --monte-carlo --scenario permian --n-runs 5000
```

输出示例：
```
Simulation Results:
  Success rate: 100.00%
  Successful runs: 5000/5000

Seawater Sr isotope evolution:
  Age (Ma)     Mean       Std
  --------------------------------
  299.0        0.71822    0.00956
  274.3        0.71822    0.00956
  252.0        0.71822    0.00956

Parameter statistics:
  F_river: 1.001e+11 ± 5.241e+10
  R_river: 7.258e-01 ± 1.384e-02
```

保存结果到CSV：

```bash
python cli.py sr --monte-carlo --n-runs 5000 --output sr_monte_carlo.csv
```

## 情景分析

Wang et al. (2021) 论文中的7种情景：

| 情景 | 固定参数 | 可变参数 |
|------|----------|----------|
| scenario1 | 成岩参数 | 其他所有 |
| scenario2 | 河流同位素 + 成岩 | 其他 |
| scenario3 | 河流通量 + 成岩 | 其他 |
| scenario4 | 高温热液通量 + 成岩 | 其他 |
| scenario5 | 低温热液通量 + 成岩 | 其他 |
| scenario6 | 所有热液通量 + 成岩 | 其他 |
| scenario7 | 河流参数 + 成岩 | 热液参数 |

运行情景分析：

```bash
python cli.py sr --scenario-analysis scenario7 --n-runs 3000
```

输出示例：
```
Results:
  Successful runs: 3000
  Success rate: 100.00%

Parameter statistics:
  F_hydrothermal_lowT:
    Mean:   2.213e+10
    Median: 2.218e+10
    Range:  [2.618e+09, 3.996e+10]
```

## 敏感性分析

分析河流同位素组成的影响：

```bash
python cli.py sr --sensitivity R_river --sens-range 0.705 0.725 --n-steps 50
```

分析高温热液通量的影响：

```bash
python cli.py sr --sensitivity F_highT --sens-range 2e9 35e9 --n-steps 50
```

输出示例：
```
Parameter range: [7.0500e-01, 7.2500e-01]
Sr ratio range:  [0.70551, 0.71930]

Key points:
  R_river=7.0500e-01: R_sw=0.70551
  R_river=7.1611e-01: R_sw=0.71317
  R_river=7.2500e-01: R_sw=0.71930
```

## 查看系统信息

```bash
python cli.py info sr
```

## 列出所有可用体系

```bash
python cli.py list
```

## 完整示例流程

### 示例1: 探索二叠纪Sr同位素驱动机制

```bash
# 1. 查看Sr系统信息
python cli.py info sr

# 2. 运行蒙特卡洛模拟探索参数空间
python cli.py sr --monte-carlo --scenario permian --n-runs 5000 \
    --time-span 299 252 --output permian_mc.csv

# 3. 运行情景7（固定河流，可变热液）
python cli.py sr --scenario-analysis scenario7 --n-runs 3000 \
    --time-span 299 252 --output scenario7.csv

# 4. 分析热液通量敏感性
python cli.py sr --sensitivity F_highT --sens-range 2e9 35e9 \
    --output sensitivity_highT.csv
```

### 示例2: 现代海洋Sr循环

```bash
# 基本计算
python cli.py sr --calculate

# 测试不同河流输入
python cli.py sr --calculate --R-river 0.708
python cli.py sr --calculate --R-river 0.714
python cli.py sr --calculate --R-river 0.720

# 敏感性分析
python cli.py sr --sensitivity R_river --sens-range 0.708 0.720
```

## 输出文件格式

CSV文件包含以下列：

**蒙特卡洛模拟输出：**
- `time_Ma`: 时间（百万年前）
- `mean_Sr_ratio`: 平均⁸⁷Sr/⁸⁶Sr
- `std_Sr_ratio`: 标准差
- `p2.5`: 2.5%分位数
- `p97.5`: 97.5%分位数

**敏感性分析输出：**
- 参数名（如`R_river`）
- `Sr_ratio`: 海水Sr同位素比值
- `sensitivity`: 敏感性（dR_sw/dparam）

## 注意事项

1. **时间单位**: 所有时间参数单位为百万年（Ma）
2. **通量单位**: 通量单位为mol/yr
3. **蒙特卡洛运行次数**: 建议至少3000次以获得稳定的统计结果
4. **过滤**: 当前CLI版本不进行观测数据过滤，所有运行都被视为成功

## 参考文献

Wang et al. (2021). Revisiting the Permian seawater Sr/86Sr record: New perspectives from brachiopod proxy data and stochastic oceanic box models. *Earth-Science Reviews*, 218, 103679.
