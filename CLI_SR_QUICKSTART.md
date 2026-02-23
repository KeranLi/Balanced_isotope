# Sr同位素CLI快速入门

## 安装与运行

无需安装，直接在项目目录下运行：

```bash
cd /Users/keranli/Desktop/Coding/Balanced_isotope
python cli.py sr [选项]
```

## 最常用命令

### 1. 基本质量平衡计算

```bash
# 现代海水Sr同位素
python cli.py sr --calculate

# 自定义参数
python cli.py sr --calculate --F-highT 20e9 --R-river 0.705
```

### 2. 蒙特卡洛随机模拟

```bash
# 二叠纪情景（5000次运行）
python cli.py sr --monte-carlo --scenario permian --n-runs 5000

# 保存结果
python cli.py sr --monte-carlo --n-runs 5000 --output results.csv
```

### 3. 情景分析（Wang et al., 2021）

```bash
# Scenario 7: 固定河流参数，探索热液影响
python cli.py sr --scenario-analysis scenario7 --n-runs 3000
```

### 4. 敏感性分析

```bash
# 河流同位素影响
python cli.py sr --sensitivity R_river --sens-range 0.705 0.725

# 高温热液通量影响  
python cli.py sr --sensitivity F_highT --sens-range 2e9 35e9
```

## 关键参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--F-river` | 47.6e9 | 河流通量 (mol/yr) |
| `--R-river` | 0.71107 | 河流⁸⁷Sr/⁸⁶Sr |
| `--F-highT` | 8.04e9 | 高温热液通量 (mol/yr) |
| `--R-highT` | 0.7037 | 高温热液⁸⁷Sr/⁸⁶Sr |
| `--F-lowT` | 10e9 | 低温热液通量 (mol/yr) |
| `--R-lowT` | 0.7084 | 低温热液⁸⁷Sr/⁸⁶Sr |

## 典型工作流程

```bash
# 1. 查看系统信息
python cli.py info sr

# 2. 计算现代海水Sr同位素
python cli.py sr --calculate

# 3. 探索热液活动影响（敏感性分析）
python cli.py sr --sensitivity F_highT --sens-range 2e9 35e9 --output sens_highT.csv

# 4. 运行蒙特卡洛模拟探索参数空间
python cli.py sr --monte-carlo --n-runs 5000 --output mc_results.csv

# 5. 运行特定情景分析
python cli.py sr --scenario-analysis scenario7 --n-runs 3000 --output scenario7.csv
```

## 输出文件

所有 `--output` 选项生成CSV文件，可用Excel或Python分析：

```python
import pandas as pd

# 读取蒙特卡洛结果
df = pd.read_csv('mc_results.csv')
print(df.head())
print(df[['mean_Sr_ratio', 'std_Sr_ratio']].describe())
```

## 获取帮助

```bash
# 主帮助
python cli.py --help

# Sr同位素子命令帮助
python cli.py sr --help

# 查看系统详细信息
python cli.py info sr
```

## 参考

Wang et al. (2021). Earth-Science Reviews 218: 103679
