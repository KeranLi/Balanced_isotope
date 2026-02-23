# 批量处理功能使用指南

本文档介绍如何使用批量处理功能从Excel文件读取同位素数据并进行分析。

## 功能概述

批量处理功能允许用户：
1. 从Excel/CSV文件读取多个样品的同位素数据
2. 自动识别数据列
3. 批量计算并输出结果
4. 支持不确定度分析
5. 将结果保存回Excel文件

## 支持的文件格式

- **Excel**: `.xlsx`, `.xls`
- **CSV**: `.csv`

## 输入文件格式要求

### U同位素数据格式

| 列名 | 说明 | 必需 |
|------|------|------|
| `sample_id` | 样品编号 | 否 |
| `elevation_m` / `depth_m` / `age_ma` | 地层信息 | 否 |
| `delta_238_u` / `delta238U` / `d238U` | δ²³⁸U值 (‰) | **是** |
| `delta_238_u_std` / `delta238U_std` | 不确定度 | 否 |

**示例Excel内容**:
```
sample_id | elevation_m | delta_238_u | delta_238_u_std | description
----------|-------------|-------------|-----------------|------------
BS-01     | 0.0         | -0.32       | 0.05            | background
BS-02     | 2.5         | -0.45       | 0.06            | anoxia
BS-03     | 5.0         | -0.55       | 0.04            | peak_anoxia
```

### Mg同位素数据格式

| 列名 | 说明 | 必需 |
|------|------|------|
| `delta_26_Mg_iso` / `delta26Mg` | δ²⁶Mg值 (‰) | **是** |
| `delta_26_Mg_iso_2sd` | 不确定度 | 否 |
| `delta_25_Mg_iso` | δ²⁵Mg值 (‰) | 否 |

### C同位素数据格式

| 列名 | 说明 | 必需 |
|------|------|------|
| `delta_13_c` / `delta13C` | δ¹³C值 (‰) | **是** |
| `delta_13_c_std` | 不确定度 | 否 |

### N同位素数据格式

| 列名 | 说明 | 必需 |
|------|------|------|
| `delta_15_n` / `delta15N` | δ¹⁵N值 (‰) | **是** |
| `delta_15_n_std` | 不确定度 | 否 |

## CLI使用

### 基本批量处理

```bash
# U同位素批量处理
python cli.py u --file data/uranium_data.xlsx --output results/u_results.xlsx

# Mg同位素批量处理
python cli.py mg --file data/magnesium_data.xlsx --output results/mg_results.xlsx

# C同位素批量处理
python cli.py c --file data/carbon_data.xlsx --output results/c_results.xlsx
```

### 带不确定度分析的批量处理

```bash
# 启用蒙特卡洛不确定度分析（默认1000次采样）
python cli.py u --file data/uranium_data.xlsx --output results/u_results.xlsx

# 禁用不确定度分析（仅计算点估计）
python cli.py u --file data/uranium_data.xlsx --output results/u_results.xlsx --no-uncertainty
```

### 指定情景

```bash
# 使用特定地质情景
python cli.py u --file data/ff_data.xlsx --scenario frasnian_famennian --output results/ff_results.xlsx
```

## Python API使用

### 基本批量处理

```python
from toolkit.io import BatchProcessor

# 创建处理器
processor = BatchProcessor(
    element='u',
    scenario='modern',
    apply_diagenetic_correction=True,
    delta_diag=0.4,
    include_uncertainty=True,
    n_monte_carlo=1000
)

# 处理文件
results_df = processor.process_file(
    'data/uranium_data.xlsx',
    output_path='results/u_results.xlsx',
    show_progress=True
)

# 查看结果
print(results_df[['sample_id', 'f_anox', 'anoxic_area_percent']])
```

### 处理结果摘要

```python
from toolkit.io import BatchProcessor

processor = BatchProcessor(element='u')
results = processor.process_file('data/input.xlsx')

# 获取摘要统计
summary = processor.get_summary(results)
print(f"Success rate: {summary['success_rate']:.1%}")
print(f"Mean f_anox: {summary['f_anox_mean']:.1%}")
```

### 数据读取和验证

```python
from toolkit.io import DataHandler

# 读取数据
handler = DataHandler('data/uranium_data.xlsx')
data = handler.read()

# 查看数据信息
info = handler.get_info()
print(info)

# 验证数据
is_valid, errors = handler.validate_data('u')
if not is_valid:
    print(f"Validation errors: {errors}")

# 获取同位素数据
u_data = handler.get_isotope_data('u')
print(u_data['delta'])  # δ²³⁸U值数组
```

### 自定义列名映射

```python
from toolkit.io import DataHandler, ColumnMapping

# 自定义列名映射
mapping = ColumnMapping(
    delta238_u=['U238', 'U_238', 'delta_U'],
    sample_id=['ID', 'Sample_Name'],
    elevation=['Height', 'Strat_Height']
)

handler = DataHandler('data/custom_format.xlsx', column_mapping=mapping)
data = handler.read()
```

## 输出文件格式

### U同位素输出列

| 列名 | 说明 |
|------|------|
| `f_anox` | 缺氧汇比例 |
| `f_oxic` | 氧化汇比例 |
| `delta238_seawater` | 海水δ²³⁸U |
| `delta238_carb_corrected` | 校正后碳酸盐δ²³⁸U |
| `anoxic_area_percent` | 估算缺氧面积(%) |
| `f_anox_std` | f_anox标准差 |
| `f_anox_ci_lower` | 95%置信区间下限 |
| `f_anox_ci_upper` | 95%置信区间上限 |
| `processing_success` | 处理是否成功 |
| `processing_error` | 错误信息 |

### Mg同位素输出列

| 列名 | 说明 |
|------|------|
| `f_carbonate` | 碳酸盐风化比例 |
| `f_silicate` | 硅酸盐风化比例 |
| `processing_success` | 处理是否成功 |

### N同位素输出列

| 列名 | 说明 |
|------|------|
| `f_assimilator` | 硝酸盐同化比例 |
| `delta15_sed_calculated` | 计算沉积物δ¹⁵N |
| `residual` | 残差 |
| `processing_success` | 处理是否成功 |

## 示例工作流

### 1. 准备数据

创建包含样品信息的Excel文件：

```python
import pandas as pd

# 创建示例数据
data = {
    'sample_id': ['S-01', 'S-02', 'S-03', 'S-04', 'S-05'],
    'elevation_m': [0, 5, 10, 15, 20],
    'delta_238_u': [-0.32, -0.45, -0.55, -0.40, -0.28],
    'delta_238_u_std': [0.05, 0.06, 0.04, 0.05, 0.04]
}

df = pd.DataFrame(data)
df.to_excel('my_data.xlsx', index=False)
```

### 2. 批量处理

```bash
python cli.py u --file my_data.xlsx --output results.xlsx
```

### 3. 分析结果

```python
import pandas as pd
import matplotlib.pyplot as plt

# 读取结果
results = pd.read_excel('results.xlsx')

# 绘制f_anox随高程变化
plt.figure(figsize=(10, 6))
plt.plot(results['elevation_m'], results['f_anox'] * 100, 'o-')
plt.xlabel('Elevation (m)')
plt.ylabel('f_anox (%)')
plt.title('Anoxic Fraction vs Stratigraphic Height')
plt.grid(True)
plt.savefig('f_anox_profile.png')

# 统计摘要
print(f"Mean f_anox: {results['f_anox'].mean():.1%}")
print(f"Range: [{results['f_anox'].min():.1%}, {results['f_anox'].max():.1%}]")
```

## 常见问题

### Q: 如何批量处理多个文件？

```bash
# 使用shell循环
for file in data/*.xlsx; do
    python cli.py u --file "$file" --output "results/$(basename $file)"
done
```

### Q: 如何处理自定义列名？

使用Python API自定义ColumnMapping（见上文示例）。

### Q: 如何处理非常大的数据文件？

```python
# 使用pandas分块读取
import pandas as pd

chunks = pd.read_excel('large_file.xlsx', chunksize=1000)
for chunk in chunks:
    # 处理每个块
    pass
```

### Q: 某些样品处理失败怎么办？

检查输出文件中的`processing_success`和`processing_error`列，可以筛选出失败的样品进行排查。

```python
results = pd.read_excel('results.xlsx')
failed = results[~results['processing_success']]
print(failed[['sample_id', 'processing_error']])
```

## 技术细节

- 使用`pandas`进行数据读写
- 使用`tqdm`显示进度条（如果已安装）
- 自动检测常见列名变体
- 支持多种不确定度分布（正态、均匀）
- 蒙特卡洛模拟支持收敛性诊断
