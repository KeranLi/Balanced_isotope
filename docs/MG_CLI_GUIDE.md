# Mg同位素CLI命令参数详解

本文档详细说明 `python cli.py mg` 命令的所有参数，基于Kasemann等(2014)论文的Mg同位素风化通量模型。

---

## 快速示例

```bash
python cli.py mg --file data/Nie_Section_B.xlsx \
  --use-raw-values \
  --delta-silicate 0.1 \
  --delta-carbonate -1.5 \
  --weathering-simulation \
  --duration 5 \
  --output results/Nie_B_full.csv
```

---

## 参数分类详解

### 1. 输入文件参数

#### `--file FILE`
- **功能**: 指定输入的Excel或CSV文件路径
- **必需**: 是（用于批量处理）
- **示例**: `--file data/Nie_Section_A.xlsx`
- **说明**: 文件必须包含 `delta_26_Mg_iso` 列，可选包含 `delta_26_Mg_iso_2sd` 不确定度列

#### `--column {delta_25_Mg_iso, delta_26_Mg_iso}`
- **功能**: 指定使用的同位素列
- **默认值**: `delta_26_Mg_iso`
- **说明**: 通常使用δ²⁶Mg进行风化分析

---

### 2. 数据解释参数

#### `--use-raw-values`
- **功能**: 直接使用原始值计算，不进行校正
- **默认值**: False（自动启用校正）
- **使用场景**: 
  - 数据已经是可直接对比端元的形式
  - 需要自定义解释数据类型时
- **示例**: `--use-raw-values`

#### `--data-type {carbonate, seawater, river}`
- **功能**: 指定输入数据的类型
- **默认值**: `carbonate`
- **选项说明**:
  - `carbonate`: 碳酸盐沉积物数据（自动应用-2.7‰校正反推海水）
  - `seawater`: 海水数据（无校正）
  - `river`: 河水数据（直接计算风化比例）
- **示例**: `--data-type river`

#### `--seawater-correction FLOAT`
- **功能**: 自定义海水校正系数
- **默认值**: -2.7‰（碳酸盐-海水分馏）
- **说明**: 仅在 `--data-type carbonate` 时生效
- **示例**: `--seawater-correction -3.0`

#### `--apply-offset FLOAT`
- **功能**: 对所有数据值应用全局偏移校正
- **默认值**: None
- **使用场景**: 不同参考标准间的转换
- **示例**: `--apply-offset -3.0`（将所有值减3‰）

---

### 3. 端元定义参数

#### `--delta-silicate FLOAT`
- **功能**: 自定义硅酸盐风化端元δ²⁶Mg值
- **默认值**: -0.3‰（上地壳平均）
- **论文参考值**: 
  - 现代硅酸盐: -0.3‰
  - Cryogenian硅酸盐: ~-0.3‰
- **使用建议**: 根据研究区域地质背景调整
- **示例**: `--delta-silicate 0.1`

#### `--delta-carbonate FLOAT`
- **功能**: 自定义碳酸盐风化端元δ²⁶Mg值
- **默认值**: -2.5‰（典型碳酸盐岩）
- **论文参考值**:
  - Cryogenian碳酸盐: -2.5‰
  - 现代碳酸盐: -4.3‰
- **使用建议**: 根据研究层位碳酸盐特征调整
- **示例**: `--delta-carbonate -1.5`

#### `--delta-seawater FLOAT`
- **功能**: 海水δ²⁶Mg值（用于风化比例计算）
- **默认值**: -0.83‰（现代海水）
- **论文参考值**:
  - 现代海水: -0.83‰
  - Cryogenian海水: 变化范围大（-0.5 to +1.0‰）
- **示例**: `--delta-seawater -0.5`

---

### 4. 风化模拟参数

#### `--weathering-simulation`
- **功能**: 运行风化体制转变时间演化模拟
- **默认值**: False
- **说明**: 基于输入数据的风化比例范围进行模拟
- **示例**: `--weathering-simulation`

#### `--f-initial FLOAT`
- **功能**: 初始硅酸盐风化比例
- **默认值**: 0.2（20%硅酸盐，80%碳酸盐）
- **范围**: 0.0 ~ 1.0
- **说明**: 仅用于 `--weathering-simulation`，不用于数据驱动的模拟
- **示例**: `--f-initial 0.3`

#### `--f-final FLOAT`
- **功能**: 最终硅酸盐风化比例
- **默认值**: 0.8（80%硅酸盐，20%碳酸盐）
- **范围**: 0.0 ~ 1.0
- **说明**: 仅用于 `--weathering-simulation`
- **示例**: `--f-final 0.9`

#### `--transition-mode {linear, exponential}`
- **功能**: 风化转变的时间模式
- **默认值**: `linear`
- **选项说明**:
  - `linear`: 线性转变
  - `exponential`: 指数衰减/增长转变
- **示例**: `--transition-mode exponential`

#### `--transition-start FLOAT`
- **功能**: 风化转变开始时间
- **默认值**: 0.0 Myr
- **单位**: 百万年（Ma）
- **示例**: `--transition-start 0.5`

#### `--transition-end FLOAT`
- **功能**: 风化转变结束时间
- **默认值**: 2.0 Myr
- **单位**: 百万年（Ma）
- **示例**: `--transition-end 3.0`

#### `--duration FLOAT`
- **功能**: 总模拟时长
- **默认值**: 5.0 Myr
- **单位**: 百万年（Ma）
- **说明**: 模拟的时间跨度
- **示例**: `--duration 10`

#### `--flux-multiplier FLOAT`
- **功能**: 风化通量放大系数（相对于现代）
- **默认值**: 1.0（现代通量）
- **论文参考值**:
  - 正常风化: 1×
  - Cryogenian高风化: 6-9×
- **示例**: `--flux-multiplier 6`

---

### 5. Cryogenian论文情景参数

#### `--cryogenian-scenario`
- **功能**: 运行Kasemann等(2014)论文中的Cryogenian冰期后情景
- **默认值**: False
- **情景设定**:
  - Phase 1 (0-0.5 Myr): 混合风化，9×现代通量
  - Phase 2 (0.5 Myr后): 硅酸盐主导，6×现代通量
- **示例**: `--cryogenian-scenario`

#### `--cryogenian-duration FLOAT`
- **功能**: Cryogenian情景模拟时长
- **默认值**: 3.0 Myr
- **单位**: 百万年（Ma）
- **示例**: `--cryogenian-duration 5`

---

### 6. 反演计算参数

#### `--inverse`
- **功能**: 从观测数据反演风化通量历史
- **默认值**: False
- **说明**: 需要配合 `--file` 使用
- **示例**: `--inverse --file data/samples.xlsx`

---

### 7. 风化比例计算参数

#### `--weathering-ratio`
- **功能**: 计算单个样品的风化端元比例
- **默认值**: False
- **配合参数**: `--delta-sample`
- **示例**: `--weathering-ratio --delta-sample -2.0`

#### `--delta-sample FLOAT`
- **功能**: 样品δ²⁶Mg值
- **默认值**: None
- **单位**: ‰
- **说明**: 用于单个样品的风化比例计算
- **示例**: `--delta-sample -1.5`

---

### 8. 输出参数

#### `--output FILE`
- **功能**: 输出CSV文件路径
- **默认值**: None（仅屏幕输出）
- **说明**: 保存分析结果到CSV文件
- **示例**: `--output results/analysis.csv`

#### `--plot`
- **功能**: 生成并保存图表
- **默认值**: False
- **说明**: 自动生成演化曲线图（待实现）
- **示例**: `--plot`

#### `--n-points INT`
- **功能**: 模拟时间点数
- **默认值**: 500
- **说明**: 时间演化的分辨率
- **示例**: `--n-points 1000`

---

## 完整参数速查表

| 参数 | 类型 | 默认值 | 说明 |
|-----|------|-------|------|
| `--file` | str | None | 输入文件路径 |
| `--column` | choice | delta_26_Mg_iso | 同位素列名 |
| `--use-raw-values` | flag | False | 使用原始值 |
| `--data-type` | choice | carbonate | 数据类型 |
| `--seawater-correction` | float | -2.7 | 海水校正系数 |
| `--apply-offset` | float | None | 全局偏移校正 |
| `--delta-silicate` | float | -0.3 | 硅酸盐端元 |
| `--delta-carbonate` | float | -2.5 | 碳酸盐端元 |
| `--delta-seawater` | float | -0.83 | 海水值 |
| `--weathering-simulation` | flag | False | 运行模拟 |
| `--f-initial` | float | 0.2 | 初始f_silicate |
| `--f-final` | float | 0.8 | 最终f_silicate |
| `--transition-mode` | choice | linear | 转变模式 |
| `--transition-start` | float | 0.0 | 转变开始(Myr) |
| `--transition-end` | float | 2.0 | 转变结束(Myr) |
| `--duration` | float | 5.0 | 总时长(Myr) |
| `--flux-multiplier` | float | 1.0 | 通量放大系数 |
| `--cryogenian-scenario` | flag | False | Cryogenian情景 |
| `--cryogenian-duration` | float | 3.0 | Cryogenian时长 |
| `--inverse` | flag | False | 反演计算 |
| `--weathering-ratio` | flag | False | 风化比例计算 |
| `--delta-sample` | float | None | 样品δ²⁶Mg |
| `--output` | str | None | 输出文件 |
| `--plot` | flag | False | 生成图表 |
| `--n-points` | int | 500 | 时间点数 |

---

## 输出文件格式

### 主分析文件 (`*_analysis.csv`)

| 列名 | 说明 |
|-----|------|
| `sample_index` | 样品序号 |
| `delta_26Mg_raw` | 原始δ²⁶Mg值 |
| `delta_26Mg_corrected` | 校正后δ²⁶Mg值 |
| `delta_26Mg_1sigma` | 不确定度(1σ) |
| `f_silicate` | 硅酸盐风化比例 |
| `f_carbonate` | 碳酸盐风化比例 |
| `delta_river_inferred` | 推断的河流δ²⁶Mg |

### 模拟文件 (`*_simulation.csv`)

| 列名 | 说明 |
|-----|------|
| `time_myr` | 时间 (Myr) |
| `M_sw_mol` | 海水Mg储库 (mol) |
| `delta_sw_permil` | 海水δ²⁶Mg (‰) |
| `f_silicate` | 硅酸盐风化比例 |
| `f_carbonate` | 碳酸盐风化比例 |
| `flux_multiplier` | 通量放大系数 |
| `delta_river_permil` | 河流δ²⁶Mg (‰) |

---

## 参考标准说明

**D3MS = DSM3**
- **全称**: Dead Sea Magnesium Standard
- **说明**: 这是同一个参考标准的两种写法
- **δ²⁶Mg值**: 定义为0‰
- **同位素比值**: ²⁶Mg/²⁴Mg = 0.13932

---

## 参考文献

1. Kasemann, S.A., et al. (2014). Continental weathering following a Cryogenian glaciation: Evidence from calcium and magnesium isotopes. *Earth and Planetary Science Letters*, 396, 66-77.

2. Tipper, E.T., et al. (2006). The magnesium isotope budget of the modern ocean: Constraints from riverine magnesium isotope ratios. *Earth and Planetary Science Letters*, 250(1-2), 241-253.

---

*文档版本: 1.0 | 更新日期: 2024*
