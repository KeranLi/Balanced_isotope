# 同位素质量平衡分馏模型框架

一个用于地球化学同位素体系质量平衡和分馏计算的可扩展框架。

## 特性

- **三层架构设计**：工具包层、同位素体系层、应用示例层
- **多元素支持**：Mg、C、S、Sr、Nd等同位素体系
- **统一接口**：所有同位素体系继承相同基类
- **可扩展**：轻松添加新的同位素体系

## 已实现的同位素体系

| 元素 | 应用场景 | 状态 |
|-----|---------|------|
| Mg  | 风化-沉积体系，碳酸盐/硅酸盐风化比例 | ✓ 已实现 |
| C   | 碳循环，DOC氧化与碳同位素负漂 | ✓ 已实现 |
| S   | 硫循环，硫酸盐还原 | ○ 模板 |
| Sr  | 风化示踪 | ○ 计划中 |
| Nd  | 洋流循环 | ○ 计划中 |

## 快速开始

### 安装依赖

```bash
pip install numpy scipy pandas matplotlib openpyxl
```

### 运行示例

```bash
# 统一入口（推荐）
python run.py all          # 运行所有演示
python run.py mg           # Mg同位素分析
python run.py c            # C同位素分析
python run.py tools        # 核心工具演示

# 具体案例
python examples/basic_usage.py

# Nie碎屑组分Mg同位素剖面（输入行按深度递增）
python cli.py mg --component-type siliciclastic \
  --file data/Nie_Section_A.xlsx \
  --output results/Nie_A_siliciclastic_flux.csv
python cli.py mg --component-type siliciclastic \
  --file data/Nie_Section_B.xlsx \
  --output results/Nie_B_siliciclastic_flux.csv

# Jupyter Notebook 交互式示例
pip install jupyter

# 基础用法
jupyter notebook examples/notebooks/01_basic_usage.ipynb

# N同位素论文复现 (Kang et al. 2023 & Ma et al. 2025)
jupyter notebook examples/notebooks/02_nitrogen_isotope_reproduction.ipynb

# CLI方式
python cli.py list
python cli.py info mg
python cli.py mg --weathering-ratio
python cli.py c --F-odoc 4e18
```

### 代码示例

```python
# Mg同位素风化分析
from systems.mg import MgIsotopeSystem

mg = MgIsotopeSystem()
ratios = mg.calculate_weathering_ratio(delta_sample=-2.5, delta_seawater=-0.83)
print(f"碳酸盐风化比例: {ratios['f_carbonate']:.1%}")

# C同位素DOC模型
from systems.c import CIsotopeSystem

c = CIsotopeSystem(scenario='dice')
result = c.solve_steady_state(F_odoc=4e18)
print(f"δ¹³C_carb: {result.get('delta13C_carb'):.2f}‰")
```

碎屑组分模型输出的绝对风化通量依赖总河流 Mg 通量先验。默认采用
Hu et al. (2023) 的长江值 `3.0e11 mol/yr`；用于其他流域或古剖面时，
应通过 `--river-mg-flux` 提供研究区约束。

在岩石供应量和原岩 Mg 含量不变时，碎屑剖面的首选相对指标是
`Mg release fraction = 1-f_Mg`。CLI 直接绘制该基线无关曲线，并自动检测
两状态变点；转变前后状态之比给出风化通量改变倍数。没有年龄数据时，
纵轴只表示从浅到深的相对地层位置。

碎屑剖面默认在每次 Monte Carlo 中为全部样品共享一个原岩端元，并从
`Uniform(-0.45, -0.25‰)` 抽样。可用 `--delta-protolith VALUE` 固定端元，
或用 `--delta-protolith-range LOW HIGH` 修改先验范围。

## 项目结构

```
Balanced_isotope/
├── toolkit/           # 🔧 工具包层（数学、物理、同位素公式）
├── systems/           # 🧪 同位素体系层（各元素实现）
├── examples/          # 📚 应用示例层
├── data/               # 本地输入；研究数据不纳入 Git
├── results/            # 本地生成结果；不纳入 Git
├── reference/          # 本地文献 PDF；不纳入 Git
├── attachments/        # 本地 Word/PDF 附件；不纳入 Git
├── tests/              # 回归测试
├── run.py              # 统一入口
└── cli.py              # 命令行接口
```

### GitHub 文件边界

`.gitignore` 会保留示例数据、源码、测试和说明文档，同时忽略研究数据、
生成结果、文献 PDF、附件稿、Python 缓存和系统临时文件。Nie、Sr 等研究
数据仍可在本机 `data/` 中使用；结果仍保存在 `results/`，但不会随项目
提交到 GitHub。若要与合作者共享数据，应使用受控的数据共享渠道，而不是
用 `git add -f` 绕过忽略规则。

详见 [ARCHITECTURE.md](docs/ARCHITECTURE.md)

## 原代码迁移

| 原文件 | 新位置 | 说明 |
|--------|--------|------|
| `main.py` + `calculation/` | `systems/mg/` | Mg同位素体系 |
| `Carbon_Modeling/Oxi_Est.py` | `examples/dice_doc_oxidation.py` | DICE事件DOC氧化示例 |
| `utils/` | `toolkit/` | 通用工具 |

迁移说明暂未作为当前项目文件保留。

## 文档

| 文档 | 内容 |
|------|------|
| [ARCHITECTURE.md](docs/ARCHITECTURE.md) | 项目架构设计 |
| [MG_CLI_GUIDE.md](docs/MG_CLI_GUIDE.md) | Mg同位素CLI使用指南 |
| [UNCERTAINTY_ANALYSIS.md](docs/UNCERTAINTY_ANALYSIS.md) | 不确定度分析方法 |
| [docs/ODESOLVER_GUIDE.md](docs/ODESOLVER_GUIDE.md) | **ODE求解器使用指南** |
| [BATCH_PROCESSING.md](docs/BATCH_PROCESSING.md) | 批量数据处理 |

## 参考

- Mg同位素模型参考：原 `calculation/` 和 `main.py`
- C同位素DOC模型：Li et al. 2020, Precambrian Research
