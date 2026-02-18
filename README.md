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
python examples/dice_doc_oxidation.py

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

## 项目结构

```
Balanced_isotope/
├── toolkit/           # 🔧 工具包层（数学、物理、同位素公式）
├── systems/           # 🧪 同位素体系层（各元素实现）
├── examples/          # 📚 应用示例层
├── run.py            # 统一入口
└── cli.py            # 命令行接口
```

详见 [ARCHITECTURE.md](ARCHITECTURE.md)

## 原代码迁移

| 原文件 | 新位置 | 说明 |
|--------|--------|------|
| `main.py` + `calculation/` | `systems/mg/` | Mg同位素体系 |
| `Carbon_Modeling/Oxi_Est.py` | `examples/dice_doc_oxidation.py` | DICE事件DOC氧化示例 |
| `utils/` | `toolkit/` | 通用工具 |

详见 [MIGRATION.md](MIGRATION.md)

## 参考

- Mg同位素模型参考：原 `calculation/` 和 `main.py`
- C同位素DOC模型：Li et al. 2020, Precambrian Research
