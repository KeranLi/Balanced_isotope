# 代码迁移说明

本文档说明原项目代码如何迁移到新架构。

## 迁移映射表

| 原文件/目录 | 新位置 | 状态 |
|------------|--------|------|
| `main.py` | `systems/mg/` | ✅ 已迁移到Mg同位素体系 |
| `calculation/` | `systems/mg/` + `systems/c/` | ✅ 已拆分迁移 |
| `Carbon_Modeling/Oxi_Est.py` | `examples/dice_doc_oxidation.py` | ✅ 已改写为示例 |
| `utils/` | `toolkit/` | ✅ 工具函数已整合 |

## 目录结构变化

### 旧结构
```
Balanced_isotope/
├── main.py                 # Mg分析入口
├── calculation/            # 计算模块（Mg专用）
│   ├── isotope.py
│   ├── mass_balance_model.py
│   └── monte_carlo.py
├── Carbon_Modeling/        # C分析（独立）
│   └── Oxi_Est.py
└── utils/                  # 工具函数
    ├── data_processor.py
    └── plotter.py
```

### 新结构
```
Balanced_isotope/
├── toolkit/                # 🔧 通用工具包（原utils + 通用化后的calculation）
│   ├── math/               #    ODE求解、插值、优化
│   ├── physics/            #    物理常数、分馏理论
│   └── isotope/            #    同位素公式（Delta转换、质量平衡等）
│
├── systems/                # 🧪 同位素体系层
│   ├── mg/                 #    Mg体系（整合原main.py + calculation）
│   ├── c/                  #    C体系（整合原Carbon_Modeling）
│   └── s/                  #    S体系（模板）
│
├── examples/               # 📚 应用示例
│   ├── basic_usage.py
│   └── dice_doc_oxidation.py  # DICE事件分析（替代原Oxi_Est.py）
│
├── run.py                  # 统一入口
└── cli.py                  # 命令行接口
```

## 使用方式对比

### Mg同位素分析

**旧写法（已废弃）：**
```python
# main.py
from calculation.mass_balance_model import MassBalanceModel
from calculation.isotope import MgIsotopeModel
from utils.data_processor import DataProcessor

# 手动加载和处理数据
data_processor = DataProcessor(file_path, ...)
age_data, depth_data, rm_data = data.data_processor.process_data()

# 创建多个模型实例
mass_balance_model = MassBalanceModel(RMg_interpolated, ...)
swpre = mass_balance_model.calculate_swpre()

model = MgIsotopeModel(...)
seawater_isotope = model.calculate_seawater_isotope(...)
```

**新写法（推荐）：**
```python
from systems.mg import MgIsotopeSystem

# 创建统一的体系实例
mg = MgIsotopeSystem()

# 风化比例计算
ratios = mg.calculate_weathering_ratio(
    delta_sample=-2.5,
    delta_seawater=-0.83
)

# 海水演化模拟
result = mg.seawater_evolution(
    time_span=(0, 100),
    flux_scenario='modern'
)
```

### C同位素DOC氧化分析

**旧写法（已废弃）：**
```python
# Carbon_Modeling/Oxi_Est.py
# 所有代码在一个文件中，参数硬编码

M_DIC = 4.0e18           # 硬编码参数
F_w = 25.0e18
delta13C_w = -4.0
# ...

def d_delta13C_org_dt(delta13C_org, t, F_odoc):
    # 直接实现微分方程
    ...

# 手动求解ODE
solution = odeint(d_delta13C_org_dt, delta_org_0, t)

# 手动绘图
plt.plot(...)
```

**新写法（推荐）：**
```python
from systems.c import CIsotopeSystem

# 使用配置化的C同位素体系
c_system = CIsotopeSystem(scenario='dice')

# 运行模型
result = c_system.doc_excursion_model(
    F_odoc_range=(0, 10e18),
    n_points=300
)

# 轻松切换不同情景
for scenario in ['dice', 'modern']:
    c = CIsotopeSystem(scenario=scenario)
    result = c.solve_steady_state(F_odoc=4e18)
    print(f"{scenario}: δ¹³C = {result.get('delta13C_carb'):.2f}‰")
```

## 主要改进

1. **命名更清晰**
   - `core/` → `toolkit/`（工具包，更直观）
   - `calculation/` → `systems/mg/`, `systems/c/`（按元素组织）
   - `Carbon_Modeling/` → `examples/`（明确是示例/应用）

2. **参数配置化**
   - 不再硬编码，通过 `parameters.py` 管理
   - 支持多情景切换（如 'dice', 'modern'）

3. **接口统一**
   - 所有同位素体系继承 `IsotopeSystem` 基类
   - 统一的方法：`mass_balance_equation()`, `fractionation_factor()`, `mixing_model()`

4. **代码复用**
   - 底层工具（ODE求解、插值等）与具体体系解耦
   - 所有体系共享 `toolkit/` 中的通用工具

5. **可扩展**
   - 添加新体系只需在 `systems/` 下创建新目录
   - 自动获得所有基础功能

## 运行方式对比

### 旧方式（已废弃）
```bash
python main.py data/Nie_Section_A.xlsx delta_26_Mg_iso
python Carbon_Modeling/Oxi_Est.py
```

### 新方式（推荐）
```bash
# 统一入口
python run.py all          # 运行所有演示
python run.py mg           # Mg同位素分析
python run.py c            # C同位素分析

# CLI方式
python cli.py list
python cli.py info mg
python cli.py mg --weathering-ratio
python cli.py c --F-odoc 4e18

# 直接运行示例
python examples/basic_usage.py
python examples/dice_doc_oxidation.py
```

## 删除的文件

以下文件已删除，功能已迁移到新位置：

- ❌ `main.py` → 功能在 `systems/mg/` 和 `run.py`
- ❌ `calculation/` → 功能拆分到 `systems/mg/` 和 `systems/c/`
- ❌ `Carbon_Modeling/` → 功能在 `examples/dice_doc_oxidation.py`
- ❌ `utils/` → 功能整合到 `toolkit/`
