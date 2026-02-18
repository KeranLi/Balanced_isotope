# 同位素质量平衡模型 - 架构设计文档

## 概述

本项目采用**三层架构**设计，将通用的数学/物理工具、同位素体系实现和应用场景分离，实现代码的高度可扩展性和可维护性。

```
┌─────────────────────────────────────────────────────────────────┐
│                        应用场景层 (Applications)                   │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────────────────┐  │
│  │  Weathering │  │    Redox    │  │     Carbon Cycle        │  │
│  │   风化分析   │  │  氧化还原   │  │       碳循环             │  │
│  └─────────────┘  └─────────────┘  └─────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
                              ▲
                              │ 调用
┌─────────────────────────────────────────────────────────────────┐
│                        同位素体系层 (Systems)                      │
│  ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌─────────┐             │
│  │   Mg    │  │    C    │  │    S    │  │   ...   │             │
│  │ 镁同位素 │  │ 碳同位素 │  │ 硫同位素 │  │ 其他体系 │             │
│  └─────────┘  └─────────┘  └─────────┘  └─────────┘             │
│                              ▲                                  │
│                              │ 继承                             │
│  ┌─────────────────────────────────────────────────────────┐   │
│  │              基类 (IsotopeSystem)                        │   │
│  │     mass_balance_equation() | fractionation_factor()    │   │
│  │     mixing_model() | time_evolution() | solve_steady_state() │ │
│  └─────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────┘
                              ▲
                              │ 调用
┌─────────────────────────────────────────────────────────────────┐
│                         核心工具层 (Core)                         │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────────────────┐  │
│  │    Math     │  │   Physics   │  │        Isotope          │  │
│  │   数学工具   │  │  物理常数   │  │      同位素公式          │  │
│  │  ODE Solver │  │  Fractionation│  │  Delta Calculator      │  │
│  │ Interpolator│  │   Constants  │  │  Mass Balance          │  │
│  │  Optimizer  │  │  Reservoir   │  │  Rayleigh Model        │  │
│  └─────────────┘  └─────────────┘  └─────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
```

---

## 第一层：核心工具层 (Core)

### 1.1 toolkit/math - 数学工具

提供与具体同位素无关的通用数学功能。

```python
from toolkit.math.numerical import ODESolver, Interpolator, Optimizer

# ODE求解
result = ODESolver.solve(
    func=dy_dt,      # 微分方程
    y0=initial,      # 初始条件
    t_span=(0, 100), # 时间范围
    args=(params,)   # 额外参数
)

# 插值
y_new = Interpolator.interpolate(x, y, x_new, method='cubic')
```

**包含模块：**
- `numerical.py`: ODE求解、插值、优化、统计工具
- 支持多种ODE求解方法（RK45, RK23, BDF, LSODA, odeint）
- 多种插值方法（线性、三次样条）

### 1.2 toolkit/physics - 物理化学常数

```python
from toolkit.physics.constants import (
    PhysicalConstants,      # 物理常数（阿伏伽德罗常数等）
    FractionationTheory,    # 分馏理论公式
    ReservoirConstants,     # 地球化学储库参数
    get_element_info        # 元素信息查询
)

# 计算温度依赖的平衡分馏
epsilon = FractionationTheory.equilibrium_fractionation(T=298, A=0, B=-0.8)
```

### 1.3 toolkit/isotope - 同位素数学公式

提供同位素地球化学中通用的计算公式。

```python
from toolkit.isotope.formulas import (
    DeltaCalculator,         # δ值转换
    MassBalance,            # 质量平衡
    RayleighFractionation,  # 瑞利分馏
    EvolutionEquations      # 演化方程
)

# Delta转换
ratio = DeltaCalculator.delta_to_ratio(delta=-5, standard_ratio=0.0112372)

# 二元混合
delta_mix = MassBalance.two_component_mixing(-4, -30, f_a=0.7)
```

---

## 第二层：同位素体系层 (Systems)

### 2.1 基类设计 (systems/base)

所有同位素体系继承 `IsotopeSystem` 基类：

```python
class IsotopeSystem(ABC):
    # 必须实现的抽象方法
    @abstractmethod
    def mass_balance_equation(self, state, fluxes, time): ...
    
    @abstractmethod
    def fractionation_factor(self, process, temperature, **kwargs): ...
    
    @abstractmethod
    def mixing_model(self, end_members, proportions): ...
    
    # 已实现的高级方法
    def time_evolution(self, initial_state, time_span, fluxes): ...
    def solve_steady_state(self, fluxes): ...
    def inverse_model(self, observations, observation_times): ...
```

### 2.2 具体体系实现

每个同位素体系是一个独立的Python包：

```
systems/
├── mg/           # Mg同位素体系
│   ├── __init__.py
│   ├── model.py      # MgIsotopeSystem类
│   └── parameters.py # Mg同位素参数
│
├── c/            # C同位素体系
│   ├── __init__.py
│   ├── model.py      # CIsotopeSystem类
│   └── parameters.py # C同位素参数
│
├── s/            # S同位素体系（模板）
├── sr/           # Sr同位素体系（计划中）
└── nd/           # Nd同位素体系（计划中）
```

#### 示例：Mg同位素体系使用

```python
from systems.mg import MgIsotopeSystem

# 创建体系实例
mg = MgIsotopeSystem()

# 风化比例计算
ratios = mg.calculate_weathering_ratio(
    delta_sample=-2.5, 
    delta_seawater=-0.83
)
print(f"碳酸盐风化: {ratios['f_carbonate']:.1%}")

# 海水演化模拟
result = mg.seawater_evolution(
    time_span=(0, 100),     # 100 Ma
    flux_scenario='modern'
)
```

#### 示例：C同位素体系使用

```python
from systems.c import CIsotopeSystem

# 创建DICE情景的C体系
c = CIsotopeSystem(scenario='dice')

# 计算特定DOC通量下的稳态同位素
result = c.solve_steady_state(F_odoc=4e18)
print(f"δ¹³C_carb: {result.get('delta13C_carb'):.2f}‰")

# 氧化剂消耗计算
oxidant = c.calculate_oxidant_consumption(
    F_odoc=4e18,
    scenario_name='modern_o2_high_sulfate'
)
```

---

## 第三层：应用场景层 (Applications)

具体的地质应用，可以组合多个同位素体系：

```
applications/
├── weathering/     # 风化通量估算
├── redox/          # 氧化还原条件
└── carbon_cycle/   # 碳循环分析
```

**示例：风化应用**

```python
# 可同时使用Mg和Sr同位素来约束风化
from systems.mg import MgIsotopeSystem
from systems.sr import SrIsotopeSystem

def comprehensive_weathering_analysis(mg_data, sr_data):
    mg_system = MgIsotopeSystem()
    sr_system = SrIsotopeSystem()
    
    # Mg约束碳酸盐/硅酸盐风化比例
    mg_ratio = mg_system.calculate_weathering_ratio(...)
    
    # Sr约束总体风化强度
    sr_flux = sr_system.calculate_weathering_flux(...)
    
    return combine_results(mg_ratio, sr_flux)
```

---

## 如何添加新的同位素体系

### 步骤1：创建参数文件

```python
# systems/x/parameters.py
from systems.base.isotope_system import IsotopeParameters

def get_x_parameters() -> IsotopeParameters:
    return IsotopeParameters(
        element='x',
        name='Element X',
        reference_standard='XXX',
        reference_ratios={'x/y': 0.123},
        fractionation_factors={...},
        end_members={...},
        reservoir_mass=1e20,
        input_fluxes={...},
        output_fluxes={...}
    )
```

### 步骤2：实现模型类

```python
# systems/x/model.py
from systems.base.isotope_system import IsotopeSystem, ModelResult
from systems.x.parameters import get_x_parameters

class XIsotopeSystem(IsotopeSystem):
    ELEMENT = 'x'
    NAME = 'Element X'
    ISOTOPES = ['x1', 'x2', 'x3']
    
    def _default_parameters(self):
        return get_x_parameters()
    
    def mass_balance_equation(self, state, fluxes, time):
        # 实现质量平衡方程
        pass
    
    def fractionation_factor(self, process, temperature, **kwargs):
        # 返回分馏系数
        pass
    
    def mixing_model(self, end_members, proportions):
        # 实现混合模型
        pass
    
    def state_dimension(self):
        return 1  # 或更多，取决于状态变量数量
```

### 步骤3：创建__init__.py

```python
# systems/x/__init__.py
from systems.x.model import XIsotopeSystem
from systems.x.parameters import get_x_parameters

__all__ = ['XIsotopeSystem', 'get_x_parameters']
```

### 步骤4：添加到CLI（可选）

在 `cli.py` 中添加新的子命令。

---

## 设计原则

1. **单一职责**：每层只负责特定功能
   - Core：通用数学/物理工具
   - Systems：特定同位素体系的物理化学模型
   - Applications：具体地质问题的解决方案

2. **依赖方向**：只能上层调用下层
   - Systems → Core
   - Applications → Systems → Core

3. **可扩展性**：添加新体系只需在systems/下创建新目录

4. **可配置性**：通过parameters配置，避免硬编码

5. **与I/O解耦**：模型只返回数据，不直接绘图或读写文件

---

## CLI使用

```bash
# 列出支持的体系
python cli.py list

# 查看体系信息
python cli.py info mg

# Mg同位素分析
python cli.py mg --weathering-ratio

# C同位素DOC模型
python cli.py c --F-odoc 4e18
python cli.py c --target-excursion -4
```

---

## 文件结构

```
Balanced_isotope/
├── toolkit/                      # 核心工具层
│   ├── math/
│   │   └── numerical.py       # 数学工具
│   ├── physics/
│   │   └── constants.py       # 物理常数
│   ├── isotope/
│   │   └── formulas.py        # 同位素公式
│   └── io/                    # IO工具（待实现）
│
├── systems/                   # 同位素体系层
│   ├── base/
│   │   └── isotope_system.py  # 抽象基类
│   ├── mg/                    # Mg体系（已实现）
│   ├── c/                     # C体系（已实现）
│   ├── s/                     # S体系（模板）
│   ├── sr/                    # Sr体系（计划中）
│   └── nd/                    # Nd体系（计划中）
│
├── applications/              # 应用场景层
│   ├── weathering/
│   ├── redox/
│   └── carbon_cycle/
│
├── examples/                  # 使用示例
│   └── basic_usage.py
│
├── cli.py                     # 命令行接口
├── data/                      # 数据文件
└── reference/                 # 参考文献
```
