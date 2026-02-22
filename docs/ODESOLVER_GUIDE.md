# ODE 求解器使用指南

本文档详细介绍同位素质量平衡框架中的常微分方程（ODE）求解工具。

## 目录

- [概述](#概述)
- [支持的求解方法](#支持的求解方法)
- [基本使用](#基本使用)
- [方法选择指南](#方法选择指南)
- [高级用法](#高级用法)
- [各同位素体系中的应用](#各同位素体系中的应用)
- [故障排除](#故障排除)

---

## 概述

`ODESolver` 是同位素质量平衡框架的核心数值计算工具，位于 `toolkit/math/numerical.py`。它封装了 SciPy 的 ODE 求解功能，为同位素演化模型提供统一、易用的接口。

### 主要特点

- **统一接口**：所有求解方法使用相同的 API
- **方法丰富**：支持 7 种不同的求解算法
- **大小写不敏感**：`'rk45'` 和 `'RK45'` 等效
- **自动错误处理**：求解失败返回友好的错误信息
- **结果封装**：返回结构化的 `ODEResult` 对象

---

## 支持的求解方法

| 方法 | 类型 | 阶数 | 适用问题 | 精度 |
|------|------|------|----------|------|
| `RK45` | 显式 Runge-Kutta | 5(4) | 非刚性 | 中等 |
| `RK23` | 显式 Runge-Kutta | 3(2) | 非刚性 | 低 |
| `DOP853` | 显式 Runge-Kutta | 8 | 非刚性 | **高** |
| `Radau` | 隐式 Runge-Kutta | 5 | 刚性 | 中等 |
| `BDF` | 隐式多步法 | 1-5 | 刚性 | 中等 |
| `LSODA` | 自动切换 | 可变 | 刚性/非刚性 | 中等 |
| `odeint` | 旧版 API | 可变 | 通用 | 中等 |

### 方法分类

**非刚性问题**（推荐 `RK45` 或 `DOP853`）：
- Mg 同位素海水演化
- 大多数同位素混合模型
- 短时间尺度演化

**刚性问题**（推荐 `BDF` 或 `Radau`）：
- 多时间尺度耦合系统
- 化学反应动力学
- 长时间尺度稳态逼近

**未知特性**（推荐 `LSODA`）：
- 自动检测刚性
- 通用求解
- 探索性计算

---

## 基本使用

### 导入

```python
from toolkit.math.numerical import ODESolver
```

### 最小示例

```python
import numpy as np
from toolkit.math.numerical import ODESolver

# 定义ODE: dy/dt = -k * y (指数衰减)
def decay_model(t, y, k):
    return -k * y

# 求解
result = ODESolver.solve(
    func=decay_model,
    y0=1.0,                    # 初始条件 y(0) = 1
    t_span=(0, 10),            # 时间范围 0 到 10
    args=(0.5,),               # 传递参数 k=0.5
    method='RK45',             # 求解方法
    n_points=100               # 输出100个点
)

if result.success:
    print(f"时间: {result.t}")
    print(f"解: {result.y}")
else:
    print(f"求解失败: {result.message}")
```

### 返回结果结构

```python
result.t      # 时间点数组 (n_points,)
result.y      # 解数组 (n_points,) 或 (n_points, n_variables)
result.success # 是否成功 (bool)
result.message # 状态信息 (str)
```

---

## 方法选择指南

### 决策流程

```
问题类型？
├── 已知非刚性
│   ├── 需要高精度？ → DOP853
│   └── 一般精度？   → RK45 (默认)
├── 已知刚性
│   ├── 需要稳定性？ → BDF
│   └── 需要效率？   → Radau
└── 不确定？
    ├── 快速尝试？   → LSODA (自动)
    └── 保守选择？   → RK45 (失败后换 BDF)
```

### 各方法详细说明

#### 1. RK45（推荐默认）

```python
result = ODESolver.solve(func, y0, t_span, method='RK45')
```

- **优点**：精度与效率平衡好，适合大多数问题
- **缺点**：对刚性问题可能极慢或失败
- **适用**：Mg同位素海水演化、U同位素稳态逼近

#### 2. DOP853（高精度）

```python
result = ODESolver.solve(func, y0, t_span, method='DOP853', n_points=2000)
```

- **优点**：8阶精度，误差控制严格
- **缺点**：计算量较大
- **适用**：需要精确捕捉小幅度同位素漂移

#### 3. BDF（刚性问题）

```python
result = ODESolver.solve(func, y0, t_span, method='BDF', max_step=1e3)
```

- **优点**：对刚性问题稳定
- **缺点**：非刚性问题上效率较低
- **适用**：多储库耦合系统、化学平衡快速反应

#### 4. LSODA（自动检测）

```python
result = ODESolver.solve(func, y0, t_span, method='LSODA')
```

- **优点**：自动切换 Adams（非刚性）和 BDF（刚性）
- **缺点**：有一定检测开销
- **适用**：问题特性未知、探索性计算

#### 5. odeint（兼容性）

```python
result = ODESolver.solve(func, y0, t_span, method='odeint')
```

- **优点**：某些问题更稳定，向后兼容
- **缺点**：功能较少，不推荐新代码使用
- **适用**：从旧代码迁移

---

## 高级用法

### 多变量系统

```python
import numpy as np
from toolkit.math.numerical import ODESolver

# Lotka-Volterra 捕食者-猎物模型
def lotka_volterra(t, y, alpha, beta, gamma, delta):
    prey, predator = y
    dprey_dt = alpha * prey - beta * prey * predator
    dpredator_dt = delta * prey * predator - gamma * predator
    return np.array([dprey_dt, dpredator_dt])

result = ODESolver.solve(
    func=lotka_volterra,
    y0=np.array([10.0, 5.0]),  # 初始种群
    t_span=(0, 50),
    args=(1.0, 0.5, 0.75, 0.25),  # alpha, beta, gamma, delta
    method='RK45',
    n_points=500
)

if result.success:
    prey = result.y[:, 0]
    predator = result.y[:, 1]
```

### 刚性问题处理

```python
# Van der Pol 振荡器（刚性参数）
def van_der_pol(t, y, mu):
    x, v = y
    dxdt = v
    dvdt = mu * (1 - x**2) * v - x
    return np.array([dxdt, dvdt])

# 当 mu 很大时，系统是刚性的
result = ODESolver.solve(
    func=van_der_pol,
    y0=np.array([2.0, 0.0]),
    t_span=(0, 100),
    args=(1000.0,),  # 大 mu → 刚性
    method='BDF'     # 必须用刚性求解器
)
```

### 事件检测（Events）

```python
from scipy.integrate import solve_ivp

# 定义事件函数（当 y[0] = 0.5 时触发）
def event_func(t, y):
    return y[0] - 0.5

event_func.terminal = False  # 不终止积分
event_func.direction = -1    # 只检测下降沿

# 通过 kwargs 传递额外参数给 solve_ivp
result = ODESolver.solve(
    func=my_ode,
    y0=y0,
    t_span=(0, 100),
    method='RK45',
    events=event_func,  # 传递给 solve_ivp
    max_step=1.0        # 限制最大步长
)
```

### 精度控制

```python
# 通过 rtol 和 atol 控制精度
result = ODESolver.solve(
    func=precise_model,
    y0=y0,
    t_span=(0, 1e6),
    method='DOP853',
    rtol=1e-10,    # 相对误差
    atol=1e-12,    # 绝对误差
    n_points=1000
)
```

---

## 各同位素体系中的应用

### Mg 同位素 - 海水演化

```python
from systems.mg import MgIsotopeSystem

system = MgIsotopeSystem()

# 模拟风化转变
result = system.simulate_weathering_transition(
    time_span=(0, 5),  # 5 Myr
    transition={
        't_start': 0, 't_end': 2,
        'f_initial': 0.2, 'f_final': 0.8,
        'mode': 'linear'
    },
    n_points=500
)

# 内部使用 ODESolver.solve(method='RK45')
```

### U 同位素 - 缺氧事件

```python
from systems.u import UIsotopeSystem

system = UIsotopeSystem(scenario='modern')

# 非稳态模拟
result = system.simulate_anoxic_event(
    event_duration=1.0,  # 1 Myr
    peak_f_anox=0.8,
    background_f_anox=0.2
)

# 内部使用 ODESolver.solve(method='RK45')
```

### 自定义模型

```python
from toolkit.math.numerical import ODESolver

def custom_isotope_model(t, state, params):
    """
    自定义同位素演化模型
    state = [M, delta]  # 储库质量和同位素值
    """
    M, delta = state
    F_in, F_out, delta_in, epsilon = params
    
    dM_dt = F_in - F_out
    ddelta_dt = (F_in * (delta_in - delta) - F_out * epsilon) / M
    
    return np.array([dM_dt, ddelta_dt])

result = ODESolver.solve(
    func=custom_isotope_model,
    y0=np.array([1e20, 0.0]),  # 初始储库和同位素
    t_span=(0, 100e6),         # 100 Myr
    args=((1e10, 1e10, -1.0, -2.7),),
    method='RK45',
    n_points=1000
)
```

---

## 故障排除

### 问题：求解失败 (`success=False`)

**可能原因与解决方案：**

1. **刚性问题使用非刚性求解器**
   ```python
   # 错误：
   result = ODESolver.solve(stiff_func, y0, t_span, method='RK45')
   
   # 正确：
   result = ODESolver.solve(stiff_func, y0, t_span, method='BDF')
   ```

2. **时间跨度过大**
   ```python
   # 增加点数或限制步长
   result = ODESolver.solve(func, y0, t_span, method='RK45', 
                            n_points=5000, max_step=1e5)
   ```

3. **初始条件不合理**
   - 检查 `y0` 是否在物理合理范围内
   - 避免除零或负数开方

4. **数值溢出**
   ```python
   # 使用更高精度或缩小时间范围
   result = ODESolver.solve(func, y0, (0, 1e3), method='DOP853',
                            rtol=1e-8, atol=1e-10)
   ```

### 问题：结果不准确

**检查清单：**

- [ ] 增加 `n_points`（输出点数）
- [ ] 降低 `rtol` 和 `atol`（误差容限）
- [ ] 尝试 `DOP853` 方法（更高阶）
- [ ] 检查 ODE 函数实现是否正确

### 问题：计算太慢

**优化建议：**

1. **降低精度要求**（如果可以接受）
   ```python
   result = ODESolver.solve(func, y0, t_span, method='RK45',
                            rtol=1e-4, atol=1e-6)  # 默认 1e-3
   ```

2. **减少输出点数**
   ```python
   result = ODESolver.solve(func, y0, t_span, n_points=100)  # 默认 1000
   ```

3. **使用更简单的方法**
   ```python
   result = ODESolver.solve(func, y0, t_span, method='RK23')  # 比 RK45 快
   ```

---

## 性能比较

基于典型同位素演化模型的基准测试：

| 方法 | 相对速度 | 精度 | 刚性处理 | 推荐指数 |
|------|----------|------|----------|----------|
| RK23 | ★★★★★ | ★★☆☆☆ | ✗ | ★★☆☆☆ |
| RK45 | ★★★★☆ | ★★★☆☆ | ✗ | ★★★★★ |
| DOP853 | ★★☆☆☆ | ★★★★★ | ✗ | ★★★★☆ |
| BDF | ★★☆☆☆ | ★★★☆☆ | ✓ | ★★★★☆ |
| Radau | ★★★☆☆ | ★★★★☆ | ✓ | ★★★☆☆ |
| LSODA | ★★★☆☆ | ★★★☆☆ | 自动 | ★★★☆☆ |
| odeint | ★★★★☆ | ★★★☆☆ | △ | ★★☆☆☆ |

**建议：**
- 默认使用 `RK45`
- 求解失败时尝试 `BDF`
- 高精度需求使用 `DOP853`
- 不确定时用 `LSODA` 自动检测

---

## 参考

- SciPy 文档：[scipy.integrate.solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
- 数值方法教材：Hairer, E., Nørsett, S. P., & Wanner, G. (1993). Solving Ordinary Differential Equations I: Nonstiff Problems.
- 刚性问题：Hairer, E., & Wanner, G. (1996). Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems.

---

**文档版本：** 1.0  
**最后更新：** 2026-02-22  
**对应代码：** `toolkit/math/numerical.py`
