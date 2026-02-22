# Jupyter Notebook 示例

本目录包含同位素质量平衡模型的 Jupyter Notebook 示例，提供交互式的学习和探索环境。

## 安装 Jupyter

### 方式1: pip 安装

```bash
# 安装 Jupyter
pip install jupyter notebook

# 或者安装 JupyterLab (推荐)
pip install jupyterlab
```

### 方式2: conda 安装

```bash
conda install -c conda-forge jupyterlab
```

## 启动 Notebook

```bash
# 进入 notebook 目录
cd examples/notebooks

# 启动 Jupyter Notebook
jupyter notebook

# 或者启动 JupyterLab
jupyter lab
```

浏览器将自动打开，显示可用的 Notebook 文件。

## 可用 Notebooks

| 文件 | 内容 | 难度 |
|------|------|------|
| `01_basic_usage.ipynb` | 核心工具、同位素公式、各体系基础用法 | ⭐ 入门 |
| `02_nitrogen_isotope_reproduction.ipynb` | **N同位素论文复现**：Kang et al. (2023) & Ma et al. (2025) | ⭐⭐ 进阶 |

## 使用建议

### 在 VS Code 中使用

VS Code 原生支持 Jupyter Notebook：
1. 安装 Python 扩展
2. 打开 `.ipynb` 文件即可直接运行

### 在 Google Colab 中使用

可以将 Notebook 上传到 Google Colab 在线运行：
1. 上传 `.ipynb` 文件到 Google Drive
2. 使用 Colab 打开
3. 安装依赖: `!pip install numpy scipy pandas`

### 导出为其他格式

```bash
# 导出为 Python 脚本
jupyter nbconvert --to script 01_basic_usage.ipynb

# 导出为 HTML
jupyter nbconvert --to html 01_basic_usage.ipynb

# 导出为 PDF
jupyter nbconvert --to pdf 01_basic_usage.ipynb
```

## Notebook 结构

每个 Notebook 通常包含：

1. **环境准备** - 导入必要的库
2. **核心工具演示** - ODE求解、插值等
3. **同位素公式** - Delta计算、混合、分馏
4. **体系应用** - Mg、C、N、U等具体示例
5. **批量计算** - 数据处理和导出

## 交互式探索

Notebook 的优势在于可以：

- **修改参数**: 直接修改代码单元格中的参数，重新运行查看结果
- **可视化**: 添加 `matplotlib` 图表（如果需要）
- **实验**: 快速测试不同输入和情景
- **文档**: 结合 Markdown 文本和代码，形成完整的工作流

## 注意事项

1. **路径问题**: Notebook 中使用相对路径 `Path.cwd().parent.parent` 指向项目根目录
2. **依赖**: 确保已安装项目依赖 `pip install numpy scipy pandas`
3. **保存**: 定期保存 Notebook (`Ctrl+S`)，避免丢失工作

## 创建新 Notebook

参考 `01_basic_usage.ipynb` 的结构，可以创建针对特定研究的 Notebook：

```python
# 标准开头
import sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd().parent.parent))

# 导入需要的模块
from systems.x import XIsotopeSystem
import numpy as np
import pandas as pd

# 你的分析代码...
```

## 贡献

如果你有好的 Notebook 示例想要分享，欢迎提交 PR！
