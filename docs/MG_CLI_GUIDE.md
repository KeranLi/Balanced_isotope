# Mg 同位素 CLI 使用指南

项目提供两个不可混用的 Mg 同位素模型：

- `carbonate`：海相碳酸盐岩和海水 Mg 箱模型，参考 Kasemann et al. (2014)。
- `siliciclastic`：陆源碎屑/黏土残余 Mg 同位素与风化通量模型，参考 Hu et al. (2023)。

Nie Section A 和 B 是碎屑组分，必须使用 `siliciclastic`。

## Nie 剖面运行命令

```bash
python cli.py mg \
  --component-type siliciclastic \
  --file data/Nie_Section_A.xlsx \
  --delta-protolith 0.40 \
  --mg-monte-carlo 2000 \
  --plot-output results/Nie_A_p040_Mg_release.png \
  --flux-plot-output results/Nie_A_p040_Changjiang_flux.png \
  --output results/Nie_A_p040_siliciclastic_flux.csv

python cli.py mg \
  --component-type siliciclastic \
  --file data/Nie_Section_B.xlsx \
  --mg-monte-carlo 2000 \
  --plot-output results/Nie_B_Mg_release.png \
  --flux-plot-output results/Nie_B_Changjiang_flux.png \
  --output results/Nie_B_siliciclastic_flux.csv
```

输入行必须按深度递增排列。没有实测深度时，程序生成：

- `sample_index`：从 1 开始的样品序号。
- `relative_depth`：从 0 到 1 的归一化剖面位置。

归一化位置只能用于剖面趋势，不能解释为米或年龄。

## 不依赖人工基线的通量代理

在岩石供应量和原岩 Mg 含量不变的假设下，单位岩石释放的 Mg 与
`1-f_Mg` 成正比。因此程序直接将其作为风化通量代理：

```text
weathering_flux_proxy(i) = Mg_release_fraction(i) = 1 - f_Mg(i)
```

剖面级 Monte Carlo 在同一次抽样中为所有样品共享分馏参数，并独立传播
每个样品的分析误差。曲线的数值范围为物理上的 `0-1`，不需要指定
哪几个样品等于 `1x`。

程序使用两段常数状态的 L1 损失自动检测一个变点。报告的变点由整条
后验中值剖面确定；逐次 Monte Carlo 变点用于估计位置不确定度。在自动
变点固定后，程序分别计算浅部段和深部段的 `1-f_Mg` 中值，并用二者之比
表示通量改变倍数：

```text
deep_to_shallow_flux_ratio =
    median(1 - f_Mg_deep) / median(1 - f_Mg_shallow)
```

这个倍数不依赖手工基线，但仍依赖“两状态”近似以及岩石供应量和原岩
Mg 含量不变的假设。自动算法总能为有足够数据的剖面选择一个最佳分割，
因此还应结合 L1 代价改善程度和倍数区间判断转变是否明显。

主要输出为：

- `Mg_release_fraction_mc_median`
- `Mg_release_fraction_mc_ci95_low/high`
- `Mg_release_fraction_mc_valid_fraction`
- `automatic_change_point_after_sample`
- `automatic_change_point_ci95_low/high`
- `shallow_Mg_release_fraction_mc_median`
- `deep_Mg_release_fraction_mc_median`
- `deep_to_shallow_weathering_flux_ratio_median`
- `deep_to_shallow_weathering_flux_ratio_ci95_low/high`
- `change_point_l1_cost_improvement`

## 竖版通量图

`--plot` 生成约 `4.1 x 8.2 inch` 的窄幅竖版图：相对地层位置位于纵轴，
`1-f_Mg` 位于横轴，并显示 95% 区间和自动变点。有效 Monte Carlo 比例
低于 50% 的样品用空心方形表示；完全没有物理解时，方形放在横轴左边界，
不赋予虚假的 `1-f_Mg=0` 数值。

`--flux-plot-output` 另外生成现今长江先验标定的条件硅酸盐溶解 Mg 通量
纵向演化图。该图使用 `F_river,Mg = 3.0e11 mol/yr`，横轴单位为
`1e11 mol/yr`。它不是总岩石风化通量；绝对尺度随河流 Mg 总通量先验
线性变化。

```bash
python cli.py mg \
  --component-type siliciclastic \
  --file data/Nie_Section_B.xlsx \
  --plot-output results/Nie_B_Mg_release.png \
  --output results/Nie_B_siliciclastic_flux.csv
```

没有年龄锚点时，纵轴为 0 到 1 的相对地层位置。若已知剖面两端年龄，
可线性映射为时间：

```bash
python cli.py mg \
  --component-type siliciclastic \
  --file data/Nie_Section_B.xlsx \
  --time-start 0.0 \
  --time-end 2.0 \
  --time-unit Myr \
  --plot \
  --output results/Nie_B_timed_flux.csv
```

两端年龄只定义线性时间轴；如果沉积速率变化明显，应先提供逐样品年龄。

## 碎屑模型

### Rayleigh 残余模型

以相对同位素比值 `R = 1 + delta/1000` 表示：

```text
R_residue = R_protolith * f_Mg^(alpha - 1)
alpha = 1 + Delta_fluid-protolith/1000
```

其中 `f_Mg` 是风化残余中保留的 Mg 比例。确定性参考值与 Monte Carlo
先验分别为：

```text
delta26Mg_protolith reference = -0.25 per mil
delta26Mg_protolith MC prior = Uniform(-0.45, -0.25) per mil
Delta_fluid-protolith = -0.25 per mil
Delta_fluid-protolith MC prior = Uniform(-0.35, -0.15) per mil
```

每次剖面 Monte Carlo 实现只抽取一个原岩值，并由全部样品共享。这样传播
了原岩端元的不确定度，同时不会把随机的逐样品端元变化误判为地层趋势。
`-0.45` 到 `-0.25‰` 是当前工作先验，不是由 Nie 数据独立反演得到的
后验约束，最终解释应报告该先验范围。

Hu et al. (2023) 的校准范围是 `f_Mg >= 0.5`，即 Mg 损失不超过 50%。
更强风化条件仍可计算，但输出标记为 `extrapolated_advanced_weathering`。

累计释放到水体的硅酸盐 Mg 端元为：

```text
R_silicate_flux = R_protolith * (1 - f_Mg^alpha) / (1 - f_Mg)
```

这一步替代了旧代码中没有质量守恒依据的
`delta_clay + fixed_offset` 方法。

### 河流 Mg 通量质量平衡

```text
F_river = F_silicate + F_carbonate

F_river * delta_river =
    F_silicate * delta_silicate_flux
  + F_carbonate * delta_carbonate
```

因此：

```text
f_silicate_flux =
    (delta_river - delta_carbonate)
  / (delta_silicate_flux - delta_carbonate)

F_silicate = F_river * f_silicate_flux
```

默认采用 Hu et al. (2023) 长江先验：

```text
delta26Mg_river = -1.14 per mil
delta26Mg_carbonate = -2.0 per mil
F_river_total = 3.0e11 mol/yr
```

绝对 `F_silicate` 与 `F_river_total` 成正比。没有研究区总河流 Mg
通量约束时，结果是基于长江先验的条件通量，不是唯一解。

## 不确定度

剖面 Monte Carlo 对分析误差逐样品抽样，对原岩和分馏参数逐次实现共享。
默认包括：

- 样品 `delta_26_Mg_iso_2sd / 2`；
- 原岩 `delta26Mg` 在 `-0.45` 到 `-0.25 per mil` 均匀分布，整条剖面共享；
- `Delta_fluid-protolith` 在 `-0.35` 到 `-0.15` 均匀分布；
- 河水 `delta26Mg = -1.14 +/- 0.15 per mil`；
- 碳酸盐端元 `delta26Mg = -2.0 +/- 0.4 per mil`，对应文献 `+/-0.8 per mil (2SD)`。

主要输出包括中值和 95% 区间：

- `f_Mg_mc_median`
- `f_Mg_mc_ci95_low`、`f_Mg_mc_ci95_high`
- `F_silicate_mc_median_mol_yr`
- `F_silicate_mc_ci95_low_mol_yr`
- `F_silicate_mc_ci95_high_mol_yr`
- `mc_valid_flux_fraction`

`mc_valid_flux_fraction` 较低表示多数随机参数组合没有物理解，此时通量约束较弱。

## 数据质量控制

当输入同时包含 delta25Mg 和 delta26Mg 时，程序计算：

```text
d25_d26_residual = delta25Mg - 0.521 * delta26Mg
```

如果两列都有 2SD，程序同时输出残差 z 分数和 `mass_dependent_qc`。
质控失败的样品不会被静默删除。

## 碎屑模型参数

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `--delta-protolith` | 无 | 固定原岩 delta26Mg，并关闭范围抽样 |
| `--delta-protolith-range LOW HIGH` | -0.45 -0.25 | 原岩 delta26Mg 均匀先验 |
| `--delta-fluid-protolith` | -0.25 | 流体-原岩分馏 |
| `--delta-river-mg` | -1.14 | 河流溶解 Mg 端元 |
| `--delta-carbonate` | -2.0 | 碳酸盐风化端元 |
| `--river-mg-flux` | 3.0e11 | 河流总溶解 Mg 通量，mol/yr |
| `--river-mg-flux-rel-std` | 0 | 河流总通量相对 1 sigma 不确定度 |
| `--mg-monte-carlo` | 2000 | 每个样品的随机抽样数 |
| `--random-seed` | 42 | 随机种子 |
| `--plot` | false | 生成竖版 Mg 释放比例图 |
| `--plot-output` | 自动 | 自定义图件路径 |
| `--flux-plot-output` | 自动 | 长江条件硅酸盐 Mg 通量图路径 |
| `--time-start` | 无 | 第一行样品的时间 |
| `--time-end` | 无 | 最后一行样品的时间 |

自定义通量先验示例：

```bash
python cli.py mg \
  --component-type siliciclastic \
  --file data/Nie_Section_A.xlsx \
  --river-mg-flux 4.0e11 \
  --river-mg-flux-rel-std 0.2 \
  --delta-protolith -0.30 \
  --output results/Nie_A_custom_flux.xlsx
```

固定原岩端元的敏感性对照：

```bash
python cli.py mg \
  --component-type siliciclastic \
  --file data/Nie_Section_B.xlsx \
  --delta-protolith -0.25 \
  --plot \
  --output results/Nie_B_fixed_protolith.csv
```

## 模型状态

| 状态 | 含义 |
|---|---|
| `within_hu_calibration` | 位于 Hu 初始至中等风化校准范围 |
| `extrapolated_advanced_weathering` | 数学上有解，但属于高级风化外推 |
| `incompatible_with_assumed_protolith` | 样品比假定原岩更轻，当前 Rayleigh 假设无物理解 |
| `incompatible_flux_endmembers` | 河水值不在两个通量端元之间 |

## 参考文献

Hu, Z. et al. (2023). Mg isotopes of siliciclastic sediments on continental
marginal sea: Insights for the potential to trace silicate weathering. Global
and Planetary Change, 231, 104307.

Kasemann, S.A. et al. (2014). Continental weathering following a Cryogenian
glaciation: Evidence from calcium and magnesium isotopes. Earth and Planetary
Science Letters, 396, 66-77.
