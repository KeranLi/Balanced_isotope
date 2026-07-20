"""
Dukou剖面Sr同位素数据读取与处理模块

处理流程:
1. 读取Sr同位素数据和Mn/Sr数据
2. 根据Mn/Sr < 1筛选数据
3. 深度-年龄转换
4. LOWESS平滑生成目标曲线

数据文件:
- data/Sr_Table_sun.xlsx: Sr同位素数据
- data/Sr_MnSr_sun.xlsx: Mn/Sr筛选数据
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Tuple, Dict, List
from dataclasses import dataclass, field
from scipy.interpolate import interp1d
from scipy import stats

# 导入statsmodels用于LOWESS平滑
try:
    from statsmodels.nonparametric.smoothers_lowess import lowess
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False


@dataclass
class DukouData:
    """Dukou剖面数据类"""
    ages: np.ndarray          # 年龄 (Ma)
    positions: np.ndarray     # 深度 (m)
    sr_ratios: np.ndarray     # ⁸⁷Sr/⁸⁶Sr
    mn_sr: np.ndarray         # Mn/Sr比值
    sample_ids: List[str]     # 样品ID
    uncertainties: np.ndarray # 2σ不确定度
    formations: Optional[np.ndarray] = None  # 地层
    stages: Optional[np.ndarray] = None      # Stage/Decline标记
    diagenesis: Optional[np.ndarray] = None  # 成岩作用标记

    # 目标曲线 (LOWESS平滑后)
    target_ages: Optional[np.ndarray] = None
    target_sr: Optional[np.ndarray] = None
    target_ci_lower: Optional[np.ndarray] = None  # 95% CI lower
    target_ci_upper: Optional[np.ndarray] = None  # 95% CI upper

    def get_interval_data(self, age_min: float, age_max: float) -> 'DukouData':
        """获取特定年龄区间数据"""
        mask = (self.ages >= age_min) & (self.ages <= age_max)
        return DukouData(
            ages=self.ages[mask],
            positions=self.positions[mask],
            sr_ratios=self.sr_ratios[mask],
            mn_sr=self.mn_sr[mask],
            sample_ids=[sid for sid, m in zip(self.sample_ids, mask) if m],
            uncertainties=self.uncertainties[mask],
            formations=self.formations[mask] if self.formations is not None else None,
            stages=self.stages[mask] if self.stages is not None else None,
            diagenesis=self.diagenesis[mask] if self.diagenesis is not None else None
        )

    def get_stage_data(self, stage: str) -> 'DukouData':
        """
        获取特定Stage的数据

        Parameters
        ----------
        stage : str
            Stage名称，如 '1', '2', '3' 或 'A', 'B', 'C' 或 'S1', 'S2', 'S3'

        Returns
        -------
        DukouData
            该Stage的数据
        """
        if self.stages is None:
            raise ValueError("数据中没有Stage标记")

        # Stage标记映射: S1->A (Roadian), S2->B (Wordian), S3->C (Capitanian)
        # 注意：根据新Excel文件，Stage标记按年龄从老到新为S1->S2->S3
        stage_map = {
            'A': ['A', 'S1', 'Decline1', 'Roadian'],
            'B': ['B', 'S2', 'Decline2', 'Wordian'],
            'C': ['C', 'S3', 'Decline3', 'Capitanian'],
        }

        # 查找匹配的stage
        valid_stages = stage_map.get(stage.upper(), [stage])

        # 处理可能的NaN值（转换为字符串比较）
        stages_str = np.array([str(s) if pd.notna(s) else '' for s in self.stages])
        mask = np.isin(stages_str, valid_stages)

        if not np.any(mask):
            raise ValueError(f"没有找到Stage '{stage}'的数据")

        return DukouData(
            ages=self.ages[mask],
            positions=self.positions[mask],
            sr_ratios=self.sr_ratios[mask],
            mn_sr=self.mn_sr[mask],
            sample_ids=[sid for sid, m in zip(self.sample_ids, mask) if m],
            uncertainties=self.uncertainties[mask],
            formations=self.formations[mask] if self.formations is not None else None,
            stages=self.stages[mask],
            diagenesis=self.diagenesis[mask] if self.diagenesis is not None else None
        )

    def get_interval_statistics(self, age_min: float, age_max: float) -> Dict:
        """获取特定区间的统计信息"""
        interval = self.get_interval_data(age_min, age_max)
        return {
            'n_samples': len(interval.sample_ids),
            'age_range': (float(interval.ages.min()), float(interval.ages.max())),
            'position_range': (float(interval.positions.min()), float(interval.positions.max())),
            'sr_range': (float(interval.sr_ratios.min()), float(interval.sr_ratios.max())),
            'sr_mean': float(interval.sr_ratios.mean()),
            'sr_std': float(interval.sr_ratios.std()),
            'sample_ids': interval.sample_ids
        }


class DukouDataLoader:
    """Dukou剖面数据加载器"""

    # 生物地层锚点 (来自manuscript Fig. 1C, 4)
    # 格式: (Sample ID, Position_m, Conodont_Zone, Age_Ma)
    BIOSTRATIGRAPHIC_ANCHORS = [
        ('DK-85', 120.10, 'Pre-Roadian', 283.0),      # Base Maokou
        ('DK-75', 144.20, 'J. nankingensis', 272.95), # Roadian base
        ('DK-53', 195.70, 'J. postserrata', 268.0),   # Wordian base
        ('DK-36', 221.30, 'J. prexuanhanensis', 265.1), # Mid-Capitanian
        ('DK-30', 232.00, 'J. xuanhanensis', 259.1),  # GLB (top Maokou)
        ('DK-18', 257.83, 'Wuchiapingian', 254.14),   # Wuchiaping base
    ]

    def __init__(self, data_dir: Optional[Path] = None):
        """
        初始化数据加载器

        Parameters
        ----------
        data_dir : Path, optional
            数据目录，默认为项目根目录下的data文件夹
        """
        if data_dir is None:
            self.data_dir = Path(__file__).parent.parent.parent / 'data'
        else:
            self.data_dir = Path(data_dir)

        # 创建年龄-深度插值函数
        self._setup_age_model()

    def _setup_age_model(self):
        """设置年龄模型 (线性插值)"""
        anchors = self.BIOSTRATIGRAPHIC_ANCHORS
        positions = np.array([a[1] for a in anchors])
        ages = np.array([a[3] for a in anchors])

        # 线性插值，超出范围用外推
        self._depth_to_age = interp1d(
            positions, ages,
            kind='linear',
            bounds_error=False,
            fill_value='extrapolate'
        )
        self._age_to_depth = interp1d(
            ages, positions,
            kind='linear',
            bounds_error=False,
            fill_value='extrapolate'
        )

    def depth_to_age(self, position: float) -> float:
        """深度转换为年龄"""
        return float(self._depth_to_age(position))

    def age_to_depth(self, age: float) -> float:
        """年龄转换为深度"""
        return float(self._age_to_depth(age))

    def load_data(
        self,
        mn_sr_threshold: float = 1.0,
        apply_lowess: bool = True,
        lowess_frac: float = 0.3,
        confidence_interval: float = 0.95
    ) -> DukouData:
        """
        加载并处理Dukou剖面数据

        Parameters
        ----------
        mn_sr_threshold : float
            Mn/Sr筛选阈值，默认<1.0
        apply_lowess : bool
            是否应用LOWESS平滑
        lowess_frac : float
            LOWESS平滑参数 (0-1)，越大越平滑
        confidence_interval : float
            置信区间 (默认95%)

        Returns
        -------
        DukouData
            处理后的数据对象
        """
        # 读取Sr同位素数据（优先使用新的Sr_sun.xlsx）
        sr_file = self.data_dir / 'Sr_sun.xlsx'
        if not sr_file.exists():
            # 回退到旧文件名
            sr_file = self.data_dir / 'Sr_Table_sun.xlsx'
        if not sr_file.exists():
            raise FileNotFoundError(f"Sr数据文件未找到: {sr_file}")

        df_sr = pd.read_excel(sr_file, skiprows=1)

        # 读取Mn/Sr数据
        mn_file = self.data_dir / 'Sr_MnSr_sun.xlsx'
        if not mn_file.exists():
            raise FileNotFoundError(f"Mn/Sr数据文件未找到: {mn_file}")

        df_mn = pd.read_excel(mn_file, skiprows=1)

        # 合并数据
        # 检查是否有stage列（用户自定义的Decline划分）
        stage_col = None
        for col in df_sr.columns:
            if col.lower() == 'stage':
                stage_col = col
                break

        # 检查是否有Diagenesis列（成岩作用标记）
        diagenesis_col = None
        for col in df_sr.columns:
            if col.lower() == 'diagenesis':
                diagenesis_col = col
                break

        # 选择需要的列
        sr_columns = ['Sample ID', 'Position', 'Formation', '87Sr/86Sr', '2σ']
        if stage_col:
            sr_columns.append(stage_col)
        if diagenesis_col:
            sr_columns.append(diagenesis_col)

        df = pd.merge(
            df_sr[sr_columns],
            df_mn[['Sample ID', 'Mn/Sr']],
            on='Sample ID',
            how='inner'
        )

        # 成岩作用筛选：排除Altered（成岩改造）的样本
        if diagenesis_col:
            n_before = len(df)
            df = df[df[diagenesis_col] != 'Altered'].copy()
            n_excluded = n_before - len(df)
            if n_excluded > 0:
                print(f"[Dukou Data] 成岩作用筛选: 排除 {n_excluded} 个Altered样本")
                print(f"  保留样本: {len(df)}/{n_before}")

        # Mn/Sr筛选
        df_filtered = df[df['Mn/Sr'] < mn_sr_threshold].copy()
        n_total = len(df)
        n_filtered = len(df_filtered)
        print(f"[Dukou Data] Mn/Sr筛选: {n_filtered}/{n_total} 样本通过 ({n_filtered/n_total*100:.1f}%)")

        # 按深度排序
        df_filtered = df_filtered.sort_values('Position')

        # 解析不确定度 (处理科学计数法如"4.25×10-5")
        uncertainties = df_filtered['2σ'].apply(self._parse_uncertainty).values

        # 深度转年龄
        ages = df_filtered['Position'].apply(self.depth_to_age).values

        # 获取stage信息（如果有）
        stages = None
        if stage_col:
            stages = df_filtered[stage_col].values

        # 获取diagenesis信息（如果有）
        diagenesis = None
        if diagenesis_col and diagenesis_col in df_filtered.columns:
            diagenesis = df_filtered[diagenesis_col].values

        # 创建数据对象
        data = DukouData(
            ages=ages,
            positions=df_filtered['Position'].values,
            sr_ratios=df_filtered['87Sr/86Sr'].values,
            mn_sr=df_filtered['Mn/Sr'].values,
            sample_ids=df_filtered['Sample ID'].tolist(),
            uncertainties=uncertainties,
            formations=df_filtered['Formation'].values,
            stages=stages,
            diagenesis=diagenesis
        )

        # 应用LOWESS平滑
        if apply_lowess:
            data = self._apply_lowess_smooth(data, lowess_frac, confidence_interval)

        return data

    def _parse_uncertainty(self, value) -> float:
        """解析不确定度字符串"""
        if isinstance(value, (int, float)):
            return float(value)

        if isinstance(value, str):
            # 处理科学计数法如 "4.25×10-5"
            value = value.replace('×10', 'e').replace('−', '-')
            try:
                return float(value)
            except ValueError:
                pass

        return 0.00001  # 默认不确定度

    def _apply_lowess_smooth(
        self,
        data: DukouData,
        frac: float = 0.3,
        confidence_interval: float = 0.95
    ) -> DukouData:
        """
        应用LOWESS平滑生成目标曲线

        Parameters
        ----------
        data : DukouData
            原始数据
        frac : float
            LOWESS参数，数据点比例用于每个局部拟合
        confidence_interval : float
            置信水平
        """
        if not HAS_STATSMODELS:
            print("[Warning] statsmodels未安装，跳过LOWESS平滑")
            print("  安装: pip install statsmodels")
            # 使用简单移动平均作为备用
            return self._apply_moving_average_smooth(data, confidence_interval)

        # 按年龄排序
        sort_idx = np.argsort(data.ages)
        ages_sorted = data.ages[sort_idx]
        sr_sorted = data.sr_ratios[sort_idx]

        # 生成均匀分布的年龄点用于平滑曲线
        age_min, age_max = ages_sorted.min(), ages_sorted.max()
        n_points = max(50, len(ages_sorted) * 2)
        target_ages = np.linspace(age_min, age_max, n_points)

        # LOWESS平滑
        # 注意: lowess需要(endog, exog)格式，即(y, x)
        lowess_result = lowess(
            sr_sorted, ages_sorted,
            frac=frac,
            it=3,  # 迭代次数
            delta=0.0,
            is_sorted=True
        )

        # 插值到目标年龄点
        interp_smooth = interp1d(
            lowess_result[:, 0],
            lowess_result[:, 1],
            kind='linear',
            bounds_error=False,
            fill_value='extrapolate'
        )
        target_sr = interp_smooth(target_ages)

        # 计算残差标准差用于置信区间
        residuals = sr_sorted - interp_smooth(ages_sorted)
        residual_std = np.std(residuals)

        # 95%置信区间
        z_score = stats.norm.ppf((1 + confidence_interval) / 2)
        ci_width = z_score * residual_std

        data.target_ages = target_ages
        data.target_sr = target_sr
        data.target_ci_lower = target_sr - ci_width
        data.target_ci_upper = target_sr + ci_width

        print(f"[Dukou Data] LOWESS平滑完成 (frac={frac})")
        print(f"  目标曲线年龄范围: {age_max:.1f} - {age_min:.1f} Ma")
        print(f"  Sr范围: {target_sr.min():.5f} - {target_sr.max():.5f}")
        print(f"  95% CI宽度: ±{ci_width:.5f}")

        return data

    def _apply_moving_average_smooth(
        self,
        data: DukouData,
        confidence_interval: float = 0.95,
        window: int = 5
    ) -> DukouData:
        """备用: 使用移动平均平滑（正确处理边界）"""
        # 按年龄排序
        sort_idx = np.argsort(data.ages)
        ages_sorted = data.ages[sort_idx]
        sr_sorted = data.sr_ratios[sort_idx]

        # 移动平均 - 使用'valid'模式避免边界效应
        # 然后使用线性插值填充边界
        half_window = window // 2

        # 计算内部点的移动平均
        ma_valid = np.convolve(sr_sorted, np.ones(window)/window, mode='valid')

        # 构建完整的平滑曲线
        target_sr = np.zeros_like(sr_sorted)

        # 边界：使用原始值（或简单平均）
        target_sr[:half_window] = sr_sorted[:half_window]  # 前half_window个保持原值
        target_sr[-half_window:] = sr_sorted[-half_window:]  # 后half_window个保持原值

        # 中间：移动平均
        target_sr[half_window:-half_window] = ma_valid

        # 计算残差
        residuals = sr_sorted - target_sr
        residual_std = np.std(residuals)

        z_score = stats.norm.ppf((1 + confidence_interval) / 2)
        ci_width = z_score * residual_std

        data.target_ages = ages_sorted
        data.target_sr = target_sr
        data.target_ci_lower = target_sr - ci_width
        data.target_ci_upper = target_sr + ci_width

        print(f"[Dukou Data] 移动平均平滑完成 (window={window})")
        print(f"  Sr范围: {target_sr.min():.5f} - {target_sr.max():.5f}")

        return data

    def get_anchors_df(self) -> pd.DataFrame:
        """获取生物地层锚点DataFrame"""
        return pd.DataFrame(
            self.BIOSTRATIGRAPHIC_ANCHORS,
            columns=['Sample_ID', 'Position_m', 'Conodont_Zone', 'Age_Ma']
        )


def load_dukou_data(
    mn_sr_threshold: float = 1.0,
    lowess_frac: float = 0.3,
    confidence_interval: float = 0.95
) -> DukouData:
    """
    便捷函数: 加载Dukou数据

    Examples
    --------
    >>> data = load_dukou_data()
    >>> print(f"样本数: {len(data.sample_ids)}")
    >>> print(f"年龄范围: {data.ages.min():.1f} - {data.ages.max():.1f} Ma")
    """
    loader = DukouDataLoader()
    return loader.load_data(
        mn_sr_threshold=mn_sr_threshold,
        apply_lowess=True,
        lowess_frac=lowess_frac,
        confidence_interval=confidence_interval
    )


if __name__ == '__main__':
    # 测试数据加载
    print("="*60)
    print("Dukou剖面数据加载测试")
    print("="*60)

    data = load_dukou_data()

    print("\n数据概览:")
    print(f"  样本数: {len(data.sample_ids)}")
    print(f"  年龄范围: {data.ages.min():.1f} - {data.ages.max():.1f} Ma")
    print(f"  深度范围: {data.positions.min():.1f} - {data.positions.max():.1f} m")
    print(f"  Sr同位素范围: {data.sr_ratios.min():.5f} - {data.sr_ratios.max():.5f}")

    print("\n生物地层锚点:")
    loader = DukouDataLoader()
    print(loader.get_anchors_df().to_string(index=False))

    # 测试区间数据提取 - 三个下降期
    print("\n" + "="*60)
    print("三个下降期数据")
    print("="*60)

    # Interval 2: Roadian (Decline 1)
    print("\nInterval 2 - Roadian (Decline 1):")
    stats2 = data.get_interval_statistics(268, 273)
    print(f"  年龄: 268-273 Ma")
    print(f"  样本: {stats2['n_samples']}个 {stats2['sample_ids'][:5]}...")
    print(f"  Sr: {stats2['sr_range'][0]:.5f} - {stats2['sr_range'][1]:.5f}")

    # Interval 3: Wordian (Decline 2)
    print("\nInterval 3 - Wordian (Decline 2):")
    stats3 = data.get_interval_statistics(265, 268)
    print(f"  年龄: 265-268 Ma")
    print(f"  样本: {stats3['n_samples']}个")
    print(f"  Sr: {stats3['sr_range'][0]:.5f} - {stats3['sr_range'][1]:.5f}")

    # Interval 4: Capitanian (Decline 3)
    print("\nInterval 4 - Capitanian (Decline 3):")
    stats4 = data.get_interval_statistics(259.1, 265)
    print(f"  年龄: 259.1-265 Ma")
    print(f"  样本: {stats4['n_samples']}个")
    print(f"  Sr: {stats4['sr_range'][0]:.5f} - {stats4['sr_range'][1]:.5f}")
