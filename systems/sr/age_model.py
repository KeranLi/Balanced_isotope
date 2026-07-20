"""
Sr同位素年龄模型模块

基于Dukou剖面生物地层锚点的深度-年龄转换
"""

import numpy as np
from scipy.interpolate import interp1d
from typing import Optional, Tuple, List
from dataclasses import dataclass


@dataclass
class BiostratigraphicAnchor:
    """生物地层锚点"""
    sample_id: str
    position_m: float
    conodont_zone: str
    age_ma: float
    uncertainty_ma: float = 0.5  # 默认±0.5 Ma不确定度


class AgeModel:
    """
    年龄模型 - 深度与年龄的双向转换

    基于生物地层锚点的线性插值模型
    """

    # 默认锚点 (Dukou剖面)
    DEFAULT_ANCHORS = [
        BiostratigraphicAnchor('DK-85', 120.10, 'Pre-Roadian', 283.0, 1.0),
        BiostratigraphicAnchor('DK-75', 144.20, 'J. nankingensis', 272.95, 0.5),
        BiostratigraphicAnchor('DK-53', 195.70, 'J. postserrata', 268.0, 0.5),
        BiostratigraphicAnchor('DK-36', 221.30, 'J. prexuanhanensis', 265.1, 0.5),
        BiostratigraphicAnchor('DK-30', 232.00, 'J. xuanhanensis', 259.1, 0.5),
        BiostratigraphicAnchor('DK-18', 257.83, 'Wuchiapingian', 254.14, 1.0),
    ]

    def __init__(self, anchors: Optional[List[BiostratigraphicAnchor]] = None):
        """
        初始化年龄模型

        Parameters
        ----------
        anchors : list, optional
            自定义生物地层锚点列表，默认使用Dukou剖面锚点
        """
        self.anchors = anchors if anchors is not None else self.DEFAULT_ANCHORS.copy()
        self._setup_interpolation()

    def _setup_interpolation(self):
        """设置插值函数"""
        positions = np.array([a.position_m for a in self.anchors])
        ages = np.array([a.age_ma for a in self.anchors])

        # 按深度排序
        sort_idx = np.argsort(positions)
        positions_sorted = positions[sort_idx]
        ages_sorted = ages[sort_idx]

        # 深度 -> 年龄
        self._depth_to_age_interp = interp1d(
            positions_sorted, ages_sorted,
            kind='linear',
            bounds_error=False,
            fill_value='extrapolate'
        )

        # 年龄 -> 深度
        self._age_to_depth_interp = interp1d(
            ages_sorted, positions_sorted,
            kind='linear',
            bounds_error=False,
            fill_value='extrapolate'
        )

        # 存储范围
        self.position_range = (positions_sorted.min(), positions_sorted.max())
        self.age_range = (ages_sorted.min(), ages_sorted.max())

    def depth_to_age(self, position: float) -> float:
        """
        深度转换为年龄

        Parameters
        ----------
        position : float
            深度 (m)

        Returns
        -------
        float
            年龄 (Ma)
        """
        return float(self._depth_to_age_interp(position))

    def age_to_depth(self, age: float) -> float:
        """
        年龄转换为深度

        Parameters
        ----------
        age : float
            年龄 (Ma)

        Returns
        -------
        float
            深度 (m)
        """
        return float(self._age_to_depth_interp(age))

    def depth_to_age_array(self, positions: np.ndarray) -> np.ndarray:
        """数组版本: 深度转年龄"""
        return self._depth_to_age_interp(positions)

    def age_to_depth_array(self, ages: np.ndarray) -> np.ndarray:
        """数组版本: 年龄转深度"""
        return self._age_to_depth_interp(ages)

    def get_sedimentation_rate(self, position: float, window: float = 10.0) -> float:
        """
        计算沉积速率 (m/Myr)

        Parameters
        ----------
        position : float
            深度位置 (m)
        window : float
            计算窗口 (m)

        Returns
        -------
        float
            沉积速率 (m/Myr)
        """
        pos_min = max(position - window/2, self.position_range[0])
        pos_max = min(position + window/2, self.position_range[1])

        age_min = self.depth_to_age(pos_min)
        age_max = self.depth_to_age(pos_max)

        return (pos_max - pos_min) / (age_max - age_min)

    def get_interval_duration(self, position_start: float, position_end: float) -> float:
        """
        计算两个深度之间的时间跨度

        Parameters
        ----------
        position_start : float
            起始深度 (m)
        position_end : float
            结束深度 (m)

        Returns
        -------
        float
            时间跨度 (Myr)
        """
        age_start = self.depth_to_age(position_start)
        age_end = self.depth_to_age(position_end)
        return abs(age_end - age_start)

    def add_anchor(self, anchor: BiostratigraphicAnchor):
        """添加新的锚点并重新计算插值"""
        self.anchors.append(anchor)
        self._setup_interpolation()

    def get_anchors_table(self) -> dict:
        """获取锚点表格"""
        return {
            'Sample ID': [a.sample_id for a in self.anchors],
            'Position (m)': [a.position_m for a in self.anchors],
            'Conodont Zone': [a.conodont_zone for a in self.anchors],
            'Age (Ma)': [a.age_ma for a in self.anchors],
            'Uncertainty (Ma)': [a.uncertainty_ma for a in self.anchors],
        }

    def print_summary(self):
        """打印模型摘要"""
        print("=" * 60)
        print("年龄模型摘要")
        print("=" * 60)
        print(f"\n深度范围: {self.position_range[0]:.1f} - {self.position_range[1]:.1f} m")
        print(f"年龄范围: {self.age_range[1]:.1f} - {self.age_range[0]:.1f} Ma")
        print(f"\n生物地层锚点:")
        for a in self.anchors:
            print(f"  {a.sample_id:8s}  {a.position_m:6.1f}m  {a.age_ma:6.2f}Ma  {a.conodont_zone}")


# 预定义的Interval时间范围 (Ma)
INTERVALS = {
    'interval_1': {  # Late Early Permian (Chihsia)
        'name': 'Late Early Permian',
        'age_range': (273.0, 283.0),
        'description': 'Chihsia Formation, Sr increase to peak'
    },
    'interval_2': {  # Roadian (Decline 1)
        'name': 'Roadian',
        'age_range': (268.0, 273.0),
        'description': 'Decline 1: P3 glaciation reduces weathering'
    },
    'interval_3': {  # Wordian-early Capitanian (Decline 2)
        'name': 'Wordian-early Capitanian',
        'age_range': (265.0, 268.0),
        'description': 'Decline 2: Neo-Tethys expansion'
    },
    'interval_4': {  # Mid-late Capitanian (Decline 3)
        'name': 'Mid-late Capitanian',
        'age_range': (259.1, 265.0),
        'description': 'Decline 3: Emeishan LIP basalt weathering'
    },
    'interval_5': {  # Wuchiapingian
        'name': 'Wuchiapingian',
        'age_range': (254.14, 259.1),
        'description': 'Post-GLB recovery'
    },
}


def get_interval_age_range(interval_name: str) -> Tuple[float, float]:
    """
    获取Interval的年龄范围

    Parameters
    ----------
    interval_name : str
        'interval_1' 到 'interval_5' 或 'roadian', 'wordian' 等

    Returns
    -------
    tuple
        (age_min, age_max) in Ma
    """
    interval_name = interval_name.lower()

    # 直接匹配
    if interval_name in INTERVALS:
        return INTERVALS[interval_name]['age_range']

    # 别名匹配
    aliases = {
        'roadian': 'interval_2',
        'decline_1': 'interval_2',
        'wordian': 'interval_3',
        'decline_2': 'interval_3',
        'capitanian': 'interval_4',
        'decline_3': 'interval_4',
        'chihsia': 'interval_1',
        'wuchiapingian': 'interval_5',
    }

    if interval_name in aliases:
        return INTERVALS[aliases[interval_name]]['age_range']

    raise ValueError(f"未知Interval: {interval_name}")


if __name__ == '__main__':
    # 测试年龄模型
    print("=" * 60)
    print("年龄模型测试")
    print("=" * 60)

    model = AgeModel()
    model.print_summary()

    # 测试转换
    print("\n" + "=" * 60)
    print("深度-年龄转换测试")
    print("=" * 60)

    test_positions = [120.1, 144.2, 195.7, 221.3, 232.0, 257.83]
    print("\n深度 -> 年龄:")
    for pos in test_positions:
        age = model.depth_to_age(pos)
        print(f"  {pos:6.1f} m -> {age:6.2f} Ma")

    # 测试反向转换
    print("\n年龄 -> 深度:")
    test_ages = [283.0, 272.95, 268.0, 265.1, 259.1, 254.14]
    for age in test_ages:
        pos = model.age_to_depth(age)
        print(f"  {age:6.2f} Ma -> {pos:6.1f} m")

    # 测试沉积速率
    print("\n沉积速率计算:")
    for pos in [130, 170, 210, 245]:
        rate = model.get_sedimentation_rate(pos)
        print(f"  {pos} m 附近: {rate:.2f} m/Myr")

    # 测试Interval
    print("\n" + "=" * 60)
    print("Interval时间范围")
    print("=" * 60)
    for key, info in INTERVALS.items():
        age_min, age_max = info['age_range']
        print(f"\n{key}: {info['name']}")
        print(f"  年龄: {age_max:.1f} - {age_min:.1f} Ma")
        print(f"  描述: {info['description']}")
