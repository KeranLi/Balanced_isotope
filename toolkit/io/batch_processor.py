"""
批量处理器模块

提供同位素数据的批量处理功能
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Union, Callable, Any
from dataclasses import dataclass, field
from tqdm import tqdm
import warnings

from toolkit.io.data_handler import DataHandler, DataFormatError


@dataclass
class ProcessingResult:
    """处理结果数据类"""
    success: bool
    index: int
    input_data: Dict[str, Any]
    output_data: Dict[str, Any] = field(default_factory=dict)
    error_message: str = ""
    processing_time: float = 0.0


class BatchProcessor:
    """
    同位素数据批量处理器
    
    支持从文件批量读取、处理和输出结果
    
    示例用法:
        processor = BatchProcessor('u')
        results = processor.process_file('data/input.xlsx')
        processor.save_results(results, 'output/results.xlsx')
    """
    
    def __init__(self, element: str, 
                 scenario: str = 'modern',
                 apply_diagenetic_correction: bool = True,
                 delta_diag: float = 0.4,
                 include_uncertainty: bool = True,
                 n_monte_carlo: int = 1000):
        """
        初始化批量处理器
        
        Parameters
        ----------
        element : str
            处理的元素，如'u', 'n', 'c', 'mg'
        scenario : str
            模型情景
        apply_diagenetic_correction : bool
            是否应用成岩校正（U同位素）
        delta_diag : float
            成岩校正值
        include_uncertainty : bool
            是否计算不确定度
        n_monte_carlo : int
            蒙特卡洛采样次数
        """
        self.element = element.lower()
        self.scenario = scenario
        self.apply_diagenetic_correction = apply_diagenetic_correction
        self.delta_diag = delta_diag
        self.include_uncertainty = include_uncertainty
        self.n_monte_carlo = n_monte_carlo
        
        # 初始化同位素体系
        self.system = self._init_system()
        self.analyzer = None
        
        if self.include_uncertainty and element.lower() == 'u':
            from systems.u import UncertaintyAnalyzer
            self.analyzer = UncertaintyAnalyzer(self.system)
    
    def _init_system(self):
        """初始化同位素体系"""
        if self.element == 'u':
            from systems.u import UIsotopeSystem
            return UIsotopeSystem(scenario=self.scenario)
        elif self.element == 'n':
            from systems.n import NIsotopeSystem
            return NIsotopeSystem(scenario=self.scenario)
        elif self.element == 'c':
            from systems.c import CIsotopeSystem
            return CIsotopeSystem(scenario=self.scenario)
        elif self.element == 'mg':
            from systems.mg import MgIsotopeSystem
            return MgIsotopeSystem()
        else:
            raise ValueError(f"Unsupported element: {self.element}")
    
    def process_file(self, 
                    input_path: Union[str, Path],
                    output_path: Optional[Union[str, Path]] = None,
                    show_progress: bool = True) -> pd.DataFrame:
        """
        处理整个数据文件
        
        Parameters
        ----------
        input_path : str or Path
            输入文件路径
        output_path : str or Path, optional
            输出文件路径，如果为None则不保存
        show_progress : bool
            是否显示进度条
            
        Returns
        -------
        pd.DataFrame
            包含原始数据和处理结果的DataFrame
        """
        # 读取数据
        handler = DataHandler(input_path)
        data = handler.read()
        
        # 验证数据
        is_valid, errors = handler.validate_data(self.element)
        if not is_valid:
            raise DataFormatError(f"Data validation failed: {', '.join(errors)}")
        
        # 获取同位素数据
        isotope_data = handler.get_isotope_data(self.element)
        
        # 处理每一行
        results = []
        n_samples = len(data)
        
        iterator = range(n_samples)
        if show_progress:
            try:
                iterator = tqdm(iterator, desc=f"Processing {self.element}")
            except ImportError:
                pass
        
        for i in iterator:
            result = self._process_single(i, isotope_data)
            results.append(result)
        
        # 合并结果到原始数据
        result_df = self._merge_results(data, results)
        
        # 保存结果
        if output_path:
            self.save_results(result_df, output_path)
            print(f"\nResults saved to: {output_path}")
        
        return result_df
    
    def _process_single(self, index: int, 
                       isotope_data: Dict[str, np.ndarray]) -> ProcessingResult:
        """
        处理单个样品
        
        Parameters
        ----------
        index : int
            样品索引
        isotope_data : dict
            同位素数据字典
            
        Returns
        -------
        ProcessingResult
            处理结果
        """
        import time
        start_time = time.time()
        
        try:
            if self.element == 'u':
                return self._process_u(index, isotope_data)
            elif self.element == 'n':
                return self._process_n(index, isotope_data)
            elif self.element == 'c':
                return self._process_c(index, isotope_data)
            elif self.element == 'mg':
                return self._process_mg(index, isotope_data)
            else:
                raise ValueError(f"Unsupported element: {self.element}")
        except Exception as e:
            return ProcessingResult(
                success=False,
                index=index,
                input_data={},
                error_message=str(e),
                processing_time=time.time() - start_time
            )
    
    def _process_u(self, index: int, 
                  isotope_data: Dict[str, np.ndarray]) -> ProcessingResult:
        """处理U同位素数据"""
        delta238 = isotope_data['delta'][index]
        
        # 获取不确定度（如果有）
        delta_std = isotope_data.get('delta_std', [None])[index] if 'delta_std' in isotope_data else None
        
        # 获取元数据
        input_data = {'delta238_carb': delta238}
        for key in ['sample_id', 'depth', 'age', 'elevation']:
            if key in isotope_data:
                input_data[key] = isotope_data[key][index]
        
        # 稳态计算
        result = self.system.calculate_f_anox_steady_state(
            delta238_carb=delta238,
            apply_diagenetic_correction=self.apply_diagenetic_correction,
            delta_diag=self.delta_diag
        )
        
        output_data = {
            'f_anox': result['f_anox'],
            'f_oxic': result['f_oxic'],
            'delta238_seawater': result['delta238_seawater'],
            'delta238_carb_corrected': result['delta238_carb_corrected'],
            'anoxic_area_percent': self.system.estimate_anoxic_area(result['f_anox'])
        }
        
        # 不确定度分析
        if self.include_uncertainty and self.analyzer and delta_std is not None:
            try:
                mc_result = self.analyzer.monte_carlo_steady_state(
                    delta238_carb=delta238,
                    n_samples=self.n_monte_carlo,
                    random_seed=42 + index  # 确保可复现但不同样品不同
                )
                output_data['f_anox_std'] = mc_result['f_anox_std']
                output_data['f_anox_ci_lower'] = mc_result['f_anox_ci'][0]
                output_data['f_anox_ci_upper'] = mc_result['f_anox_ci'][1]
            except Exception:
                pass
        
        return ProcessingResult(
            success=True,
            index=index,
            input_data=input_data,
            output_data=output_data
        )
    
    def _process_n(self, index: int,
                  isotope_data: Dict[str, np.ndarray]) -> ProcessingResult:
        """处理N同位素数据"""
        delta15 = isotope_data['delta'][index]
        
        input_data = {'delta15_sed': delta15}
        for key in ['sample_id', 'depth', 'age', 'elevation']:
            if key in isotope_data:
                input_data[key] = isotope_data[key][index]
        
        # 反向模型：从沉积物δ¹⁵N计算f_assimilator
        result = self.system.inverse_model(delta15_sed=delta15)
        
        output_data = {
            'f_assimilator': result['f_assimilator'],
            'delta15_sed_calculated': result['delta15N_sed_calculated'],
            'residual': result['residual']
        }
        
        return ProcessingResult(
            success=True,
            index=index,
            input_data=input_data,
            output_data=output_data
        )
    
    def _process_c(self, index: int,
                  isotope_data: Dict[str, np.ndarray]) -> ProcessingResult:
        """处理C同位素数据"""
        delta13 = isotope_data['delta'][index]
        
        input_data = {'delta13_carb': delta13}
        for key in ['sample_id', 'depth', 'age', 'elevation']:
            if key in isotope_data:
                input_data[key] = isotope_data[key][index]
        
        # C同位素模型通常用于稳态计算
        # 这里简化为返回输入值
        output_data = {
            'delta13_carb': delta13
        }
        
        return ProcessingResult(
            success=True,
            index=index,
            input_data=input_data,
            output_data=output_data
        )
    
    def _process_mg(self, index: int,
                   isotope_data: Dict[str, np.ndarray]) -> ProcessingResult:
        """处理Mg同位素数据"""
        input_data = {}
        
        if 'delta26' in isotope_data:
            input_data['delta26'] = isotope_data['delta26'][index]
        if 'delta25' in isotope_data:
            input_data['delta25'] = isotope_data['delta25'][index]
        
        for key in ['sample_id', 'depth', 'age', 'elevation']:
            if key in isotope_data:
                input_data[key] = isotope_data[key][index]
        
        # Mg同位素风化比例计算
        output_data = {}
        if 'delta26' in input_data:
            try:
                ratios = self.system.calculate_weathering_ratio(
                    delta_sample=input_data['delta26'],
                    delta_seawater=-0.83
                )
                output_data['f_carbonate'] = ratios['f_carbonate']
                output_data['f_silicate'] = ratios['f_silicate']
            except Exception:
                pass
        
        return ProcessingResult(
            success=True,
            index=index,
            input_data=input_data,
            output_data=output_data
        )
    
    def _merge_results(self, original_data: pd.DataFrame,
                      results: List[ProcessingResult]) -> pd.DataFrame:
        """
        合并处理结果到原始数据
        
        Parameters
        ----------
        original_data : pd.DataFrame
            原始数据
        results : list
            处理结果列表
            
        Returns
        -------
        pd.DataFrame
            合并后的DataFrame
        """
        df = original_data.copy()
        
        # 提取输出数据
        output_keys = set()
        for result in results:
            output_keys.update(result.output_data.keys())
        
        # 添加新列
        for key in output_keys:
            values = []
            for result in results:
                values.append(result.output_data.get(key, np.nan))
            df[key] = values
        
        # 添加处理状态列
        df['processing_success'] = [r.success for r in results]
        df['processing_error'] = [r.error_message if r.error_message else '' for r in results]
        
        return df
    
    def save_results(self, 
                    results: Union[pd.DataFrame, List[ProcessingResult]],
                    output_path: Union[str, Path],
                    sheet_name: str = 'Results'):
        """
        保存结果到文件
        
        Parameters
        ----------
        results : DataFrame or list
            处理结果
        output_path : str or Path
            输出文件路径
        sheet_name : str
            Excel工作表名称
        """
        if isinstance(results, list):
            # 转换为DataFrame
            data = []
            for result in results:
                row = {**result.input_data, **result.output_data}
                row['success'] = result.success
                row['error'] = result.error_message
                data.append(row)
            df = pd.DataFrame(data)
        else:
            df = results
        
        handler = DataHandler()
        handler.write(df, output_path, sheet_name=sheet_name)
    
    def get_summary(self, results: List[ProcessingResult]) -> Dict[str, Any]:
        """
        获取处理结果摘要
        
        Parameters
        ----------
        results : list
            处理结果列表
            
        Returns
        -------
        dict
            摘要统计信息
        """
        n_total = len(results)
        n_success = sum(1 for r in results if r.success)
        n_failed = n_total - n_success
        
        summary = {
            'total_samples': n_total,
            'successful': n_success,
            'failed': n_failed,
            'success_rate': n_success / n_total if n_total > 0 else 0
        }
        
        # 统计输出变量
        if results and results[0].success:
            for key in results[0].output_data.keys():
                values = [r.output_data.get(key, np.nan) 
                         for r in results if r.success and key in r.output_data]
                if values:
                    summary[f'{key}_mean'] = np.nanmean(values)
                    summary[f'{key}_std'] = np.nanstd(values)
                    summary[f'{key}_min'] = np.nanmin(values)
                    summary[f'{key}_max'] = np.nanmax(values)
        
        return summary
