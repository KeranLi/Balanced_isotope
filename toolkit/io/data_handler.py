"""
数据处理器模块

提供同位素数据的读取、验证和写入功能
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple, Any
from dataclasses import dataclass
import warnings


class DataFormatError(Exception):
    """数据格式错误"""
    pass


@dataclass
class ColumnMapping:
    """列名映射配置"""
    # 铀同位素
    delta238_u: List[str] = None
    delta238_u_std: List[str] = None
    
    # 氮同位素
    delta15_n: List[str] = None
    delta15_n_std: List[str] = None
    
    # 碳同位素
    delta13_c: List[str] = None
    delta13_c_std: List[str] = None
    
    # 镁同位素
    delta26_mg: List[str] = None
    delta26_mg_std: List[str] = None
    delta25_mg: List[str] = None
    delta25_mg_std: List[str] = None
    
    # 通用列
    sample_id: List[str] = None
    depth: List[str] = None
    age: List[str] = None
    elevation: List[str] = None
    
    def __post_init__(self):
        # 设置默认列名映射
        if self.delta238_u is None:
            self.delta238_u = ['delta_238_u', 'delta238U', 'd238U', 'δ238U']
        if self.delta238_u_std is None:
            self.delta238_u_std = ['delta_238_u_std', 'delta238U_std', 'd238U_2sd', 'δ238U_err']
            
        if self.delta15_n is None:
            self.delta15_n = ['delta_15_n', 'delta15N', 'd15N', 'δ15N']
        if self.delta15_n_std is None:
            self.delta15_n_std = ['delta_15_n_std', 'delta15N_std', 'd15N_2sd', 'δ15N_err']
            
        if self.delta13_c is None:
            self.delta13_c = ['delta_13_c', 'delta13C', 'd13C', 'δ13C']
        if self.delta13_c_std is None:
            self.delta13_c_std = ['delta_13_c_std', 'delta13C_std', 'd13C_2sd', 'δ13C_err']
            
        if self.delta26_mg is None:
            self.delta26_mg = ['delta_26_mg', 'delta26Mg', 'd26Mg', 'δ26Mg', 'delta_26_Mg_iso']
        if self.delta26_mg_std is None:
            self.delta26_mg_std = ['delta_26_mg_std', 'delta26Mg_std', 'd26Mg_2sd', 'δ26Mg_err', 'delta_26_Mg_iso_2sd']
        if self.delta25_mg is None:
            self.delta25_mg = ['delta_25_mg', 'delta25Mg', 'd25Mg', 'δ25Mg', 'delta_25_Mg_iso']
        if self.delta25_mg_std is None:
            self.delta25_mg_std = ['delta_25_mg_std', 'delta25Mg_std', 'd25Mg_2sd', 'δ25Mg_err', 'delta_25_Mg_iso_2sd']
            
        if self.sample_id is None:
            self.sample_id = ['sample_id', 'sample', 'id', 'Sample', 'Sample_ID']
        if self.depth is None:
            self.depth = ['depth', 'Depth', 'depth_m', 'stratigraphic_depth']
        if self.age is None:
            self.age = ['age', 'Age', 'age_ma', 'Age_Ma']
        if self.elevation is None:
            self.elevation = ['elevation', 'Elevation', 'elev', 'elev_m', 'stratigraphic_elevation']


class DataHandler:
    """
    同位素数据处理器
    
    负责读取、验证和写入同位素数据文件
    """
    
    # 支持的文件格式
    SUPPORTED_FORMATS = ['.xlsx', '.xls', '.csv']
    
    def __init__(self, file_path: Optional[Union[str, Path]] = None,
                 column_mapping: Optional[ColumnMapping] = None):
        """
        初始化数据处理器
        
        Parameters
        ----------
        file_path : str or Path, optional
            数据文件路径
        column_mapping : ColumnMapping, optional
            自定义列名映射
        """
        self.file_path = Path(file_path) if file_path else None
        self.column_mapping = column_mapping or ColumnMapping()
        self._data: Optional[pd.DataFrame] = None
        self._detected_columns: Dict[str, str] = {}
    
    def read(self, file_path: Optional[Union[str, Path]] = None,
             sheet_name: Union[str, int] = 0,
             **kwargs) -> pd.DataFrame:
        """
        读取数据文件
        
        Parameters
        ----------
        file_path : str or Path, optional
            文件路径，如果初始化时已提供则可为None
        sheet_name : str or int
            Excel工作表名称或索引
        **kwargs
            传递给pandas读取函数的额外参数
            
        Returns
        -------
        pd.DataFrame
            读取的数据
            
        Raises
        ------
        DataFormatError
            文件格式不支持或读取失败
        """
        path = Path(file_path) if file_path else self.file_path
        if path is None:
            raise DataFormatError("No file path provided")
        
        if not path.exists():
            raise DataFormatError(f"File not found: {path}")
        
        suffix = path.suffix.lower()
        
        try:
            if suffix in ['.xlsx', '.xls']:
                self._data = pd.read_excel(path, sheet_name=sheet_name, **kwargs)
            elif suffix == '.csv':
                self._data = pd.read_csv(path, **kwargs)
            else:
                raise DataFormatError(
                    f"Unsupported file format: {suffix}. "
                    f"Supported: {', '.join(self.SUPPORTED_FORMATS)}"
                )
        except Exception as e:
            raise DataFormatError(f"Failed to read file: {e}")
        
        # 自动检测列
        self._detect_columns()
        
        return self._data
    
    def _detect_columns(self):
        """自动检测数据列"""
        if self._data is None:
            return
        
        columns = self._data.columns.tolist()
        self._detected_columns = {}
        
        # 检测各同位素列
        for col_type, col_names in [
            ('delta238_u', self.column_mapping.delta238_u),
            ('delta238_u_std', self.column_mapping.delta238_u_std),
            ('delta15_n', self.column_mapping.delta15_n),
            ('delta15_n_std', self.column_mapping.delta15_n_std),
            ('delta13_c', self.column_mapping.delta13_c),
            ('delta13_c_std', self.column_mapping.delta13_c_std),
            ('delta26_mg', self.column_mapping.delta26_mg),
            ('delta26_mg_std', self.column_mapping.delta26_mg_std),
            ('delta25_mg', self.column_mapping.delta25_mg),
            ('delta25_mg_std', self.column_mapping.delta25_mg_std),
            ('sample_id', self.column_mapping.sample_id),
            ('depth', self.column_mapping.depth),
            ('age', self.column_mapping.age),
            ('elevation', self.column_mapping.elevation),
        ]:
            for col in columns:
                if col in col_names:
                    self._detected_columns[col_type] = col
                    break
    
    def get_column(self, col_type: str) -> Optional[str]:
        """
        获取检测到的列名
        
        Parameters
        ----------
        col_type : str
            列类型，如'delta238_u', 'sample_id'等
            
        Returns
        -------
        str or None
            检测到的列名，未检测到则返回None
        """
        return self._detected_columns.get(col_type)
    
    def has_isotope_data(self, element: str) -> bool:
        """
        检查是否包含指定元素的同位素数据
        
        Parameters
        ----------
        element : str
            元素符号，如'u', 'n', 'c', 'mg'
            
        Returns
        -------
        bool
            是否包含该元素的同位素数据
        """
        element = element.lower()
        mapping = {
            'u': ['delta238_u'],
            'n': ['delta15_n'],
            'c': ['delta13_c'],
            'mg': ['delta26_mg', 'delta25_mg']
        }
        
        required_cols = mapping.get(element, [])
        return any(col in self._detected_columns for col in required_cols)
    
    def get_isotope_data(self, element: str, 
                        include_uncertainty: bool = True) -> Dict[str, np.ndarray]:
        """
        获取指定元素的同位素数据
        
        Parameters
        ----------
        element : str
            元素符号
        include_uncertainty : bool
            是否包含不确定度数据
            
        Returns
        -------
        dict
            同位素数据字典
        """
        if self._data is None:
            raise DataFormatError("No data loaded. Call read() first.")
        
        element = element.lower()
        result = {}
        
        # 获取数据列
        if element == 'u':
            col = self.get_column('delta238_u')
            if col:
                result['delta'] = self._data[col].values
            if include_uncertainty:
                std_col = self.get_column('delta238_u_std')
                if std_col:
                    result['delta_std'] = self._data[std_col].values
                    
        elif element == 'n':
            col = self.get_column('delta15_n')
            if col:
                result['delta'] = self._data[col].values
            if include_uncertainty:
                std_col = self.get_column('delta15_n_std')
                if std_col:
                    result['delta_std'] = self._data[std_col].values
                    
        elif element == 'c':
            col = self.get_column('delta13_c')
            if col:
                result['delta'] = self._data[col].values
            if include_uncertainty:
                std_col = self.get_column('delta13_c_std')
                if std_col:
                    result['delta_std'] = self._data[std_col].values
                    
        elif element == 'mg':
            col_26 = self.get_column('delta26_mg')
            col_25 = self.get_column('delta25_mg')
            if col_26:
                result['delta26'] = self._data[col_26].values
            if col_25:
                result['delta25'] = self._data[col_25].values
            if include_uncertainty:
                std_26 = self.get_column('delta26_mg_std')
                std_25 = self.get_column('delta25_mg_std')
                if std_26:
                    result['delta26_std'] = self._data[std_26].values
                if std_25:
                    result['delta25_std'] = self._data[std_25].values
        
        # 获取元数据列
        for meta_col in ['sample_id', 'depth', 'age', 'elevation']:
            col = self.get_column(meta_col)
            if col:
                result[meta_col] = self._data[col].values
        
        return result
    
    def validate_data(self, element: str) -> Tuple[bool, List[str]]:
        """
        验证数据完整性
        
        Parameters
        ----------
        element : str
            要验证的元素
            
        Returns
        -------
        (is_valid, errors) : tuple
            验证结果和错误信息列表
        """
        errors = []
        
        if self._data is None:
            return False, ["No data loaded"]
        
        if not self.has_isotope_data(element):
            errors.append(f"No isotope data found for element '{element}'")
        
        # 检查缺失值
        data = self.get_isotope_data(element, include_uncertainty=False)
        if 'delta' in data:
            n_missing = np.sum(np.isnan(data['delta']))
            if n_missing > 0:
                errors.append(f"Found {n_missing} missing values in isotope data")
        
        return len(errors) == 0, errors
    
    def write(self, data: Union[pd.DataFrame, Dict, List[Dict]],
              output_path: Union[str, Path],
              sheet_name: str = 'Results',
              **kwargs):
        """
        写入数据到文件
        
        Parameters
        ----------
        data : DataFrame, dict, or list of dict
            要写入的数据
        output_path : str or Path
            输出文件路径
        sheet_name : str
            Excel工作表名称
        **kwargs
            传递给pandas写入函数的额外参数
        """
        path = Path(output_path)
        
        # 转换为DataFrame
        if isinstance(data, dict):
            df = pd.DataFrame([data])
        elif isinstance(data, list) and len(data) > 0 and isinstance(data[0], dict):
            df = pd.DataFrame(data)
        elif isinstance(data, pd.DataFrame):
            df = data
        else:
            raise DataFormatError(f"Unsupported data type: {type(data)}")
        
        suffix = path.suffix.lower()
        
        try:
            if suffix in ['.xlsx', '.xls']:
                df.to_excel(path, sheet_name=sheet_name, index=False, **kwargs)
            elif suffix == '.csv':
                df.to_csv(path, index=False, **kwargs)
            else:
                # 默认使用Excel
                df.to_excel(path, sheet_name=sheet_name, index=False)
        except Exception as e:
            raise DataFormatError(f"Failed to write file: {e}")
    
    def add_results(self, results: Dict[str, np.ndarray]):
        """
        将计算结果添加到数据中
        
        Parameters
        ----------
        results : dict
            结果字典，键为新列名，值为数据数组
        """
        if self._data is None:
            raise DataFormatError("No data loaded")
        
        for col_name, values in results.items():
            self._data[col_name] = values
    
    def get_dataframe(self) -> pd.DataFrame:
        """获取当前数据DataFrame"""
        if self._data is None:
            raise DataFormatError("No data loaded")
        return self._data.copy()
    
    def get_info(self) -> Dict[str, Any]:
        """
        获取数据信息
        
        Returns
        -------
        dict
            数据信息字典
        """
        if self._data is None:
            return {"status": "No data loaded"}
        
        return {
            "n_rows": len(self._data),
            "n_columns": len(self._data.columns),
            "columns": self._data.columns.tolist(),
            "detected_columns": self._detected_columns,
            "detected_isotopes": [
                elem for elem in ['u', 'n', 'c', 'mg']
                if self.has_isotope_data(elem)
            ]
        }
