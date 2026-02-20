"""
数据I/O工具模块

提供同位素数据的读写、批量处理和结果导出功能。

支持格式:
- Excel (.xlsx, .xls)
- CSV (.csv)
- 未来可扩展: JSON, NetCDF等

示例用法:
    from toolkit.io import DataHandler, BatchProcessor
    
    # 读取数据
    handler = DataHandler('data/input.xlsx')
    data = handler.read()
    
    # 批量处理
    processor = BatchProcessor('u')
    results = processor.process_file('data/input.xlsx')
    
    # 保存结果
    handler.write(results, 'output/results.xlsx')
"""

from toolkit.io.data_handler import DataHandler, DataFormatError
from toolkit.io.batch_processor import BatchProcessor, ProcessingResult

__all__ = [
    'DataHandler',
    'DataFormatError',
    'BatchProcessor',
    'ProcessingResult',
]
