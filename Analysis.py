#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
基因组分析主调度程序
使用方式: python Analysis.py <工具名> [工具参数]
"""

import argparse
import sys
from tools import CompletenessAndContamination, QSDistribution, SpeciesAbundance, ComparisonOfKeyParameters

def main():
    parser = argparse.ArgumentParser(
        description='河海微生物基因组MAGs分析工具集',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    subparsers = parser.add_subparsers(dest='tool', title='可用工具', help='选择分析工具')
    
    # 1. 完整性与污染分析工具
    cple_parser = subparsers.add_parser(
        'cple-ctam', 
        help='分析中高质量MAGs的完整性和污染分布',
        description='分析中高质量MAGs的完整性和污染分布\n使用示例: python Analysis.py cple-ctam -c db.ini -t HeHaiQSData -p out.png'
    )
    cple_parser.add_argument('-c', '--config', required=True, help='数据库配置文件路径 (db_config.ini)')
    cple_parser.add_argument('-t', '--table', required=True, help='需要分析的表名 (例如: HeHaiQSData)')
    cple_parser.add_argument('-p', '--plot', required=True, help='分布图保存路径 (例如: output.png 或 output.tiff)')
    cple_parser.add_argument('--dpi', type=int, default=300, help='图像分辨率DPI (默认: 300)')
    cple_parser.set_defaults(func=CompletenessAndContamination.run)
    
    # 2. 质量分布分析工具
    qs_parser = subparsers.add_parser(
        'qs-dist', 
        help='质量分类分布分析',
        description='质量分类分布分析\n使用示例: python Analysis.py qs-dist -c db.ini -t HeHaiQSData -p pie.png'
    )
    qs_parser.add_argument('-c', '--config', required=True, help='数据库配置文件路径 (db_config.ini)')
    qs_parser.add_argument('-t', '--table', required=True, help='需要分析的表名 (例如: HeHaiQSData)')
    qs_parser.add_argument('-p', '--plot', required=True, help='分布图保存路径 (例如: pie.png 或 pie.tiff)')
    qs_parser.add_argument('--dpi', type=int, default=300, help='图像分辨率DPI (默认: 300)')
    qs_parser.set_defaults(func=QSDistribution.run)
    
    # 3. 物种丰度分析工具
    tax_parser = subparsers.add_parser(
        'tax-abundance', 
        help='物种分类丰度分析',
        description='物种分类丰度分析(界，门，纲，目，科，属，种)\n使用示例: python Analysis.py tax-abundance -c db.ini -t table1,table2 -p abundance.png'
    )
    tax_parser.add_argument('-c', '--config', required=True, help='数据库配置文件路径 (db_config.ini)')
    tax_parser.add_argument('-t', '--tables', required=True, help='需要分析的表名，多个表用逗号分隔 (例如: table1,table2)')
    tax_parser.add_argument('-p', '--plot', required=True, help='丰度图保存路径 (例如: abundance.png 或 abundance.tiff)')
    tax_parser.add_argument('--dpi', type=int, default=300, help='图像分辨率DPI (默认: 300)')
    tax_parser.set_defaults(func=SpeciesAbundance.run)
    
    # 4. 关键参数比较工具
    param_parser = subparsers.add_parser(
        'key-params', 
        help='基因组关键参数比较',
        description='基因组组装关键参数比较\n使用示例: python Analysis.py key-params -c db.ini -t QS表.外键:参数表.外键 -o comparison.png'
    )
    param_parser.add_argument('-c', '--config', required=True, help='数据库配置文件路径')
    param_parser.add_argument('-t', '--tables', required=True, 
                             help='表对列表，格式: QS表.外键:关键参数表.外键 (多个用逗号分隔)')
    param_parser.add_argument('-p', '--plot', required=True, 
                             help='输出图像路径 (如: output.png, output.tiff)')
    param_parser.add_argument('--dpi', type=int, default=300, 
                             help='图像DPI (默认: 300)')
    param_parser.add_argument('--gomc-sql', default='./GOMCSQL/GOMCTables.sql', 
                             help='基准比较对象GOMC表SQL文件路径 (默认: ./GOMCSQL/GOMCTables.sql)')
    param_parser.set_defaults(func=ComparisonOfKeyParameters.run)
    
    # 如果没有参数，显示帮助信息
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()