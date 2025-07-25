#!/usr/bin/env python3
import os
import sys
import gzip
import glob
import argparse
import hashlib
from datetime import datetime
from collections import defaultdict

def main():
    # 设置默认值
    default_input_dir = "../Bins"
    default_output_file = f"Genome_Statistics_{datetime.now().strftime('%Y%m%d')}.csv"
    
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="基因组组装统计工具")
    parser.add_argument("-i", "--input-dir", default=default_input_dir,
                        help="输入目录路径（默认：../Bins）")
    parser.add_argument("-o", "--output-file", default=default_output_file,
                        help="输出CSV文件路径（默认：带日期的文件名）")
    parser.add_argument("-t", "--file-type", choices=["auto", "flat", "nested"], default="auto",
                        help="文件组织结构：'flat'(扁平), 'nested'(嵌套), 'auto'(自动检测)")
    args = parser.parse_args()

    # 验证输入目录
    if not os.path.isdir(args.input_dir):
        sys.exit(f"错误：输入目录不存在 - {args.input_dir}")
    
    # 创建输出目录（如果需要）
    output_dir = os.path.dirname(args.output_file) or "."
    os.makedirs(output_dir, exist_ok=True)
    
    # 处理文件并生成统计结果
    process_files(args.input_dir, args.output_file, args.file_type)
    print(f"处理完成! 结果已保存至：{args.output_file}")

def detect_file_type(input_dir):
    """自动检测文件组织结构（扁平或嵌套）"""
    for entry in os.scandir(input_dir):
        if entry.is_dir() and not entry.name.startswith('.'):
            return "nested"
    return "flat"

def calculate_md5(file_path):
    """计算文件的MD5校验值"""
    hash_md5 = hashlib.md5()
    try:
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    except Exception as e:
        return f"MD5_ERROR: {str(e)}"

def calculate_n50_l50(lengths, total_bp):
    """计算N50和L50指标"""
    if not lengths:
        return 0, 0
    
    sorted_lengths = sorted(lengths, reverse=True)
    half_total = total_bp / 2
    cumulative = 0
    for i, seq_len in enumerate(sorted_lengths, 1):
        cumulative += seq_len
        if cumulative >= half_total:
            return seq_len, i
    return sorted_lengths[0], 1

def calc_stats(file_handle):
    """计算FASTA文件的组装指标"""
    total_bp = 0
    seq_lengths = []
    current_length = 0
    seq_count = 0

    for line in file_handle:
        if isinstance(line, bytes):
            line = line.decode('utf-8', errors='ignore')
        line = line.strip()
        
        if line.startswith('>'):
            if current_length > 0:
                seq_lengths.append(current_length)
                total_bp += current_length
                seq_count += 1
            current_length = 0
        elif line:  # 忽略空行
            current_length += len(line)
    
    if current_length > 0:
        seq_lengths.append(current_length)
        total_bp += current_length
        seq_count += 1
    
    # 空文件处理
    if seq_count == 0:
        return (0, 0, 0, 0, 0, 0)
    
    max_len = max(seq_lengths) if seq_lengths else 0
    min_len = min(seq_lengths) if seq_lengths else 0
    n50, l50 = calculate_n50_l50(seq_lengths, total_bp)
    
    return (total_bp, seq_count, max_len, min_len, n50, l50)

def process_files(input_dir, output_file, file_type):
    """处理目录中的所有FASTA文件"""
    patterns = ["*.fa", "*.fa.gz"]
    all_files = []
    use_type = file_type  # 中间变量
    
    # 类型检测（只执行一次）▼
    if use_type == "auto":
        use_type = detect_file_type(input_dir)
        print(f"检测到文件结构类型: {use_type} (自动检测模式)")
    elif use_type == "nested":
        print(f"文件结构类型: {use_type}")
    elif use_type == "flat":
        print(f"文件结构类型: {use_type}")
    
    # 文件搜索（无打印语句）▼
    for pattern in patterns:
        if use_type == "nested":
            full_pattern = os.path.join(input_dir, "**", pattern)
            found_files = glob.glob(full_pattern, recursive=True)
        elif use_type == "flat":
            found_files = glob.glob(os.path.join(input_dir, pattern))
        all_files.extend([f for f in found_files if os.path.isfile(f)])

    print(f"找到 {len(all_files)} 个FASTA文件：")
    for f in all_files[:3]:
        if use_type == "nested":
            dir_name = os.path.basename(os.path.dirname(f))
            print(f" - {dir_name}_{os.path.basename(f)}")
        else:
            print(f" - {os.path.basename(f)}")
    
    if not all_files:
        print(f"警告：在 {input_dir} 中未找到任何FASTA文件")
        return

    with open(output_file, 'w') as out_f:
        out_f.write("fasta_file_name,fasta_file_md5,total_size(bp),sequences,largest_seq(bp),smallest_seq(bp),N50(bp),L50\n")
        
        for file_path in all_files:
            if use_type == "nested":  # 使用检测后的类型
                dir_name = os.path.basename(os.path.dirname(file_path))
                file_name = f"{dir_name}_{os.path.basename(file_path)}"
            else:
                file_name = os.path.basename(file_path)
            
            try:
                md5_val = calculate_md5(file_path)
                
                if file_path.endswith('.gz'):
                    with gzip.open(file_path, 'rt') as fh:
                        stats = calc_stats(fh)
                else:
                    with open(file_path, 'r') as fh:
                        stats = calc_stats(fh)
                
                out_f.write(f"{file_name},{md5_val},{','.join(map(str, stats))}\n")
                
            except Exception as e:
                print(f"处理文件 {file_path} 时出错：{str(e)}")
                error_stats = (0, 0, 0, 0, 0, 0)
                out_f.write(f"{file_name},ERROR,{','.join(map(str, error_stats))}\n")

if __name__ == "__main__":
    main()
