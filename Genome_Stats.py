#!/usr/bin/env python3
import os
import sys
import gzip
import glob
import argparse
import hashlib
import shutil  # 主要使用shutil进行文件复制
import csv
import random
import subprocess
import tempfile
from datetime import datetime
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

def main():
    # 参数解析扩展
    parser = argparse.ArgumentParser(description="基因组统计与GTDB-Tk并行分类工具")
    parser.add_argument("-i", "--input-dir", required=True, help="输入目录路径")
    parser.add_argument("-o", "--output-file", default=f"Genome_Stats_{datetime.now().strftime('%Y%m%d')}.csv", help="输出CSV路径")
    parser.add_argument("-j", "--threads", type=int, default=1, help="并行任务数(CPU核数)")
    parser.add_argument("-x", "--extension", default="fa", help="文件扩展名(默认:fa)")
    parser.add_argument("-k", "--keep-temp", default="False", help="保留临时文件(默认:False)")
    parser.add_argument("-r", "--temp-dir", default="/tmp", help="临时目录路径(默认:/tmp)")
    parser.add_argument("-t", "--file-type", choices=["auto", "flat", "nested"], default="auto", help="文件组织结构")
    args = parser.parse_args()

    # 执行核心流程
    process_files(args)

def run_gtdbtk(task_dir, task_out, file_ext):
    """执行GTDB-Tk分类任务"""
    try:
        cmd = f"gtdbtk classify_wf --genome_dir {task_dir} --out_dir {task_out} -x {file_ext} --cpus 1 --skip_ani_screen"
        subprocess.run(cmd, shell=True, check=True)
        return True, ""
    except subprocess.CalledProcessError as e:
        return False, f"GTDB-Tk失败: {str(e)}"

def parse_gtdb_summary(summary_path, args):
    """解析GTDB分类结果文件"""
    classifications = {}
    if not os.path.exists(summary_path):
        return classifications
        
    with open(summary_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            genome_id = row['user_genome'] + "." + args.extension
            classifications[genome_id] = row['classification']
    return classifications

def find_fasta_files(input_dir, file_type):
    """查找所有FASTA文件（新增函数）"""
    patterns = ["*.fa", "*.fa.gz"]
    all_files = []
    use_type = file_type
    
    # 类型检测
    if use_type == "auto":
        use_type = detect_file_type(input_dir)
        print(f"检测到文件结构类型: {use_type} (自动检测模式)")
    elif use_type == "nested":
        print(f"文件结构类型: {use_type}")
    elif use_type == "flat":
        print(f"文件结构类型: {use_type}")
    
    # 文件搜索
    for pattern in patterns:
        if use_type == "nested":
            full_pattern = os.path.join(input_dir, "**", pattern)
            found_files = glob.glob(full_pattern, recursive=True)
        elif use_type == "flat":
            found_files = glob.glob(os.path.join(input_dir, pattern))
        all_files.extend([f for f in found_files if os.path.isfile(f)])
    
    return all_files, use_type

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

def calculate_file_stats(file_path):
    """封装函数：计算文件统计信息"""
    md5_val = calculate_md5(file_path)
    
    if file_path.endswith('.gz'):
        with gzip.open(file_path, 'rt') as fh:
            stats = calc_stats(fh)
    else:
        with open(file_path, 'r') as fh:
            stats = calc_stats(fh)
    
    return (md5_val,) + stats

def process_files(args):
    """处理主流程"""
    # 文件搜索与类型检测
    all_files, use_type = find_fasta_files(args.input_dir, args.file_type)
    if int(args.threads) > len(all_files):
        print(f"实际文件数{len(all_files)}小于线程数{args.threads}")
        args.threads = len(all_files)
        print(f"调整线程数为实际文件数{args.threads}")
    
    # 创建临时根目录
    temp_root = f"{args.temp_dir}/gtdbtk_parallel_{random.randint(1000,9999)}"
    os.makedirs(temp_root, exist_ok=True)
    
    # 分桶创建任务目录
    task_dirs = []
    for i in range(args.threads):
        task_dir = os.path.join(temp_root, f"task_{i}")
        os.makedirs(task_dir, exist_ok=True)
        task_dirs.append(task_dir)
    
    # 文件分配到任务目录（使用文件复制而非软链接）
    for idx, file_path in enumerate(all_files):
        task_idx = idx % args.threads
        if use_type == "nested":
            # 文件路径+文件名
            dest_name = os.path.basename(os.path.dirname(file_path)) + "_" +os.path.basename(file_path)
        else:
            dest_name = os.path.basename(file_path)
        dest_path = os.path.join(task_dirs[task_idx], dest_name)
        
        try:
            # 使用shutil.copy2进行文件复制（保留元数据）
            shutil.copy2(os.path.abspath(file_path), dest_path)
        except Exception as e:
            print(f"错误：复制文件失败 {file_path} -> {dest_path}: {str(e)}")
            sys.exit(1)
    
    # 并行执行GTDB-Tk
    task_outputs = [os.path.join(temp_root, f"output_{i}") for i in range(args.threads)]
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for task_dir, task_out in zip(task_dirs, task_outputs):
            os.makedirs(task_out, exist_ok=True)
            futures.append(executor.submit(run_gtdbtk, task_dir, task_out, args.extension))
        
        # 检查任务状态
        for future in futures:
            success, error = future.result()
            if not success:
                print(f"警告: {error}")

    # 合并分类结果
    all_classifications = {}
    for task_out in task_outputs:
        bac_file = os.path.join(task_out, "gtdbtk.bac120.summary.tsv")
        arc_file = os.path.join(task_out, "gtdbtk.ar53.summary.tsv")
        
        if os.path.exists(bac_file):
            all_classifications.update(parse_gtdb_summary(bac_file, args))
        if os.path.exists(arc_file):
            all_classifications.update(parse_gtdb_summary(bac_file, args))
    
    # 生成最终CSV (扩展分类列)
    generate_final_csv(all_files, args.output_file, all_classifications, args.file_type)
    
    # 清理临时文件
    if args.keep_temp != False :
        print(f"清理临时目录: {temp_root}")
        shutil.rmtree(temp_root)

def generate_final_csv(file_list, output_path, classifications, file_type):
    """生成带分类结果的CSV文件"""
    with open(output_path, 'w') as f:
        f.write("fasta_file_name,fasta_file_md5,total_size(bp),sequences,largest_seq(bp),smallest_seq(bp),N50(bp),L50,classification\n")
        
        for file_path in file_list:
            # 文件名处理
            if file_type == "nested":
                dir_name = os.path.basename(os.path.dirname(file_path))
                file_name = f"{dir_name}_{os.path.basename(file_path)}"
            else:
                file_name = os.path.basename(file_path)
            
            # 获取分类结果
            classification = classifications.get(file_name, "Not Found")
            
            try:
                # 计算文件统计信息
                stats = calculate_file_stats(file_path)
                # 写入CSV行
                f.write(f"{file_name},{','.join(map(str, stats))},{classification}\n")
            except Exception as e:
                print(f"处理文件 {file_path} 时出错：{str(e)}")
                # 写入错误占位符
                error_stats = (f"ERROR-{str(e)[:30]}", 0, 0, 0, 0, 0, 0)
                f.write(f"{file_name},{','.join(map(str, error_stats))},ERROR\n")

if __name__ == "__main__":
    main()
