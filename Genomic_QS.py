#!/usr/bin/env python3
"""
Genome Quality Assessment Tool
Version: 1.1.3
Author: FaZhang
Date: 2025-07-10
"""
import os
import sys
import argparse
import hashlib
import subprocess
import tempfile
import shutil
import csv
import glob
from collections import OrderedDict

def parse_arguments():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(
        description="Genome Quality Assessment with CheckM Integration",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input directory containing genome files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output CSV file path")
    parser.add_argument("-t", "--type", default="auto", choices=["auto", "flat", "nested"],
                        help="Directory structure type: auto/flat/nested")
    parser.add_argument("-c", "--checkm-data", 
                        default=os.getenv("CHECKM_DATA_PATH", "/opt/checkmData/"),
                        help="CheckM database path (默认: /opt/checkmData/)")
    parser.add_argument("-x", "--extension", default="fa",
                        help="Genome file extension for CheckM to recognize (default: fa)")
    parser.add_argument("-j", "--threads", type=int, default=os.cpu_count(),
                        help="Number of CPU threads for CheckM")
    parser.add_argument("-k", "--keep-temp", action="store_true",
                        help="Keep temporary files after execution")
    return parser.parse_args()

def detect_structure(input_dir, structure_type):
    """Detect directory structure (flat or nested)"""
    if structure_type != "auto":
        return structure_type
        
    subdirs = [d for d in os.listdir(input_dir) 
              if os.path.isdir(os.path.join(input_dir, d))]
    return "nested" if subdirs else "flat"

def calculate_md5(file_path):
    """Calculate MD5 checksum for a file (supports gzip)"""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def classify_quality(completeness, contamination):
    """Classify genome quality based on MIMAG standards"""
    qs = completeness - 5 * contamination
    if completeness >= 90 and contamination <= 5:
        return "near-complete", qs
    elif completeness >= 70 and contamination <= 10:
        return "high-quality", qs
    elif completeness >= 50 and contamination <= 10:
        return "medium-quality", qs
    else:
        return "low-quality", qs

def run_checkm(input_dir, output_dir, checkm_db, threads, extension):
    """Run CheckM analysis on genome collection"""
    # 通过环境变量设置数据库路径
    env = os.environ.copy()
    env["CHECKM_DATA_PATH"] = checkm_db
    
    lineage_cmd = [
        "checkm", "lineage_wf",
        "-x", extension,
        "-t", str(threads),
        "--pplacer_threads", str(max(1, threads//4)),
        input_dir,
        output_dir
    ]
    
    print("执行CheckM命令:")
    print(" ".join(lineage_cmd))
    print(f"使用数据库路径: {checkm_db}")
    print(f"使用文件扩展名: {extension}")
    
    try:
        process = subprocess.Popen(
            lineage_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )
        
        # 实时输出进度
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                print(output.strip())
                
        # 检查返回码
        return_code = process.wait()
        if return_code != 0:
            error_output = process.stderr.read()
            print(f"CheckM错误输出:\n{error_output}")
            raise subprocess.CalledProcessError(return_code, lineage_cmd, output=error_output)
            
    except subprocess.CalledProcessError as e:
        print(f"CheckM命令执行失败: {e}")
        print(f"目录内容: {os.listdir(input_dir)}")
        # 检查文件权限
        for file in os.listdir(input_dir):
            file_path = os.path.join(input_dir, file)
            print(f"{file}: {'可读' if os.access(file_path, os.R_OK) else '不可读'}")
        raise
    
    # 质量评估
    qa_cmd = [
        "checkm", "qa",
        os.path.join(output_dir, "lineage.ms"),
        output_dir,
        "-o", "2",
        "--tab_table"
    ]

    result_file = os.path.join(output_dir, "checkm_results.tsv")
    with open(result_file, "w") as fh:
        subprocess.run(qa_cmd, stdout=fh, check=True, env=env)
    
    return result_file

def parse_checkm_results(result_file):
    """Parse CheckM results into structured data"""
    results = {}
    with open(result_file, "r") as f:
        lines = f.readlines()
    
    # Find the start of the tabular data
    start_index = 0
    for i, line in enumerate(lines):
        if line.startswith("Bin Id"):
            start_index = i
            break
    
    if start_index == 0:
        raise ValueError("CheckM results format error: No header found")
    
    # Parse header and data
    header = lines[start_index].strip().split("\t")
    completeness_idx = header.index("Completeness")
    contamination_idx = header.index("Contamination")
    
    for line in lines[start_index+1:]:
        if not line.strip():
            continue
        fields = line.strip().split("\t")
        bin_id = fields[0]
        try:
            completeness = float(fields[completeness_idx].rstrip('%'))
            contamination = float(fields[contamination_idx].rstrip('%'))
            results[bin_id] = (completeness, contamination)
        except (ValueError, IndexError):
            continue
    
    return results

def process_genomes(args):
    """Main processing workflow"""
    # 验证数据库路径
    print(f"使用CheckM数据库: {args.checkm_data}")
    if not os.path.exists(args.checkm_data):
        raise FileNotFoundError(f"数据库路径不存在: {args.checkm_data}")
    
    # 检测目录结构
    dir_structure = detect_structure(args.input, args.type)
    
    # 创建临时工作空间
    temp_dir = tempfile.mkdtemp(prefix="genomic_qs_")
    genome_dir = os.path.join(temp_dir, "genomes")
    os.makedirs(genome_dir, exist_ok=True)
    
    # 收集基因组文件 - 使用文件复制代替符号链接
    print(f"使用文件扩展名: {args.extension}")
    genomes = []
    file_counter = 1
    found_files = 0
    
    if dir_structure == "flat":
        pattern = os.path.join(args.input, f"*.{args.extension}")
        files = glob.glob(pattern)
        print(f"找到 {len(files)} 个匹配 {args.extension} 扩展名的文件")
        
        for file_path in files:
            if not os.path.isfile(file_path):
                continue
            file_name = os.path.basename(file_path)
            bin_id = f"GENOME_{file_counter:04d}"
            temp_path = os.path.join(genome_dir, file_name)
            
            # 复制文件而不是符号链接
            try:
                shutil.copy2(file_path, temp_path)
                print(f"复制: {file_path} → {temp_path}")
                # 验证文件大小
                src_size = os.path.getsize(file_path)
                dest_size = os.path.getsize(temp_path)
                if src_size != dest_size:
                    print(f"警告: 文件大小不匹配! 源: {src_size}字节, 目标: {dest_size}字节")
            except Exception as e:
                print(f"复制错误: {e}")
                continue
            
            genomes.append({
                "src_path": file_path,
                "temp_path": temp_path,
                "bin_id": bin_id,
                "fasta_name": file_name,
                "md5": calculate_md5(file_path)
            })
            file_counter += 1
            found_files += 1
    
    else:  # 嵌套结构
        for sample_dir in os.listdir(args.input):
            sample_path = os.path.join(args.input, sample_dir)
            if not os.path.isdir(sample_path):
                continue
            
            pattern = os.path.join(sample_path, f"*.{args.extension}")
            files = glob.glob(pattern)
            print(f"样本 '{sample_dir}' 找到 {len(files)} 个匹配 {args.extension} 扩展名的文件")
            
            for file_path in files:
                if not os.path.isfile(file_path):
                    continue
                file_name = f"{sample_dir}_{os.path.basename(file_path)}"
                bin_id = f"GENOME_{file_counter:04d}"
                temp_path = os.path.join(genome_dir, file_name)
                
                # 复制文件而不是符号链接
                try:
                    shutil.copy2(file_path, temp_path)
                    print(f"复制: {file_path} → {temp_path}")
                    # 验证文件大小
                    src_size = os.path.getsize(file_path)
                    dest_size = os.path.getsize(temp_path)
                    if src_size != dest_size:
                        print(f"警告: 文件大小不匹配! 源: {src_size}字节, 目标: {dest_size}字节")
                except Exception as e:
                    print(f"复制错误: {e}")
                    continue
                
                genomes.append({
                    "src_path": file_path,
                    "temp_path": temp_path,
                    "bin_id": bin_id,
                    "fasta_name": file_name,
                    "md5": calculate_md5(file_path)
                })
                file_counter += 1
                found_files += 1
    
    if found_files == 0:
        raise FileNotFoundError(
            f"在 {args.input} 中未找到基因组文件\n"
            f"使用的扩展名: {args.extension}\n"
            "可能的解决方案:\n"
            "1. 使用正确的文件扩展名 (通过 -x 参数指定)\n"
            "2. 确认文件格式正确 (FASTA格式)"
        )
    
    print(f"共复制 {len(genomes)} 个基因组文件")
    
    # 运行CheckM分析 - 使用命令行指定的扩展名
    checkm_output = os.path.join(temp_dir, "checkm_output")
    os.makedirs(checkm_output, exist_ok=True)
    result_file = run_checkm(genome_dir, checkm_output, args.checkm_data, args.threads, args.extension)
    
    # 解析CheckM结果
    checkm_results = parse_checkm_results(result_file)
    
    # 准备最终结果
    final_results = []
    for genome in genomes:
        # 从文件名中去掉扩展名作为bin_id
        base_name = os.path.splitext(os.path.basename(genome["temp_path"]))[0]
        
        if base_name in checkm_results:
            completeness, contamination = checkm_results[base_name]
            quality_class, qs = classify_quality(completeness, contamination)
        else:
            completeness, contamination, quality_class, qs = "NA", "NA", "NA", "NA"
        
        final_results.append(OrderedDict([
            ("fasta_file_name", genome["fasta_name"]),
            ("fasta_file_md5", genome["md5"]),
            ("completeness(%)", completeness),
            ("contamination(%)", contamination),
            ("QS", qs),
            ("quality_class", quality_class)
        ]))

    output_dir = os.path.dirname(args.output) or '.'  # 空路径默认为当前目录
    os.makedirs(output_dir, exist_ok=True)  # 安全创建目录
    
    # 添加路径验证
    if not os.access(output_dir, os.W_OK):
        raise PermissionError(f"输出目录不可写: {output_dir}")
    
    print(f"最终输出路径: {os.path.abspath(args.output)}")
    
    # 写入文件
    with open(args.output, "w") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=final_results[0].keys())
        writer.writeheader()
        writer.writerows(final_results)
    
    # 清理临时文件
    if not args.keep_temp:
        shutil.rmtree(temp_dir)
    
    return final_results

def main():
    try:
        args = parse_arguments()
        results = process_genomes(args)
        print(f"成功处理 {len(results)} 个基因组")
        print(f"结果保存至: {args.output}")
    except Exception as e:
        print(f"错误: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()