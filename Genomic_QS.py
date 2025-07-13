#!/usr/bin/env python3
"""
Genome Quality Assessment Tool
Version: 1.1.5
Author: FaZhang
Date: 2025-07-13
修复嵌套结构文件名前缀问题
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
import pandas as pd

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
    parser.add_argument("-k", "--keep-temp", default="False", action="store_true",
                        help="Keep temporary files after execution")
    parser.add_argument("-n", "--batch-size", type=int, default=0,
                        help="Number of files to process per batch (0=process all at once)")
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
    
    # 确保HMMER可执行文件在PATH中
    env["PATH"] = os.path.dirname(shutil.which("hmmsearch")) + ":" + env["PATH"]
    
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

def process_batch(genomes_batch, batch_idx, temp_dir, args, dir_structure):
    """Process a batch of genomes"""
    # 创建当前批次的临时目录
    batch_genome_dir = os.path.join(temp_dir, f"genomes_batch_{batch_idx}")
    os.makedirs(batch_genome_dir, exist_ok=True)
    
    batch_results = []
    for genome in genomes_batch:
        # 关键修改：根据目录结构生成新文件名
        if dir_structure == "nested":
            # 使用文件夹名_文件名格式
            new_file_name = f"{genome['parent_dir']}_{os.path.basename(genome['src_path'])}"
        else:
            new_file_name = os.path.basename(genome['src_path'])
        
        temp_path = os.path.join(batch_genome_dir, new_file_name)
        try:
            shutil.copy2(genome["src_path"], temp_path)
            print(f"复制批次 {batch_idx}: {genome['src_path']} → {temp_path}")
        except Exception as e:
            print(f"复制错误: {e}")
            continue
        
        # 保存带前缀的文件名用于后续结果关联
        genome["prefixed_name"] = new_file_name
    
    # 运行CheckM分析
    checkm_output = os.path.join(temp_dir, f"checkm_output_batch_{batch_idx}")
    os.makedirs(checkm_output, exist_ok=True)
    result_file = run_checkm(batch_genome_dir, checkm_output, args.checkm_data, args.threads, args.extension)
    
    # 解析CheckM结果
    checkm_results = parse_checkm_results(result_file)
    
    # 准备当前批次结果
    for genome in genomes_batch:
        # 使用带前缀的文件名（去掉扩展名）作为匹配键
        base_name = os.path.splitext(genome["prefixed_name"])[0]
        
        if base_name in checkm_results:
            completeness, contamination = checkm_results[base_name]
            quality_class, qs = classify_quality(completeness, contamination)
        else:
            completeness, contamination, quality_class, qs = "NA", "NA", "NA", "NA"
        
        batch_results.append(OrderedDict([
            ("batch", batch_idx),
            ("fasta_file_name", genome["prefixed_name"]),  # 使用带前缀的文件名
            ("fasta_file_md5", genome["md5"]),
            ("completeness(%)", completeness),
            ("contamination(%)", contamination),
            ("QS", qs),
            ("quality_class", quality_class)
        ]))
    
    # 保存当前批次结果到临时文件
    batch_output = os.path.join(temp_dir, f"batch_{batch_idx}_results.csv")
    with open(batch_output, "w") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=batch_results[0].keys())
        writer.writeheader()
        writer.writerows(batch_results)
    
    print(f"批次 {batch_idx} 完成，处理了 {len(genomes_batch)} 个基因组")
    return batch_output

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
    
    # 收集所有基因组文件
    print(f"使用文件扩展名: {args.extension}")
    all_genomes = []
    found_files = 0
    
    if dir_structure == "flat":
        pattern = os.path.join(args.input, f"*.{args.extension}")
        files = glob.glob(pattern)
        print(f"找到 {len(files)} 个匹配 {args.extension} 扩展名的文件")
        
        for file_path in files:
            if not os.path.isfile(file_path):
                continue
            all_genomes.append({
                "src_path": file_path,
                "md5": calculate_md5(file_path),
                "parent_dir": None  # 扁平结构没有父目录
            })
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
                all_genomes.append({
                    "src_path": file_path,
                    "md5": calculate_md5(file_path),
                    "parent_dir": sample_dir  # 记录父目录名
                })
                found_files += 1
    
    if found_files == 0:
        raise FileNotFoundError(
            f"在 {args.input} 中未找到基因组文件\n"
            f"使用的扩展名: {args.extension}\n"
            "可能的解决方案:\n"
            "1. 使用正确的文件扩展名 (通过 -x 参数指定)\n"
            "2. 确认文件格式正确 (FASTA格式)"
        )
    
    print(f"共找到 {len(all_genomes)} 个基因组文件")
    
    # 批次处理逻辑
    batch_results_files = []
    batch_size = args.batch_size if args.batch_size > 0 else len(all_genomes)
    
    if batch_size == len(all_genomes):
        print(f"一次性处理所有 {len(all_genomes)} 个文件")
        batch_output = process_batch(all_genomes, 0, temp_dir, args, dir_structure)
        batch_results_files.append(batch_output)
    else:
        print(f"将 {len(all_genomes)} 个文件分成 {len(all_genomes)//batch_size + 1} 批次，每批 {batch_size} 个")
        
        for i in range(0, len(all_genomes), batch_size):
            batch_idx = i // batch_size + 1
            batch_end = min(i + batch_size, len(all_genomes))
            genomes_batch = all_genomes[i:batch_end]
            
            print(f"\n=== 开始处理批次 {batch_idx} ({len(genomes_batch)} 个文件) ===")
            batch_output = process_batch(genomes_batch, batch_idx, temp_dir, args, dir_structure)
            batch_results_files.append(batch_output)
    
    # 合并所有批次结果
    final_results = []
    for batch_file in batch_results_files:
        try:
            # 使用pandas合并CSV文件更高效
            df = pd.read_csv(batch_file)
            final_results.append(df)
        except Exception as e:
            print(f"加载批次结果文件 {batch_file} 时出错: {e}")
    
    if final_results:
        # 使用pandas合并所有批次
        final_df = pd.concat(final_results, ignore_index=True)
        
        # 删除批次列（如果不需要）
        if "batch" in final_df.columns:
            final_df = final_df.drop(columns=["batch"])
        
        # 保存最终结果
        output_dir = os.path.dirname(args.output) or '.' 
        os.makedirs(output_dir, exist_ok=True)
        
        if not os.access(output_dir, os.W_OK):
            raise PermissionError(f"输出目录不可写: {output_dir}")
        
        # 使用pandas保存为CSV
        final_df.to_csv(args.output, index=False)
        print(f"合并结果保存至: {os.path.abspath(args.output)}")
    else:
        print("警告: 没有生成有效结果")
    
    # 清理临时文件
    if not args.keep_temp:
        shutil.rmtree(temp_dir)
    
    return final_df.to_dict("records") if not final_df.empty else []

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
