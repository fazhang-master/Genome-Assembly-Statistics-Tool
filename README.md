# 基因组组装统计工具

## 概述
该脚本用于计算FASTA文件(.fa/.fa.gz)的基因组统计信息，包括文件大小、序列数量、N50等关键指标。支持两种文件组织方式并自动检测目录结构。

## 安装依赖
```bash
sudo apt install gawk # 必需依赖
```

## 使用方式

```bash
./Genome_Stats.sh [-i INPUT_DIR] [-o OUTPUT_FILE] [-t FILE_TYPE] [-h]
```

### 参数说明

| 参数 |                             描述                             |              默认值              |
| :--: | :----------------------------------------------------------: | :------------------------------: |
| `-i` |                         输入目录路径                         |            `../Bins`             |
| `-o` |                       输出CSV文件路径                        | `Genome_Statistics_YYYYMMDD.csv` |
| `-t` | 文件组织类型：`flat`(扁平)、`nested`(嵌套)或`auto`(自动检测) |              `auto`              |
| `-h` |                         显示帮助信息                         |                -                 |

### 文件目录格式

**1. 扁平结构**：

```bash
input_dir/
├── file1.fa
├── file2.fa.gz
└── file3.fa
```

- 保存文件名时直接使用原文件名，例如：`file1.fa`

**2. 嵌套结构**：

input_dir/

```bash
└── sample_folder/
    ├── bin1.fa
    ├── bin2.fa.gz
    └── bin3.fa
```

- 文件名保存格式：`sample_folder_bin1.fa`

### 输出字段说明

|        字段        |     描述      |
| :----------------: | :-----------: |
| `fasta_file_name`  |  FASTA文件名  |
|  `fasta_file_md5`  | 文件MD5校验值 |
|  `total_size(bp)`  |   总碱基数    |
|    `sequences`     |   序列总数    |
| `largest_seq(bp)`  | 最长序列长度  |
| `smallest_seq(bp)` | 最短序列长度  |
|     `N50(bp)`      |     N50值     |
|       `L50`        |     L50值     |

### 使用示例

**1. 基本用法**：

```bash
./Genome_Stats.sh
```

**2. 指定输入输出**（自动检测目录结构）：

```bash
./Genome_Stats.sh -i /data/genomes -o results.csv
```

**3. 强制扁平结构**：

```bash
./Genome_Stats.sh -t flat
```

**4. 强制嵌套结构**：

```bash
./Genome_Stats.sh -t nested
```

## 注意事项

1. 当使用`-t auto`(默认)时，脚本自动检测目录结构：
   - 如果输入目录包含子目录 → 按嵌套结构处理
   - 如果输入目录只有文件 → 按扁平结构处理
2. 支持递归处理所有子目录中的`.fa`和`.fa.gz`文件

## 脚本输出效果：

```markdown
fasta_file_name,fasta_file_md5,total_size(bp),sequences,largest_seq(bp),smallest_seq(bp),N50(bp),L50
sample_folder_bin.10.fa,4262639e5730dd421297675550c8e174,5002757,873,42836,1504,6898,234
direct_file.fa,b2960bdf4ec93e4089be887dbeffaac4,2184416,500,36272,1503,5386,123
```

# 基因组质量检测工具

## 概述

该工具自动化执行以下流程：
1. 扫描输入目录中的基因组文件（支持FASTA格式，包括压缩文件）
2. 计算原始文件的MD5校验值
3. 使用CheckM评估基因组的完整度(completeness)和污染度(contamination)
4. 计算质量分数(QS = completeness - 5×contamination)
5. 根据MIMAG标准对基因组质量进行分类
6. 生成包含所有统计指标的CSV报告

## 功能特点
- 支持扁平结构和嵌套结构目录
- 自动检测目录结构类型（auto模式）
- 保留原始文件完整性（通过文件复制处理）
- 兼容压缩基因组文件（.fa.gz, .fasta.gz）
- 生成符合MIMAG标准的分类结果

## 安装要求

### 系统依赖
- Python 3.9
- CheckM v1.2.3 ([安装指南](https://anaconda.org/bioconda/checkm-genome))
- CheckM数据库 ([数据下载](https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz))
- hmmer 3.4 ([安装指南](https://anaconda.org/bioconda/hmmer))
- pplacer 1.1.alpha19 ([安装指南](https://anaconda.org/bioconda/pplacer))  

### 安装CheckM
```bash
# 创建conda环境
conda create -n GenomicQS python=3.9
conda activate GenomicQS

# 安装CheckM及其依赖
conda install -c bioconda checkm-genome==1.2.3
conda install -c bioconda hmmer==3.4
conda install -c bioconda pplacer==1.1.alpha19
```

### 下载CheckM数据库

```bash
# 下载数据库
cd /opt
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -zxvf checkm_data_2015_01_16.tar.gz
mv checkm_data_2015_01_16 checkmData && cd checkmData

# 设置环境变量
export CHECKM_DATA_PATH=$(pwd)
```

## 使用说明

```python
python Genomic_QS.py -i <input_dir> -o <output.csv> [-t STRUCTURE_TYPE]
```

### 参数说明

|         参数          |                  描述                  |      默认值       |
| :-------------------: | :------------------------------------: | :---------------: |
|    `-i`, `--input`    |              输入目录路径              |     **必需**      |
|   `-o`, `--output`    |            输出CSV文件路径             |     **必需**      |
|    `-t`, `--type`     | 目录结构类型：`auto`, `flat`, `nested` |      `auto`       |
| `-c`, `--checkm-data` |            CheckM数据库路径            | `/opt/checkmData` |
|  `-x`, `--extension`  |            基因组文件扩展名            |       `fa`        |
|   `-j`, `--threads`   |         CheckM使用的CPU线程数          |   系统CPU核心数   |
|  `-k`, `--keep-temp`  |              保留临时文件              |       False       |
|    `-h`, `--help`     |              显示帮助信息              |                   |

### 目录结构示例

见《基因组组装统计工具》文件目录格式章节

### 输出文件格式

CSV文件包含以下列：

- `fasta_file_name`: 文件名（嵌套结构为`样本名_文件名`）
- `fasta_file_md5`: 原始文件的MD5校验值
- `completeness(%)`: 基因组完整度百分比
- `contamination(%)`: 基因组污染度百分比
- `QS`: 质量分数 (completeness - 5×contamination)
- `quality_class`: 质量分类（基于MIMAG标准）

### 质量分类标准

|           分类            |            标准            |
| :-----------------------: | :------------------------: |
| 接近完整 (near-complete)  | 完整度 ≥90% 且 污染度 ≤5%  |
|   高质量 (high-quality)   | 完整度 ≥70% 且 污染度 ≤10% |
| 中等质量 (medium-quality) | 完整度 ≥50% 且 污染度 ≤10% |
|   低质量 (low-quality)    |          其他情况          |

### 使用示例

**扁平结构处理**

```python
python Genomic_QS.py -i flat_genomes/ -o flat_results.csv -t flat
```

**嵌套结构处理**

```python
python Genomic_QS.py -i nested_samples/ -o nested_results.csv -t nested
```

**自动检测结构**

```python
python Genomic_QS.py -i mixed_samples/ -o results.csv
```

**指定数据库**

```python
python Genomic_QS.py -i mixed_samples/ -o results.csv -c /opt/checkmData
```

## 注意事项

1. 首次运行前需设置`CHECKM_DATA_PATH`环境变量
2. 处理大量基因组时建议使用高线程数（如`-j 32`）
3. 临时文件默认自动删除（保留请使用`-k`参数）
4. 输出CSV可直接导入Excel进行进一步分析

## 错误处理

- **未找到CheckM数据库**：通过`-c`参数指定或设置环境变量
- **无基因组文件**：检查目录结构和文件扩展名
- **CheckM执行失败**：检查日志中的具体错误信息

### 使用示例输出

```csv
fasta_file_name,fasta_file_md5,completeness(%),contamination(%),QS,quality_class
Sample1_bin1.fa,4262639e5730dd421297675550c8e174,97.21,1.67,88.86,near-complete
Sample1_bin2.fa,b2960bdf4ec93e4089be887dbeffaac4,92.96,2.85,78.71,high-quality
Sample2_bin1.fa,cfc07118448e1e1bcc6f475a19bb9fb3,16.61,1.72,7.81,low-quality
