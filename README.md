# 基因组组装统计工具

## 概述
该脚本用于计算FASTA文件(.fa/.fa.gz)的基因组统计信息，包括文件大小、序列数量、N50等关键指标。支持两种文件组织方式并自动检测目录结构。

## 安装依赖
```bash
sudo apt install gawk # 必需依赖
```

## 使用方式

```bash
./genome_GOMC_stats.sh [-i INPUT_DIR] [-o OUTPUT_FILE] [-t FILE_TYPE] [-h]
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
./genome_GOMC_stats.sh
```

**2. 指定输入输出**（自动检测目录结构）：

```bash
./genome_GOMC_stats.sh -i /data/genomes -o results.csv
```

**3. 强制扁平结构**：

```bash
./genome_GOMC_stats.sh -t flat
```

**4. 强制嵌套结构**：

```bash
./genome_GOMC_stats.sh -t nested
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

