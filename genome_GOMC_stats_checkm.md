### 主要修改内容

1. **新增命令行选项**：

   - 添加了 `-c` 参数用于指定 CheckM 结果文件
   - 添加了 CheckM 文件存在性验证

2. **添加了 4 个新列**：

   ```csv
   completeness(%),contamination(%),QS,quality_class
   ```

   其中：

   - `completeness(%)`: 基因组完整性百分比
   - `contamination(%)`: 基因组污染率百分比
   - `QS`: 质量分数 (Completeness - 5 × Contamination)
   - `quality_class`: 质量分类 (near-complete/high-quality/medium-quality/low-quality)

3. **新增质量分类函数**：

- `classify_quality()` 根据 MIMAG 标准进行基因组质量分类
- 使用 `bc` 工具进行精确的浮点数计算

4. **标识符匹配优化**：

- 创建专门的 `generate_bin_id()` 函数处理标识符生成
- 完全移除文件扩展名（兼容 .fa, .fna, .gz 等格式）
- 统一嵌套结构和平坦结构的标识符生成逻辑

### CheckM 文件存放建议

#### 对于扁平结构：

```markdown
project_root/
├── Bins/                  # 基因组文件
│   ├── bin1.fa
│   ├── bin2.fa
│   └── bin3.fa
├── checkm_results/        # CheckM 结果目录
│   ├── lineage.ms
│   └── checkm_results.tsv  # 质量评估结果
└── genome_GOMC_stats.sh   # 统计脚本
```

运行命令：

```bash
# 生成CheckM结果
checkm lineage_wf -x fa ./Bins ./checkm_results
checkm qa ./checkm_results/lineage.ms ./checkm_results -o 2 --tab_table > ./checkm_results/checkm_results.tsv

# 运行统计脚本
./genome_GOMC_stats.sh -c ./checkm_results/checkm_results.tsv
```

#### 对于嵌套结构：

```markdown
project_root/
├── Samples/               # 样本目录
│   ├── Sample1/           # 样本1
│   │   ├── bin1.fa
│   │   └── bin2.fa
│   └── Sample2/           # 样本2
│       ├── bin1.fa
│       └── bin2.fa
├── checkm_results/        # CheckM 结果目录
│   ├── Sample1_results.tsv
│   └── Sample2_results.tsv
└── genome_GOMC_stats.sh   # 统计脚本
```

运行命令：

```
# 为每个样本单独生成CheckM结果
for sample in Samples/*; do
    sample_name=$(basename "$sample")
    checkm lineage_wf -x fa "$sample" "./checkm_results/${sample_name}_checkm"
    checkm qa "./checkm_results/${sample_name}_checkm/lineage.ms" "./checkm_results/${sample_name}_checkm" -o 2 --tab_table > "./checkm_results/${sample_name}_results.tsv"
done

# 合并所有结果
cat ./checkm_results/*_results.tsv > ./checkm_results/combined_checkm_results.tsv

# 运行统计脚本
./genome_GOMC_stats.sh -t nested -c ./checkm_results/combined_checkm_results.tsv
```

### 使用说明

1. **生成 CheckM 文件**：

   ```bash
   # 基本命令
   checkm lineage_wf -x fa ./bins checkm_output
   checkm qa checkm_output/lineage.ms checkm_output -o 2 --tab_table > checkm_results.tsv
   ```

2. **运行统计脚本**：

   ```bash
   # 扁平结构
   ./genome_GOMC_stats.sh -c /path/to/checkm_results.tsv
   
   # 嵌套结构
   ./genome_GOMC_stats.sh -t nested -c /path/to/checkm_results.tsv
   ```

3. **输出示例**：

   ```csv
   fasta_file_name,fasta_file_md5,...,completeness(%),contamination(%),QS,quality_class
   Sample1_bin1.fa,4262639e5...,...,98.5,0.5,96.0,near-complete
   Sample1_bin2.fa,b2960bdf4...,...,87.3,1.2,81.3,high-quality
   Sample2_bin1.fa,c3d259427...,...,45.2,8.3,3.7,low-quality
   ```

### 注意事项

**本代码使用checkm版本：**

```bash
CheckM v1.2.3 conda_env python3.9
```

**依赖安装**：

```bash
sudo apt-get install bc  # 浮点数计算依赖
```

**标识符匹配规则**：

| 文件结构 |     文件名示例     |  生成的标识符   |
| :------: | :----------------: | :-------------: |
| 扁平结构 |   `bin.1.fa.gz`    |     `bin.1`     |
| 嵌套结构 | `Sample1/bin.1.fa` | `Sample1_bin.1` |