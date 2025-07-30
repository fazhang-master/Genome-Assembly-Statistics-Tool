import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import configparser
import mysql.connector
import os
import seaborn as sns
from matplotlib.ticker import MaxNLocator

def generate_taxonomy_query(table_name):
    """生成分类统计查询语句"""
    return f"""
    SELECT 
        COUNT(DISTINCT NULLIF(TRIM(SUBSTRING_INDEX(SUBSTRING_INDEX(classification, ';', 1), 'd__', -1)), '')) AS domain,
        COUNT(DISTINCT NULLIF(TRIM(SUBSTRING_INDEX(SUBSTRING_INDEX(classification, ';', 2), 'p__', -1)), '')) AS phylum,
        COUNT(DISTINCT NULLIF(TRIM(SUBSTRING_INDEX(SUBSTRING_INDEX(classification, ';', 3), 'c__', -1)), '')) AS class,
        COUNT(DISTINCT NULLIF(TRIM(SUBSTRING_INDEX(SUBSTRING_INDEX(classification, ';', 4), 'o__', -1)), '')) AS order_name,
        COUNT(DISTINCT NULLIF(TRIM(SUBSTRING_INDEX(SUBSTRING_INDEX(classification, ';', 5), 'f__', -1)), '')) AS family,
        COUNT(DISTINCT NULLIF(TRIM(SUBSTRING_INDEX(SUBSTRING_INDEX(classification, ';', 6), 'g__', -1)), '')) AS genus,
        COUNT(DISTINCT NULLIF(TRIM(SUBSTRING_INDEX(SUBSTRING_INDEX(classification, ';', 7), 's__', -1)), '')) AS species
    FROM {table_name}
    WHERE classification IS NOT NULL
    """

def get_taxonomy_data(db_config, table_list):
    """从数据库获取多个表的分类统计数据"""
    data = {}
    
    try:
        conn = mysql.connector.connect(**db_config)
        cursor = conn.cursor()
        
        for table in table_list:
            query = generate_taxonomy_query(table)
            cursor.execute(query)
            results = cursor.fetchone()
            
            if results:
                # 转换为整数并存储
                data[table] = [int(x) for x in results]
        
        cursor.close()
        conn.close()
        return data
    
    except mysql.connector.Error as err:
        raise RuntimeError(f"数据库错误: {err}")
    except Exception as e:
        raise RuntimeError(f"数据处理错误: {e}")

def plot_taxonomy_comparison(data, output_path, dpi=300):
    """绘制物种丰度比较图"""
    # 设置全局样式
    sns.set_style("whitegrid")
    mpl.rcParams.update({
        'font.family': 'DejaVu Sans',
        'axes.edgecolor': '#2c3e50',
        'axes.linewidth': 1.5,
        'grid.color': '#ecf0f1',
        'grid.linestyle': '--'
    })
    
    # 准备数据
    tables = list(data.keys())
    categories = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    
    # 创建画布
    fig, ax = plt.subplots(figsize=(14, 8), dpi=dpi)
    
    # 设置颜色方案
    colors = plt.cm.viridis(np.linspace(0, 1, len(tables)))
    
    # 设置柱状图位置和宽度
    x = np.arange(len(categories))
    width = 0.8 / len(tables)  # 动态宽度调整
    
    # 绘制每个表的柱状图
    for i, table in enumerate(tables):
        values = data[table]
        position = x + i * width
        ax.bar(position, values, width, label=table, color=colors[i], 
               edgecolor='white', linewidth=0.8, alpha=0.85)
        
        # 添加数据标签
        for j, v in enumerate(values):
            ax.text(position[j], v + max(values)*0.02, str(v), 
                    ha='center', va='bottom', fontsize=9, rotation=90)
    
    # 设置图表元素
    ax.set_title('Taxonomic Diversity Comparison', fontsize=16, pad=15, fontweight='bold')
    ax.set_ylabel('Unique Count', fontsize=12, labelpad=10)
    ax.set_xlabel('Taxonomic Level', fontsize=12, labelpad=10)
    ax.set_xticks(x + width * (len(tables) - 1) / 2)
    ax.set_xticklabels(categories, fontsize=11)
    
    # 设置Y轴为整数刻度
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    
    # 添加图例和网格
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), 
              ncol=min(3, len(tables)), frameon=True, framealpha=0.9)
    ax.grid(True, axis='y', alpha=0.4)
    
    # 优化布局并保存
    plt.tight_layout()
    
    # 根据扩展名确定格式
    output_format = os.path.splitext(output_path)[1][1:].lower()
    if output_format not in ['png', 'tiff', 'jpg', 'jpeg', 'svg']:
        output_format = 'png'
    
    plt.savefig(output_path, format=output_format, dpi=dpi, bbox_inches='tight')
    plt.close()
    return output_format

def run(args):
    """模块主入口，由Analysis.py调用"""
    # 创建日志文件路径（与图片同路径、同文件名，扩展名为.txt）
    log_path = os.path.splitext(args.plot)[0] + '.txt'
    
    # 打开日志文件用于写入
    with open(log_path, 'w', encoding='utf-8') as log_file:
        def log(message):
            """同时输出到控制台和日志文件"""
            print(message)
            log_file.write(message + '\n')
    
        log(f"\n{'='*50}")
        log(f"执行物种丰度分析:")
        log(f"配置文件: {args.config}")
        log(f"分析表名: {args.tables}")
        log(f"输出路径: {args.plot}")
        log(f"图像DPI: {args.dpi}")
        log(f"{'='*50}")
    
        try:
            # 读取数据库配置
            config = configparser.ConfigParser()
            config.read(args.config)

            if 'database' not in config:
                raise ValueError("配置文件中缺少[database]部分")

            db_config = {
                'host': config['database']['host'],
                'user': config['database']['user'],
                'password': config['database']['password'],
                'database': config['database']['database']
            }

            # 处理表名列表
            table_list = [t.strip() for t in args.tables.split(',')]

            # 获取分类数据
            taxonomy_data = get_taxonomy_data(db_config, table_list)

            if not taxonomy_data:
                raise RuntimeError("未获取到有效数据，请检查表名和数据库连接")

            # 打印数据统计
            log("\n物种分类统计:")
            for table, counts in taxonomy_data.items():
                log(f"  - {table}:")
                log(f"    Domain: {counts[0]}, Phylum: {counts[1]}, Class: {counts[2]}")
                log(f"    Order: {counts[3]}, Family: {counts[4]}, Genus: {counts[5]}, Species: {counts[6]}")

            # 生成图表
            img_format = plot_taxonomy_comparison(taxonomy_data, args.plot, args.dpi)
            log(f"\n物种丰度比较图已保存为 '{args.plot}' (格式: {img_format}, DPI: {args.dpi})")
            log(f"分析完成!")

        except Exception as e:
            log(f"\n错误: {str(e)}")
            log("可能原因:")
            log("1. 配置文件路径错误或格式不正确")
            log("2. 数据库服务不可用")
            log("3. 指定表不存在或未包含classification字段")
            log("4. 数据库连接权限不足")
            log("\n配置文件示例格式:")
            log("[database]")
            log("host = localhost")
            log("user = username")
            log("password = your_password")
            log("database = db_name")
            raise  # 向上抛出异常供主程序处理
    print(f"分析日志已保存至: {log_path}")