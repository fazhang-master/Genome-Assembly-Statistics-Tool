import configparser
import matplotlib.pyplot as plt
import mysql.connector
import pandas as pd
import seaborn as sns
import os

# 质量分类映射字典
QUALITY_MAP = {
    'low-quality': 'low',
    'medium-quality': 'medium',
    'high-quality': 'high',
    'near-complete': 'near complete'
}

def get_db_config(config_path):
    """从配置文件读取数据库连接参数"""
    config = configparser.ConfigParser()
    config.read(config_path)
    
    if 'database' not in config:
        raise ValueError("配置文件缺少[database]部分")
    
    return {
        'host': config['database']['host'],
        'user': config['database']['user'],
        'password': config['database']['password'],
        'database': config['database']['database']
    }

def get_quality_class_stats(db_config, table_name):
    """从数据库获取质量分类统计数据"""
    try:
        conn = mysql.connector.connect(**db_config)
        query = f"SELECT quality_class FROM {table_name} WHERE quality_class IS NOT NULL"
        df = pd.read_sql(query, conn)
        conn.close()
        
        # 应用质量映射
        df['quality_class'] = df['quality_class'].map(QUALITY_MAP)
        
        # 统计类别分布
        class_counts = df['quality_class'].value_counts()
        total = class_counts.sum()
        percentages = (class_counts / total * 100).round(1)
        
        return class_counts, percentages, total
    
    except Exception as e:
        raise RuntimeError(f"数据库查询失败: {str(e)}")

def plot_quality_class_pie(class_counts, percentages, table_name, output_file, dpi=300):
    """生成专业质量分布饼图"""
    # 设置全局样式
    sns.set_style("whitegrid")
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'axes.edgecolor': '#2c3e50',
        'axes.linewidth': 1.5,
    })
    
    # 创建画布
    fig = plt.figure(figsize=(10, 8), dpi=dpi)
    
    # 颜色方案
    palette = sns.color_palette("pastel", len(class_counts))
    explode = [0.05] * len(class_counts)  # 轻微分离所有扇区
    
    # 创建饼图
    wedges, texts, autotexts = plt.pie(
        class_counts,
        labels=[f"{label}\n({percent}%)" for label, percent in zip(class_counts.index, percentages)],
        colors=palette,
        autopct=lambda p: f'{p:.1f}%',
        startangle=90,
        explode=explode,
        wedgeprops={'linewidth': 1.5, 'edgecolor': 'white'},
        textprops={'fontsize': 10}
    )
    
    # 设置文本样式
    plt.setp(autotexts, size=10, weight="bold", color='black')
    
    # 添加标题和图例
    plt.title(f'Genome Quality Distribution: {table_name}', fontsize=14, pad=20, weight='bold')
    plt.legend(wedges, class_counts.index,
               title="Quality Classes",
               loc='upper left',
               bbox_to_anchor=(1, 0.5),
               frameon=True,
               edgecolor='#ddd',
               fancybox=True)
    
    # 添加中心统计信息
    total = class_counts.sum()
    plt.text(0, 0, f"Total Genomes: {total}", 
             ha='center', va='center', fontsize=12, weight='bold')
    
    # 保存图像
    plt.axis('equal')
    plt.tight_layout()
    
    # 根据扩展名确定格式
    output_format = os.path.splitext(output_file)[1][1:].lower()
    if output_format not in ['png', 'tiff', 'jpg', 'jpeg', 'svg']:
        output_format = 'png'
    
    plt.savefig(output_file, format=output_format, dpi=dpi, bbox_inches='tight')
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
        log(f"执行质量分布分析:")
        log(f"配置文件: {args.config}")
        log(f"分析表名: {args.table}")
        log(f"输出路径: {args.plot}")
        log(f"图像DPI: {args.dpi}")
        log(f"{'='*50}")

        try:
            # 获取数据库配置
            db_config = get_db_config(args.config)

            # 获取质量统计数据
            class_counts, percentages, total = get_quality_class_stats(db_config, args.table)
            if total == 0:
                raise RuntimeError("未获取到有效数据，请检查表名和数据库连接")

            # 打印统计信息
            log("\n质量类别分布统计:")
            log("-" * 50)
            for cls, count in class_counts.items():
                log(f"  - {cls}: {count} genomes ({percentages[cls]}%)")
            log(f"\n总计: {total} 个基因组")

            # 生成可视化图表
            img_format = plot_quality_class_pie(class_counts, percentages, args.table, args.plot, args.dpi)
            log(f"\n质量分布图已保存为 '{args.plot}' (格式: {img_format}, DPI: {args.dpi})")
            log(f"分析完成!")

        except Exception as e:
            log(f"\n错误: {str(e)}")
            log("可能原因:")
            log("1. 配置文件路径错误或格式不正确")
            log("2. 数据库服务不可用")
            log("3. 指定表不存在或未包含quality_class字段")
            log("4. 数据库连接权限不足")
            log("\n配置文件示例格式:")
            log("[database]")
            log("host = localhost")
            log("user = username")
            log("password = your_password")
            log("database = db_name")
            raise  # 向上抛出异常供主程序处理
    print(f"分析日志已保存至: {log_path}")