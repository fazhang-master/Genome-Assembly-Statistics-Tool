import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sqlalchemy import create_engine
from urllib.parse import quote_plus
import matplotlib as mpl
from matplotlib.patches import Patch
import configparser
import os

def get_mags_data(config_path, table_name):
    """从MySQL数据库获取中高质量MAGs的完整性和污染数据"""
    config = configparser.ConfigParser()
    config.read(config_path)
    
    # 验证配置有效性
    if 'database' not in config:
        raise ValueError("配置文件缺少[database]部分")
    
    db_config = config['database']
    required_keys = ['host', 'user', 'password', 'database']
    for key in required_keys:
        if key not in db_config:
            raise ValueError(f"配置文件缺少必需的数据库参数: {key}")
    
    # 安全处理密码特殊字符
    password = db_config['password']
    encoded_pwd = quote_plus(password)
    
    # 创建数据库连接
    engine = create_engine(
        f"mysql+mysqlconnector://{db_config['user']}:{encoded_pwd}@{db_config['host']}/{db_config['database']}"
    )
    
    # 查询中高质量MAGs
    query = f"""
    SELECT completeness, contamination 
    FROM {table_name} 
    WHERE completeness >= 50 
      AND contamination <= 10
    """
    df = pd.read_sql(query, engine)
    
    # 质量分类（按MIMAG标准）
    def classify_quality(row):
        comp = row['completeness']
        contam = row['contamination']
        
        if comp >= 90 and contam < 5:
            return 'Near complete'
        elif (comp >= 70 and comp < 90 and contam < 10) or (comp >= 90 and contam >= 5 and contam <= 10):
            return 'High quality'
        elif comp >= 50 and comp < 70 and contam < 10:
            return 'Medium quality'
        else:
            return 'Other'
    
    df['quality'] = df.apply(classify_quality, axis=1)
    return df

def plot_mags_distribution(df, table_name, output_path, dpi=300):
    """可视化完整性与污染分布"""
    plt.style.use('default')
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['axes.edgecolor'] = '#333F4B'
    mpl.rcParams['axes.linewidth'] = 1.2
    
    fig = plt.figure(figsize=(14, 10), dpi=dpi)
    grid = plt.GridSpec(4, 4, hspace=0.5, wspace=0.5)
    
    # 主散点图
    ax_main = fig.add_subplot(grid[1:4, :3])
    
    # 质量分类颜色编码
    colors = {
        'Near complete': '#e74c3c',
        'High quality': '#3498db',
        'Medium quality': '#2ecc71',
        'Other': '#95a5a6'
    }
    
    # 绘制各质量等级散点
    for quality, color in colors.items():
        subset = df[df['quality'] == quality]
        if not subset.empty:
            ax_main.scatter(
                x=subset['completeness'],
                y=subset['contamination'],
                s=45,
                alpha=0.7 if quality != 'Other' else 0.4,
                color=color,
                edgecolor='w',
                linewidth=0.5,
                label=f'{quality} ({len(subset)} MAGs, {len(subset)/len(df)*100:.1f}%)'
            )
    
    # 坐标轴设置
    ax_main.set_xticks([50, 60, 70, 80, 90, 100])
    ax_main.set_yticks([0, 2.5, 5, 7.5, 10])
    ax_main.set_xlim(49, 101)
    ax_main.set_ylim(-0.3, 10.3)
    ax_main.set_xlabel('Completeness (%)', fontsize=14, labelpad=10, fontweight='bold')
    ax_main.set_ylabel('Contamination (%)', fontsize=14, labelpad=10, fontweight='bold')
    ax_main.grid(True, linestyle='--', alpha=0.3)
    
    # 质量阈值参考线
    ax_main.axvline(90, color='#E74C3C', linestyle='--', alpha=0.8, linewidth=1.5)
    ax_main.axhline(5, color='#3498DB', linestyle='--', alpha=0.8, linewidth=1.5)
    ax_main.text(91, 9.5, "High-Quality Threshold", color='#E74C3C', fontsize=11, ha='left')
    ax_main.text(52, 5.3, "Low-Contamination Threshold", color='#3498DB', fontsize=11)
    
    # 完整性分布直方图（顶部）
    ax_histx = fig.add_subplot(grid[0, :3])
    bins = np.arange(50, 101, 2)
    for quality, color in colors.items():
        if quality != 'Other':
            subset = df[df['quality'] == quality]
            ax_histx.hist(
                subset['completeness'], bins=bins, 
                color=color, alpha=0.85, edgecolor='white',
                label=quality, stacked=True
            )
    ax_histx.set_title('Completeness Distribution', fontsize=12, pad=10)
    ax_histx.set_ylabel('Count', fontsize=10)
    ax_histx.tick_params(axis='x', labelbottom=False)
    ax_histx.spines[['top', 'right']].set_visible(False)
    
    # 污染分布直方图（右侧）
    ax_histy = fig.add_subplot(grid[1:4, 3])
    bins_y = np.arange(0, 10.5, 0.5)
    for quality, color in colors.items():
        if quality != 'Other':
            subset = df[df['quality'] == quality]
            ax_histy.hist(
                subset['contamination'], bins=bins_y, 
                orientation='horizontal', color=color, alpha=0.85, 
                edgecolor='white', stacked=True
            )
    ax_histy.set_title('Contamination Distribution', fontsize=12, pad=10)
    ax_histy.set_xlabel('Count', fontsize=10)
    ax_histy.tick_params(axis='y', labelleft=False)
    ax_histy.spines[['top', 'right']].set_visible(False)
    
    # 质量统计信息框
    near_complete = df[df['quality'] == 'Near complete']
    high_quality = df[df['quality'] == 'High quality']
    medium_quality = df[df['quality'] == 'Medium quality']
    
    stats_text = (
        f"Near complete: {len(near_complete)} MAGs, {len(near_complete)/len(df)*100:.1f}%\n"
        f"High quality: {len(high_quality)} MAGs, {len(high_quality)/len(df)*100:.1f}%\n"
        f"Medium quality: {len(medium_quality)} MAGs, {len(medium_quality)/len(df)*100:.1f}%"
    )
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='#BDC3C7', linewidth=1.2)
    ax_main.text(0.95, 0.95, stats_text, transform=ax_main.transAxes,
                 fontsize=12, verticalalignment='top', horizontalalignment='right',
                 bbox=props)
    
    # 图例
    legend_elements = [
        Patch(facecolor=colors['Near complete'], label='Near complete'),
        Patch(facecolor=colors['High quality'], label='High quality'),
        Patch(facecolor=colors['Medium quality'], label='Medium quality'),
        Patch(facecolor=colors['Other'], label='Other')
    ]
    ax_main.legend(handles=legend_elements, loc='upper left', frameon=True, 
                   edgecolor='#333', fontsize=11, bbox_to_anchor=(0.01, 0.99),
                   title="Quality Classes", title_fontsize=12)
    
    # 标题和输出
    plt.suptitle(f'Distribution of Completeness and Contamination in {table_name}', 
                 fontsize=18, y=0.97, fontweight='bold')
    
    output_format = os.path.splitext(output_path)[1][1:].lower()
    if output_format not in ['png', 'tiff', 'jpg', 'jpeg', 'svg']:
        output_format = 'png'
    
    plt.savefig(output_path, format=output_format, dpi=dpi, bbox_inches='tight')
    plt.close()
    return len(df)

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
        log(f"执行完整性与污染分析:")
        log(f"配置文件: {args.config}")
        log(f"分析表名: {args.table}")
        log(f"输出路径: {args.plot}")
        log(f"日志文件: {log_path}")
        log(f"图像DPI: {args.dpi}")
        log(f"{'='*50}")
        
        try:
            # 获取数据
            mags_df = get_mags_data(args.config, args.table)
            log(f"成功获取 {len(mags_df)} 个中高质量MAGs数据")
            
            # 质量统计
            log("\n质量分类统计:")
            for quality in ['Near complete', 'High quality', 'Medium quality', 'Other']:
                count = len(mags_df[mags_df['quality'] == quality])
                percentage = count / len(mags_df) * 100
                log(f"  - {quality}: {count} MAGs ({percentage:.1f}%)")
            
            # 生成可视化
            mags_count = plot_mags_distribution(mags_df, args.table, args.plot, args.dpi)
            img_format = os.path.splitext(args.plot)[1][1:]
            log(f"\n图表已保存: '{args.plot}' (格式: {img_format}, DPI: {args.dpi})")
            log(f"分析完成! 共处理 {mags_count} 个中高质量MAGs")
            
        except Exception as e:
            import traceback
            error_msg = f"\n错误: {str(e)}\n{traceback.format_exc()}"
            log(error_msg)
            log("可能原因:")
            log("1. 配置文件路径错误或格式不正确")
            log("2. 数据库服务不可用")
            log("3. 指定表不存在")
            log("4. 数据库连接权限不足")
            log("\n配置文件示例格式:")
            log("[database]")
            log("host = localhost")
            log("user = username")
            log("password = your_password")
            log("database = db_name")
            raise  # 向上抛出异常供主程序处理
    
    print(f"分析日志已保存至: {log_path}")