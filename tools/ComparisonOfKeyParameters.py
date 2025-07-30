#!/usr/bin/env python
# -*- coding: utf-8 -*-

import configparser
import pymysql
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

def check_and_import_gomc_tables(conn, sql_file_path):
    """检查并导入GOMC表（如果不存在）"""
    try:
        with conn.cursor() as cursor:
            # 检查GOMCBasicData表是否存在
            cursor.execute("""
                SELECT COUNT(*)
                FROM information_schema.tables
                WHERE table_schema = DATABASE()
                AND table_name IN ('GOMCBasicData', 'GOMCSuppleData')
            """)
            table_count = cursor.fetchone()[0]
            
            # 如果两个表都不存在，则导入SQL文件
            if table_count < 2:
                print(f"检测到缺少GOMC表（找到{table_count}/2），正在从 {sql_file_path} 导入...")
                
                # 检查SQL文件是否存在
                if not os.path.exists(sql_file_path):
                    raise FileNotFoundError(f"SQL文件不存在: {sql_file_path}")
                
                # 读取并执行SQL文件
                with open(sql_file_path, 'r', encoding='utf-8') as f:
                    sql_script = f.read()
                    
                # 分割SQL语句并执行
                for statement in sql_script.split(';'):
                    if statement.strip():
                        cursor.execute(statement)
                conn.commit()
                print("GOMC表导入成功")
            else:
                print("GOMC表已存在，跳过导入")
    except Exception as e:
        raise RuntimeError(f"导入GOMC表时出错: {str(e)}")

def extract_gomc_data(conn):
    """从GOMCBasicData表提取数据"""
    try:
        # 执行SQL查询
        query = """
            SELECT 
                'GOMC' AS source,
                total_sizebp AS `total_size(bp)`,
                sequences,
                largest_seqbp AS `largest_seq(bp)`,
                smaller_seqbp AS `smallest_seq(bp)`,
                n50bp AS `N50(bp)`,
                l50 AS `L50`
            FROM GOMCBasicData Limit 1000
        """
        return pd.read_sql(query, conn)
    except Exception as e:
        raise RuntimeError(f"提取GOMC数据时出错: {str(e)}")

def extract_user_data(conn, table_pairs):
    """从用户指定的表对中提取数据"""
    all_data = pd.DataFrame()
    
    for pair in table_pairs:
        try:
            # 解析表对参数
            parts = pair.split(':')
            if len(parts) != 2:
                raise ValueError(f"无效的表对格式: {pair}. 格式应为 QS表.外键:关键参数表.外键")
            
            qs_table, param_table = parts
            qs_table = qs_table.split('.')[0]  # 提取表名
            param_table = param_table.split('.')[0]  # 提取表名
            
            print(f"处理表对: QS表={qs_table}, 关键参数表={param_table}")
            
            # 执行SQL查询获取中高质量数据
            query = f"""
                SELECT 
                    '{param_table}' AS source,
                    b.total_size_bp AS `total_size(bp)`,
                    b.sequences,
                    b.largest_seq_bp AS `largest_seq(bp)`,
                    b.smallest_seq_bp AS `smallest_seq(bp)`,
                    b.n50_bp AS `N50(bp)`,
                    b.l50 AS `L50`
                FROM {param_table} b
                JOIN {qs_table} q ON b.fasta_file_md5 = q.fasta_file_md5
                WHERE q.quality_class != 'low'  -- 仅中高质量数据
            """
            
            # 将查询结果添加到总数据框
            table_data = pd.read_sql(query, conn)
            all_data = pd.concat([all_data, table_data], ignore_index=True)
            
        except Exception as e:
            raise RuntimeError(f"处理表对 {pair} 时出错: {str(e)}")
    
    return all_data

def generate_boxplot(df, output_path, dpi=300, log=print):
    """生成箱线图并保存"""
    # 参数列表
    params = ['total_size(bp)', 'sequences', 'largest_seq(bp)', 
              'smallest_seq(bp)', 'N50(bp)', 'L50']
    
    # 创建多子图布局
    plt.figure(figsize=(18, 12))
    plt.suptitle('基因组组装关键参数比较', fontsize=24, y=0.98)
    
    # 统一设置字号
    AXIS_LABEL_FONTSIZE = 16
    TICK_LABEL_FONTSIZE = 12
    
    # 为每个参数创建子图
    for i, param in enumerate(params, 1):
        ax = plt.subplot(2, 3, i)
        
        # 对大范围参数使用对数尺度
        use_log = param in ['total_size(bp)', 'largest_seq(bp)', 'N50(bp)']
        
        # 创建箱线图
        sns.boxplot(
            x='source', 
            y=param, 
            data=df,
            palette='viridis',
            hue='source',
            legend=False,
            showfliers=True,
            width=0.6
        )
        
        # 添加抖动点显示数据分布
        sns.stripplot(
            x='source', 
            y=param, 
            data=df,
            color='black',
            alpha=0.5,
            jitter=0.2,
            size=4
        )
        
        # 设置Y轴标签
        if use_log:
            ax.set_yscale('log')
            ax.set_ylabel(f'{param} (对数尺度)', fontsize=AXIS_LABEL_FONTSIZE)
        else:
            ax.set_ylabel(param, fontsize=AXIS_LABEL_FONTSIZE)
        
        # 设置X轴标签
        ax.set_xlabel('数据来源', fontsize=AXIS_LABEL_FONTSIZE)
        ticks = ax.get_xticks()
        ax.set_xticks(ticks) 
        ax.set_xticklabels(
            ax.get_xticklabels(), 
            rotation=30,  # 关键修改：旋转30度
            ha='right',   # 右对齐防止标签重叠
            fontsize=TICK_LABEL_FONTSIZE
        )
        
        # 设置刻度标签
        ax.tick_params(axis='both', labelsize=TICK_LABEL_FONTSIZE)
        ax.grid(True, linestyle='--', alpha=0.3)
    
    # 调整布局并保存图像
    plt.tight_layout()
    
    # 根据扩展名确定格式
    output_format = os.path.splitext(output_path)[1][1:].lower()
    if output_format not in ['png', 'tiff', 'jpg', 'jpeg', 'svg']:
        output_format = 'png'
    
    plt.savefig(output_path, format=output_format, dpi=dpi, bbox_inches='tight')
    print(f"结果已保存至: {output_path}")
    
    # 输出基本统计信息
    log("\n基本统计信息:")
    for param in params:
        log(f"\n{param}统计:")
        # 将统计结果转为字符串再记录
        stats_table = df.groupby('source')[param].describe().to_string()
        log(stats_table)
    
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
        log(f"执行基因组关键参数比较分析:")
        log(f"配置文件: {args.config}")
        log(f"表对参数: {args.tables}")
        log(f"输出路径: {args.plot}")
        log(f"图像DPI: {args.dpi}")
        log(f"GOMC SQL文件: {args.gomc_sql}")
        log(f"{'='*50}")

        try:
            # 读取数据库配置
            config = configparser.ConfigParser()
            config.read(args.config)

            if 'database' not in config:
                raise ValueError("配置文件中缺少[database]部分")

            db_config = config['database']
            db_params = {
                'host': db_config.get('host', 'localhost'),
                'user': db_config.get('user', 'root'),
                'password': db_config.get('password', ''),
                'database': db_config.get('database', ''),
                'port': db_config.getint('port', 3306),
                'charset': 'utf8mb4'
            }

            # 连接数据库
            conn = pymysql.connect(**db_params)
            log("数据库连接成功")

            # 检查并导入GOMC表
            check_and_import_gomc_tables(conn, args.gomc_sql)

            # 提取GOMC数据
            gomc_df = extract_gomc_data(conn)
            log(f"成功提取 {len(gomc_df)} 条GOMC数据")

            # 提取用户表对数据
            table_pairs = [pair.strip() for pair in args.tables.split(',')]
            user_df = extract_user_data(conn, table_pairs)
            log(f"成功提取 {len(user_df)} 条用户数据")

            # 合并数据
            combined_df = pd.concat([gomc_df, user_df], ignore_index=True)
            log(f"合并后数据集大小: {len(combined_df)} 条记录")

            # 生成箱线图
            img_format = generate_boxplot(combined_df, args.plot, args.dpi, log)
            log(f"\n关键参数比较图已保存为 '{args.plot}' (格式: {img_format}, DPI: {args.dpi})")
            log(f"分析完成!")

            # 关闭数据库连接
            conn.close()

        except Exception as e:
            log(f"\n错误: {str(e)}")
            log("可能原因:")
            log("1. 配置文件路径错误或格式不正确")
            log("2. 数据库服务不可用")
            log("3. 指定表不存在或外键关联错误")
            log("4. GOMC SQL文件路径错误")
            log("5. 数据库连接权限不足")
            log("\n配置文件示例格式:")
            log("[database]")
            log("host = localhost")
            log("user = username")
            log("password = your_password")
            log("database = db_name")
            log("port = 3306")
            raise  # 向上抛出异常供主程序处理
    print(f"分析日志已保存至: {log_path}") 