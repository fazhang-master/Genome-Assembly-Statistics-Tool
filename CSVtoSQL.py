import argparse
import configparser
import pandas as pd
import numpy as np
import mysql.connector
from mysql.connector import Error

def sanitize_column_names(df):
    df.columns = [col.replace('(%)', '').replace('%', '_percent')
                 .replace(' ', '_').replace('-', '_')
                 .lower() for col in df.columns]
    return df

def create_tables(cursor, basic_table, supply_table, uk_column):
    """创建基础表和扩增表，并设置外键关联"""
    # 基础表结构
    basic_table_sql = f"""
    CREATE TABLE IF NOT EXISTS `{basic_table}` (
        `fasta_file_name` VARCHAR(255) NOT NULL,
        `fasta_file_md5` VARCHAR(32) NOT NULL,
        `total_size_bp` BIGINT,
        `sequences` INT,
        `largest_seq_bp` INT,
        `smallest_seq_bp` INT,
        `n50_bp` INT,
        `l50` INT,
        `classification` TEXT,
        PRIMARY KEY (`{uk_column}`),
        UNIQUE KEY `idx_fasta_file_name` (`fasta_file_name`),
        UNIQUE KEY `idx_fasta_file_md5` (`fasta_file_md5`)
    ) ENGINE=InnoDB;
    """
    
    # 扩增表结构
    supply_table_sql = f"""
    CREATE TABLE IF NOT EXISTS `{supply_table}` (
        `id` INT AUTO_INCREMENT PRIMARY KEY,
        `fasta_file_name` VARCHAR(255) NOT NULL,
        `fasta_file_md5` VARCHAR(32) NOT NULL,
        `completeness` FLOAT,
        `contamination` FLOAT,
        `qs` FLOAT,
        `quality_class` VARCHAR(50),
        CONSTRAINT `fk_{supply_table}_md5` 
            FOREIGN KEY (`fasta_file_md5`) 
            REFERENCES `{basic_table}` (`fasta_file_md5`)
            ON DELETE CASCADE
    ) ENGINE=InnoDB;
    """
    
    cursor.execute(basic_table_sql)
    cursor.execute(supply_table_sql)
    print(f"Tables created: {basic_table}, {supply_table}")

def import_df_to_mysql(cursor, table_name, df, batch_size=500):
    """核心修正：处理数据类型转换"""
    df = sanitize_column_names(df)
    # 关键修复：转换numpy类型为Python原生类型
    df = df.replace({np.nan: None})
    df = df.astype(object)
    
    # 生成插入语句
    columns = ', '.join([f'`{col}`' for col in df.columns])
    placeholders = ', '.join(['%s'] * len(df.columns))
    insert_sql = f"INSERT IGNORE INTO `{table_name}` ({columns}) VALUES ({placeholders})"
    
    # 分批插入
    total_rows = len(df)
    for i in range(0, total_rows, batch_size):
        batch = df.iloc[i:i+batch_size]
        data_tuples = [tuple(x) for x in batch.to_records(index=False)]
        cursor.executemany(insert_sql, data_tuples)
        print(f"Inserted {min(i+batch_size, total_rows)}/{total_rows} rows")
    
    print(f"✅ Imported {total_rows} rows into {table_name}")

def main():
    # 命令行参数解析
    parser = argparse.ArgumentParser(description='Import CSV to MySQL with foreign key relation')
    parser.add_argument('-bd', '--basicdata', required=True, help='Path to basic data CSV')
    parser.add_argument('-sd', '--supplydata', required=True, help='Path to supply data CSV')
    parser.add_argument('-uk', '--uniquekey', default='fasta_file_md5', 
                        choices=['fasta_file_name', 'fasta_file_md5'],
                        help='Unique key for table relation')
    parser.add_argument('-obd', '--outputbasicdata', required=True, help='Output table name for basic data')
    parser.add_argument('-osd', '--outputsupplydata', required=True, help='Output table name for supply data')
    parser.add_argument('-c', '--config', default='db_config.ini', help='Database config file path')
    args = parser.parse_args()

    # 读取数据库配置
    config = configparser.ConfigParser()
    config.read(args.config)
    db_config = {
        'host': config['database']['host'],
        'user': config['database']['user'],
        'password': config['database']['password'],
        'database': config['database']['database']
    }

    try:
        # 连接数据库
        connection = mysql.connector.connect(**db_config)
        cursor = connection.cursor()
        
        # 读取CSV文件
        df_basic = pd.read_csv(args.basicdata, dtype={
                'total_size_bp': 'Int64', 
                'sequences': 'Int64',
                'n50_bp': 'Int64'
            })
        df_supply = pd.read_csv(args.supplydata)
        
        # 创建表结构
        create_tables(cursor, args.outputbasicdata, args.outputsupplydata, args.uniquekey)
        
        # 导入数据（先基础表后扩增表）
        import_df_to_mysql(cursor, args.outputbasicdata, df_basic)
        import_df_to_mysql(cursor, args.outputsupplydata, df_supply)
        
        connection.commit()
        print("✅ Data imported successfully!")
        
    except Error as e:
        print(f"❌ Database error: {e}")
        connection.rollback()
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()
            print("MySQL connection closed")

if __name__ == "__main__":
    main()