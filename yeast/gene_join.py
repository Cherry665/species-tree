import os
import pandas as pd
from tqdm import tqdm # 用于显示进度条，如果没有可以 pip install tqdm

# --- 配置区 ---
base_dir = "./"  # 存放 440 个文件夹的根目录
file_suffix = "_copy_num.tsv"    # 文件的后缀名
output_file = "all_strains_copy_number_matrix.tsv"
# --------------

all_data = []

# 1. 遍历所有文件夹
folders = [f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f))]

print(f"开始处理 {len(folders)} 个文件夹...")

for strain_name in tqdm(folders):
    # 构建文件的完整路径，例如: ./strains_folders/AAA/AAA_copy_num.tsv
    file_path = os.path.join(base_dir, strain_name, f"{strain_name}{file_suffix}")
    
    if os.path.exists(file_path):
        # 读取单件文件，根据你提供的示例，第一列是值，第二列是基因名
        try:
            df = pd.read_csv(file_path, sep='\t', header=None, names=['copy_num', 'gene_id'])
            
            # 转换为长格式：增加一列菌株名
            df['strain'] = strain_name
            all_data.append(df)
        except Exception as e:
            print(f"读取 {strain_name} 失败: {e}")
    else:
        # 如果有些文件夹里没这个文件，跳过或记录
        pass

# 2. 合并所有数据
print("正在合并数据并转置矩阵...")
combined_df = pd.concat(all_data, ignore_index=True)

# 3. 使用数据透视表 (Pivot) 将其转换为矩阵
# 行 (Index) 为基因名，列 (Columns) 为菌株名，值为复制数
matrix_df = combined_df.pivot(index='gene_id', columns='strain', values='copy_num')

# 4. 保存结果
matrix_df.to_csv(output_file, sep='\t')

print(f"✅ 整合完成！大表已保存至: {output_file}")
print(f"结果预览 (前5行5列):")
print(matrix_df.iloc[:5, :5])