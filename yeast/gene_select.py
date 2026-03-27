import pandas as pd

# 1. 读取整合后的大表 (假设你已经生成了 tsv)
input_file = "all_strains_copy_number_matrix.tsv"
df = pd.read_csv(input_file, sep='\t', index_col='gene_id')

# 2. 定义“存在”的筛选逻辑
# 逻辑 A：如果矩阵中没有 NaN，且每个菌株的值都 > 0
# 我们先将 NaN 填充为 0 (代表不存在)
df_filled = df.fillna(0)

# 3. 执行筛选
# 统计每一行（每个基因）中，值大于 1 的列数是否等于总列数（菌株数）
total_strains = df_filled.shape[1]
core_genes_df = df_filled[(df_filled > 1).all(axis=1)]

# 4. 保存结果
output_file = "core_genes_matrix.tsv"
core_genes_df.to_csv(output_file, sep='\t')

# 5. 输出统计结果
print(f"📊 统计信息：")
print(f"  - 总菌株数: {total_strains}")
print(f"  - 原始总基因数: {df.shape[0]}")
print(f"  - 筛选出的核心基因数 (每个菌株都存在): {core_genes_df.shape[0]}")
print(f"✅ 核心基因列表已保存至: {output_file}")