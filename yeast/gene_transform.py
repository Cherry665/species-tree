import pandas as pd

# --- 文件路径配置 ---
matrix_file = "all_strains_copy_number_matrix.tsv"
sgd_file = "SGD_features.tab"
output_file = "final_copy_number_matrix_named.tsv"

# --- 1. 解析 SGD_features.tab 建立映射字典 ---
# 目标：{ 'YIL064W': 'EFM4', ... }
gene_map = {}
with open(sgd_file, 'r', encoding='utf-8') as f:
    for line in f:
        # 跳过以 ! 开头的注释行
        if line.startswith('!'):
            continue
        parts = line.strip().split('\t')
        # 根据文件结构：
        # parts[3] 是系统编号 (如 YIL064W)
        # parts[4] 是标准基因名 (如 EFM4)
        if len(parts) > 4:
            sys_name = parts[3].strip()
            std_name = parts[4].strip()
            # 只有当标准名不为空时才加入映射，否则后面保留原名
            if std_name:
                gene_map[sys_name] = std_name

print(f"✅ 已从 SGD 文件提取了 {len(gene_map)} 个基因名映射。")

# --- 2. 处理大表 ---
# 读取原始矩阵，假设第一列是基因编号 (index_col=0)
df = pd.read_csv(matrix_file, sep='\t', index_col=0)

# 使用 map 函数转换行索引（基因编号 -> 基因名）
# 如果在字典里找不到（比如一些未命名的 ORF），则保留原始编号
df.index = df.index.map(lambda x: gene_map.get(x, x))

# --- 3. 转置矩阵并让菌株名在第一列 ---
# 目前：行=基因，列=菌株
# 转置后：行=菌株，列=基因
df_t = df.T

# 将索引（即菌株名）重命名为 'Strain'，并变成普通列
df_t.index.name = 'Strain'
df_t = df_t.reset_index()

# --- 4. 保存结果 ---
df_t.to_csv(output_file, sep='\t', index=False)

print(f"✅ 转换与转置完成！结果已保存至: {output_file}")
print(f"💡 现在的格式：第一列是菌株名，表头是基因名。")




cols=$(head -n 1 final_copy_number_matrix_named.tsv | \
      tr '\t' '\n' | \
      grep -nE "^ID$|CUP1-1|CUP1-2" | \
      cut -d: -f1 | \
      paste -sd,)

cut -f"$cols" final_copy_number_matrix_named.tsv > strains_CUP1_only.tsv