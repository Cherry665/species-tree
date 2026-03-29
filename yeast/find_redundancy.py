import os

# --- 阈值设置 (针对 859个菌株, 61个基因) ---
MIN_STRAIN_PRESENCE = 817  # 基因门槛：至少在 817 个菌株中存在 (约 95%)
MAX_STRAIN_MULTI = 86      # 基因门槛：最多允许在 86 个菌株中是多拷贝 (约 5%)
MAX_GENE_MULTI = 10        # 菌株门槛：如果一个菌株有超过 10 个 marker 都是多拷贝，说明该菌株质量极差

matrix_file = "copy_number_matrix.tsv" # 确保这个文件在当前目录下
gene_omit_file = "Protein/Fungi61_marker.omit.lst"
strain_omit_file = "Protein/Fungi61_strain.omit.lst"

# 确保输出目录存在
os.makedirs("Protein", exist_ok=True)

gene_stats = {}   # 记录每个基因的状态
strain_stats = {} # 记录每个菌株的状态

print("1. 正在读取拷贝数矩阵并计算双向统计...")
with open(matrix_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 3: continue
        strain, gene, count = parts[0], parts[1], int(parts[2])
        
        # 初始化字典
        if gene not in gene_stats: gene_stats[gene] = {'present':0, 'multi':0}
        if strain not in strain_stats: strain_stats[strain] = {'present':0, 'multi':0}
        
        # 统计频次
        gene_stats[gene]['present'] += 1
        strain_stats[strain]['present'] += 1
        if count > 1:
            gene_stats[gene]['multi'] += 1
            strain_stats[strain]['multi'] += 1

print("2. 正在生成剔除名单...")
bad_genes = 0
bad_strains = 0

# 筛基因
with open(gene_omit_file, 'w') as f_gene:
    for gene, stats in gene_stats.items():
        if stats['present'] < MIN_STRAIN_PRESENCE or stats['multi'] > MAX_STRAIN_MULTI:
            f_gene.write(f"{gene}\n")
            bad_genes += 1

# 筛菌株
with open(strain_omit_file, 'w') as f_strain:
    for strain, stats in strain_stats.items():
        if stats['multi'] > MAX_GENE_MULTI:
            f_strain.write(f"{strain}\n")
            bad_strains += 1

print(f"执行完毕！")
print(f" -> 剔除了 {bad_genes} 个高冗余/高缺失基因，名单保存在 {gene_omit_file}")
print(f" -> 剔除了 {bad_strains} 个高冗余烂菌株，名单保存在 {strain_omit_file}")