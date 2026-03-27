import re

# --- 文件路径 ---
tree_id_file = "NR.lst"                
old_annotation = "itol_annotation.txt" 
new_annotation = "itol_perfect_branches.txt"

# 1. 提取颜色和环境，准备图例
color_map = {}
legend_colors = []
legend_labels = []

with open(old_annotation, 'r', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith(("DATASET", "SEPARATOR", "COLOR", "STRIP", "DATA")):
            continue
        parts = line.split('\t')
        if len(parts) >= 3:
            short_id = parts[0].strip()
            color = parts[1].strip()
            origin = parts[2].strip()
            
            color_map[short_id] = color
            
            # 收集图例信息 (去重)
            if color not in legend_colors:
                legend_colors.append(color)
                legend_labels.append(origin)

# 2. 生成完全符合 iTOL 标准的 6 列文件
with open(new_annotation, 'w', encoding='utf-8', newline='\n') as f:
    f.write("DATASET_STYLE\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tBranch_Colors_Perfect\n")
    f.write("COLOR\t#000000\n")
    
    # 写入图例
    f.write("LEGEND_TITLE\tEcological_Origins\n")
    f.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(legend_colors)) + "\n")
    f.write("LEGEND_COLORS\t" + "\t".join(legend_colors) + "\n")
    f.write("LEGEND_LABELS\t" + "\t".join(legend_labels) + "\n")
    
    f.write("DATA\n")
    
    with open(tree_id_file, 'r', encoding='utf-8') as ids:
        for full_id in ids:
            full_id = full_id.strip()
            if not full_id: continue
            
            parts = re.split(r'[_.]', full_id)
            clean_id = "_".join(parts[:2]) if full_id.startswith("SACE") else parts[0]
            
            if clean_id in color_map:
                color = color_map[clean_id]
                # ！！！真正的修复：加入了 'node' 作为第 3 列 ！！！
                # 格式：ID \t TYPE \t WHAT \t COLOR \t WIDTH \t STYLE
                f.write(f"{full_id}\tbranch\tnode\t{color}\t2\tnormal\n")

print(f"✅ 完美修复！已生成: {new_annotation}")