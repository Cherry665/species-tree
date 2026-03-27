import pandas as pd

# 1. 定义颜色字典
color_dict = {
    "Wine": "#e41a1c", "Beer": "#e41a1c", "Sake": "#e41a1c", "Distillery": "#e41a1c",
    "Human, clinical": "#984ea3", "Probiotic": "#984ea3",
    "Soil": "#4daf4a", "Tree": "#4daf4a", "Nature": "#4daf4a", "Flower": "#4daf4a", "Fruit": "#4daf4a",
    "Industrial": "#ff7f00", "Bioethanol": "#ff7f00",
    "Dairy": "#377eb8",
    "Water": "#a6cee3",
    "Palm wine": "#fb9a99",
    "Fermentation": "#cab2d6",
    "Unknown": "#d9d9d9"
}

# 2. 读取你的数据文件
# 请确保该目录下确实存在 "Ecological origins.txt" 这个文件
df = pd.read_csv("Ecological origins.txt", sep="\t")

# 3. 准备生成 iTOL 格式
output_file = "itol_annotation.txt"
with open(output_file, "w") as f:
    f.write("DATASET_COLORSTRIP\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tEcological_Origins\n")
    f.write("COLOR\t#ff0000\n")
    f.write("STRIP_WIDTH\t30\n")
    f.write("DATA\n")
    
    for _, row in df.iterrows():
        strain = str(row['Strain'])
        origin = str(row['Ecological origins'])
        
        # 匹配颜色
        color = color_dict.get(origin, "#d9d9d9")
        
        # 写入 iTOL 格式
        f.write(f"{strain}\t{color}\t{origin}\n")

print(f"✅ 已生成 iTOL 注释文件: {output_file}")
