import pandas as pd

# --- 1. 设置文件路径 ---
excel_file = "菌株具体信息.xlsx"
my_list_file = "isolation.txt"
output_file = "Ecological origins.txt"

# --- 2. 读取数据 ---
# header=1 表示从 Excel 的第 2 行（索引为 1）开始读取作为表头
df_excel = pd.read_excel(excel_file, header=1)

# 读取你现有的菌株列表
my_strains = pd.read_csv(my_list_file, header=None, names=['Strain'])

# --- 3. 核心匹配逻辑 ---
# 根据截图，环境在 'Isolation'
excel_strain_col = 'Standardized name'
excel_env_col = 'Ecological origins'

# 提取需要的两列
df_mapping = df_excel[[excel_strain_col, excel_env_col]]

# 匹配逻辑
result = pd.merge(my_strains, df_mapping, 
                  left_on='Strain', 
                  right_on=excel_strain_col, 
                  how='left')

# --- 4. 清理并保存 ---
# 填充缺失值为 Unknown
result['Final_Isolation'] = result[excel_env_col].fillna('Unknown')

# 保存为两列：原始名 和 匹配到的环境
result[['Strain', 'Final_Isolation']].to_csv(output_file, sep='\t', index=False, header=False)

print(f"✅ 匹配成功！结果已保存至: {output_file}")