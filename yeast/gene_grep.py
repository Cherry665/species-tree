import os
import shutil

# --- 配置区 ---
list_file = "NR.lst"          # 你的清单文件
base_dir = "./"               # 菌株文件夹所在的根目录（"." 表示当前目录）
# --------------

# 1. 读取清单，存入 set 提高查询效率
with open(list_file, 'r', encoding='utf-8') as f:
    keep_names = {line.strip() for line in f if line.strip()}

print(f"📦 清单中共有 {len(keep_names)} 个待保留菌株。")

# 2. 遍历目录下的所有项
all_items = os.listdir(base_dir)
removed_count = 0
keep_count = 0

for item in all_items:
    item_path = os.path.join(base_dir, item)
    
    # 只处理文件夹
    if os.path.isdir(item_path):
        if item in keep_names:
            keep_count += 1
            # print(f"保持: {item}") # 如果想看保留了哪些，取消注释
        else:
            # --- 危险操作：执行删除 ---
            shutil.rmtree(item_path) 
            removed_count += 1
            print(f"已清理: {item}")

print("-" * 30)
print(f"✅ 处理完成！")
print(f"保留文件夹: {keep_count} 个")
print(f"删除文件夹: {removed_count} 个")