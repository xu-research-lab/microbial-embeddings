import os
import math

# 配置区（必须修改！）
REPO_ROOT = "/home/cjj/projects/microbial-embeddings"  # Git仓库根目录（绝对路径）
TARGET_DIR = "analysis/Metabolic/Data/OTU_bigg_gene"           # 要分批上传的目录（相对于仓库根目录）
OUTPUT_DIR = "./"         # 脚本输出目录
MAX_BATCH_SIZE_GB = 1.8              # 每批最大大小（GB）
BRANCH_NAME = "main"                 # 目标分支

def calculate_batches():
    """计算文件分批"""
    target_path = os.path.join(REPO_ROOT, TARGET_DIR)
    if not os.path.exists(target_path):
        raise FileNotFoundError(f"目标目录不存在: {target_path}")

    batches = []
    current_batch = []
    current_size = 0
    max_size = MAX_BATCH_SIZE_GB * 1024**3

    for root, _, files in os.walk(target_path):
        for file in files:
            file_path = os.path.join(root, file)
            rel_path = os.path.relpath(file_path, REPO_ROOT)
            file_size = os.path.getsize(file_path)

            # GitHub文件大小限制检查
            if file_size > 100 * 1024**2:
                print(f"⚠️ 注意: {rel_path} 超过GitHub 100MB限制，需要Git LFS")
                continue

            if current_size + file_size > max_size:
                batches.append(current_batch)
                current_batch = []
                current_size = 0

            current_batch.append(rel_path)
            current_size += file_size

    if current_batch:
        batches.append(current_batch)
    
    return batches

def generate_scripts(batches):
    """生成可执行的批处理脚本"""
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    print(f"\n🗂 将在 {os.path.abspath(OUTPUT_DIR)} 生成 {len(batches)} 个脚本：")

    for i, batch in enumerate(batches, 1):
        script_path = os.path.join(OUTPUT_DIR, f"batch_{i}.sh")
        batch_size_gb = sum(os.path.getsize(os.path.join(REPO_ROOT, f)) for f in batch) / 1024**3

        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(f"# 批次 {i} - {len(batch)} 个文件, {batch_size_gb:.2f}GB\n")
            f.write(f"cd {REPO_ROOT}\n\n")  # 切换到仓库目录
            
            # 添加文件
            for file in batch:
                f.write(f"git add '{file}'\n")
            
            # 提交并推送
            f.write(f'\ngit commit -m "添加文件批次 {i} ({len(batch)}个)"\n')
            f.write(f"git push origin {BRANCH_NAME}\n")
            f.write('echo "✅ 批次 ${i} 完成"\n')

        os.chmod(script_path, 0o755)  # 添加执行权限
        print(f"✅ 生成: batch_{i}.sh ({len(batch)} files, {batch_size_gb:.2f}GB)")

if __name__ == "__main__":
    print("🔍 正在分析文件结构...")
    try:
        batches = calculate_batches()
        
        if not batches:
            print("❌ 没有找到可处理的文件（注意：已跳过>100MB文件）")
        else:
            total_files = sum(len(b) for b in batches)
            total_size = sum(sum(os.path.getsize(os.path.join(REPO_ROOT, f)) for f in b) for b in batches) / 1024**3
            print(f"\n📦 总计: {len(batches)} 个批次, {total_files} 个文件, {total_size:.2f}GB")
            generate_scripts(batches)
            
            print("\n🛠 使用说明:")
            print(f"1. 确保已初始化Git仓库: cd {REPO_ROOT} && git init")
            print(f"2. 对大文件执行: git lfs track '*.大文件后缀'")
            print(f"3. 按顺序运行脚本:")
            print(f"   cd {OUTPUT_DIR}")
            print(f"   for script in batch_*.sh; do ./$script; done")
            
    except Exception as e:
        print(f"❌ 错误: {str(e)}")