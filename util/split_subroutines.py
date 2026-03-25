import re
import os

source_file = "velest.f"
output_dir = "subroutines"
os.makedirs(output_dir, exist_ok=True)

# 读取源文件
with open(source_file, "r") as f:
    lines = f.readlines()

# 1️⃣ 找到主逻辑结束的第一个 "      END" 或 "      end"
# 找到主逻辑结束的第一个真正的 END
main_end_idx = None
for idx, line in enumerate(lines):
    if re.match(r'^\s{6}(END|end)\s*$', line):
        main_end_idx = idx
        break

if main_end_idx is None:
    raise ValueError("没有找到主逻辑结束的 END 行")

# 2️⃣ 找到所有子程序起始行
sub_start_indices = []
sub_names = []
for idx, line in enumerate(lines[main_end_idx+1:], start=main_end_idx+1):
    if line.startswith("      subroutine"):
        sub_start_indices.append(idx)
        # 提取子程序名字
        match = re.match(r" {6}(?i:subroutine)\s+(\w+)", line)
        if match:
            sub_names.append(match.group(1))
        else:
            sub_names.append(f"sub_{len(sub_start_indices)}")  # 防止匹配失败

# 3️⃣ 生成每个子程序文件
sub_start_indices.append(len(lines))  # 方便计算最后一个子程序的结束行
for i in range(len(sub_start_indices)-1):
    start = sub_start_indices[i]
    end_search_range = lines[start:sub_start_indices[i+1]]  # 下一sub开始之前
    # 找最后一个 "      END" 或 "      end" 作为子程序结尾
    end = None
    for idx_offset, line in enumerate(end_search_range):
        if line.startswith("      END") or line.startswith("      end"):
            end = start + idx_offset + 1  # +1 包括结束行
    if end is None:
        end = sub_start_indices[i+1]  # 万一没找到 END，就取下一sub开始前
    sub_file = os.path.join(output_dir, f"{sub_names[i].lower()}.f")
    with open(sub_file, "w") as f:
        f.writelines(["c\n"])
        f.writelines(["c\n"])
        f.writelines(["c\n"])
        f.writelines(["c\n"])
        f.writelines(["c\n"])
        f.writelines(["c\n"])
        f.writelines(lines[start:end])
    print(f"Saved subroutine '{sub_names[i]}' to '{sub_file}'")

# 4️⃣ 更新原文件，只保留主逻辑
with open(source_file, "w") as f:
    f.writelines(lines[:main_end_idx+1])

print(f"Updated '{source_file}' to keep main logic only.")