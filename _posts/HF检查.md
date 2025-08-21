```
import os
from huggingface_hub import constants

# 1. 检查 Python 进程是否能看到 HF_HOME 环境变量
hf_home_from_env = os.getenv('HF_HOME')
print(f"环境变量 'HF_HOME' 的值是: {hf_home_from_env}")

# 2. 获取 huggingface_hub 库实际使用的缓存根目录
#    这个函数会自动遵循 HF_HOME > XDG_CACHE_HOME > ~/.cache/huggingface 的优先级
actual_cache_root = constants.HUGGINGFACE_HUB_CACHE
print(f"Hugging Face 库实际使用的缓存目录: {actual_cache_root}")

# 3. 构造并检查特定模型的完整路径
#    注意：Hugging Face 会将模型名中的 '/' 替换为 '--' 来创建文件夹名
model_folder_name = "models--jinaai--jina-embeddings-v2-base-zh"
model_full_path = os.path.join(actual_cache_root, model_folder_name)
print(f"正在检查模型是否存在于: {model_full_path}")

if os.path.exists(model_full_path):
    print("\n[成功] 模型已在上述路径中找到！")
else:
    print("\n[失败] 在预期的缓存路径中未找到模型。")
    
```