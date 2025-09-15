---
title: "技术杂记：常见问题与解决方案"
date: "2025-09-01"
tags: [technical-notes, troubleshooting, rdkit, nvidia, cuda, linux, system-administration, development-environment]
---

# 技术杂记

## RDKit安装错误解决

```
  File "/home/gxf1212/data/local-programs/.conda/envs/pytorch/lib/python3.12/site-packages/rdkit/Chem/__init__.py", line 16, in <module>
    from rdkit.Chem import rdchem
ImportError: /lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.31' not found (required by /home/gxf1212/data/local-programs/.conda/envs/pytorch/lib/python3.12/site-packages/rdkit/Chem/../../../../libboost_serialization.so.1.86.0)
```

根据错误信息，您的系统缺少 `libstdc++.so.6` 中的 `GLIBCXX_3.4.31` 版本，而该版本是 RDKit 的依赖项。以下是逐步解决方案：

---

### ✅ **1. 检查当前系统 `libstdc++` 版本**
运行以下命令查看系统中已安装的 `libstdc++` 版本：

```bash
strings /usr/lib/x86_64-linux-gnu/libstdc++.so.6 | grep GLIBCXX
```

如果输出中没有 `GLIBCXX_3.4.31`，说明需要升级 `libstdc++`。

---

### ✅ **2. 安装更高版本的 `libstdc++`**

#### **方法 A: 使用 Conda 安装**
如果您的环境是 Conda 环境，可以直接安装更新的 `libstdcxx-ng`：

```bash
conda install -c conda-forge libstdcxx-ng
```

这会安装包含 `GLIBCXX_3.4.31` 的版本。

---

#### **方法 B: 手动下载并替换库文件**
如果 Conda 无法解决，手动下载并替换 `libstdc++.so.6`：

1. **下载 GCC 12 或更高版本**  
   `GLIBCXX_3.4.31` 需要 GCC 13.1.0 或更高版本。您可以通过以下方式安装 GCC：

   ```bash
   sudo apt update
   sudo apt install g++-12  # 或 g++-13，根据需求选择
   ```

2. **查找新版本的 `libstdc++.so.6` 文件**  
   安装完成后，查找新版本的库文件路径：

   ```bash
   find /usr -name "libstdc++.so*"
   ```

   例如：
   ```
   /usr/lib/gcc/x86_64-linux-gnu/12/libstdc++.so
   ```

3. **备份并替换系统库文件**  
   备份旧文件并创建软链接到新版本：

   ```bash
   sudo mv /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/x86_64-linux-gnu/libstdc++.so.6.bak
   sudo ln -s /usr/lib/gcc/x86_64-linux-gnu/12/libstdc++.so /usr/lib/x86_64-linux-gnu/libstdc++.so.6
   ```

4. **更新动态链接库缓存**  
   ```bash
   sudo ldconfig
   ```

---

### ✅ **3. 验证修复**
运行以下命令检查 `GLIBCXX_3.4.31` 是否已存在：

```bash
strings /usr/lib/x86_64-linux-gnu/libstdc++.so.6 | grep GLIBCXX_3.4.31
```

如果输出中包含 `GLIBCXX_3.4.31`，则修复成功。

---

### ✅ **4. 测试 RDKit 导入**
在 Python 中尝试导入 RDKit 的 Chem 模块：

```python
from rdkit import Chem
print("RDKit Chem imported successfully!")
```

如果不再报错，说明问题已解决。

---

### 📌 **注意事项**
- **备份重要数据**：替换系统库文件前，务必备份关键数据。
- **环境变量**：确保 `LD_LIBRARY_PATH` 指向正确的库路径（例如 Conda 环境的 `lib` 目录）。
- **兼容性**：升级 `libstdc++` 可能影响其他依赖旧版本的应用程序，建议在测试环境中验证。

---

### 🛠️ **附加建议**
如果问题依旧存在，可能是 Conda 环境的路径冲突。尝试以下操作：
1. **强制使用 Conda 的库路径**：
   ```bash
   export LD_LIBRARY_PATH=/home/gxf1212/data/local-programs/.conda/envs/pytorch/lib:$LD_LIBRARY_PATH
   ```
2. **重建 Conda 环境**：
   ```bash
   conda create -n pytorch_new python=3.12
   conda activate pytorch_new
   conda install -c conda-forge rdkit
   ```

---

如果需要进一步帮助，请提供以下信息：
1. `strings /usr/lib/x86_64-linux-gnu/libstdc++.so.6 | grep GLIBCXX` 的输出。
2. 当前系统 GCC 版本：`gcc --version`。
3. Conda 环境的 `LD_LIBRARY_PATH` 内容：`echo $LD_LIBRARY_PATH`。

## Ubuntu系统NVIDIA驱动和CUDA安装指南

在 Ubuntu 上安装 NVIDIA 驱动（570 系列）和 CUDA 的最简便方法如下，假设你只需要命令行环境且无需禁用原有驱动。以下是详细步骤：

---

### **1. 安装 NVIDIA 驱动（570 系列）**
#### **方法一：通过 APT 仓库安装（推荐）**
1. **添加官方驱动仓库**：
   ```bash
   sudo add-apt-repository ppa:graphics-drivers/ppa
   sudo apt update
   ```

2. **安装推荐的驱动版本**：
   ```bash
   sudo ubuntu-drivers autoinstall
   ```
   - 如果需要指定版本为 570，可直接安装：
     ```bash
     sudo apt install nvidia-driver-570
     ```

3. **重启系统**：
   ```bash
   sudo reboot
   ```

4. **验证驱动安装**：
   ```bash
   nvidia-smi
   ```
   - 如果看到 GPU 信息，说明驱动已成功安装。

---

#### **方法二：通过 `.run` 文件安装（手动）**
如果 APT 仓库中没有 570 驱动，可手动下载 `.run` 文件安装：
1. **下载驱动**：
   - 访问 [NVIDIA 驱动下载页面](https://www.nvidia.cn/Download/index.aspx)：https://www.nvidia.cn/Download/index.aspx。
   - 选择你的显卡型号和系统版本，下载对应版本的 `.run` 文件（例如 `NVIDIA-Linux-x86_64-570.xx.run`）。

2. **赋予执行权限**：
   ```bash
   chmod +x NVIDIA-Linux-x86_64-570.xx.run
   ```

3. **安装驱动**：
   ```bash
   sudo ./NVIDIA-Linux-x86_64-570.xx.run --no-opengl-files
   ```
   - 安装过程中，按提示选择选项（通常默认即可）。
   - 如果提示已安装旧驱动，可选择 **继续安装** 或 **卸载旧驱动**。

4. **重启系统**：
   ```bash
   sudo reboot
   ```

5. **验证驱动安装**：
   ```bash
   nvidia-smi
   ```

---

### **2. 安装 CUDA 工具包**
#### **方法一：通过 APT 仓库安装**
1. **添加 NVIDIA CUDA 仓库**：
   ```bash
   sudo apt install curl
   curl -fsSL https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/7fa2af80.pub | sudo gpg --dearmor -o /usr/share/keyrings/cuda-archive-keyring.gpg
   echo "deb [arch=amd64 signed-by=/usr/share/keyrings/cuda-archive-keyring.gpg] https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/ /" | sudo tee /etc/apt/sources.list.d/cuda.list
   sudo apt update
   ```

2. **安装 CUDA 工具包**：
   - 根据驱动版本选择对应的 CUDA 版本（570 驱动通常支持 CUDA 12.5 或 12.6）：
     ```bash
     sudo apt install cuda-toolkit-12-5  # CUDA 12.5
     # 或
     sudo apt install cuda-toolkit-12-6  # CUDA 12.6
     ```

3. **配置环境变量**：
   - 编辑 `~/.bashrc`：
     ```bash
     echo 'export PATH=/usr/local/cuda/bin:$PATH' >> ~/.bashrc
     echo 'export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
     source ~/.bashrc
     ```

4. **验证 CUDA 安装**：
   ```bash
   nvcc --version
   ```

---

#### **方法二：通过 `.run` 文件安装**
1. **下载 CUDA 安装包**：
   - 访问 [CUDA Toolkit 下载页面](https://developer.nvidia.com/cuda-downloads)：https://developer.nvidia.com/cuda-downloads，选择与驱动 570 兼容的版本（如 CUDA 12.5）。
   - 下载 `.run` 文件（例如 `cuda_12.5.0_570.53.1_linux.run`）。

2. **安装 CUDA**：
   ```bash
   chmod +x cuda_12.5.0_570.53.1_linux.run
   sudo ./cuda_12.5.0_570.53.1_linux.run --toolkit --samples --silent --override
   ```
   - 参数说明：
     - `--toolkit`: 安装 CUDA 工具包。
     - `--samples`: 安装示例代码。
     - `--silent`: 静默安装（无需交互）。
     - `--override`: 覆盖旧版本。

3. **配置环境变量**（同上）。

4. **验证 CUDA 安装**：
   ```bash
   nvcc --version
   ```

---

### **3. 注意事项**
- **驱动与 CUDA 版本匹配**：  
  确保安装的 CUDA 版本与 NVIDIA 驱动兼容。例如：
  - CUDA 12.5 → 驱动 570.53.05+
  - CUDA 12.6 → 驱动 570.100.09+

- **仅命令行环境**：  
  如果系统没有 GUI，安装驱动时选择 **不安装 X Server 组件**（如 `.run` 安装时取消勾选）。

- **无需禁用 Nouveau**：  
  如果系统已禁用 Nouveau（默认情况下可能未启用），无需额外操作。若需禁用，可参考知识库中的步骤。

---

### **4. 常见问题**
#### **Q: 安装后黑屏或无法登录？**
- 可能是驱动与系统内核冲突，尝试以下步骤：
  ```bash
  sudo apt purge nvidia*
  sudo reboot
  ```
  然后重新安装驱动。

#### **Q: CUDA 环境变量未生效？**
- 检查 `~/.bashrc` 是否正确配置，并执行：
  ```bash
  source ~/.bashrc
  ```

---

远程连接PyCharm时报错：
```
==== FAILURES ====

The following exception failed the deployment

com.jetbrains.gateway.ssh.deploy.DeployException: 

Details:

An error occurred while executing command: 'host-status --ide-path=/public/home/gxf1212/programs/pycharm-2025.1.2 --project-path=/public/home/gxf1212/data/work/NanoMedicine/pynanomed'

Exit code: 1

	at com.jetbrains.gateway.ssh.DeployFlowUtil$fullDeployCycleImpl$2.invokeSuspend(DeployFlowUtil.kt:308)

	at kotlin.coroutines.jvm.internal.BaseContinuationImpl.resumeWith(ContinuationImpl.kt:33)

	at kotlinx.coroutines.DispatchedTask.run(DispatchedTask.kt:108)

	at kotlinx.coroutines.internal.LimitedDispatcher$Worker.run(LimitedDispatcher.kt:115)

	at kotlinx.coroutines.scheduling.TaskImpl.run(Tasks.kt:103)

	at kotlinx.coroutines.scheduling.CoroutineScheduler.runSafely(CoroutineScheduler.kt:584)

	at kotlinx.coroutines.scheduling.CoroutineScheduler$Worker.executeTask(CoroutineScheduler.kt:793)

	at kotlinx.coroutines.scheduling.CoroutineScheduler$Worker.runWorker(CoroutineScheduler.kt:697)

	at kotlinx.coroutines.scheduling.CoroutineScheduler$Worker.run(CoroutineScheduler.kt:684)

Caused by: com.jetbrains.gateway.ssh.deploy.DeployException: 

Details:

An error occurred while executing command: 'host-status --ide-path=/public/home/gxf1212/programs/pycharm-2025.1.2 --project-path=/public/home/gxf1212/data/work/NanoMedicine/pynanomed'

Exit code: 1

	at com.jetbrains.gateway.ssh.DeployFlowUtil$fullDeployCycleImpl$2.invokeSuspend(DeployFlowUtil.kt:303)

... 8 more

Caused by: com.jetbrains.gateway.ssh.deploy.DeployException: 

Details:

An error occurred while executing command: 'host-status --ide-path=/public/home/gxf1212/programs/pycharm-2025.1.2 --project-path=/public/home/gxf1212/data/work/NanoMedicine/pynanomed'

Exit code: 1

	at com.jetbrains.gateway.ssh.DeployFlowUtil$fullDeployCycleImpl$2.invokeSuspend(DeployFlowUtil.kt:301)

... 8 more

Caused by: com.jetbrains.gateway.ssh.RemoteCommandException: 

Details:

An error occurred while executing command: 'host-status --ide-path=/public/home/gxf1212/programs/pycharm-2025.1.2 --project-path=/public/home/gxf1212/data/work/NanoMedicine/pynanomed'

Exit code: 1

	at com.jetbrains.gateway.ssh.GoHighLevelHostAccessor.createException(GoHighLevelHostAccessor.kt:272)

	at com.jetbrains.gateway.ssh.GoHighLevelHostAccessor.callAndGetError(GoHighLevelHostAccessor.kt:218)

	at com.jetbrains.gateway.ssh.GoHighLevelHostAccessor.access$callAndGetError(GoHighLevelHostAccessor.kt:38)

	at com.jetbrains.gateway.ssh.GoHighLevelHostAccessor$callAndGetError$2.invokeSuspend(GoHighLevelHostAccessor.kt)

... 8 more
```

2. 清理远程IDE后端残留进程

```bash
# 查找PyCharm相关进程

ps aux | grep -E 'pycharm|idea|jetbrains'

# 强制终止（替换PID）

kill -9 <PID>
```






