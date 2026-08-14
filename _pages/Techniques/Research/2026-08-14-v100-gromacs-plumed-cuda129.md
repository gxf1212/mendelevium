---
title: "V100 上跑 GROMACS 2026.1 加 PLUMED 遇到的 CUDA 和 MPI 的坑"
date: "2026-08-14"
last_modified_at: 2026-08-14
tags: [gromacs, plumed, cuda, v100, mpi, gpu-computing, openmpi, molecular-dynamics]
description: "在 8 张 V100 的机器上跑 GROMACS 2026.1 加 PLUMED 时踩的 CUDA 和 MPI 的坑以及最终解决方案"
image: "/assets/img/thumbnail_mine/wh-vg7vk8.jpg"
thumbnail: "/assets/img/thumbnail_mine/wh-vg7vk8.jpg"
author: Xufan Gao
lang: zh-CN
---

# V100 上跑 GROMACS 2026.1 加 PLUMED 遇到的 CUDA 和 MPI 的坑

手头这台机器配置不错：2 颗 Intel Xeon Gold 6132（每颗 28 物理核，双路共 56 逻辑核），8 张 Tesla V100-SXM2-32GB，NVIDIA 580.95.05 驱动。

目标是跑通 GROMACS 2026.1 加 PLUMED metadynamics（增强采样，用于自由能计算），让 8 张 V100 真正派上用场。

> take away message: V100不能和cuda13.x一起工作；安装辅助cuda编译Gromacs的方法。

## 1. 现状：纯 CPU 跑得太慢，MPI 又崩

非 MPI 版 GROMACS 2026.1，16 核纯 CPU，速度约 **2.6 ns/day**（含 PLUMED 开销）。跑了半小时到 53 ps，10 ns 需要约 3.8 天。有 8 张 V100 空在那儿不用，太浪费了。

于是尝试走 MPI。用 OpenMPI 5.0.5 编译的 GROMACS 2026.1 加了 PLUMED 之后，测试结果：

| 配置 | 结果 | 速度 |
| ----- | ----- | ----- |
| `mpirun -np 2` 纯 GROMACS 无 PLUMED | 正常 | 约 2.5 ns/day |
| `mpirun -np 2` 加 PLUMED | 崩溃 | `MPI_ERR_REVOKED` |
| 单 rank MPI 加 PLUMED | 正常 | 约 2.6 ns/day |

可以看到，`mpirun -np 2` 的纯 GROMACS 并没有比单 rank 更快，反而略慢一点。这不是 MPI 的问题，而是编译时指定的 CUDA 架构 `75;80;86;89;90` 不包含 V100 的 `sm_70`，GPU 根本没有参与计算（详见第 3 节）。加 PLUMED 的多 rank 直接崩溃：

```
Using 2 MPI processes
Using 8 OpenMP threads per MPI process
+++ Loading the PLUMED kernel runtime +++
+++ PLUMED_KERNEL="/opt/plumed/2.10.0/lib/libplumedKernel.so" +++
[zju-X795-G30:00000] *** An error occurred in MPI_Comm_dup
[zju-X795-G30:00000] *** MPI_ERR_REVOKED: Communication Object Revoked
[zju-X795-G30:00000] *** MPI_ERRORS_ARE_FATAL
```

崩溃发生在 PLUMED kernel 加载后第一次 `MPI_Comm_dup`。`PLUMED_MPI_CTX_DUP=0` 无效，`mpirun --mca mpi_ft_detector false` 也试过，没用。

## 2. 根本原因排查，一度被误判成 PLUMED 已知 bug

最初一度以为这是 PLUMED 的已知问题，跟编译无关。回头看，虽然 PLUMED GitHub 上确实有相关 issue，但那不是唯一的可能。

查 GROMACS 官方文档才发现，**PLUMED 不支持的是多个 thread-MPI rank**，并没有说 external MPI 不能用。`mpirun -np 2 gmx_mpi ... -plumed` 本来就是应该能工作的。

参考：[GROMACS 2026.1 PLUMED 文档](https://manual.gromacs.org/documentation/2026.1/reference-manual/special/plumed.html)：https://manual.gromacs.org/documentation/2026.1/reference-manual/special/plumed.html

**PLUMED 和 GROMACS 必须用同一套 MPI 库**。检查方法是分别用 `ldd` 查看两边链接了哪套 `libmpi`：

```bash
ldd /opt/gromacs/2026.1-mpi/bin/gmx_mpi | grep -E 'libmpi|pmix'
ldd /opt/plumed/2.10.0/lib/libplumedKernel.so | grep -E 'libmpi|pmix'
```

如果 PLUMED 不是用 `/opt/openmpi/5.0.5/` 这套编译的，需要重新编译。原来的 PLUMED 安装在 `/opt/plumed/2.10.0/`，重新编译时改路径加 `-ompi505` 后缀，避免覆盖原安装：

```bash
./configure --prefix=/opt/plumed/2.10.0-ompi505 \
  CXX=/opt/openmpi/5.0.5/bin/mpicxx
make -j16
make install
```

## 3. CUDA 13 不支持 V100

除了 MPI 之外，CUDA 版本也有问题。CUDA 13 已经删除了对 `compute_70` 架构（Volta）的支持，V100 就是 `sm_70`。用 CUDA 13 的 nvcc 加 `70` 编译时直接报错：

```
nvcc fatal: Unsupported gpu architecture 'compute_70'
```

`nvidia-smi` 显示 `CUDA Version: 13.0` 只是说明 580.95.05 驱动最高支持 CUDA 13 运行时，不是说 V100 能用 CUDA 13 编译。CUDA 12.9 仍明确支持 `sm_70`。

参考 [CUDA 13 发布说明](https://docs.nvidia.com/cuda/archive/13.0.3/cuda-toolkit-release-notes/index.html)：https://docs.nvidia.com/cuda/archive/13.0.3/cuda-toolkit-release-notes/index.html

## 4. 装 CUDA 12.9 踩坑，cuda-keyring 不生成源文件

按 NVIDIA 官方文档执行以下命令：

```bash
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt update
apt-cache policy cuda-toolkit-12-9
```

结果只有 `Candidate: (none)`，什么都没有。

原因是 `cuda-keyring` 装上去之后并没有生成 CUDA 仓库的 source 文件，`/etc/apt/sources.list.d/` 目录里根本没有 CUDA 源。

解决办法是手工加源：

```bash
echo "deb [signed-by=/usr/share/keyrings/cuda-archive-keyring.gpg] https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/ /" \
  > /etc/apt/sources.list.d/cuda-ubuntu2204-x86_64.list

apt update
```

国内网络环境下，可以把仓库地址里的 `nvidia.com` 改成 `nvidia.cn`：

```bash
sed -i 's#developer.download.nvidia.com#developer.download.nvidia.cn#' \
  /etc/apt/sources.list.d/cuda-ubuntu2204-x86_64.list
```

加完源之后，装 `cuda-toolkit-12-9` 仍然失败。原因是系统上已经装了 CUDA 13 的 `cuda-toolkit-config-common`，跟 12.9 对应的 config 包冲突，APT 报 `held broken packages`。这其实不是真有包被 `apt-mark hold`，而是依赖解析根本解不出来。

## 5. 最终方案，用 NVIDIA runfile 装 12.9

与其继续跟 APT 的依赖解析纠缠，不如绕过它。直接用 NVIDIA 官方 12.9.1 runfile：

```bash
cd /root
wget https://developer.download.nvidia.com/compute/cuda/12.9.1/local_installers/cuda_12.9.1_575.57.08_linux.run
chmod +x cuda_12.9.1_575.57.08_linux.run

./cuda_12.9.1_575.57.08_linux.run \
  --silent \
  --toolkit \
  --toolkitpath=/usr/local/cuda-12.9
```

`--toolkit` 参数**只装 Toolkit 部分，不动 580.95.05 驱动**。不加会连 575.57.08 驱动一起装。两套并存：

```
/usr/local/cuda-12.9/
/usr/local/cuda-13.0/
```

验证：

```bash
/usr/local/cuda-12.9/bin/nvcc --version
/usr/local/cuda-12.9/bin/nvcc --list-gpu-arch | grep 70
```

## 6. 编译 GROMACS，干净目录加 `sm_70`

```bash
cd ~/gromacs-2026.1-mpi
rm -rf build
mkdir build && cd build
```

```bash
export PATH=/opt/openmpi/5.0.5/bin:$PATH
export LD_LIBRARY_PATH=/opt/openmpi/5.0.5/lib:$LD_LIBRARY_PATH

/opt/cmake/4.2.1/bin/cmake .. \
  -DGMX_FFT_LIBRARY=fftw3 \
  -DFFTWF_LIBRARY=/opt/fftw3/3.3.10/float-sse2-avx2/lib/libfftw3f.so \
  -DFFTWF_INCLUDE_DIR=/opt/fftw3/3.3.10/float-sse2-avx2/include \
  -DGMX_GPU=CUDA \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.9/bin/nvcc \
  -DCUDAToolkit_ROOT=/usr/local/cuda-12.9 \
  -DCMAKE_CUDA_ARCHITECTURES=70 \
  -DGMX_MPI=ON \
  -DGMX_USE_PLUMED=ON \
  -DGMX_USE_COLVARS=internal \
  -DCMAKE_INSTALL_PREFIX=/opt/gromacs/2026.1-mpi \
  -DCMAKE_CXX_COMPILER=g++ \
  -DCMAKE_C_COMPILER=gcc \
  -DMPI_CXX_COMPILER=/opt/openmpi/5.0.5/bin/mpicxx \
  -DMPI_C_COMPILER=/opt/openmpi/5.0.5/bin/mpicc
```

编译时有几个注意点：

- **只编译 `sm_70` 架构就够**：没必要编译 `75;80;86;89;90`，以前加这些是为了兼容其他 GPU，这里没有。
- `--toolkitpath` 和 `-DCMAKE_CUDA_COMPILER` 的**路径必须一致**：防止 nvcc 来自 12.9 而头文件和库却从 CUDA 13 拿这种混装情况。
- **PLUMED 必须和 GROMACS 用同一套 OpenMPI**：两者编译时必须指向同一个 `mpicxx`，否则多 rank 时 `MPI_Comm_dup` 会崩溃。
- 不要碰 CUDA 软链接，也不要改 580.95.05 驱动，**让两套 CUDA 并存才安全**。

## 经验总结

> 不要拿“已知 bug”搪塞，**先做基础检查**。出问题时先问：PLUMED 和 GROMACS 用的是同一套 MPI 库吗？编译用的 nvcc 真的支持这台 GPU 吗？
>
> V100 上跑 GROMACS 2026.1，CUDA 13 已经不支持 `sm_70`。装 CUDA 12.9 最简单的方式是用 NVIDIA runfile 的 `--toolkit --toolkitpath`，不动现有驱动和软链接。
>
> `cuda-keyring` 不会自动添加 APT source，必须手工 echo 到 `/etc/apt/sources.list.d/`。

修好之后，下一步就是 benchmark。56 核配 8 张 V100，可以先试 `2 rank × 14 OMP` 和 `4 rank × 7 OMP` 这两组配置，看哪个效率更高。

## 参考链接

- [GROMACS 2026.1 PLUMED 文档](https://manual.gromacs.org/documentation/2026.1/reference-manual/special/plumed.html)：https://manual.gromacs.org/documentation/2026.1/reference-manual/special/plumed.html
- [CUDA 13.0.3 发布说明](https://docs.nvidia.com/cuda/archive/13.0.3/cuda-toolkit-release-notes/index.html)：https://docs.nvidia.com/cuda/archive/13.0.3/cuda-toolkit-release-notes/index.html
- [CUDA 版本兼容性说明](https://docs.nvidia.com/deploy/cuda-compatibility/latest/why-cuda-compatibility.html)：https://docs.nvidia.com/deploy/cuda-compatibility/latest/why-cuda-compatibility.html
- [PLUMED 2.10 安装文档](https://www.plumed.org/doc-v2.10/user-doc/html/_installation.html)：https://www.plumed.org/doc-v2.10/user-doc/html/_installation.html
- [CUDA 12.9.1 runfile 下载](https://developer.download.nvidia.com/compute/cuda/12.9.1/local_installers/)：https://developer.download.nvidia.com/compute/cuda/12.9.1/local_installers/
- [OpenMPI ULFM 文档](https://docs.open-mpi.org/en/v5.0.x/features/ulfm.html)（`MPI_ERR_REVOKED` 含义）：https://docs.open-mpi.org/en/v5.0.x/features/ulfm.html