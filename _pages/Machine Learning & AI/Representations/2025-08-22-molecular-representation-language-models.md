---
title: "Comprehensive Guide to Molecular Representation Language Models: From Proteins to Small Molecules"
date: "2025-06-15"
tags: [molecular-representation, language-models, protein-modeling, small-molecules, transformers, deep-learning, ai-drug-discovery]
---
# 分子表示学习模型全览：从蛋白质到小分子的语言模型

**分子表示学习已成为计算化学和生物信息学的核心技术**。随着Transformer架构在自然语言处理中的成功，研究者们将其应用到分子数据的表示学习中，取得了显著进展。本文全面介绍从蛋白质到小分子的各种语言模型，为读者提供完整的技术栈和实用代码。

## 环境配置

### 基础依赖安装

```bash
# PyTorch安装（根据CUDA版本调整）
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu126

# HuggingFace Transformers
pip install transformers

# 检查GPU可用性
python -c "import torch; print(f'CUDA Available: {torch.cuda.is_available()}'); print(f'GPU Count: {torch.cuda.device_count()}')"
```

### 可选：设置模型缓存路径

```python
import os
os.environ['TORCH_HOME'] = '/your/path/to/model'
os.environ['HF_HOME'] = '/your/path/to/hf_model'
```

## 一、蛋白质语言模型

### 1.1 ESM-2系列

#### 模型简介
ESM-2（Evolutionary Scale Modeling）是Meta开发的大规模蛋白质语言模型[1]，**在进化规模的蛋白质序列数据上进行预训练，能够捕获蛋白质的进化和结构信息**。

#### 可用模型规模

| 模型名称 | 层数 | 参数量 | 模型大小 |
|----------|------|--------|----------|
| esm2_t48_15B_UR50D | 48 | 15B | ~60GB |
| esm2_t36_3B_UR50D | 36 | 3B | ~12GB |
| esm2_t33_650M_UR50D | 33 | 650M | 2.5GB |
| esm2_t30_150M_UR50D | 30 | 150M | ~600MB |
| esm2_t12_35M_UR50D | 12 | 35M | ~140MB |
| esm2_t6_8M_UR50D | 6 | 8M | ~32MB |

#### 安装和使用

```bash
pip install fair-esm
```

```python
import torch
import esm

# 检查GPU
print("Number of GPUs:", torch.cuda.device_count())

# 加载模型（选择合适的规模）
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # 禁用dropout以获得确定性结果

# 如果有GPU，移动到GPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = model.to(device)

# 准备序列数据
data = [
    ("protein1", "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"),
    ("protein2", "KALTARQQEVFDLIRDHISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSGASRGIRLLQEE"),
]

# 批量转换
batch_labels, batch_strs, batch_tokens = batch_converter(data)
batch_tokens = batch_tokens.to(device)
batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

# 提取表示
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=True)

# 获取token表示（每个氨基酸的embedding）
token_representations = results["representations"][33]

# 获取序列级表示（整个蛋白质的embedding）
sequence_representations = []
for i, tokens_len in enumerate(batch_lens):
    # 移除特殊token（开始和结束）
    seq_repr = token_representations[i, 1 : tokens_len - 1].mean(0)
    sequence_representations.append(seq_repr)

print(f"Token representation shape: {token_representations.shape}")
print(f"Sequence representation shape: {sequence_representations[0].shape}")
```

#### 高级用法：注意力权重和接触预测

```python
# 获取注意力权重和接触预测
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=True)

# 接触预测（用于蛋白质结构预测）
contacts = results["contacts"]
print(f"Contacts shape: {contacts.shape}")

# 注意力权重
attentions = results["attentions"]
print(f"Attention shape: {attentions.shape}")
```

### 1.2 ESM-C (ESM Cambrian)

#### 模型简介
ESM-C是ESM3模型家族中专注于表示学习的平行模型[2]，**相比ESM-2在相同参数量下提供更高效的性能和更低的内存消耗**。ESM-C设计为ESM-2的直接替代品，具有重大性能优势。

#### 性能对比

| ESM-C参数量 | 对应ESM-2参数量 | ESM-C优势 |
|-------------|------------------|-----------|
| 300M | 650M | 更低内存消耗，更快推理 |
| 600M | 3B | 高效达到甚至超越更大规模ESM-2性能 |
| 6B | - | 性能远超最佳ESM-2模型 |

#### 安装和使用

```bash
pip install esm
```

#### 方法一：使用ESM SDK API（推荐）

```python
from esm.sdk.api import ESMProtein, LogitsConfig
from esm.models.esmc import ESMC

# 创建蛋白质对象
protein = ESMProtein(sequence="MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG")

# 加载模型（如果遇到tokenizer错误，使用方法二）
try:
    client = ESMC.from_pretrained("esmc_600m").to("cuda")  # 或 "cpu"
    
    # 编码蛋白质
    protein_tensor = client.encode(protein)
    
    # 获取logits和embeddings
    logits_output = client.logits(
        protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
    )
    
    print(f"Logits shape: {logits_output.logits.sequence.shape}")
    print(f"Embeddings shape: {logits_output.embeddings.shape}")
    
    # 提取序列级表示
    sequence_embedding = logits_output.embeddings.mean(dim=1)  # 平均池化
    print(f"Sequence embedding shape: {sequence_embedding.shape}")
    
except AttributeError as e:
    print(f"ESM-C错误: {e}")
    print("请使用方法二或方法三")
```

If you see

```
ESM-C错误: property 'cls_token' of 'EsmSequenceTokenizer' object has no setter
```

please do this according to https://github.com/evolutionaryscale/esm/issues/214

```bash
pip install esm==3.1.1
```

The output is like

```
Logits shape: torch.Size([1, 67, 64])
Embeddings shape: torch.Size([1, 67, 1152])
Sequence embedding shape: torch.Size([1, 1152])
```

#### 方法二：使用远程API（需要注册）

```python
from esm.sdk.forge import ESM3ForgeInferenceClient
from esm.sdk.api import ESMProtein, LogitsConfig

# 需要先在 https://forge.evolutionaryscale.ai 注册获取token
forge_client = ESM3ForgeInferenceClient(
    model="esmc-6b-2024-12", 
    url="https://forge.evolutionaryscale.ai", 
    token="<your_forge_token>"
)

protein = ESMProtein(sequence="MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG")
protein_tensor = forge_client.encode(protein)
logits_output = forge_client.logits(
    protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
)

print(f"Remote embeddings shape: {logits_output.embeddings.shape}")
```

### 1.3 CARP

#### 模型简介
CARP（Contrastive Autoregressive Protein model）是微软开发的蛋白质语言模型[3]，**采用对比学习和自回归训练目标，在蛋白质序列建模方面表现优异**。

#### 安装和使用

**在线安装**：
```bash
pip install git+https://github.com/microsoft/protein-sequence-models.git
```

**离线安装**：
1. 下载仓库：https://github.com/microsoft/protein-sequence-models
2. 解压并安装：
```bash
cd /path/to/protein-sequence-models
pip install .
```

#### 代码实现

```python
from sequence_models.pretrained import load_model_and_alphabet

# 加载模型和序列处理器
model, collater = load_model_and_alphabet('carp_640M')

# 准备序列数据（注意：需要嵌套列表格式）
seqs = [['MDREQ'], ['MGTRRLLP'], ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG']]

# 将序列转换为模型输入格式
x = collater(seqs)[0]  # (n, max_len)

# 获取表示（第56层的表示）
with torch.no_grad():
    rep = model(x)['representations'][56]  # (n, max_len, d_model)

print(f"Input shape: {x.shape}")
print(f"Representation shape: {rep.shape}")

# 获取序列级表示（平均池化）
sequence_repr = rep.mean(dim=1)
print(f"Sequence representation shape: {sequence_repr.shape}")
```

### 1.4 ProtT5

#### 模型简介
ProtT5是基于T5架构的蛋白质语言模型[4]，**采用编码器-解码器结构，在大规模蛋白质数据上预训练，支持多种下游任务**。

#### 从本地路径加载模型

```python
import torch
import re
from transformers import T5Tokenizer, T5EncoderModel

# 设备配置
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

# 本地模型路径（如果已下载）
tokenizer_path = '/your/path/to/prot_t5_xl_half_uniref50-enc/'

# 加载tokenizer和模型
try:
    tokenizer = T5Tokenizer.from_pretrained(tokenizer_path, do_lower_case=False)
    print(f"Tokenizer loaded from local path: {tokenizer_path}")
except OSError:
    # 如果本地路径不存在，从HuggingFace下载
    tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_half_uniref50-enc', do_lower_case=False)
    print("Tokenizer loaded from HuggingFace")

# 加载模型
model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc").to(device)

# 示例蛋白质序列
sequence_examples = ["PRTEINO", "SEQWENCE"]

# 预处理：替换稀有氨基酸，添加空格
sequence_examples = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in sequence_examples]

# Tokenization
ids = tokenizer(sequence_examples, add_special_tokens=True, padding="longest", return_tensors="pt")
input_ids = ids['input_ids'].to(device)
attention_mask = ids['attention_mask'].to(device)

# 生成embeddings
with torch.no_grad():
    embedding_repr = model(input_ids=input_ids, attention_mask=attention_mask)

# 提取每个序列的残基embeddings
emb_0 = embedding_repr.last_hidden_state[0, :7]  # 第一个序列
emb_1 = embedding_repr.last_hidden_state[1, :8]  # 第二个序列

print("Shape of embedding for sequence 1:", emb_0.shape)
print("Shape of embedding for sequence 2:", emb_1.shape)
print("Protein embeddings generated successfully!")
```

### 1.5 Ankh

#### 模型简介
Ankh是专门为阿拉伯语蛋白质序列优化的多语言蛋白质模型[5]，**基于T5架构，支持多种语言和蛋白质表示任务**。

#### 实现代码

```python
from transformers import AutoTokenizer, AutoModelForSeq2SeqLM
import torch

# 本地模型路径
local_model_path = "/your/path/to/ankh-large/"

# 加载tokenizer和模型
tokenizer = AutoTokenizer.from_pretrained(local_model_path)
model = AutoModelForSeq2SeqLM.from_pretrained(local_model_path)

# 示例序列
sequence_examples = ["MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"]
inputs = tokenizer(sequence_examples, return_tensors="pt", padding=True)

# 设备配置
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = model.to(device)
inputs = {key: value.to(device) for key, value in inputs.items()}

# 生成编码器embeddings
with torch.no_grad():
    encoder_outputs = model.encoder(**inputs)
    embeddings = encoder_outputs.last_hidden_state

# 提取有效序列的embeddings（移除padding）
emb_0 = embeddings[0, :inputs['attention_mask'][0].sum()]

print("Shape of encoder embeddings for sequence 1:", emb_0.shape)
print("Model loaded successfully from:", local_model_path)
```

## 二、肽语言模型

### 2.1 PepBERT

#### 模型简介
PepBERT是专门为肽序列设计的BERT模型[6]，**针对短肽序列进行优化，在肽-蛋白质相互作用预测等任务中表现优异**。

#### 模型特点
- 专门针对肽序列（通常长度较短）
- 基于BERT架构，采用掩码语言建模
- 在UniParc数据库的大规模肽序列上预训练
- 输出维度：320

#### 安装和使用

```python
import os
import torch
import importlib.util
from tokenizers import Tokenizer

# 设置环境变量
os.environ['TORCH_HOME'] = '/home/gxf1212/data/local-programs/model' 
os.environ['HF_HOME'] = '/home/gxf1212/data/local-programs/hf_model'

# 本地模型路径
snapshot_path = "/home/gxf1212/data/local-programs/hf_model/hub/models--dzjxzyd--PepBERT-large-UniParc/snapshots/7b0cbb2f925d05c9fca42c63c1712f94200fdb41" 

def load_module_from_local(file_path):
    """从本地文件加载Python模块"""
    module_name = os.path.splitext(os.path.basename(file_path))[0]
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

# 1) 动态加载模型配置
model_module = load_module_from_local(os.path.join(snapshot_path, "model.py"))
config_module = load_module_from_local(os.path.join(snapshot_path, "config.py"))
build_transformer = model_module.build_transformer
get_config = config_module.get_config

# 2) 加载tokenizer
tokenizer_path = os.path.join(snapshot_path, "tokenizer.json")
tokenizer = Tokenizer.from_file(tokenizer_path)

# 3) 加载模型权重
weights_path = os.path.join(snapshot_path, "tmodel_17.pt")

# 4) 初始化模型
device = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_built() else "cpu"
config = get_config()
model = build_transformer(
    src_vocab_size=tokenizer.get_vocab_size(),
    src_seq_len=config["seq_len"],
    d_model=config["d_model"]
)

# 加载预训练权重
state = torch.load(weights_path, map_location=torch.device(device))
model.load_state_dict(state["model_state_dict"])
model.eval()

# 5) 生成embeddings
def get_peptide_embedding(sequence):
    """生成肽序列的embedding"""
    # 添加特殊token [SOS] 和 [EOS]
    encoded_ids = (
        [tokenizer.token_to_id("[SOS]")]
        + tokenizer.encode(sequence).ids
        + [tokenizer.token_to_id("[EOS]")]
    )
    
    input_ids = torch.tensor([encoded_ids], dtype=torch.int64)
    
    with torch.no_grad():
        # 创建注意力掩码
        encoder_mask = torch.ones((1, 1, 1, input_ids.size(1)), dtype=torch.int64)
        
        # 前向传播获取token embeddings
        emb = model.encode(input_ids, encoder_mask)
        
        # 移除特殊token的embeddings
        emb_no_special = emb[:, 1:-1, :]
        
        # 平均池化获取序列级表示
        emb_avg = emb_no_special.mean(dim=1)
    
    return emb_avg

# 使用示例
sequence = "KRKGFLGI"
embedding = get_peptide_embedding(sequence)
print("Shape of peptide embedding:", embedding.shape)  # (1, 320)
print("Peptide embedding generated successfully!")
```

## 三、小分子语言模型

### 3.1 ChemBERTa系列

#### 模型简介
ChemBERTa是首个大规模的分子BERT模型[7]，**在7700万PubChem分子上预训练，采用掩码语言建模目标，为分子性质预测提供强大的预训练表示**。

#### 主要版本
- **ChemBERTa-77M-MLM**: 在77M分子上用掩码语言建模预训练
- **ChemBERTa-2**: 改进版本，支持多任务预训练
- **参数量**: 约12M-77M参数

#### 安装和使用

```bash
# 安装依赖
pip install transformers torch rdkit
```

```python
from transformers import AutoTokenizer, AutoModel
import torch
from rdkit import Chem

# 加载预训练模型
model_name = "DeepChem/ChemBERTa-77M-MLM"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModel.from_pretrained(model_name)

# 设备配置
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = model.to(device)

def get_molecular_embedding(smiles_list):
    """获取分子的ChemBERTa embedding"""
    # Tokenization
    inputs = tokenizer(smiles_list, return_tensors="pt", padding=True, truncation=True, max_length=512)
    inputs = {key: value.to(device) for key, value in inputs.items()}
    
    # 前向传播
    with torch.no_grad():
        outputs = model(**inputs)
        # 使用[CLS] token的表示作为分子级表示
        molecular_embeddings = outputs.last_hidden_state[:, 0, :]  # [CLS] token
        # 或者使用平均池化
        # molecular_embeddings = outputs.last_hidden_state.mean(dim=1)
    
    return molecular_embeddings

# 使用示例
smiles_examples = [
    "CCO",  # 乙醇
    "CC(=O)O",  # 乙酸
    "c1ccccc1",  # 苯
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # 咖啡因
]

# 验证SMILES有效性
valid_smiles = []
for smi in smiles_examples:
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        valid_smiles.append(smi)
    else:
        print(f"Invalid SMILES: {smi}")

# 生成embeddings
embeddings = get_molecular_embedding(valid_smiles)
print(f"Generated embeddings shape: {embeddings.shape}")
print(f"Embedding dimension: {embeddings.shape[1]}")

# 单个分子的embedding
single_embedding = get_molecular_embedding(["CCO"])
print(f"Single molecule embedding shape: {single_embedding.shape}")
```

#### 高级用法：微调ChemBERTa

```python
from transformers import AutoTokenizer, AutoModelForSequenceClassification, TrainingArguments, Trainer
import torch.nn as nn

# 加载用于分类任务的模型
model = AutoModelForSequenceClassification.from_pretrained(
    "DeepChem/ChemBERTa-77M-MLM",
    num_labels=2  # 二分类任务
)

# 准备数据集和训练参数
class MolecularDataset(torch.utils.data.Dataset):
    def __init__(self, smiles_list, labels, tokenizer, max_length=512):
        self.smiles_list = smiles_list
        self.labels = labels
        self.tokenizer = tokenizer
        self.max_length = max_length
    
    def __len__(self):
        return len(self.smiles_list)
    
    def __getitem__(self, idx):
        smiles = self.smiles_list[idx]
        label = self.labels[idx]
        
        encoding = self.tokenizer(
            smiles,
            truncation=True,
            padding='max_length',
            max_length=self.max_length,
            return_tensors='pt'
        )
        
        return {
            'input_ids': encoding['input_ids'].flatten(),
            'attention_mask': encoding['attention_mask'].flatten(),
            'labels': torch.tensor(label, dtype=torch.long)
        }

# 微调代码示例（需要准备训练数据）
# training_args = TrainingArguments(
#     output_dir='./results',
#     num_train_epochs=3,
#     per_device_train_batch_size=16,
#     per_device_eval_batch_size=64,
#     warmup_steps=500,
#     weight_decay=0.01,
#     logging_dir='./logs',
# )
```

### 3.2 MolFormer系列

#### 模型简介
MolFormer是IBM开发的大规模化学语言模型[8]，**在11亿分子上预训练，采用线性注意力机制和旋转位置编码，在多个分子性质预测任务上达到SOTA性能**。

#### 模型特点
- **预训练数据**: 11亿分子（PubChem + ZINC）
- **架构**: 线性注意力Transformer + 旋转位置编码
- **高效性**: 线性时间复杂度，支持长序列
- **性能**: 在多个基准数据集上超越GNN模型

#### 安装和使用

```bash
git clone https://github.com/IBM/molformer.git
cd molformer
pip install -e .
```

```python
import torch
from molformer.models import MolFormer
from molformer.tokenizer import MolTranBertTokenizer

# 加载预训练模型和tokenizer
model_path = "ibm/MoLFormer-XL-both-10pct"  # HuggingFace模型路径
tokenizer = MolTranBertTokenizer.from_pretrained(model_path)
model = MolFormer.from_pretrained(model_path)

# 设备配置
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = model.to(device)
model.eval()

def get_molformer_embedding(smiles_list, max_length=512):
    """获取MolFormer分子embedding"""
    # Tokenization
    encoded = tokenizer(
        smiles_list,
        padding=True,
        truncation=True,
        max_length=max_length,
        return_tensors="pt"
    )
    
    # 移动到设备
    input_ids = encoded['input_ids'].to(device)
    attention_mask = encoded['attention_mask'].to(device)
    
    # 前向传播
    with torch.no_grad():
        outputs = model(input_ids=input_ids, attention_mask=attention_mask)
        # 使用最后一层的隐藏状态
        hidden_states = outputs.last_hidden_state
        
        # 计算分子级表示（掩码平均池化）
        mask_expanded = attention_mask.unsqueeze(-1).expand(hidden_states.size()).float()
        sum_embeddings = torch.sum(hidden_states * mask_expanded, 1)
        sum_mask = torch.clamp(mask_expanded.sum(1), min=1e-9)
        molecular_embeddings = sum_embeddings / sum_mask
    
    return molecular_embeddings

# 使用示例
smiles_examples = [
    "CCO",
    "CC(=O)O", 
    "c1ccccc1",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
]

embeddings = get_molformer_embedding(smiles_examples)
print(f"MolFormer embeddings shape: {embeddings.shape}")
print(f"Embedding dimension: {embeddings.shape[1]}")
```

#### MolFormer-XL超大规模版本

```python
# 对于MolFormer-XL（需要更多内存）
model_xl_path = "ibm/MoLFormer-XL-both-10pct"
tokenizer_xl = MolTranBertTokenizer.from_pretrained(model_xl_path)
model_xl = MolFormer.from_pretrained(model_xl_path)

# 使用混合精度以节省内存
model_xl = model_xl.half().to(device)  # 使用半精度

# 对于大批量处理，建议分批处理
def batch_process_molecules(smiles_list, batch_size=32):
    """分批处理大量分子"""
    all_embeddings = []
    
    for i in range(0, len(smiles_list), batch_size):
        batch = smiles_list[i:i+batch_size]
        embeddings = get_molformer_embedding(batch)
        all_embeddings.append(embeddings.cpu())
        
        # 清理GPU缓存
        torch.cuda.empty_cache()
    
    return torch.cat(all_embeddings, dim=0)
```

### 3.3 SMILES Transformer

#### 模型简介
SMILES Transformer是首个专门为SMILES序列设计的Transformer模型[9]，**采用自编码任务进行预训练，学习分子的潜在表示，适用于低数据量的药物发现任务**。

#### 特点
- **预训练任务**: 自编码（去噪自编码器）
- **数据**: 170万ChEMBL分子（不超过100字符）
- **SMILES增强**: 使用SMILES枚举增加数据多样性
- **应用**: 低数据药物发现

#### 安装和使用

```bash
git clone https://github.com/DSPsleeporg/smiles-transformer.git
cd smiles-transformer
pip install -r requirements.txt
```

```python
import torch
import torch.nn as nn
from torch.nn import Transformer
import numpy as np
from rdkit import Chem

class SMILESTransformer(nn.Module):
    """SMILES Transformer模型"""
    def __init__(self, vocab_size, d_model=512, nhead=8, num_layers=6, max_seq_len=100):
        super(SMILESTransformer, self).__init__()
        self.d_model = d_model
        self.embedding = nn.Embedding(vocab_size, d_model)
        self.pos_encoder = PositionalEncoding(d_model, max_seq_len)
        self.transformer = Transformer(
            d_model=d_model,
            nhead=nhead,
            num_encoder_layers=num_layers,
            num_decoder_layers=num_layers,
            dim_feedforward=2048,
            dropout=0.1
        )
        self.fc_out = nn.Linear(d_model, vocab_size)
        
    def forward(self, src, tgt=None, src_mask=None, tgt_mask=None):
        # 编码器
        src_emb = self.pos_encoder(self.embedding(src) * np.sqrt(self.d_model))
        
        if tgt is not None:
            # 训练模式（编码器-解码器）
            tgt_emb = self.pos_encoder(self.embedding(tgt) * np.sqrt(self.d_model))
            output = self.transformer(src_emb, tgt_emb, src_mask=src_mask, tgt_mask=tgt_mask)
            return self.fc_out(output)
        else:
            # 推理模式（仅编码器）
            memory = self.transformer.encoder(src_emb, src_mask)
            return memory

class PositionalEncoding(nn.Module):
    """位置编码"""
    def __init__(self, d_model, max_len=100):
        super(PositionalEncoding, self).__init__()
        
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * 
                           (-np.log(10000.0) / d_model))
        
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0).transpose(0, 1)
        
        self.register_buffer('pe', pe)
        
    def forward(self, x):
        return x + self.pe[:x.size(0), :]

class SMILESTokenizer:
    """SMILES分词器"""
    def __init__(self):
        # 基础SMILES字符集
        self.chars = ['<PAD>', '<SOS>', '<EOS>', '<UNK>'] + list("()[]1234567890=+-#@CNOSPFIBrClcnos")
        self.char_to_idx = {char: idx for idx, char in enumerate(self.chars)}
        self.idx_to_char = {idx: char for char, idx in self.char_to_idx.items()}
        self.vocab_size = len(self.chars)
    
    def encode(self, smiles, max_length=100):
        """编码SMILES字符串"""
        tokens = ['<SOS>'] + list(smiles) + ['<EOS>']
        indices = [self.char_to_idx.get(token, self.char_to_idx['<UNK>']) for token in tokens]
        
        # 填充或截断
        if len(indices) < max_length:
            indices += [self.char_to_idx['<PAD>']] * (max_length - len(indices))
        else:
            indices = indices[:max_length]
            
        return torch.tensor(indices, dtype=torch.long)
    
    def decode(self, indices):
        """解码回SMILES字符串"""
        chars = [self.idx_to_char[idx.item()] for idx in indices]
        # 移除特殊token
        chars = [c for c in chars if c not in ['<PAD>', '<SOS>', '<EOS>', '<UNK>']]
        return ''.join(chars)

def get_smiles_embedding(smiles_list, model, tokenizer, device):
    """获取SMILES的分子embedding"""
    model.eval()
    embeddings = []
    
    with torch.no_grad():
        for smiles in smiles_list:
            # 编码SMILES
            encoded = tokenizer.encode(smiles).unsqueeze(0).to(device)
            
            # 获取编码器输出
            encoder_output = model(encoded)
            
            # 平均池化获取分子级表示
            # 忽略padding token
            mask = (encoded != tokenizer.char_to_idx['<PAD>']).float()
            pooled = (encoder_output * mask.unsqueeze(-1)).sum(dim=1) / mask.sum(dim=1, keepdim=True)
            
            embeddings.append(pooled)
    
    return torch.cat(embeddings, dim=0)

# 使用示例
def demo_smiles_transformer():
    """演示SMILES Transformer的使用"""
    # 初始化模型和分词器
    tokenizer = SMILESTokenizer()
    model = SMILESTransformer(vocab_size=tokenizer.vocab_size)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    
    # 示例SMILES
    smiles_examples = [
        "CCO",
        "CC(=O)O",
        "c1ccccc1",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    ]
    
    # 验证SMILES
    valid_smiles = []
    for smi in smiles_examples:
        if Chem.MolFromSmiles(smi) is not None:
            valid_smiles.append(smi)
    
    # 获取embeddings（注意：这里使用的是未训练的模型，仅用于演示）
    embeddings = get_smiles_embedding(valid_smiles, model, tokenizer, device)
    print(f"SMILES embeddings shape: {embeddings.shape}")
    
    return embeddings

# 运行演示
# embeddings = demo_smiles_transformer()
```

### 3.4 SMILES-BERT

#### 模型简介
SMILES-BERT是Wang等人开发的基于BERT的分子语言模型[10]，**专门设计用于处理SMILES序列，采用掩码SMILES恢复任务进行大规模无监督预训练**。该模型使用基于注意力机制的Transformer层，能够有效捕获分子序列中的长程依赖关系。

#### 模型特点
- **半监督学习**: 结合大规模无标签数据预训练和下游任务微调
- **注意力机制**: 基于Transformer的注意力机制捕获分子内原子关系
- **可迁移性**: 预训练模型可轻松迁移到不同的分子性质预测任务

#### 使用示例

```python
# SMILES-BERT通常需要从源码安装或使用类似的实现
from transformers import AutoTokenizer, AutoModel
import torch
from rdkit import Chem

def create_smiles_bert_embedding(smiles_list, model_name="DeepChem/ChemBERTa-77M-MLM"):
    """
    使用BERT-like模型生成SMILES embedding
    注：这里使用ChemBERTa作为SMILES-BERT的替代实现
    """
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    
    # 验证SMILES
    valid_smiles = [smi for smi in smiles_list if Chem.MolFromSmiles(smi) is not None]
    
    # Tokenization和编码
    inputs = tokenizer(valid_smiles, return_tensors="pt", padding=True, truncation=True, max_length=512)
    inputs = {key: value.to(device) for key, value in inputs.items()}
    
    # 生成embeddings
    with torch.no_grad():
        outputs = model(**inputs)
        # 使用[CLS] token表示或平均池化
        embeddings = outputs.last_hidden_state.mean(dim=1)  # 平均池化
    
    return embeddings

# 使用示例
smiles_examples = ["CCO", "CC(=O)O", "c1ccccc1", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]
embeddings = create_smiles_bert_embedding(smiles_examples)
print(f"SMILES-BERT embeddings shape: {embeddings.shape}")
```

### 3.5 Smile-to-Bert

#### 模型简介
Smile-to-Bert是最新发布的BERT架构模型[11]，**专门预训练用于从SMILES表示预测113个分子描述符，将分子结构和理化性质信息整合到embeddings中**。该模型在22个分子性质预测数据集上进行了评估，表现优异。

#### 模型特点
- **多任务预训练**: 同时预测113个RDKit计算的分子描述符
- **理化性质感知**: embeddings包含分子结构和理化性质信息
- **最新技术**: 2024年发布，代表最新的分子BERT技术

#### 使用示例

```python
# Smile-to-Bert的概念实现
from transformers import BertModel, BertTokenizer
import torch
from rdkit import Chem

class SmileToBert:
    """Smile-to-Bert模型的概念实现"""
    
    def __init__(self, model_path="smile-to-bert"):
        """
        初始化Smile-to-Bert模型
        注：实际使用需要从官方仓库获取预训练权重
        """
        # 这里使用通用BERT作为示例，实际应使用预训练的Smile-to-Bert权重
        self.tokenizer = BertTokenizer.from_pretrained('bert-base-uncased')
        self.model = BertModel.from_pretrained('bert-base-uncased')
        
        # 添加分子特定的特殊token
        special_tokens = ['[MOL]', '[BOND]', '[RING]']
        self.tokenizer.add_tokens(special_tokens)
        self.model.resize_token_embeddings(len(self.tokenizer))
        
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.model.to(self.device)
    
    def preprocess_smiles(self, smiles):
        """预处理SMILES字符串"""
        # 在SMILES中添加空格以便tokenization
        processed = ' '.join(list(smiles))
        return processed
    
    def get_molecular_embedding(self, smiles_list):
        """获取分子的embedding"""
        # 预处理SMILES
        processed_smiles = [self.preprocess_smiles(smi) for smi in smiles_list]
        
        # Tokenization
        inputs = self.tokenizer(
            processed_smiles,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=512
        )
        inputs = {key: value.to(self.device) for key, value in inputs.items()}
        
        # 获取embeddings
        with torch.no_grad():
            outputs = self.model(**inputs)
            # 使用[CLS] token或平均池化
            embeddings = outputs.last_hidden_state[:, 0, :]  # [CLS] token
        
        return embeddings

# 使用示例
def demo_smile_to_bert():
    """演示Smile-to-Bert使用"""
    # 初始化模型
    smile_bert = SmileToBert()
    
    # 示例SMILES
    smiles_examples = [
        "CCO",  # 乙醇
        "CC(=O)O",  # 乙酸
        "c1ccccc1",  # 苯
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # 咖啡因
    ]
    
    # 验证SMILES有效性
    valid_smiles = []
    for smi in smiles_examples:
        if Chem.MolFromSmiles(smi) is not None:
            valid_smiles.append(smi)
    
    # 生成embeddings
    embeddings = smile_bert.get_molecular_embedding(valid_smiles)
    print(f"Smile-to-Bert embeddings shape: {embeddings.shape}")
    print("Note: 这是概念实现，实际使用需要官方预训练权重")
    
    return embeddings

# 运行演示
# embeddings = demo_smile_to_bert()
```

### 3.6 MolBERT

#### 模型简介
MolBERT是专门为化学领域定制的BERT模型[12]，**针对处理SMILES字符串进行了优化，能够提取丰富的上下文分子表示**。该模型在大规模化学语料库上预训练，特别适合分子相似性搜索和药物发现任务。

#### 模型特点
- **化学特异性**: 专门为化学SMILES数据定制
- **双向上下文**: 利用BERT的双向注意力机制
- **迁移学习**: 在小数据集上表现优异

#### 使用示例

```python
import os
import torch
import yaml
from typing import Sequence, Tuple, Union
import numpy as np

# 这里需要根据实际情况修改类的定义，为了代码完整，从原始文件中提取相关部分
class MolBertFeaturizer:
    def __init__(
        self,
        checkpoint_path: str,
        device: str = None,
        embedding_type: str = 'pooled',
        max_seq_len: int = None,
        permute: bool = False,
    ) -> None:
        super().__init__()
        self.checkpoint_path = checkpoint_path
        self.model_dir = os.path.dirname(os.path.dirname(checkpoint_path))
        self.hparams_path = os.path.join(self.model_dir, 'hparams.yaml')
        self.device = device or 'cuda' if torch.cuda.is_available() else 'cpu'
        self.embedding_type = embedding_type
        self.output_all = False if self.embedding_type in ['pooled'] else True
        self.max_seq_len = max_seq_len
        self.permute = permute

        # load config
        with open(self.hparams_path) as yaml_file:
            config_dict = yaml.load(yaml_file, Loader=yaml.FullLoader)

        # 假设这里有一个简单的 logger 实现，实际使用时需要导入 logging 模块
        class SimpleLogger:
            def debug(self, msg):
                print(msg)
        logger = SimpleLogger()
        logger.debug('loaded model trained with hparams:')
        logger.debug(config_dict)

        # 这里假设 SmilesIndexFeaturizer 已经定义，为了简化，省略其实现
        class SmilesIndexFeaturizer:
            @staticmethod
            def bert_smiles_index_featurizer(max_seq_len, permute):
                return None

        # load smiles index featurizer
        self.featurizer = self.load_featurizer(config_dict)

        # 这里假设 SmilesMolbertModel 已经定义，为了简化，省略其实现
        class SmilesMolbertModel:
            def __init__(self, config):
                self.config = config
            def load_from_checkpoint(self, checkpoint_path, hparam_overrides):
                pass
            def load_state_dict(self, state_dict):
                pass
            def eval(self):
                pass
            def freeze(self):
                pass
            def to(self, device):
                return self

        # load model
        from types import SimpleNamespace
        self.config = SimpleNamespace(**config_dict)
        self.model = SmilesMolbertModel(self.config)
        self.model.load_from_checkpoint(self.checkpoint_path, hparam_overrides=self.model.__dict__)

        # HACK: manually load model weights since they don't seem to load from checkpoint (PL v.0.8.5)
        checkpoint = torch.load(self.checkpoint_path, map_location=lambda storage, loc: storage)
        self.model.load_state_dict(checkpoint['state_dict'])

        self.model.eval()
        self.model.freeze()

        self.model = self.model.to(self.device)

        if self.output_all:
            self.model.model.config.output_hidden_states = True

    def load_featurizer(self, config_dict):
        # load smiles index featurizer
        if self.max_seq_len is None:
            max_seq_len = config_dict.get('max_seq_length')
            # 假设这里有一个简单的 logger 实现，实际使用时需要导入 logging 模块
            class SimpleLogger:
                def debug(self, msg):
                    print(msg)
            logger = SimpleLogger()
            logger.debug('getting smiles index featurizer of length: ', max_seq_len)
        else:
            max_seq_len = self.max_seq_len
        return SmilesIndexFeaturizer.bert_smiles_index_featurizer(max_seq_len, permute=self.permute)

    @staticmethod
    def trim_batch(input_ids, valid):
        # trim input horizontally if there is at least 1 valid data point
        if any(valid):
            _, cols = np.where(input_ids[valid] != 0)
        # else trim input down to 1 column (avoids empty batch error)
        else:
            cols = np.array([0])

        max_idx: int = int(cols.max().item() + 1)

        input_ids = input_ids[:, :max_idx]

        return input_ids

    def transform(self, molecules: Sequence[Any]) -> Tuple[Union[Dict, np.ndarray], np.ndarray]:
        # 这里假设 self.featurizer.transform 已经实现
        input_ids, valid = self.featurizer.transform(molecules)

        input_ids = self.trim_batch(input_ids, valid)

        token_type_ids = np.zeros_like(input_ids, dtype=np.long)
        attention_mask = np.zeros_like(input_ids, dtype=np.long)

        attention_mask[input_ids != 0] = 1

        input_ids = torch.tensor(input_ids, dtype=torch.long, device=self.device)
        token_type_ids = torch.tensor(token_type_ids, dtype=torch.long, device=self.device)
        attention_mask = torch.tensor(attention_mask, dtype=torch.long, device=self.device)

        with torch.no_grad():
            # 这里假设 self.model.model.bert 已经实现
            outputs = self.model.model.bert(
                input_ids=input_ids, token_type_ids=token_type_ids, attention_mask=attention_mask
            )

        if self.output_all:
            sequence_output, pooled_output, hidden = outputs
        else:
            sequence_output, pooled_output = outputs

        # set invalid outputs to 0s
        valid_tensor = torch.tensor(
            valid, dtype=sequence_output.dtype, device=sequence_output.device, requires_grad=False
        )

        pooled_output = pooled_output * valid_tensor[:, None]

        # concatenate and sum last 4 layers
        if self.embedding_type == 'average-sum-4':
            sequence_out = torch.sum(torch.stack(hidden[-4:]), dim=0)  # B x L x H
        # concatenate and sum last 2 layers
        elif self.embedding_type == 'average-sum-2':
            sequence_out = torch.sum(torch.stack(hidden[-2:]), dim=0)  # B x L x H
        # concatenate last four hidden layer
        elif self.embedding_type == 'average-cat-4':
            sequence_out = torch.cat(hidden[-4:], dim=-1)  # B x L x 4*H
        # concatenate last two hidden layer
        elif self.embedding_type == 'average-cat-2':
            sequence_out = torch.cat(hidden[-2:], dim=-1)  # B x L x 2*H
        # only last layer - same as default sequence output
        elif self.embedding_type == 'average-1':
            sequence_out = hidden[-1]  # B x L x H
        # only penultimate layer
        elif self.embedding_type == 'average-2':
            sequence_out = hidden[-2]  # B x L x H
        # only 3rd to last layer
        elif self.embedding_type == 'average-3':
            sequence_out = hidden[-3]  # B x L x H
        # only 4th to last layer
        elif self.embedding_type == 'average-4':
            sequence_out = hidden[-4]  # B x L x H
        # defaults to last hidden layer
        else:
            sequence_out = sequence_output  # B x L x H

        sequence_out = sequence_out * valid_tensor[:, None, None]

        sequence_out = sequence_out.detach().cpu().numpy()
        pooled_output = pooled_output.detach().cpu().numpy()

        if self.embedding_type == 'pooled':
            out = pooled_output
        elif self.embedding_type == 'average-1-cat-pooled':
            sequence_out = np.mean(sequence_out, axis=1)
            out = np.concatenate([sequence_out, pooled_output], axis=-1)
        elif self.embedding_type.startswith('average'):
            out = np.mean(sequence_out, axis=1)
        else:
            out = dict(sequence_output=sequence_out, pooled_output=pooled_output)

        return out, valid

# 示例使用
if __name__ == "__main__":
    # 从 README 中获取预训练模型的下载链接
    checkpoint_path = 'path/to/your/downloaded/checkpoint.ckpt'
    featurizer = MolBertFeaturizer(checkpoint_path=checkpoint_path)

    # 示例分子的 SMILES 字符串
    smiles_list = ['CCO', 'CCN']
    features, valid = featurizer.transform(smiles_list)
    print("Features:", features)
    print("Valid:", valid)
```

### 3.7 通用大语言模型在分子数据上的应用

#### LLaMA和GPT在SMILES上的应用

最近的研究表明，**通用大语言模型如LLaMA和GPT在处理SMILES字符串方面表现出了惊人的能力**[13]。这些模型虽然没有专门为化学领域设计，但其强大的语言理解能力使其能够有效处理分子表示。

#### 性能对比
- **LLaMA**: 在分子性质预测和药物-药物相互作用预测中表现优于GPT
- **GPT**: 虽然性能略逊于LLaMA，但仍能产生有意义的分子表示
- **与专用模型对比**: LLaMA在某些任务上可与专门的分子预训练模型相媲美

#### 使用示例

```python
# 使用HuggingFace接口调用通用大语言模型
from transformers import LlamaTokenizer, LlamaModel, GPT2Tokenizer, GPT2Model
import torch
from rdkit import Chem

class UniversalLLMForMolecules:
    """通用大语言模型用于分子表示学习"""
    
    def __init__(self, model_type='llama', model_name=None):
        """
        初始化通用LLM
        
        参数:
            model_type: 'llama' 或 'gpt2'
            model_name: 具体模型名称
        """
        if model_type == 'llama':
            # 注意：需要申请LLaMA访问权限
            model_name = model_name or "meta-llama/Llama-2-7b-hf"
            self.tokenizer = LlamaTokenizer.from_pretrained(model_name)
            self.model = LlamaModel.from_pretrained(model_name)
        elif model_type == 'gpt2':
            model_name = model_name or "gpt2"
            self.tokenizer = GPT2Tokenizer.from_pretrained(model_name)
            self.model = GPT2Model.from_pretrained(model_name)
            # GPT2需要设置pad_token
            self.tokenizer.pad_token = self.tokenizer.eos_token
        else:
            raise ValueError(f"Unsupported model type: {model_type}")
        
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.model.to(self.device)
        self.model.eval()
    
    def get_molecular_embeddings(self, smiles_list):
        """使用通用LLM获取分子embeddings"""
        # 验证SMILES
        valid_smiles = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                valid_smiles.append(smi)
        
        # 为SMILES添加描述性前缀以提高理解
        prompted_smiles = [f"Molecule with SMILES: {smi}" for smi in valid_smiles]
        
        # Tokenization
        inputs = self.tokenizer(
            prompted_smiles,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=512
        )
        inputs = {key: value.to(self.device) for key, value in inputs.items()}
        
        # 生成embeddings
        with torch.no_grad():
            outputs = self.model(**inputs)
            hidden_states = outputs.last_hidden_state
            
            # 使用平均池化获取序列级表示
            attention_mask = inputs['attention_mask'].unsqueeze(-1)
            masked_embeddings = hidden_states * attention_mask
            embeddings = masked_embeddings.sum(dim=1) / attention_mask.sum(dim=1)
        
        return embeddings

# 使用示例（需要相应的模型访问权限）
def demo_universal_llm():
    """演示通用LLM在分子数据上的应用"""
    try:
        # 使用GPT-2（更容易获取）
        llm = UniversalLLMForMolecules(model_type='gpt2', model_name='gpt2')
        
        smiles_examples = ["CCO", "CC(=O)O", "c1ccccc1"]
        embeddings = llm.get_molecular_embeddings(smiles_examples)
        
        print(f"Universal LLM embeddings shape: {embeddings.shape}")
        print("注意：通用LLM可能需要更多的提示工程以获得最佳性能")
        
    except Exception as e:
        print(f"Error loading universal LLM: {e}")
        print("请确保已安装相应的模型和权限")

# demo_universal_llm()
```

## 四、模型对比与选择指南

### 4.1 主要模型对比表

| 类别 | 模型 | 参数量 | 输出维度 | 预训练数据规模 | 主要优势 | 适用场景 |
|------|------|--------|----------|----------------|----------|----------|
| **蛋白质** | ESM-2 | 8M-15B | 320-5120 | 250M序列 | 进化信息丰富，多规模选择 | 蛋白质结构预测、功能注释 |
| | ESM-C | 300M-6B | 1152 | >1B序列 | 更高效率，更强性能 | 大规模蛋白质分析 |
| | CARP | 640M | 1280 | ~1.7M序列 | 对比学习，自回归建模 | 蛋白质生成、设计 |
| | ProtT5 | ~3B | 1024 | 45M序列 | T5架构，编码器-解码器 | 多任务蛋白质预测 |
| | Ankh | ~3B | 1536 | 多语言数据 | 多语言支持 | 跨语言蛋白质研究 |
| **肽** | PepBERT | ~300M | 320 | UniParc肽序列 | 专门优化短肽 | 肽-蛋白质相互作用 |
| **小分子** | ChemBERTa | 12M-77M | 384-768 | 77M分子 | 首个分子BERT，成熟生态 | 分子性质预测 |
| | MolFormer | 47M | 512-768 | 1.1B分子 | 线性注意力，处理长序列 | 大规模分子筛选 |
| | SMILES Transformer | ~10M | 512 | 1.7M分子 | 自编码，低数据优化 | 小数据集药物发现 |
| | SMILES-BERT | ~12M | 768 | 大规模SMILES | 掩码语言建模，半监督 | 分子性质预测 |
| | Smile-to-Bert | ~110M | 768 | PubChem+113描述符 | 多任务预训练，理化性质感知 | 综合分子性质预测 |
| | MolBERT | ~12M | 768 | 化学语料库 | 化学特异性，双向上下文 | 分子相似性搜索 |
| | LLaMA (分子) | 7B+ | 4096+ | 通用+SMILES | 强大语言理解，泛化能力 | 复杂分子推理任务 |
| | GPT (分子) | 175B+ | 12288+ | 通用+SMILES | 生成能力强，对话式交互 | 分子生成和解释 |

### 4.2 性能与效率对比

#### 计算资源需求

| 模型类别 | 内存需求 | 推理速度 | 训练复杂度 | GPU要求 |
|----------|----------|----------|------------|---------|
| ESM-2 (650M) | ~3GB | 中等 | 高 | V100/A100推荐 |
| ESM-C (600M) | ~2.5GB | 快 | 中等 | GTX 1080Ti可用 |
| ChemBERTa | ~500MB | 快 | 低 | GTX 1060可用 |
| MolFormer | ~1GB | 快 | 中等 | RTX 2080可用 |
| SMILES-BERT | ~500MB | 快 | 中等 | GTX 1060可用 |
| Smile-to-Bert | ~1GB | 中等 | 中等 | RTX 2080可用 |
| MolBERT | ~500MB | 快 | 低 | GTX 1060可用 |
| LLaMA (7B) | ~14GB | 慢 | 极高 | A100推荐 |
| GPT (175B) | >350GB | 极慢 | 极高 | 多卡A100 |

#### 准确性表现

1. **蛋白质任务**
   - **结构预测**: ESM-2 > ESM-C > ProtT5
   - **功能预测**: ESM-C ≥ ESM-2 > CARP
   - **肽相互作用**: PepBERT > 通用蛋白质模型

2. **分子性质预测**
   - **通用性能**: MolFormer > Smile-to-Bert > ChemBERTa-2 > ChemBERTa
   - **小数据集**: SMILES Transformer > SMILES-BERT > 大模型
   - **多任务学习**: Smile-to-Bert > MolBERT > ChemBERTa
   - **理化性质**: Smile-to-Bert > 传统描述符方法
   - **通用推理**: LLaMA > GPT > 专用模型（在某些复杂任务上）

### 4.3 选择建议

#### 根据应用场景选择

**蛋白质研究**
- **结构生物学**: ESM-2 (t33或更大)
- **大规模分析**: ESM-C (600M)
- **蛋白质设计**: CARP
- **多任务预测**: ProtT5

**小分子研究**
- **药物发现**: MolFormer或Smile-to-Bert
- **新药研发**: ChemBERTa-2或MolBERT
- **分子生成**: 结合GPT/LLaMA的方法
- **概念验证**: ChemBERTa或SMILES Transformer
- **理化性质预测**: Smile-to-Bert（专门优化）

**肽研究**
- **肽-蛋白质相互作用**: PepBERT
- **抗菌肽设计**: PepBERT + 微调

#### 根据资源条件选择

**高性能计算环境**
- 推荐: ESM-2大模型、MolFormer-XL、LLaMA/GPT分子应用
- 优势: 最佳性能，支持复杂推理

**标准工作站**
- 推荐: ESM-C、ChemBERTa、MolFormer标准版、Smile-to-Bert
- 平衡性能与资源需求

**资源受限环境**
- 推荐: ESM-2小模型、SMILES Transformer、SMILES-BERT
- 确保基本功能

#### 根据数据特点选择

**大规模数据**
- 使用预训练大模型: MolFormer、ESM-C、LLaMA/GPT
- 利用规模优势

**小规模数据**
- 使用专门优化的模型: SMILES Transformer、PepBERT、SMILES-BERT
- 或使用预训练+微调

**特定领域**
- 理化性质预测: Smile-to-Bert
- 短肽: PepBERT
- 分子生成: GPT/LLaMA方法
- 化学推理: 通用大语言模型

## 五、最佳实践与技巧

### 5.1 模型选择策略

1. **原型阶段**: 使用小模型快速验证想法
2. **性能优化**: 逐步升级到大模型
3. **生产部署**: 平衡性能与资源需求
4. **特殊需求**: 选择专门优化的模型

### 5.2 优化技巧

#### 内存优化
```python
# 使用混合精度
model = model.half()

# 梯度检查点
model.gradient_checkpointing_enable()

# 批处理优化
def batch_inference(data, model, batch_size=32):
    results = []
    for i in range(0, len(data), batch_size):
        batch = data[i:i+batch_size]
        with torch.no_grad():
            result = model(batch)
        results.append(result.cpu())
        torch.cuda.empty_cache()
    return torch.cat(results)
```

#### 速度优化
```python
# 模型编译（PyTorch 2.0+）
model = torch.compile(model)

# TensorRT优化（NVIDIA GPU）
import torch_tensorrt
optimized_model = torch_tensorrt.compile(model)
```

### 5.3 实用工具函数

```python
def standardize_molecular_input(smiles_list):
    """标准化分子输入"""
    from rdkit import Chem
    standardized = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            # 标准化SMILES
            canonical_smi = Chem.MolToSmiles(mol, canonical=True)
            standardized.append(canonical_smi)
        else:
            print(f"Invalid SMILES: {smi}")
    return standardized

def validate_protein_sequence(sequence):
    """验证蛋白质序列"""
    valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    return all(aa in valid_amino_acids for aa in sequence.upper())

def estimate_memory_usage(model_name, batch_size, sequence_length):
    """估算内存使用量"""
    memory_map = {
        'esm2_t33_650M': lambda b, l: b * l * 1280 * 4 * 1e-9 + 2.5,
        'chemberta': lambda b, l: b * l * 768 * 4 * 1e-9 + 0.5,
        'molformer': lambda b, l: b * l * 768 * 4 * 1e-9 + 1.0,
    }
    
    if model_name in memory_map:
        estimated_gb = memory_map[model_name](batch_size, sequence_length)
        return f"Estimated memory usage: {estimated_gb:.2f} GB"
    else:
        return "Memory estimation not available for this model"
```

## 参考文献

[1] Lin Z, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science. 2023;379(6637):1123-1130.

[2] EvolutionaryScale. ESM Cambrian: Focused on creating representations of proteins. 2024. Available: https://github.com/evolutionaryscale/esm

[3] Rao R, et al. MSA Transformer. In: International Conference on Machine Learning. 2021:8844-8856.

[4] Elnaggar A, et al. ProtTrans: towards cracking the language of Life's code through self-supervised deep learning and high performance computing. IEEE Transactions on Pattern Analysis and Machine Intelligence. 2021;44(10):7112-7127.

[5] ElNaggar A, et al. Ankh: Optimized protein language model unlocks general-purpose modelling. 2023. Available: https://huggingface.co/ElnaggarLab/ankh-large

[6] Zhang H, et al. PepBERT: A BERT-based model for peptide representation learning. 2023. Available: https://github.com/dzjxzyd/PepBERT-large

[7] Chithrananda S, Grand G, Ramsundar B. ChemBERTa: Large-scale self-supervised pretraining for molecular property prediction. arXiv preprint arXiv:2010.09885. 2020.

[8] Ross J, et al. Large-scale chemical language representations capture molecular structure and properties. Nature Machine Intelligence. 2022;4(12):1256-1264.

[9] Honda S, Shi S, Ueda HR. SMILES transformer: Pre-trained molecular fingerprint for low data drug discovery. 2019. Available: https://github.com/DSPsleeporg/smiles-transformer

[10] Wang S, Guo Y, Wang Y, Sun H, Huang J. SMILES-BERT: Large scale unsupervised pre-training for molecular property prediction. Proceedings of the 10th ACM International Conference on Bioinformatics, Computational Biology and Health Informatics. 2019:429-436.

[11] Barranco-Altirriba M, Würf V, Manzini E, Pauling JK, Perera-Lluna A. Smile-to-Bert: A BERT architecture trained for physicochemical properties prediction and SMILES embeddings generation. bioRxiv. 2024. doi:10.1101/2024.10.31.621293.

[12] MolBERT: A BERT-based model for molecular representation learning. GitHub. Available: https://github.com/BenevolentAI/MolBERT

[13] Al-Ghamdi A, et al. Can large language models understand molecules? BMC Bioinformatics. 2024;25:347.

[14] Molecular Transformer. Schwaller P, et al. Molecular transformer: a model for uncertainty-calibrated chemical reaction prediction. ACS Central Science. 2019;5(9):1572-1583.

[15] ST-KD. Li S, et al. Stepping back to SMILES transformers for fast molecular representation inference. 2021. Available: https://openreview.net/forum?id=CyKQiiCPBEv
