# xVERSE 上一版本：代码技术记录

本文件只记录 `/Users/xiaohui/LocalFiles/Codes/xverse-code` 中实际存在的代码设计，不补充论文中的未实现描述。

## 1. 模型结构

核心模型是 `main/utils_model.py` 中的 `XVerseModel`：

- 输入是固定全局 gene universe 中的 dense count matrix，以及每个 cell 的 `observed_mask`。
- `FiLMMaskEncoder` 分别编码 `log1p(counts)` 和 observed-gene mask。
- mask encoder 产生 FiLM 的 gamma/beta，调制 expression representation。
- encoder 输出 `z_bio_raw`，再根据配置决定是否进行 L2 normalization，得到公开的 `z_bio`。
- `DenseExpressionDecoder` 将 `z_bio` 解码为 gene logits；`library_head` 预测 library size；softmax 后得到非负的 `mu_bio`。
- 如果模型启用 sample embedding，decoder 可进一步生成 sample-conditioned 的 `mu`。

## 2. 观测 mask 与随机 mask

`FiLMMaskEncoder._apply_random_mask()` 在训练时对已观测基因随机隐藏：

- 观测基因较少时，随机隐藏一部分基因，至少保留一个观测基因。
- 观测基因较多时，按不同概率随机隐藏 0--70% 的基因。
- 对足够大的 panel，部分情况下与预配置 panel mask 相交，用于模拟空间 panel。

模型输入中的 `observed_mask` 与 encoder 使用的随机 mask 是两个概念：前者决定哪些基因可以进入 likelihood，后者决定本次 view 的 encoder 输入。

## 3. 计数似然

代码实现了三类 reconstruction loss：

- Poisson negative log-likelihood。
- Negative Binomial negative log-likelihood。
- Zero-Inflated Negative Binomial negative log-likelihood。

NB dispersion 支持 `global`、`gene`、`cell`、`factorized`、`cell_gene` 和 `lowrank` 六种 parameterization。数值计算在 fp32 中进行，并对 `mu`、dispersion 和 count 做有限值处理与边界裁剪。

## 4. 辅助目标与训练流程

`pretrain_one_epoch()` 对每个 batch 做两次独立 forward，形成两个随机 mask view。训练目标是：

$$
\mathcal{L} =
\lambda_{sample}\mathcal{L}_{recon}(\mu)
+ \lambda_{bio}\mathcal{L}_{recon}(\mu_{bio})
+ \lambda_{celltype}\mathcal{L}_{celltype}
+ \lambda_{contrast}\mathcal{L}_{contrast}.
$$

- `mu` 和 `mu_bio` 都只在 observed genes 上计算重构损失。
- cell type 使用带 label smoothing 的 cross-entropy；无标签 cell 被忽略。
- 两个 view 的 projection 使用双向 InfoNCE 对比损失。
- 优化器是 Adam，使用 AMP、gradient clipping 和 `ReduceLROnPlateau`。
- checkpoint 主要依据 validation `loss_nb_bio` 选择。

代码支持将 cell-type label 替换为预计算文本 embedding 的 cosine-style matching，这属于辅助预测头。

## 5. 数据加载

`main/data.py` 有两套数据路径：

### 旧式 NPZ block

`FastXVerseBatchDataset` 延迟加载稀疏 CSR block，通过局部 gene list 映射到全局 gene order，并用 LRU cache 控制内存。`SparseBatchCollator` 将稀疏 row 组织成 batch payload。

### compiled `xverse_train_v1`

`CompiledShardDataset` 使用 memory-mapped shard 和 manifest。它保存 global cell range、sample ID、cell type ID 及稀疏 gene/value arrays。若 manifest 中的 `source_pair_id` 和 source metadata 可访问，loader 会恢复完整 measured panel，从而把“测量为零”和“未测量”区分开；否则只能用非零 gene index 作为 fallback，并打印警告。

## 6. 采样策略

`BalancedSampleSampler` 和 `CompiledBalancedSampler` 按 sample ID 平衡抽样，并支持：

- 每个 sample 的固定抽样数量。
- 可复现的 sample 内排序。
- 固定训练步数下的数据 fraction ablation。
- shard locality、active shard 和 reorder window。

这个设计使训练数据量实验改变 cell diversity，同时尽量保持 optimizer steps 和 validation split 不变。

## 7. 推理与 fine-tuning

`main/cli_xverse.py` 提供 embedding 和 generation 两类任务：

- embedding 将 `z_bio` 写入 `adata.obsm["xVerse"]`。
- generation 将 `mu_bio` 及 NB samples 写入新的 `.h5ad` 文件。

`XVerseFineTuneModel` 复用 base model 的 encoder、decoder、library head 和 dispersion，并可以增加 `sample_emb_ft`。它输出 `mu_bio`，在 sample ID 有效时额外输出 sample-conditioned `mu`。fine-tuning loop 使用 reconstruction、contrastive loss、validation split、scheduler 和 early stopping。

代码中虽然定义了 `freeze_base_model()`，但 CLI 的 fine-tuning loop 没有自动调用它；因此默认会更新复用的 base encoder/decoder 参数以及 sample embedding（若启用）。

## 8. 代码层面的技术精神

1. 用显式 observed mask 处理不同 gene panel，而不是把缺失 panel 当作真实零值。
2. 用随机 mask 的两个 view 学习 panel-robust representation。
3. 用 `z_bio` 生成 biological mean，用 sample embedding 表达 sample-specific variation。
4. 用 count likelihood 而非普通 MSE 约束表达生成。
5. 通过稀疏加载、memory mapping 和 sample-balanced sampling 支撑大规模训练。
6. 让同一模型同时服务 embedding、imputation、virtual cell generation 和 downstream augmentation。
