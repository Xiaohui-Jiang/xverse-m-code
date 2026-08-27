# xVERSE 上一版本：文章故事记录

本文件只记录 `xVERSE_manu/manu/V1/main/main.tex` 及其 `sections/introduction.tex`、`results.tex`、`discussion.tex`、`methods.tex` 中的论文叙事，不对代码实现做判断。

## 1. 核心主张

文章把 xVERSE 定义为一个 transcriptomics-native foundation model：它不仅学习 universal cell and gene representations，还直接学习整个 transcriptome 的表达分布，因此能够生成 virtual cells，并把生成能力用于 imputation、small-data analysis 和 cross-modality prediction。

## 2. 问题起点

文章首先建立三重限制：

1. 现有 single-cell foundation models 多借用 BERT/GPT 等语言模型，把转录组转换成 token 序列，可能忽视基因表达的无序性、高维性、稀疏性和计数分布。
2. 许多模型主要解决 representation learning，并不真正建模完整 transcriptomic probability distribution，因此难以生成高保真的 virtual cell。
3. targeted spatial panels、稀有细胞和小规模临床样本都受到实验成本与观测能力限制，真实生物信号可能无法被充分测量或统计确认。

文章据此提出：foundation model 不应只帮助研究者表示已有数据，也应帮助研究者扩展可观测和可分析的生物数据空间。

## 3. 文章提出的解决方向

文章将 xVERSE 的设计概括为四个创新点：

- 直接面向 transcriptomic data 建模，而不是强加人工 sequential structure。
- 用 panel-aware stochastic gene masking 适应不同 gene panels。
- 用 GRL 将 biological variation 与 technical confounders 分开。
- 用 distribution reconstruction loss 建模 cell- and gene-specific expression distributions。

论文将最终能力归纳为三个领域：universal representation、cell profile synthesis 和 biological discovery。

## 4. Results 的叙事顺序

### 第一部分：Universal representation

文章先在独立的 human pediatric liver 和 ALS motor cortex 数据上进行 zero-shot benchmark，比较 xVERSE、scGPT、Nicheformer、Geneformer 和 Harmony。重点展示：

- 保留 cell-type heterogeneity。
- 减少 batch effect。
- 在 whole transcriptome、Xenium Prime panel 和 tissue-specific panel 下保持表现。
- 同时具备较高 inference efficiency。

这一部分要证明 xVERSE 学到的是可迁移的 biological representation，而非某个训练数据集的特征。

### 第二部分：Gene2Cell interpretability

文章随后引入 Gene2Cell score，将模型的细胞表示解释为基因对 cellular identity 的贡献。高分基因被用于：

- 检验哪些基因最能保留 cell-type heterogeneity。
- 指导 spatial transcriptomics 的 gene-panel design。
- 在减少测量基因数量的同时尽量保留表示与插补能力。

这一部分把 foundation model 从黑箱 representation 工具推进到可用于 gene prioritization 的工具。

### 第三部分：Virtual cell synthesis

文章接着证明 xVERSE 能生成高保真的 virtual cells。证据包括：

- 重现真实细胞的 UMI count distributions。
- 保留 HVG expression rankings。
- 在联合 UMAP 中与 biological cells 融合，同时维持 cell-type clusters。
- biological-versus-virtual classifier 的 AUROC 接近 0.5。

这里的叙事重点是“indistinguishable from biological data”，而不只是“生成结果看起来相似”。

### 第四部分：Spatial imputation

文章把 targeted spatial panel 描述为一个硬件和实验效率造成的观测瓶颈。xVERSE 的 whole-transcriptome prior 被用于从 partial panel 推断 unmeasured genes，并与 SpaGE、gimVI 比较。

叙事强调两点：xVERSE 可以 zero-shot 工作，并且不依赖 external single-cell reference；同时在不同 reference 条件变化时比专门方法更稳定。

### 第五部分：Small-data biological discovery

文章把 virtual cell synthesis 从生成任务推进到 biological analysis：

- 在 minor population 只有 4--10 个细胞时，增强数据帮助 Leiden 找到 rare cell types。
- 在极小样本 DEG 分析中，xVERSE augmentation 比简单复制细胞更接近 full-data ground truth。

这里的中心概念是：virtual cells 能放大已经存在但统计功效不足的 biological signal。

### 第六部分：Cross-modality generalization

最后，文章使用 heart-transplant CITE-seq/VDJ-seq 数据，训练集只包含 normal controls，测试集包含 NGD 和 CAV pathological states。xVERSE virtual cells 被用于增强 ADT prediction、B-cell heavy-chain isotype classification 和 T-cell lineage prediction。

这一部分把模型价值提升到 out-of-distribution generalization：合成数据帮助下游模型面对训练期间未见过的 pathological state。

## 5. 文章的完整证据链

文章的逻辑推进是：

$$
\text{universal representation}
\rightarrow \text{expression distribution}
\rightarrow \text{virtual cells}
\rightarrow \text{imputation and augmentation}
\rightarrow \text{biological discovery and OOD generalization}.
$$

大规模跨组织数据提供 universal prior 的基础；zero-shot benchmark 证明 representation 泛化；Gene2Cell 提供可解释性；virtual-cell fidelity 实验验证生成质量；imputation、rare-cell、DEG 和 cross-modality 实验说明生成结果能解决实际生物学问题。

## 6. 文章中的核心对比

文章通过几组对比建立必要性：

- language-model adaptation 对比 transcriptomics-native modeling；
- foundation model 的普适性对比 specialized method 的任务精度；
- biological cells 对比 computationally synthesized virtual cells；
- external-reference-dependent imputation 对比 internalized universal prior；
- 实验数据稀缺对比计算生成带来的数据扩展。

这些对比共同服务于文章的中心叙事：xVERSE 不只是提高已有分析的性能，而是在计算层面扩展实验和生物发现的边界。

## 7. Discussion 中的最终定位

文章最终将 xVERSE 定位为：

- 可以跨任务服务 single-cell machine learning 的通用生成引擎；
- 可以通过 virtual cells 缓解 rare population 和 low-n study 的统计限制；
- 可以帮助 targeted spatial technology 获得更广的 transcriptomic insight；
- 可以为未来整合 chromatin accessibility 和 protein abundance 的多组学 latent space 提供基础。

文章同时承认其适用边界：模型主要面向 UMI-based single-cell 和 imaging-based spatial count distributions，对 Smart-seq 等非 UMI protocol 的适应性有限。
