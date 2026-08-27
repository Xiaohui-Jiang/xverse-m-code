# Previous xVERSE Version: Story Summary

This document records only the narrative presented in `xVERSE_manu/manu/V1/main/main.tex` and its `introduction.tex`, `results.tex`, `discussion.tex`, and `methods.tex` sections.

## 1. Central Claim

The paper presents xVERSE as a transcriptomics-native foundation model that learns universal cell and gene representations while directly learning the expression distribution of the transcriptome. This allows it to generate virtual cells and apply that generative capability to imputation, small-data analysis, and cross-modality prediction.

## 2. Starting Problem

The paper establishes three related limitations:

1. Many single-cell foundation models adapt BERT- or GPT-like language models and represent transcriptomes as token sequences, potentially overlooking the unordered, high-dimensional, sparse, and count-based nature of gene expression.
2. Many models primarily address representation learning rather than modeling the full transcriptomic probability distribution, making high-fidelity virtual-cell synthesis difficult.
3. Targeted spatial panels, rare cell populations, and small clinical datasets are constrained by experimental cost and measurement capacity, so important biological signals may not be sufficiently observed or statistically confirmed.

The paper therefore argues that a foundation model should not only represent existing data, but also expand the biological data space that researchers can observe and analyze.

## 3. Proposed Direction

The paper describes four architectural innovations:

- Direct transcriptomic modeling instead of imposing an artificial sequential structure.
- Panel-aware stochastic gene masking for varying gene panels.
- Gradient Reversal Layer training to separate biological variation from technical confounders.
- Distribution reconstruction loss for cell- and gene-specific expression distributions.

The paper organizes the resulting capabilities into three domains: universal representation, cell profile synthesis, and biological discovery.

## 4. Results Narrative

### Part I: Universal Representation

The paper first evaluates zero-shot representations on independent human pediatric liver and ALS motor cortex datasets, comparing xVERSE with scGPT, Nicheformer, Geneformer, and Harmony. It emphasizes:

- Preservation of cell-type heterogeneity.
- Reduction of batch effects.
- Consistent performance across whole-transcriptome, Xenium Prime, and tissue-specific panels.
- High inference efficiency.

This section is intended to show that xVERSE learns transferable biological representations rather than features tied to one training dataset.

### Part II: Gene2Cell Interpretability

The paper then introduces the Gene2Cell score, interpreting cell representations through the contribution of individual genes to cellular identity. High-scoring genes are used to:

- Test which genes preserve cell-type heterogeneity.
- Guide spatial transcriptomics gene-panel design.
- Retain representation and imputation performance with fewer measured genes.

This section moves the foundation model from a black-box representation tool toward a tool for gene prioritization.

### Part III: Virtual Cell Synthesis

The paper next demonstrates high-fidelity virtual-cell generation. The evidence includes:

- Reproduction of real-cell UMI count distributions.
- Preservation of highly variable gene expression rankings.
- Integration with biological cells in joint UMAPs while maintaining cell-type clusters.
- Biological-versus-virtual classifiers with AUROC near 0.5.

The narrative emphasis is that virtual cells are not merely visually similar; they are difficult to distinguish from biological data in a statistical classification test.

### Part IV: Spatial Imputation

The paper frames targeted spatial panels as an observation bottleneck caused by hardware and experimental efficiency constraints. xVERSE uses a whole-transcriptome prior to infer unmeasured genes from partial panels and is compared with SpaGE and gimVI.

The narrative highlights that xVERSE can work in a zero-shot setting without an external single-cell reference, and that its performance is more stable across reference conditions than the specialized methods.

### Part V: Small-Data Biological Discovery

The paper moves virtual-cell synthesis into biological analysis:

- Data augmentation helps Leiden identify rare cell types when a minor population contains only 4--10 cells.
- In extremely small-sample DEG analysis, xVERSE augmentation performs better than simply copying cells and more closely approaches the full-data ground truth.

The central concept is that virtual cells amplify biological signals that are present but statistically underpowered.

### Part VI: Cross-Modality Generalization

Finally, the paper uses heart-transplant CITE-seq/VDJ-seq data. Models are trained on normal controls and tested on NGD and CAV pathological states. xVERSE-generated virtual cells augment prediction of ADT levels, B-cell heavy-chain isotypes, and T-cell lineages.

This section elevates the model's value to out-of-distribution generalization: synthetic data help downstream models handle pathological states that were not present during training.

## 5. Complete Evidence Chain

The paper's progression is:

$$
\text{universal representation}
\rightarrow \text{expression distribution}
\rightarrow \text{virtual cells}
\rightarrow \text{imputation and augmentation}
\rightarrow \text{biological discovery and OOD generalization}.
$$

Large cross-tissue data establish the basis for a universal prior; zero-shot benchmarks test representation transfer; Gene2Cell provides interpretability; virtual-cell fidelity experiments test generation quality; and imputation, rare-cell, DEG, and cross-modality experiments show that generated data can address practical biological problems.

The story is therefore not simply that strong benchmark scores imply usefulness. It is:

> Learn a universal biological representation, learn the full expression distribution, generate credible cells, extend missing or scarce observations, and improve biological discovery and cross-state generalization.

## 6. Roles in the Paper's Narrative

- `z_bio` represents a biological space shared across experimental conditions.
- `mu_bio` represents the underlying biological expression profile after removing technical factors.
- The sample-conditioned output represents a generative profile that retains the technical characteristics of the target data.
- Gene2Cell scores connect the internal model to an interpretation of biological identity.
- Virtual cells connect the foundation model to downstream machine learning and biological discovery.

This role assignment allows the paper to cover representation learning, generative modeling, spatial imputation, small-data analysis, and OOD generalization while maintaining one central narrative.

## 7. Key Contrasts

The paper builds its case through several contrasts:

- Language-model adaptation versus transcriptomics-native modeling.
- Foundation-model universality versus specialized-method task precision.
- Experimentally observed biological cells versus computationally synthesized virtual cells.
- External-reference-dependent imputation versus an internalized universal prior.
- Experimental data scarcity versus computational expansion of the data space.

Together, these contrasts support the paper's central narrative: xVERSE does not merely improve existing analyses; it expands the boundary of experimentation and biological discovery computationally.

## 8. Final Positioning in the Discussion

The paper positions xVERSE as:

- A general generative engine that can support many single-cell machine-learning tasks.
- A way to mitigate the statistical limitations of rare populations and low-n studies through virtual cells.
- A means of extracting broader transcriptomic insight from targeted spatial technologies.
- A foundation for future latent spaces that integrate chromatin accessibility and protein abundance.

The paper also states a limitation: the model is designed primarily for UMI-based single-cell and imaging-based spatial count distributions, and its adaptability to non-UMI protocols such as Smart-seq remains limited.
