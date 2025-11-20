# TRADE vs Alternative Perturb-seq Analysis Methods

**Comprehensive Comparison for Agent-Based Evaluation**

---

## Methods Landscape Overview

### Primary Analysis Categories

**1. Gene-Level Differential Expression**
- DESeq2, edgeR, limma (pseudo-bulk)
- MAST, Seurat-Wilcox/NB, scVI (single-cell)

**2. Distribution-Level Analysis**
- **TRADE** (this method)

**3. Topic Modeling / Dimensionality Reduction**
- MUSIC, Mixscape, scRNA-CRISPR
- PCA/UMAP-based clustering

**4. Causal Inference**
- CINEMA-OT
- Propensity score matching

**5. Deep Learning**
- scVI, scGEN, CPA
- VAE-based perturbation prediction

---

## Head-to-Head Comparisons

### TRADE vs DESeq2

| Aspect | TRADE | DESeq2 |
|--------|-------|--------|
| **Analysis Type** | Distribution estimation | Gene-level testing |
| **Output** | TI, πDEG, enrichments | p-values, log2FC per gene |
| **Power-dependence** | ✅ Independent | ❌ Dependent |
| **Speed** | Fast (seconds) | Fast (seconds-minutes) |
| **Best for** | Quantitative comparison | Gene discovery |
| **Statistical model** | Empirical Bayes mixture | Negative binomial GLM |
| **Pseudo-bulk/single-cell** | Pseudo-bulk only | Pseudo-bulk only |
| **Strengths** | Cross-study comparison | Gene lists for validation |
| **Weaknesses** | No gene-level output | Cannot compare across power |

**When to use each:**
- **DESeq2:** Identify which genes to validate
- **TRADE:** Quantify overall impact, compare across conditions
- **Both:** Comprehensive analysis

**Relationship:** TRADE operates on DESeq2 output (complementary, not competitive)

**Benchmarking:**
- Replicability: TRADE TI r=0.90, DESeq2 gene counts r<0.5
- Cross-dataset: TRADE TI R²=59.7%, DESeq2 counts R²=28.4%

---

### TRADE vs MAST

| Aspect | TRADE | MAST |
|--------|-------|------|
| **Analysis Type** | Distribution estimation | Gene-level testing (single-cell) |
| **Output** | TI, πDEG, enrichments | p-values, log2FC per gene |
| **Resolution** | Population-level | Cell-level |
| **Power-dependence** | ✅ Independent | ❌ Dependent |
| **Speed** | Fast | Moderate (slower than TRADE) |
| **Best for** | Cross-study comparison | Cell-type-specific DE |
| **Statistical model** | Empirical Bayes mixture | Hurdle model (binomial + Gaussian) |
| **Handles zeros** | Pseudo-bulk (fewer zeros) | Explicitly models zero-inflation |
| **Strengths** | Power-independence | Single-cell resolution |
| **Weaknesses** | Pseudo-bulk only | Inflated p-values (some datasets) |

**When to use each:**
- **MAST:** Detect subpopulation-specific effects
- **TRADE:** Compare perturbations across contexts
- **Both:** MAST for heterogeneity, TRADE for quantification

**Relationship:** Competitive for different goals

**Benchmarking (2024 Genome Biology study):**
- MAST showed inflated p-values on some Perturb-seq datasets
- False rejection of ~2000 null perturbation-gene pairs
- TRADE not tested in this benchmark (different analysis type)

---

### TRADE vs Seurat

| Aspect | TRADE | Seurat |
|--------|-------|--------|
| **Analysis Type** | Distribution estimation | Gene-level testing (single-cell) |
| **Output** | TI, πDEG, enrichments | p-values, log2FC per gene |
| **Methods** | ash/mashr | Wilcoxon, NB, MAST, others |
| **Power-dependence** | ✅ Independent | ❌ Dependent |
| **Speed** | Fast | Fast (Wilcox), slower (NB) |
| **Best for** | Quantitative comparison | Integrated workflow |
| **Ecosystem** | Standalone | Full scRNA-seq pipeline |
| **Strengths** | Power-independence | Workflow integration, visualization |
| **Weaknesses** | No gene-level output | Cannot compare across power |

**When to use each:**
- **Seurat:** Full scRNA-seq analysis pipeline
- **TRADE:** Quantitative perturbation comparison
- **Integration:** Can use Seurat → DESeq2 → TRADE

**Relationship:** Seurat is broader ecosystem, TRADE is specialized tool

**Benchmarking (2024 study):**
- Seurat-Wilcox and Seurat-NB best performers for CRISPR screens
- Still showed miscalibration
- TRADE orthogonal (different goal)

---

### TRADE vs scVI

| Aspect | TRADE | scVI |
|--------|-------|------|
| **Analysis Type** | Distribution estimation | Deep learning (VAE) + DE |
| **Output** | TI, πDEG, enrichments | Latent embeddings, DE p-values |
| **Power-dependence** | ✅ Independent | ⚠️ Partial (learns from data) |
| **Speed** | Fast (CPU) | Slow (GPU recommended) |
| **Best for** | Cross-study comparison | Large datasets, batch correction |
| **Statistical model** | Empirical Bayes | Variational autoencoder |
| **Interpretability** | ✅ High (clear parameters) | ⚠️ Moderate (black-box embeddings) |
| **Strengths** | Theoretical foundation | Handles complex batch effects |
| **Weaknesses** | Pseudo-bulk only | Slow, GPU needed, less interpretable |

**When to use each:**
- **scVI:** Large datasets with complex batch effects
- **TRADE:** Interpretable quantitative comparison
- **Both:** scVI for dimensionality reduction → TRADE for comparison

**Relationship:** Different philosophies (deep learning vs statistical modeling)

**Benchmarking (2024 study):**
- scVI performed well (AUPRC ~0.25-0.28)
- Comparable to DESeq2, better than some single-cell methods
- TRADE not in this benchmark (different analysis type)

---

### TRADE vs MUSIC

| Aspect | TRADE | MUSIC |
|--------|-------|-------|
| **Analysis Type** | Distribution estimation | Topic modeling |
| **Output** | TI, πDEG, enrichments | Topic probabilities per cell |
| **Handles imbalance** | ✅ Via power-independence | ✅ Via topic modeling |
| **Speed** | Fast | Slow (iterative topic fitting) |
| **Best for** | Quantitative comparison | Imbalanced samples, clustering |
| **Statistical model** | Empirical Bayes | Latent Dirichlet Allocation (LDA) |
| **Resolution** | Population-level | Cell-level (soft clustering) |
| **Strengths** | Power-independence, interpretable | Handles imbalance, soft clustering |
| **Weaknesses** | Pseudo-bulk only | Slow, requires topic number selection |

**When to use each:**
- **MUSIC:** Highly imbalanced samples, want cell-level topics
- **TRADE:** Quantitative effect magnitude
- **Both:** MUSIC for clustering → TRADE for quantification

**Relationship:** Complementary (different goals)

**Benchmarking:**
- MUSIC not in standard benchmarks (specialized tool)
- TRADE validated on standard datasets

---

### TRADE vs CINEMA-OT

| Aspect | TRADE | CINEMA-OT |
|--------|-------|-----------|
| **Analysis Type** | Distribution estimation | Causal inference (optimal transport) |
| **Output** | TI, πDEG, enrichments | Causal effects, matched cell pairs |
| **Handles confounding** | ⚠️ Assumes clean data | ✅ Explicitly models confounders |
| **Speed** | Fast | Moderate |
| **Best for** | Quantitative comparison | Confounding control |
| **Statistical model** | Empirical Bayes | ICA + optimal transport |
| **Causal interpretation** | ⚠️ Association only | ✅ Causal claims |
| **Strengths** | Power-independence | Confounder separation |
| **Weaknesses** | No causal claims | More complex, requires careful setup |

**When to use each:**
- **CINEMA-OT:** Confounding present, need causal interpretation
- **TRADE:** Clean experimental design, quantitative comparison
- **Both:** CINEMA-OT for causal effects → TRADE for quantification

**Relationship:** Complementary (address different concerns)

**Benchmarking:**
- CINEMA-OT specialized for causal inference (different goal)
- TRADE focused on power-independent estimation

---

## Comparison Matrix: All Methods

| Feature | TRADE | DESeq2 | MAST | Seurat | scVI | MUSIC | CINEMA-OT |
|---------|-------|--------|------|--------|------|-------|-----------|
| **Power-independent** | ✅ | ❌ | ❌ | ❌ | ⚠️ | ⚠️ | ⚠️ |
| **Gene lists** | ❌ | ✅ | ✅ | ✅ | ✅ | ⚠️ | ✅ |
| **Single-cell resolution** | ❌ | ❌ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Quantitative magnitude** | ✅ | ⚠️ | ⚠️ | ⚠️ | ⚠️ | ❌ | ✅ |
| **Speed** | Fast | Fast | Moderate | Fast | Slow | Slow | Moderate |
| **Memory** | Moderate | Low | Moderate | Moderate | High | High | Moderate |
| **GPU required** | ❌ | ❌ | ❌ | ❌ | ⚠️ | ❌ | ❌ |
| **Statistical rigor** | ✅ | ✅ | ✅ | ⚠️ | ⚠️ | ⚠️ | ✅ |
| **Interpretability** | ✅ | ✅ | ✅ | ✅ | ⚠️ | ⚠️ | ⚠️ |
| **Handles imbalance** | ✅ | ❌ | ❌ | ❌ | ⚠️ | ✅ | ⚠️ |
| **Confounding control** | ❌ | ⚠️ | ⚠️ | ⚠️ | ⚠️ | ❌ | ✅ |
| **Batch correction** | ⚠️ | ✅ | ⚠️ | ✅ | ✅ | ⚠️ | ⚠️ |
| **Subpopulation effects** | ❌ | ❌ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Cross-study comparison** | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| **Typical use case** | Meta-analysis | Gene discovery | Cell-type DE | Exploratory | Large data | Imbalanced | Causal |

**Legend:** ✅ Excellent, ⚠️ Moderate/Conditional, ❌ Poor/Not applicable

---

## Method Selection Guide

### Decision Tree

```
START: What is your primary goal?

├─ Identify specific genes to validate
│  └─ Use: DESeq2/MAST/Seurat (gene-level methods)
│
├─ Compare perturbations quantitatively
│  ├─ Same experiment, balanced power → DESeq2 sufficient
│  └─ Different contexts, imbalanced power → TRADE ★
│
├─ Detect cell-type-specific responses
│  └─ Use: MAST/scVI (single-cell methods)
│
├─ Handle severe sample imbalance
│  └─ Use: MUSIC or TRADE ★
│
├─ Control for confounding
│  └─ Use: CINEMA-OT
│
├─ Integrate large multi-batch datasets
│  └─ Use: scVI
│
└─ Comprehensive analysis
   └─ Use: Multiple methods (recommended)
       ├─ DESeq2 → gene lists
       ├─ TRADE → quantitative comparison ★
       └─ MAST/scVI → heterogeneity
```

---

## Performance Benchmarks (Literature Review)

### 2024 Genome Biology Study (CRISPR screens)

**Methods tested:** Seurat-Wilcox, Seurat-NB, MAST, DESeq2, edgeR, limma, t-test

**Key findings:**
- **Best performers:** Seurat-Wilcox, Seurat-NB
- **Still miscalibrated:** All methods showed signs of p-value inflation
- **MAST issues:** Inflated p-values, ~2000 false rejections
- **Pseudo-bulk advantage:** DESeq2/edgeR/limma outperformed many single-cell methods

**TRADE position:** Not tested (different analysis type: distribution vs gene-level)

### 2024 Nested Single-Cell Study

**Methods tested:** t-test, MAST, scVI, DESeq2, DREAM, Seurat, MAST, others

**Key findings:**
- **Best performers:** t-test (AUPRC=0.302), MAST (0.282)
- **Close behind:** scVI, DESeq2, DREAM
- **No advantage:** Single-cell methods didn't outperform pseudo-bulk
- **Recommendation:** Pseudo-bulk + DESeq2/edgeR/limma if sensitivity important

**TRADE position:** Not tested, but uses pseudo-bulk (consistent with recommendation)

### Nadig et al. 2025 (TRADE validation)

**TRADE-specific benchmarks:**
- **Replicability:** Between-replicate TI r=0.90 (vs gene count r=0.16)
- **Cross-dataset:** K562 genome-wide vs essential R²=59.7% (vs 28.4% for counts)
- **Downsampling:** TI stable at 50% power reduction
- **Simulations:** TI unbiased in 100 replicates

**Competitive advantage:** TRADE far superior for cross-study comparison

---

## Unique Value Propositions

### TRADE
**Unique:** Power-independent metrics for fair comparison
**Value:** Only method enabling cross-study meta-analysis
**Niche:** Quantitative comparison across imbalanced datasets

### DESeq2
**Unique:** Most mature, validated pseudo-bulk method
**Value:** Gold standard for gene discovery
**Niche:** Standard DE analysis with well-calibrated p-values

### MAST
**Unique:** Explicitly models zero-inflation (hurdle model)
**Value:** Cell-type-specific DE detection
**Niche:** Single-cell DE with covariate adjustment

### scVI
**Unique:** Deep learning for complex batch effects
**Value:** Latent space useful for visualization
**Niche:** Large datasets with severe batch effects

### MUSIC
**Unique:** Topic modeling for imbalanced samples
**Value:** Soft clustering of perturbations
**Niche:** Extremely imbalanced Perturb-seq data

### CINEMA-OT
**Unique:** Causal inference with optimal transport
**Value:** Separates treatment from confounding
**Niche:** Observational studies with confounding

---

## Complementarity Analysis

### Natural Combinations

**1. DESeq2 + TRADE (Recommended)**
- DESeq2: Gene discovery
- TRADE: Quantitative comparison
- **Synergy:** Gene lists + effect magnitude

**2. MAST + TRADE**
- MAST: Cell-type-specific effects
- TRADE: Population-level quantification
- **Synergy:** Heterogeneity + overall impact

**3. scVI + TRADE**
- scVI: Batch correction, dimensionality reduction
- TRADE: Power-independent comparison
- **Synergy:** Clean embeddings + quantitative metrics

**4. CINEMA-OT + TRADE**
- CINEMA-OT: Causal effects
- TRADE: Effect magnitude
- **Synergy:** Causality + quantification

### Redundant Combinations

**DESeq2 + edgeR + limma:** All pseudo-bulk gene-level (pick one)

**MAST + Seurat-NB:** Both single-cell gene-level (pick one)

**scVI + other deep learning:** Diminishing returns

---

## Method Maturity Assessment

| Method | Maturity | Evidence | Adoption |
|--------|----------|----------|----------|
| **DESeq2** | ★★★★★ | 10+ years, extensive validation | Very high |
| **edgeR** | ★★★★★ | 10+ years, extensive validation | High |
| **MAST** | ★★★★☆ | 5+ years, well-validated | Moderate |
| **Seurat** | ★★★★★ | 5+ years, ecosystem | Very high |
| **scVI** | ★★★☆☆ | 3+ years, growing | Moderate |
| **TRADE** | ★★★☆☆ | <1 year (Nat Genet 2025) | Low (new) |
| **MUSIC** | ★★★☆☆ | 2+ years | Low |
| **CINEMA-OT** | ★★★☆☆ | <2 years (Nat Methods 2023) | Low |

**Maturity factors:**
- Years in use
- Number of citations
- Independent validation studies
- Community adoption

**TRADE status:** New method (2025) with strong initial validation, but limited field testing

---

## Computational Comparison

### Speed (1000 perturbations, 20,000 genes)

| Method | Time | Bottleneck |
|--------|------|------------|
| TRADE | ~1 hour | DESeq2 + ash/mashr |
| DESeq2 | ~1 hour | GLM fitting |
| MAST | ~3-5 hours | Hurdle model per gene |
| Seurat-Wilcox | ~30 min | Wilcoxon tests |
| scVI | ~5-10 hours | VAE training (GPU) |
| MUSIC | ~2-4 hours | Topic fitting |

**Note:** Times approximate, depend on hardware, parameters

### Memory (1000 perturbations, 20,000 genes, 100,000 cells)

| Method | RAM | Notes |
|--------|-----|-------|
| TRADE | ~3-5 GB | Pseudo-bulk + mixture models |
| DESeq2 | ~2-3 GB | Pseudo-bulk only |
| MAST | ~10-20 GB | All cells in memory |
| Seurat | ~10-15 GB | All cells + metadata |
| scVI | ~20-30 GB | Neural network + all cells |
| MUSIC | ~15-25 GB | Topic modeling on all cells |

### Scalability

| Method | Scales with | Bottleneck |
|--------|-------------|------------|
| TRADE | # perturbations (linear) | Parallelizable |
| DESeq2 | # perturbations (linear) | Parallelizable |
| MAST | # cells × # genes | Memory |
| Seurat | # cells | Memory |
| scVI | # cells × # latent dims | GPU memory |
| MUSIC | # cells × # topics | Computation |

**Best scalability:** TRADE, DESeq2 (pseudo-bulk reduces dimensionality)

---

## Recommendations for Agent

### For Comprehensive Evaluation

**Consider TRADE as:**
1. ✅ Unique for power-independent comparison
2. ✅ Complementary to gene-level methods (not competitive)
3. ✅ Best-in-class for cross-study meta-analysis
4. ⚠️ Not suitable for gene discovery
5. ⚠️ New method (limited field testing)

**Include in comparison because:**
- Fills unique niche (only power-independent method)
- Strong theoretical foundation (Empirical Bayes)
- Exceptional validation (Nature Genetics 2025)
- Addresses real problem (power heterogeneity)

**Evaluate alongside:**
- DESeq2 (gene-level standard)
- MAST (single-cell representative)
- scVI (deep learning representative)
- MUSIC (topic modeling alternative)

### Rating Summary

| Criterion | TRADE Score | Rationale |
|-----------|-------------|-----------|
| **Novelty** | ★★★★★ | Unique power-independence |
| **Rigor** | ★★★★★ | Strong theory |
| **Utility** | ★★★★☆ | High for specific use cases |
| **Flexibility** | ★★★☆☆ | Limited scope |
| **Adoption** | ★★☆☆☆ | Too new |
| **Documentation** | ★★★★☆ | Good vignette + Nature Genetics paper |

**Overall:** ★★★★☆ (4/5)

**Best for:** Power-independent comparison, cross-study meta-analysis
**Not for:** Gene discovery, single-cell heterogeneity

---

**Agent instruction:** Evaluate TRADE as complementary to gene-level methods, not as alternative. Its unique value is power-independent comparison, which no other method provides.
