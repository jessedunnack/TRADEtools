# TRADEtools: Executive Summary for Comparative Analysis

**Report Date:** 2025-11-18
**Package Version:** 0.99.0
**For:** Agent-based comparative analysis of Perturb-seq methods

---

## One-Sentence Summary

TRADEtools estimates the transcriptome-wide distribution of true differential expression effects using Empirical Bayes methods (adaptive shrinkage), producing power-independent metrics for fair comparison of genetic perturbations across experimental contexts.

---

## Quick Facts

| Aspect | Details |
|--------|---------|
| **Method Type** | Distribution estimation (not gene-level testing) |
| **Statistical Approach** | Empirical Bayes (ash/mashr packages) |
| **Input** | Pseudo-bulk DE summary statistics (log2FC, SE, p-value) |
| **Output** | Transcriptome-wide impact, effective DEG count, enrichments, correlations |
| **Language** | R |
| **Dependencies** | ashr, mashr, DESeq2 (recommended), ggplot2, doBy |
| **Computational Speed** | Fast (~seconds per perturbation for univariate) |
| **Memory Requirements** | Moderate (~GB scale for large screens) |
| **Unique Innovation** | Power-independent metrics enabling cross-study comparison |

---

## Core Innovation

**Problem:** Traditional significance testing produces biased results when perturbations have different statistical power (cell counts, sequencing depth).

**Solution:** Estimate the full distribution of effects using adaptive shrinkage, then compute distributional properties (variance, kurtosis) that are robust to power differences.

**Analogy:** Like GWAS (genome-wide association studies) methods that estimate heritability rather than listing significant SNPs.

---

## Key Metrics

### 1. Transcriptome-wide Impact (TI)
- **Definition:** Variance of effect size distribution
- **Units:** (log2FC)²
- **Interpretation:** Total transcriptional effect magnitude
- **Power-dependence:** ✅ Power-independent

### 2. Effective Number of DEGs (Me/πDEG)
- **Definition:** 3M/κ where κ = normalized kurtosis
- **Interpretation:** How many genes are "effectively" affected
- **Power-dependence:** ✅ Power-independent (but requires non-zero TI)

### 3. Gene Set Enrichment
- **Two types:** Response (how gene set responds) vs Impact (effect of perturbing gene set)
- **Power-dependence:** ✅ Power-independent

### 4. TI Correlation (bivariate)
- **Definition:** Correlation of true effects between two perturbations
- **Interpretation:** Transcriptional similarity
- **Power-dependence:** ✅ Power-independent

---

## Methodological Classification

### Workflow Position
**Stage:** Post-differential-expression analysis (operates on summary statistics)

```
Single-cell data → Pseudo-bulk aggregation → DESeq2 → TRADE → Distribution parameters
                                                    ↓
                                            Other DE tools
                                            (must validate SEs)
```

### Analysis Philosophy

| Traditional Approach | TRADE Approach |
|---------------------|----------------|
| Which genes changed? | What is the distribution of changes? |
| Gene-level significance | Distribution-level parameters |
| Threshold-dependent | Threshold-free |
| Power-sensitive | Power-independent |
| Compare gene lists | Compare effect size distributions |

---

## Competitive Positioning

### Unique Strengths
1. ✅ **Only method** with demonstrated power-independence for cross-study comparison
2. ✅ **Recovers hidden signal** below significance thresholds (64-87% of total signal)
3. ✅ **Stable metrics** robust to downsampling
4. ✅ **Theory-grounded** (Empirical Bayes with formal statistical properties)
5. ✅ **Fast computation** (no MCMC, uses efficient EM algorithm)

### Limitations
1. ❌ **Pseudo-bulk only** (misses cell-type heterogeneity within perturbations)
2. ❌ **Cannot identify specific genes** (outputs distributions, not gene lists)
3. ❌ **Requires calibrated SEs** (validated for DESeq2, others need validation)
4. ❌ **πDEG unstable** when TI ≈ 0 (low signal cases)

---

## Evidence Base

### Publication
- **Journal:** Nature Genetics 57(5):1228-1237 (May 2025)
- **Authors:** Nadig, Replogle, Pogson, McCarroll, Weissman, Robinson, O'Connor
- **DOI:** 10.1038/s41588-025-02169-3

### Validation
- ✅ Extensive simulations (100 replicates each, multiple scenarios)
- ✅ Negative controls (non-targeting guides, replicates)
- ✅ Cross-dataset consistency (5 large Perturb-seq datasets)
- ✅ Independent replication (between-replicate r=0.90)
- ✅ Cross-platform validation (microarray vs RNA-seq: r=0.78-0.96)

### Data Scale
- **Analyzed:** 9,866 genetic perturbations across K562, RPE1, Jurkat, HepG2
- **Largest dataset:** K562 genome-wide screen with ~2,000 perturbations
- **Cell counts:** 10-1000+ cells per perturbation

---

## When to Use TRADE vs Alternatives

### Use TRADE when:
1. **Comparing perturbations with different power** (different cell counts, depths)
2. **Cross-study comparison** (same gene, different contexts)
3. **Low-power pilot screens** (want to extract maximal signal from noisy data)
4. **Quantifying overall effect** (how much does this perturbation affect transcriptome?)
5. **Gene set analysis** (pathway enrichments)

### Use alternatives when:
1. **Need specific gene lists** for follow-up experiments → Traditional DE
2. **Cell-type-specific effects** within perturbations → Single-cell DE (MAST, scVI)
3. **Very sparse effects** (<10 genes) → Gene-level testing may be more powerful
4. **Network inference** or causal relationships → Specialized tools (CINEMA-OT, MUSIC)

---

## Comparison to Other Methods (Summary)

| Method | Type | Power-Independent? | Speed | Best Use Case |
|--------|------|-------------------|-------|---------------|
| **TRADE** | Distribution estimation | ✅ Yes | Fast | Cross-study comparison |
| **DESeq2** | Gene-level testing | ❌ No | Fast | Gene discovery |
| **MAST** | Gene-level testing (single-cell) | ❌ No | Moderate | Cell-type-specific DE |
| **Seurat** | Gene-level testing | ❌ No | Fast | General scRNA-seq |
| **MUSIC** | Topic modeling | ⚠️ Partial | Slow | Imbalanced samples |
| **CINEMA-OT** | Causal matching | ⚠️ Partial | Moderate | Confounding control |
| **scVI** | Deep learning DE | ❌ No | Slow (GPU) | Large datasets |

**Key distinction:** TRADE is the only method focused on power-independent distribution estimation rather than gene-level testing.

---

## Major Findings from Publication

1. **Hidden signal:** Only 13-36% of transcriptome-wide impact appears in FDR-significant genes
2. **Essential genes:** Affect >500 genes on average (vs ~45 for typical genes)
3. **Cell-type specificity:** 44% of perturbations show cell-type-dependent effects
4. **Dosage nonlinearity:** Weak and strong perturbations have qualitatively different effects
5. **Disease convergence:** Neuropsychiatric disorders share downstream transcriptional programs

---

## Computational Characteristics

### Runtime
- **Univariate:** ~1-5 seconds per perturbation (1,000-20,000 genes)
- **Bivariate:** ~10-60 seconds per pair (depends on covariance matrix strategy)
- **Bottleneck:** DESeq2 pseudo-bulk analysis (minutes per perturbation)

### Memory
- **Univariate:** <1 GB per perturbation
- **Bivariate (combined strategy):** 2-5 GB (many covariance matrices)
- **Scalability:** Linear with number of perturbations (embarrassingly parallel)

### Dependencies
- **R version:** ≥3.6
- **Critical packages:** ashr (≥2.2), mashr (≥0.2.79)
- **System requirement:** GNU Scientific Library (GSL) for mashr

---

## Critical Assessment for Comparative Analysis

### Where TRADE Excels
1. **Cross-context comparison** is impossible with other methods (power confounding)
2. **Signal recovery** from low-power experiments outperforms alternatives
3. **Theoretical foundation** provides interpretable, unbiased estimates
4. **Replicability** (r=0.90) far exceeds gene-level methods (r=0.16)

### Where TRADE Struggles
1. **Gene-level resolution:** Cannot provide specific gene lists
2. **Heterogeneity:** Misses cell-type-specific responses within perturbations
3. **Sparse effects:** πDEG estimation unstable when few genes affected
4. **DE tool dependence:** Requires well-calibrated standard errors

### Complementary vs Competitive
**TRADE is complementary to gene-level methods**, not competitive:
- Use DESeq2/MAST/etc. to identify specific genes → Follow up experimentally
- Use TRADE to quantify overall impact → Compare across studies, contexts, dosages

**Analogy:** Gene-level DE is like listing individual SNPs, TRADE is like estimating heritability.

---

## Bottom Line for Agent

**TRADE fills a unique niche:** It's the only method enabling fair comparison of perturbations with different statistical power. This is critical for:
- Cross-cell-type studies
- Dosage-response experiments
- Meta-analysis of multiple screens
- Low-power pilot studies

**Trade-offs:**
- ✅ Gain: Power-independent, stable metrics
- ❌ Lose: Gene-level specificity

**Recommendation:** TRADE should be part of a comprehensive Perturb-seq analysis pipeline alongside traditional DE, not a replacement for it.

---

## Files in This Report

1. **00_EXECUTIVE_SUMMARY.md** (this file) - Quick overview
2. **01_BIOLOGICAL_CONTEXT.md** - Perturb-seq biology
3. **02_MATHEMATICAL_FOUNDATIONS.md** - Complete mathematical theory
4. **03_STATISTICAL_METHODS.md** - ash/mashr details
5. **04_IMPLEMENTATION_DETAILS.md** - Code architecture
6. **05_STRENGTHS_AND_WEAKNESSES.md** - Critical analysis
7. **06_METHOD_COMPARISON.md** - vs other Perturb-seq methods
8. **07_USE_CASES.md** - Application scenarios
9. **08_PERFORMANCE_BENCHMARKS.md** - Speed/memory/accuracy
10. **09_COMPLETE_FUNCTION_REFERENCE.md** - All functions
11. **10_INPUT_OUTPUT_SPECIFICATIONS.md** - Data formats
12. **11_REFERENCES_AND_CITATIONS.md** - Bibliography
13. **12_FUTURE_DIRECTIONS.md** - Limitations and extensions

**Read these files for complete understanding of TRADEtools in context of Perturb-seq method landscape.**

---

**Report Compiled By:** Claude (Sonnet 4.5)
**Source Materials:**
- Nadig et al. Nature Genetics 2025
- TRADEtools R package (v0.99.0)
- Comparative benchmarking studies (2024)
- 676 lines of source code analysis
