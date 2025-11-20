# Biological Context: Perturb-seq and the Need for TRADE

---

## What is Perturb-seq?

### Technology Overview

**Perturb-seq** (Perturbation sequencing) combines CRISPR genetic perturbations with single-cell RNA sequencing to measure transcriptional responses to thousands of genetic perturbations simultaneously.

**Key components:**
1. **CRISPR system:** CRISPRi (interference) or CRISPRko (knockout)
2. **sgRNA library:** Each sgRNA targets a specific gene
3. **Cellular barcoding:** Each sgRNA has unique barcode
4. **Single-cell RNA-seq:** 10x Genomics or similar platform
5. **Barcode readout:** Identify which perturbation in each cell

### Experimental Workflow

```
1. Engineer cells with CRISPR machinery (Cas9/dCas9-KRAB)
   ↓
2. Transduce with pooled sgRNA library
   ↓
3. Perturbations take effect (days to weeks)
   ↓
4. Single-cell RNA-seq + Feature Barcoding
   ↓
5. Bioinformatics: Assign cells to perturbations
   ↓
6. Differential expression analysis per perturbation
```

### Data Structure

**Output:** Matrix of cells × genes + perturbation assignments

| Cell | Gene1 | Gene2 | ... | GeneN | Perturbation | Batch |
|------|-------|-------|-----|-------|--------------|-------|
| Cell_1 | 0 | 3 | ... | 1 | GATA1_sgRNA1 | Batch1 |
| Cell_2 | 1 | 0 | ... | 0 | Control | Batch1 |
| Cell_3 | 2 | 5 | ... | 2 | MED12_sgRNA2 | Batch2 |
| ... | ... | ... | ... | ... | ... | ... |

**Typical scale:**
- **Cells:** 10,000-500,000 per experiment
- **Genes:** 15,000-25,000 (depending on filter cutoffs)
- **Perturbations:** 100-5,000 (depends on library design)
- **Cells per perturbation:** 10-1,000+ (highly variable!)

---

## The Central Challenge: Power Heterogeneity

### Problem Statement

**Not all perturbations are equal in statistical power:**

| Perturbation | Cells | Power | FDR-sig genes |
|--------------|-------|-------|---------------|
| Gene A (popular) | 856 | High | 234 |
| Gene B (moderate) | 124 | Medium | 67 |
| Gene C (rare) | 23 | Low | 8 |

**Question:** Does Gene A truly affect more genes than Gene C, or does it just have more cells (better power)?

**Traditional answer:** Inconclusive. Can't separate true biology from statistical power.

**TRADE's answer:** Estimate the underlying effect size distribution for each. TI and πDEG adjust for power differences.

### Sources of Power Heterogeneity

**1. Biological factors:**
- **Toxicity:** Perturbations causing cell death → fewer cells captured
- **Growth advantage:** Some perturbations increase proliferation → more cells
- **Example:** Essential gene knockout → cells die → low cell count

**2. Technical factors:**
- **sgRNA transduction efficiency:** Variable across library
- **Batch effects:** Different sequencing runs have different depths
- **Dropout:** Random sampling of cells in scRNA-seq
- **Library design:** Some sgRNAs included more frequently

**3. Experimental design:**
- **Targeted screens:** Some genes deliberately over-represented
- **Time-series:** Early timepoints have more cells (before toxic effects)
- **Multi-condition:** Some conditions have deeper sequencing

### Consequences of Power Heterogeneity

**Biased comparisons:**
- High-power perturbations appear to have "more" effects
- Low-power perturbations appear "uninteresting"
- Cannot compare across studies with different cell counts

**Missed discoveries:**
- True effects below detection threshold
- Particularly problematic for subtle but biologically important perturbations

**Statistical paradoxes:**
- Same gene, different experiments → different # significant genes
- Doesn't mean biology changed—just power changed

---

## Biological Questions Perturb-seq Addresses

### 1. Gene Function Discovery

**Question:** What does gene X do?

**Traditional approach:**
- Knockout gene X → List differentially expressed genes
- Interpret via pathway enrichment

**TRADE addition:**
- TI quantifies overall transcriptional impact
- Me estimates scope of effects (45 genes vs 500 genes)
- Enrichment analysis: Which pathways respond?

**Example from Nadig et al.:**
- **GATA1** (erythroid TF): Large TI in K562 (erythroid-like), small in RPE1
- **Interpretation:** Lineage-specific function confirmed

### 2. Essential vs Non-Essential Genes

**Question:** What distinguishes essential genes transcriptionally?

**Finding:** Essential genes affect >500 genes on average (vs ~45 for typical genes)

**Biological interpretation:**
- Essential genes are in central pathways (replication, transcription, translation)
- Perturbation cascades broadly through cellular machinery
- Explains why they're essential—too many dependencies

**TRADE advantage:** Can compare essential vs non-essential fairly despite power differences (essential genes may have fewer cells due to toxicity)

### 3. Cell-Type-Specific Gene Functions

**Question:** Does gene X have same function in all cell types?

**Approach:** Perturb X in multiple cell types, compare TI correlations

**Finding (Nadig et al.):**
- 56% of genes: Universal (high correlation across cell types)
- 44% of genes: Cell-type-specific (low correlation across types)

**Example:**
- **GATA1:** K562 (erythroid) TI = 0.45, RPE1 (epithelial) TI = 0.02
- **RPS3:** Ribosomal protein, TI ≈ 0.3 in all cell types (universal)

**Biological insight:** Can systematically classify genes as universal vs context-dependent functions.

### 4. Dosage-Response Relationships

**Question:** Does stronger perturbation just amplify effects, or change biology?

**Experimental design:** Titration series (weak CRISPRi guides, dTAG degrons)

**Finding (Nadig et al.):** Three patterns:
1. **Constant kinetics** (e.g., BCR): Linear dose-response
2. **Gradient** (e.g., ATP5E): Smooth transcriptional shift
3. **Threshold** (e.g., RAN, Polycomb): Abrupt change at critical dosage

**Biological interpretation:**
- Some genes: Proportional responses (homeostatic buffering)
- Others: Bistable switches or cooperative binding
- Haploinsufficient genes show thresholds

**TRADE advantage:** Can fairly compare weak vs strong perturbations (power differs due to magnitude of effect)

### 5. Genetic Interaction Networks

**Question:** Which perturbations have similar downstream effects?

**Approach:** Compute TI correlation matrix across all perturbations

**Use:**
- Cluster perturbations by transcriptional similarity
- Infer functional modules
- Predict genetic interactions

**Example:** Polycomb complex members (EED, SUZ12, EZH2) show high TI correlations (r>0.8) → confirming they function together.

### 6. Disease-Relevant Pathways

**Question:** Do disease risk genes converge on shared pathways?

**Application:** PsychENCODE reanalysis (Nadig et al.)

**Finding:**
- Autism, bipolar, schizophrenia, depression: High TI correlations
- Despite different GWAS hits (low genetic correlation)
- Suggests convergence on downstream transcriptional programs

**Biological implication:** "Final common pathway" hypothesis for psychiatric disorders.

---

## Why Traditional DE Analysis is Insufficient

### Limitation 1: Gene Lists are Power-Dependent

**Problem:** Number of significant genes is confounded by power.

**Example:**
```
Gene A knockout (n=500 cells): 234 FDR-sig genes
Gene B knockout (n=50 cells):  23 FDR-sig genes

Cannot conclude: A affects 10× more genes than B
Could be: A has 10× more power, true effects similar
```

**Traditional workaround:** Only compare within same experiment (same power)
**TRADE solution:** TI and πDEG are power-independent

### Limitation 2: Thresholds are Arbitrary

**Problem:** FDR 5% cutoff is arbitrary. Effect at p=0.049 vs p=0.051?

**Consequence:**
- Small changes in power → genes cross threshold
- Gene list instability
- Difficult to interpret "borderline" genes

**TRADE solution:** Uses full distribution (no threshold)

### Limitation 3: Hidden Signal Below Threshold

**Finding (Nadig et al.):** 64-87% of transcriptome-wide variance is in non-significant genes

**Biological reality:**
- Many genes have small but real effects (e.g., 0.1 log2FC)
- Individually non-significant
- Collectively substantial signal

**Example:**
```
Gene A: β = 0.5, p = 1e-10 → Significant
Gene B: β = 0.1, p = 0.3   → Non-significant (but true effect!)
Gene C: β = 0.1, p = 0.4   → Non-significant (but true effect!)
...
Gene Z: β = 0.1, p = 0.2   → Non-significant (but true effect!)

Traditional: Focuses on Gene A only
TRADE: Captures collective signal from B-Z
```

**TRADE advantage:** ash shrinkage borrows strength across genes to detect collective subtle effects

### Limitation 4: Cannot Compare Across Studies

**Problem:** Same gene, different experiments → incomparable results

**Example:**
```
Study 1 (2018, n=100 cells/pert): GATA1 → 45 sig genes
Study 2 (2023, n=500 cells/pert): GATA1 → 178 sig genes

Question: Did biology change? Or just power?
```

**Traditional approach:** Meta-analysis (combine p-values)—but assumes independence
**TRADE solution:** Compare TI directly (power-independent)

---

## What Perturb-seq Can vs Cannot Answer

### Perturb-seq Strengths

✅ **Systematic functional genomics**
- Screen thousands of genes in parallel
- Cost-effective vs one-by-one knockouts

✅ **Transcriptome-wide readout**
- Not limited to candidate genes
- Discover unexpected pathways

✅ **Quantitative dose-response**
- Graded perturbations (CRISPRi, degrons)
- Measure response curves

✅ **Cell-type resolution**
- Can profile heterogeneous cell populations
- Identify cell-type-specific functions (with appropriate analysis)

✅ **Temporal dynamics**
- Time-course experiments
- Measure immediate vs delayed responses

### Perturb-seq Limitations

❌ **Pseudo-bulk averages**
- Most analysis (including TRADE) uses pseudo-bulk
- Misses cell-to-cell variability within perturbations
- Cannot detect subpopulation-specific responses without specialized methods

❌ **Transcriptional proxy**
- Measures mRNA, not protein or function
- Post-transcriptional regulation not captured
- Phenotypic outcomes require separate assays

❌ **Acute perturbations**
- Typically days-weeks duration
- Chronic/developmental effects may differ
- Adaptation over longer timescales not measured

❌ **Artificial system**
- CRISPR perturbations ≠ natural genetic variation
- Complete knockouts may not reflect hypomorphic alleles
- CRISPRi may have off-target effects

❌ **Scalability trade-offs**
- More perturbations → fewer cells per perturbation → less power
- Deep profiling (many cells) limits # perturbations

---

## Why TRADE Matters for Biology

### Enables New Types of Questions

**1. Cross-context comparison:**
- "Does this gene have same function in disease vs healthy?"
- "Is dosage-response linear or threshold?"
- "Which genes are universally important vs cell-type-specific?"

**Previously impossible due to power confounding. TRADE enables fair comparison.**

**2. Quantitative perturbation characterization:**
- "How much does this perturbation affect the transcriptome?"
- "Is this a subtle or dramatic effect?"
- "Do these two perturbations have similar or distinct effects?"

**TI provides quantitative answer in interpretable units (log2FC²).**

**3. Low-power signal extraction:**
- "This pilot screen has only 30 cells/perturbation. Can we still learn something?"
- "Can we identify high-impact genes even without full significance?"

**TRADE recovers hidden signal that traditional methods miss.**

### Changes Experimental Design Considerations

**Before TRADE:**
- Needed balanced cell counts across all perturbations
- Low-power perturbations essentially uninformative
- Cross-study comparison nearly impossible

**With TRADE:**
- Imbalanced designs OK (TI adjusts for power)
- Can analyze low-power perturbations meaningfully
- Can integrate data across studies/labs

**Practical impact:**
- More cost-effective screens (don't need deep per-perturbation sampling)
- Can analyze "failed" experiments with imbalanced data
- Enables meta-analysis across Perturb-seq studies

---

## Biological Validation of TRADE Metrics

### TI Correlates with Known Biology

**Essential genes:**
- High TI (affect many processes)
- 4.22× impact enrichment
- ✅ Consistent with central cellular roles

**Tissue-specific transcription factors:**
- High TI in relevant cell type, low in others
- Example: GATA1 in erythroid cells
- ✅ Validates lineage-specific functions

**Housekeeping genes:**
- Moderate-high TI across all cell types (universal effects)
- ✅ Consistent with conserved functions

### πDEG Reflects Perturbation Scope

**Typical gene:** πDEG ≈ 45
- Affects specific pathway
- Local transcriptional changes

**Essential gene:** πDEG > 500
- Affects >10 pathways
- Broad transcriptional reprogramming

**Validates biological intuition:** Essential genes have broader impact.

### Enrichments Match Known Pathways

**Perturbing erythroid TFs → enriched response in erythroid genes**
- TRADE enrichment: 6.8×
- ✅ Validates pathway-specific effects

**Perturbing ribosomal proteins → enriched response in stress pathways**
- TRADE enrichment: 3.2×
- ✅ Validates ribosome stress response

### TI Correlation Predicts Genetic Interactions

**Polycomb complex members (EED, SUZ12, EZH2):**
- TI correlation: r > 0.8
- ✅ Known to function together

**Unrelated genes:**
- TI correlation: r ≈ 0.1-0.3
- ✅ Expected for independent functions

---

## Integration with Biological Knowledge

### Complementary to Other Assays

**TRADE + Fitness screens:**
- Fitness: Which genes affect growth?
- TRADE: What transcriptional changes mediate fitness effects?
- **Combined insight:** Mechanism of essentiality

**TRADE + Protein measurements:**
- Protein: Post-transcriptional effects
- TRADE: Transcriptional responses
- **Combined insight:** Transcriptional vs translational control

**TRADE + Phenotypic screens:**
- Phenotype: Cellular outcomes (morphology, function)
- TRADE: Transcriptional predictors of phenotype
- **Combined insight:** Transcriptional signatures of phenotypes

### Future: Multi-omic Integration

**Perturb-seq + ATAC-seq (chromatin):**
- TRADE for transcription
- Chromatin accessibility changes
- **Insight:** Regulatory mechanisms

**Perturb-seq + CITE-seq (surface proteins):**
- TRADE for transcriptome
- Surface protein changes
- **Insight:** Cell state transitions

---

## Summary: Why Biologists Should Care

**Traditional question:** "Which genes changed?"
- Gives list for follow-up
- Power-dependent
- Threshold-sensitive

**TRADE question:** "What is the distribution of changes?"
- Gives quantitative effect magnitude
- Power-independent
- Enables cross-study comparison

**Key insight:** Both questions are important and complementary!

**Use traditional DE for:**
- Identifying specific genes for validation
- Pathway enrichment (gene-level)
- Immediate experimental follow-up

**Use TRADE for:**
- Comparing perturbations fairly
- Quantifying overall transcriptional impact
- Cross-study/cross-context comparison
- Low-power signal extraction
- Systematic perturbation classification

**Bottom line:** TRADE enables biologists to ask quantitative questions about perturbations that were previously impossible to answer rigorously.

---

**Next:** See `02_MATHEMATICAL_FOUNDATIONS.md` for the statistical theory underlying these biological applications.
