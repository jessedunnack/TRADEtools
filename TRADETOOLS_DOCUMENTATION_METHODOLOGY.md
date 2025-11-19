# Comprehensive Documentation Task: TRADEtools R Package
## Systematic Analysis Following Proven Methodology

**Repository:** TRADEtools (Nadig et al. Nature Genetics 2025)
**Target:** R package for transcriptome-wide analysis of differential expression
**Date:** 2025-11-18
**Branch:** claude/perturbseq-library-documentation-01SUE1wkcvpPxXo6Ag1m5gVm
**Methodology:** 6-phase systematic documentation approach (proven on Mixscale R package)

---

## Executive Summary

**Objective:** Create comprehensive, verified technical documentation for the TRADEtools R package following the systematic methodology that achieved 100% accuracy on Mixscale.

**Scope:**
- R package for estimating distribution of differential expression effects
- Located in `/home/user/TRADEtools/R/`
- Initial file inventory: 3 R source files (~676 total lines)
- Associated with Nadig et al. Nature Genetics 2025 publication

**Publication Reference:**
- **Title:** "Transcriptome-wide characterization of genetic perturbations"
- **Authors:** Ajay Nadig, Joseph M. Replogle, Angela N. Pogson, Steven A McCarroll, Jonathan S. Weissman, Elise B. Robinson, Luke J. O'Connor
- **bioRxiv:** 10.1101/2024.07.03.601903 (July 2024)
- **Published:** Nature Genetics 57(5):1228-1237 (May 2025), doi: 10.1038/s41588-025-02169-3

**Expected Deliverables:**
1. Comprehensive technical documentation (function-by-function)
2. Function inventory with full specifications
3. Algorithm descriptions with code verification
4. Workflow examples from vignette
5. Comprehensive audit report with certification
6. Supplemental documentation for any undocumented components

---

## Phase 0: Initial Discovery & Context Gathering

### Tasks

#### 0.1 Read Core Package Files (COMPLETED)

**Already Read:**
- ✅ README.md - Package overview
- ✅ DESCRIPTION - Package metadata
- ✅ NAMESPACE - Exported functions
- ✅ vignettes/TRADEtools-intro.Rmd - Tutorial and examples
- ✅ R/TRADE.R - Main wrapper function

**Questions Answered:**
- **What does TRADEtools do?** Estimates the transcriptome-wide distribution of true differential expression effects, accounting for estimation error. Produces metrics like "transcriptome-wide impact" and "effective number of DEGs".
- **Core capabilities:**
  - Univariate analysis: Estimate distribution of DE effects for single perturbation
  - Bivariate analysis: Estimate joint distribution and correlation of DE effects between two perturbations
  - Gene set enrichment analysis
  - Significance testing (Bonferroni and FDR)
- **Dependencies:** ashr, mashr, doBy, ggplot2
- **Integration:** Works with DESeq2 output (or similar DE tools), integrates with Perturb-seq workflows

#### 0.2 Documentation Resources (COMPLETED)

**Identified Resources:**

1. ✅ **Publication:** Nadig et al. Nature Genetics 2025
   - bioRxiv preprint accessible
   - Contains detailed methods section
   - Statistical models and algorithms described

2. ✅ **Tutorial:** vignettes/TRADEtools-intro.Rmd
   - Complete quick-start guide
   - Univariate and bivariate examples
   - Detailed parameter descriptions
   - Output interpretation guide

3. ✅ **Existing Documentation:**
   - man/TRADE.rd file exists
   - Roxygen2 comments in R/TRADE.R
   - Need to check other R files for documentation

4. ✅ **Examples:**
   - example/ directory exists
   - Real data from Replogle et al. 2020 (K562 Perturb-seq)
   - Usage patterns in vignette

#### 0.3 File Structure Inventory

**R Source Files:**
```
R/
├── TRADE.R              (116 lines)   - Main wrapper function, input validation
├── TRADE_bivariate.R    (270 lines)   - Bivariate TRADE analysis (mashr-based)
└── TRADE_univariate.R   (290 lines)   - Univariate TRADE analysis (ashr-based)
                         ─────────────
                         Total: 676 lines
```

**Function Inventory (Preliminary):**
1. **TRADE()** - Main exported wrapper function (TRADE.R:30)
2. **TRADE_bivariate()** - Bivariate analysis engine (TRADE_bivariate.R:1)
3. **TRADE_univariate()** - Univariate analysis engine (TRADE_univariate.R:1)
4. **fit_ash()** - Helper: Fit adaptive shrinkage model (TRADE_univariate.R:49)
5. **get_distribution_output()** - Helper: Process distribution estimates (TRADE_univariate.R:76)
6. **frac_subset()** - Helper: Calculate variance fractions (TRADE_univariate.R:245)

**Total Functions:** 6 (1 exported, 5 internal)

**Supporting Files:**
```
man/
└── TRADE.rd             (4.2K)        - Documentation for TRADE() function

vignettes/
├── TRADEtools-intro.Rmd (19K)        - Main tutorial
└── TRADEtools-intro.html (660K)      - Compiled tutorial

example/
└── [Example data files]
```

---

## Phase 1: Systematic Code Reading & Initial Documentation

### Module Structure

Given the package structure, organize documentation by functional modules:

**Module 1: Main Interface**
- TRADE.R
- Primary user-facing function
- Input validation and routing

**Module 2: Univariate Analysis Engine**
- TRADE_univariate.R
- Core univariate statistical methods
- Helper functions: fit_ash(), get_distribution_output(), frac_subset()

**Module 3: Bivariate Analysis Engine**
- TRADE_bivariate.R
- Core bivariate statistical methods
- Joint distribution estimation

### Documentation Template (R Package Specific)

**For each function:**

```markdown
### Function: function_name()

#### Classification
- Module: [file.R]
- Export Status: [EXPORTED (in NAMESPACE) / INTERNAL]
- Type: [Main API / Analysis Engine / Helper Function]
- File Location: [file:line]

#### Purpose
[1-2 sentence description of what the function does]

#### When to Use
[Bulleted list of use cases]

#### Algorithm
[Step-by-step breakdown of computational approach]
[Mathematical formulas if applicable]
[Reference to Nadig et al. methods section]

#### Function Signature
```r
function_name <- function(
  param1 = default1,    # Description
  param2 = default2,    # Description
  ...
)
```

#### Parameters
- **param1** (type): Description
  - Default: value
  - Constraints: [if any]
  - Required/Optional: [specify]

#### Return Value
**Type:** List

**Structure:**
- `element1`: Description
- `element2`: Description
- ...

#### Dependencies
- **Imports:** [packages used]
- **Calls:** [other functions called]
- **Called by:** [functions that call this]

#### Usage Examples
```r
# Example 1: Basic usage from vignette
result <- function_name(...)

# Example 2: Advanced usage
...
```

#### Statistical Methods
[For statistical functions:]
- **Model:** [Statistical model used]
- **Estimation method:** [E.g., Empirical Bayes, EM algorithm]
- **Reference:** [Nadig et al. section/equation]

#### Code Reference
- Implementation: [file:line-range]
- Tests: [if any]
- Vignette usage: [section of vignette]
```

### Systematic Reading Process

**For each R source file:**

1. **Read entire file** using Read tool
2. **Identify all functions:**
   ```bash
   grep -n "^[[:space:]]*[a-zA-Z_][a-zA-Z0-9_]*[[:space:]]*<-[[:space:]]*function" R/[file].R
   ```

3. **Classify components:**
   - Exported (in NAMESPACE) vs internal
   - Main API vs helper functions
   - Statistical methods vs utility functions

4. **Document each function** using template above

5. **Extract algorithms** from implementation
   - Identify statistical models (ashr, mashr calls)
   - Document preprocessing steps
   - Identify quality control filters
   - Cross-reference with Nadig et al. methods

6. **Extract Roxygen documentation** if present
   - @param descriptions
   - @return specifications
   - @examples
   - Verify accuracy against actual code

7. **Work incrementally** - Report progress after each module

**Reading Order:**
1. TRADE.R (main interface - understand API first)
2. TRADE_univariate.R (core method 1)
3. TRADE_bivariate.R (core method 2)

---

## Phase 2: Cross-References & Integration Documentation

### Tasks

#### 2.1 Workflow Diagrams

**Create flowcharts showing:**

1. **Typical univariate analysis pipeline**
```
DESeq2/other DE tool → Summary statistics → TRADE(mode="univariate") →
  → Transcriptome-wide impact + Me + Enrichments + Significant gene analysis
```

2. **Typical bivariate analysis pipeline**
```
DE results 1 + DE results 2 → TRADE(mode="bivariate") →
  → TI correlation + Effect size correlation matrix
```

3. **Internal function flow**
```
TRADE() → Input validation →
  ├─ mode="univariate" → TRADE_univariate() → fit_ash() + get_distribution_output() + frac_subset()
  └─ mode="bivariate"  → TRADE_bivariate() → mashr workflow
```

#### 2.2 Cross-Reference Tables

**Functions → Purpose:**
| Function | Type | Purpose | File |
|----------|------|---------|------|
| TRADE | Exported | Main user interface | TRADE.R |
| TRADE_univariate | Internal | Univariate analysis engine | TRADE_univariate.R |
| TRADE_bivariate | Internal | Bivariate analysis engine | TRADE_bivariate.R |
| fit_ash | Internal | Adaptive shrinkage fitting | TRADE_univariate.R |
| get_distribution_output | Internal | Extract distribution summaries | TRADE_univariate.R |
| frac_subset | Internal | Variance fraction calculations | TRADE_univariate.R |

**Parameters → Functions:**
| Parameter | Used In | Purpose |
|-----------|---------|---------|
| mode | TRADE | Select univariate vs bivariate |
| results1 | TRADE, TRADE_univariate, TRADE_bivariate | First DE results |
| results2 | TRADE, TRADE_bivariate | Second DE results (bivariate only) |
| annot_table | TRADE, TRADE_univariate | Gene set annotations |
| estimate_sampling_covariance | TRADE, TRADE_bivariate | Account for shared controls |
| ... | ... | ... |

**Statistical Methods → Publication:**
| Method/Algorithm | Code Location | Nadig et al. Reference | Purpose |
|------------------|---------------|----------------------|---------|
| Adaptive shrinkage (ash) | fit_ash() | Methods: "Statistical model" | Univariate distribution estimation |
| Multivariate adaptive shrinkage (mash) | TRADE_bivariate() | Methods: "Bivariate analysis" | Joint distribution estimation |
| Transcriptome-wide impact | get_distribution_output() | Methods: Definition | Variance of effect size distribution |
| Effective number of DEGs (Me) | get_distribution_output() | Methods: Definition | Kurtosis-based DEG count |
| ... | ... | ... | ... |

#### 2.3 Output Structure Documentation

**Univariate Output:**
```r
list(
  distribution_summary = list(
    transcriptome_wide_impact,
    Me,
    mean
  ),
  significant_genes_Bonferroni = list(...),
  significant_genes_FDR = list(...),
  annot_output = data.frame(...),
  fit = [ashr object],
  qc = list(...)
)
```

**Bivariate Output:**
```r
list(
  TI_correlation,
  correlation_matrix,
  covariance_matrix,
  cor_raw,
  fitted = [mash object],
  samples,
  V,
  loglik,
  runtime
)
```

#### 2.4 End-to-End Examples

**Provide complete working examples:**

1. **Minimal univariate example**
2. **Realistic univariate with enrichments** (from vignette)
3. **Bivariate correlation analysis** (from vignette)
4. **Interpreting output** - What each metric means

---

## Phase 3: Comprehensive Audit & Verification

### Tasks

#### 3.1 File Coverage Verification

**Execute:**
```bash
ls R/*.R
# Verify all .R files documented
```

**Deliverable:** ✅ Confirmation that 100% of source files examined (3/3)

#### 3.2 Component Inventory

**Extract all functions:**
```bash
# Find all function definitions
grep -n "^[[:space:]]*[a-zA-Z_][a-zA-Z0-9_]*[[:space:]]*<-[[:space:]]*function" R/*.R

# Check NAMESPACE exports
cat NAMESPACE
```

**Create inventory:**
1. Total functions: 6
2. Exported functions (public API): 1 (TRADE)
3. Internal functions: 5
4. Functions with Roxygen docs: [count after reading]
5. **GAP ANALYSIS:** [N-X] undocumented components

**Deliverable:** Complete inventory table

#### 3.3 Algorithm Cross-Check

**For each function:**

1. **Read actual implementation** (specific line ranges)
2. **Identify statistical method:**
   - ashr::ash() calls → Document parameters
   - mashr::mash() calls → Document covariance matrices
   - Custom calculations → Extract formulas
3. **Compare to Nadig et al. methods section**
4. **Verify parameter transformations** (e.g., log2FC filtering, SE handling)

**Verification Template:**
```markdown
### Verification: [function_name]

**Documented Algorithm:**
[What documentation says]

**Actual Code (file:lines):**
```r
[Actual code snippet]
```

**Statistical Method:**
- Package/function used: [e.g., ashr::ash]
- Key parameters: [...]
- Nadig et al. reference: [Methods section X]

**Status:** ✅ MATCH / ❌ DISCREPANCY

**Notes:** [Any differences]
```

**Priority verification targets:**
- fit_ash() - Core univariate method
- TRADE_bivariate() - Core bivariate method
- get_distribution_output() - Metric calculations
- Transcriptome-wide impact calculation
- Me (effective DEG) calculation
- Enrichment calculations

**Target:** 100% algorithmic accuracy

#### 3.4 Parameter Specification Verification

**For each documented function:**

1. **Read function signature** from source
2. **Extract parameter names, defaults**
3. **Compare to Roxygen @param tags** (if present)
4. **Compare to vignette descriptions**
5. **Check:**
   - All parameters documented? ✓
   - Correct default values? ✓
   - Type specifications accurate? ✓
   - Descriptions match behavior? ✓
   - Required vs optional correctly indicated? ✓

**Create verification table:**

| Function | Parameter | Roxygen Default | Code Default | Vignette Description | Match? |
|----------|-----------|----------------|--------------|---------------------|--------|
| TRADE | mode | NULL | NULL | Required: "univariate" or "bivariate" | ✅ |
| TRADE | covariance_matrix_set | "combined" | "combined" | Default "combined" | ✅ |
| ... | ... | ... | ... | ... | ... |

**Target:** 100% parameter accuracy

#### 3.5 Publication Concordance Check

**Read Nadig et al. Nature Genetics 2025 (or bioRxiv preprint):**

1. **Methods section** - Extract algorithmic descriptions
   - Statistical model formulation
   - Empirical Bayes approach
   - ashr and mash usage
   - Transcriptome-wide impact definition
   - Me (effective DEG) definition
   - Gene set enrichment approach

2. **Compare to code implementation**
   - Does code match described methods?
   - Are formulas implemented correctly?
   - Are preprocessing steps consistent?

3. **Verify documentation matches both code AND paper**

**Specific items to verify:**
- [ ] Transcriptome-wide impact = variance of effect size distribution
- [ ] Me calculation from kurtosis
- [ ] Gene set enrichment = var_annotation / var_total
- [ ] Bivariate correlation estimation method
- [ ] Quality control filters (|log2FC| > 10, NA removal)
- [ ] Bonferroni and FDR correction procedures
- [ ] Sampling covariance estimation (when enabled)

**Target:** 100% concordance with publication

---

## Phase 4: Gap Filling & Supplemental Documentation

### Tasks

#### 4.1 Prioritize Missing Components

**Classify by importance:**

**Priority 1 - CRITICAL:**
- TRADE() function (exported, main API) - **MUST** be fully documented
- Any undocumented parameters
- Any undocumented return values

**Priority 2 - HIGH:**
- TRADE_univariate() and TRADE_bivariate() (core engines)
- Statistical method implementations
- Output interpretation guidance

**Priority 3 - MEDIUM:**
- Helper functions (fit_ash, get_distribution_output, frac_subset)
- Internal data transformations
- QC procedures

#### 4.2 Document Missing Components

**For each gap identified:**
1. Read implementation carefully
2. Understand statistical method
3. Document using full template
4. Include code references
5. Add usage examples from vignette
6. Cross-reference to publication methods

**Output:** TRADETOOLS_DOCUMENTATION_SUPPLEMENT.md (if needed)

#### 4.3 Document Internal Helpers

**Even if not exported, document:**
- **fit_ash()**: How does it call ashr::ash? What parameters?
- **get_distribution_output()**: How are TI and Me calculated?
- **frac_subset()**: How are variance fractions computed?

---

## Phase 5: Audit Report & Certification

### Tasks

#### 5.1 Create Audit Report

**Structure:**

```markdown
# TRADEtools R Package Documentation - Comprehensive Audit Report

## Executive Summary
- Audit outcome: PASS/FAIL
- Key metrics
- Summary of findings

## Audit Methodology
- Phase 0: Context gathering (publication, vignette, package structure)
- Phase 1: Systematic documentation (3 R files, 6 functions)
- Phase 2: Integration analysis (workflows, cross-references)
- Phase 3: Verification (code, parameters, algorithms, publication)
- Phase 4: Gap filling
- Phase 5: Certification

## Findings
- Files examined: X/3
- Functions documented: X/6
- Exported functions: 1/1
- Internal functions: X/5
- Missing components discovered: N
- Errors found: M
- Discrepancies: K

## Verification Summary by Module

### Module 1: Main Interface (TRADE.R)
- Functions documented: X/1
- Algorithm accuracy: %
- Parameter accuracy: %
- Roxygen coverage: ✅/❌
- Status: ✅/⚠️/❌

### Module 2: Univariate Engine (TRADE_univariate.R)
- Functions documented: X/4
- Statistical methods verified: ✅/❌
- Publication concordance: %
- Status: ✅/⚠️/❌

### Module 3: Bivariate Engine (TRADE_bivariate.R)
- Functions documented: X/1
- Statistical methods verified: ✅/❌
- Publication concordance: %
- Status: ✅/⚠️/❌

## Publication Concordance
- Methods from Nadig et al. verified: X/Y
- Algorithms matching publication: %
- Statistical models verified: ✅/❌
- Status: ✅/⚠️

## Vignette Integration
- Examples working: ✅/❌
- Parameters match vignette: ✅/❌
- Output interpretation verified: ✅/❌

## Recommendations
1. Priority 1 items (critical)
2. Priority 2 items (important)
3. Future improvements

## Audit Certification
I certify that:
✅ All R source files examined (3/3)
✅ All functions inventoried (6/6)
✅ Algorithms verified against code
✅ Parameters cross-checked
✅ Public API fully documented (TRADE function)
✅ Statistical methods verified against Nadig et al.
✅ Vignette examples verified

**Auditor:** Claude (Sonnet 4.5)
**Date:** 2025-11-18
**Repository:** TRADEtools
**Branch:** claude/perturbseq-library-documentation-01SUE1wkcvpPxXo6Ag1m5gVm
**Commit:** [Git hash]
```

#### 5.2 Metrics Dashboard

| Metric | Value |
|--------|-------|
| Total R files | 3 |
| Total functions | 6 |
| Exported functions | 1 |
| Internal functions | 5 |
| Documentation coverage | X% |
| Roxygen documentation present | X/6 |
| Algorithmic accuracy | Y% |
| Parameter accuracy | Z% |
| Publication concordance | W% |
| Code lines audited | 676 |
| Vignette integration verified | ✅/❌ |

#### 5.3 Completion Checklist

```markdown
## Verification Checklist

### File Coverage
- [ ] R/TRADE.R examined and documented
- [ ] R/TRADE_univariate.R examined and documented
- [ ] R/TRADE_bivariate.R examined and documented
- [ ] NAMESPACE exports verified
- [ ] DESCRIPTION metadata reviewed

### Component Documentation
- [ ] TRADE() exported function fully documented
- [ ] TRADE_univariate() internal function documented
- [ ] TRADE_bivariate() internal function documented
- [ ] fit_ash() helper documented
- [ ] get_distribution_output() helper documented
- [ ] frac_subset() helper documented
- [ ] All parameters specified
- [ ] All return values documented

### Accuracy Verification
- [ ] Function signatures verified
- [ ] Algorithms cross-checked against code
- [ ] Parameters validated
- [ ] Default values confirmed
- [ ] Roxygen docs verified (if present)

### Statistical Methods
- [ ] ashr usage documented
- [ ] mashr usage documented
- [ ] Transcriptome-wide impact calculation verified
- [ ] Me (effective DEG) calculation verified
- [ ] Enrichment calculations verified
- [ ] Correlation estimation verified

### Publication Integration
- [ ] Nadig et al. methods cross-referenced
- [ ] Statistical models verified
- [ ] Formulas confirmed
- [ ] Quality control steps matched
- [ ] Citations included

### Vignette Integration
- [ ] Univariate example verified
- [ ] Bivariate example verified
- [ ] Parameter descriptions matched
- [ ] Output interpretation verified
- [ ] Usage patterns documented

### Deliverables
- [ ] Main documentation complete
- [ ] Audit report created
- [ ] Supplemental documentation (if needed)
- [ ] All committed to git
- [ ] All pushed to remote
```

---

## Phase 6: Deliverables & Handoff

### Final Outputs

#### 6.1 Documentation Files

**Primary Documentation:**
```
TRADETOOLS_COMPREHENSIVE_DOCUMENTATION.md
- Complete technical documentation
- Organized by 3 modules (Main Interface, Univariate, Bivariate)
- Function-by-function breakdown (6 functions total)
- Algorithm descriptions with statistical methods
- Usage examples from vignette
- Cross-reference tables
- Output structure documentation
```

**Audit Report:**
```
TRADETOOLS_DOCUMENTATION_AUDIT_REPORT.md
- Verification methodology
- Findings and gap analysis
- Metrics dashboard
- Publication concordance analysis
- Certification
```

**Supplemental Documentation** (if needed):
```
TRADETOOLS_DOCUMENTATION_SUPPLEMENT.md
- Advanced statistical details
- Internal implementation notes
- Extended examples
```

#### 6.2 Git Workflow

```bash
# All work on designated branch
git status

# Stage documentation files
git add TRADETOOLS_COMPREHENSIVE_DOCUMENTATION.md
git add TRADETOOLS_DOCUMENTATION_AUDIT_REPORT.md
git add TRADETOOLS_DOCUMENTATION_SUPPLEMENT.md  # if created
git add TRADETOOLS_DOCUMENTATION_METHODOLOGY.md  # this file

# Commit with descriptive message
git commit -m "$(cat <<'EOF'
docs: comprehensive TRADEtools R package documentation and audit

Complete systematic documentation of TRADEtools R package following
proven 6-phase methodology (previously applied to Mixscale R package).

Deliverables:
- Comprehensive technical documentation (3 modules, 6 functions)
- Full function inventory and specifications
- Algorithm verification against Nadig et al. Nature Genetics 2025
- Statistical methods verified (ashr, mashr)
- Vignette integration verified
- Complete audit report with certification

Metrics:
- Files documented: 3/3
- Functions documented: 6/6
- Algorithmic accuracy: 100%
- Parameter accuracy: 100%
- Publication concordance: 100%
- Audit status: PASS

References:
- Nadig et al. Nature Genetics 2025 (doi: 10.1038/s41588-025-02169-3)
- bioRxiv preprint: 10.1101/2024.07.03.601903
EOF
)"

# Push to remote
git push -u origin claude/perturbseq-library-documentation-01SUE1wkcvpPxXo6Ag1m5gVm
```

#### 6.3 Handoff Checklist

**Before declaring completion:**

- [ ] All documentation files created
- [ ] All files committed to git
- [ ] All commits pushed to remote
- [ ] Audit report complete with certification
- [ ] Metrics dashboard provided
- [ ] User informed of completion
- [ ] Location of all deliverables specified

---

## Success Metrics

### Quantitative Targets

✅ **100% file coverage** - All 3 R files documented
✅ **100% function coverage** - All 6 functions documented
✅ **100% exported function coverage** - TRADE() fully documented
✅ **100% algorithmic accuracy** - All algorithms verified against code
✅ **100% parameter accuracy** - All parameter specs verified
✅ **100% publication concordance** - Methods match Nadig et al.

### Qualitative Targets

✅ **Comprehensiveness** - User can understand package from docs alone
✅ **Accuracy** - No misleading or incorrect information
✅ **Usability** - Clear examples from vignette
✅ **Traceability** - Code references enable easy source lookup
✅ **Reproducibility** - Another analyst can verify the audit
✅ **Statistical rigor** - Methods properly explained with formulas

---

## R Package Specific Considerations

### R Package Elements

| R Package Component | Location | Documentation Need |
|---------------------|----------|-------------------|
| DESCRIPTION | Root | Metadata (already documented) |
| NAMESPACE | Root | Exports (simple - exportPattern) |
| R/ directory | R/*.R | **PRIMARY FOCUS** - all functions |
| man/ directory | man/*.rd | Existing Roxygen docs to verify |
| vignettes/ | vignettes/*.Rmd | Usage examples to reference |
| example/ | example/ | Test data (document usage) |

### Roxygen2 Documentation

**Check for and verify:**
```r
#' Title
#'
#' Description
#'
#' @param param1 description
#' @return return value
#' @export
#' @examples
```

**Compare Roxygen tags to:**
- Actual function signature
- Actual behavior in code
- Vignette descriptions
- User expectations

### Statistical Package Specifics

**Document statistical methods thoroughly:**
1. **Model specification** - What statistical model?
2. **Estimation method** - Maximum likelihood? Empirical Bayes? EM algorithm?
3. **Prior specification** - For Bayesian methods (ashr/mashr)
4. **Inference** - How are estimates derived?
5. **Uncertainty quantification** - Standard errors? Posterior distributions?

**For TRADEtools specifically:**
- ashr (Adaptive SHrinkage in R) - Empirical Bayes method
- mashr (Multivariate Adaptive SHrinkage in R) - Multivariate extension
- Transcriptome-wide impact - Variance estimate
- Me (effective DEG) - Kurtosis-based metric

---

## Tools & Commands

### R Package Analysis

```bash
# List R source files
ls -lh R/*.R

# Count lines of code
wc -l R/*.R

# Find all function definitions
grep -n "^[[:space:]]*[a-zA-Z_][a-zA-Z0-9_]*[[:space:]]*<-[[:space:]]*function" R/*.R

# Check NAMESPACE exports
cat NAMESPACE

# Check Roxygen documentation
grep -n "#'" R/*.R

# View existing man pages
ls -lh man/

# Check vignettes
ls -lh vignettes/
```

### Statistical Method Identification

```bash
# Find ashr calls
grep -n "ashr::" R/*.R
grep -n "ash(" R/*.R

# Find mashr calls
grep -n "mashr::" R/*.R
grep -n "mash(" R/*.R

# Find statistical calculations
grep -n "var\|sd\|mean\|median" R/*.R
grep -n "sum\|prod" R/*.R
```

---

## Execution Instructions

### Start Here

**Phase 0: COMPLETED ✅**
- Repository explored
- Publication identified
- File inventory complete
- Vignette reviewed

**Phase 1: Systematic documentation - BEGIN HERE**

**Reading order:**
1. **TRADE.R** (116 lines)
   - Main wrapper function
   - Understand user-facing API
   - Document parameter validation
   - Document routing logic

2. **TRADE_univariate.R** (290 lines)
   - Core univariate analysis
   - Document all 4 functions:
     - TRADE_univariate()
     - fit_ash()
     - get_distribution_output()
     - frac_subset()
   - Verify ashr usage
   - Verify TI and Me calculations

3. **TRADE_bivariate.R** (270 lines)
   - Core bivariate analysis
   - Document TRADE_bivariate()
   - Verify mashr usage
   - Verify correlation estimation

**Phase 2-6: Follow methodology exactly as written**

### Time Estimate

Based on Mixscale results and TRADEtools scope:

**TRADEtools package:**
- File size: 676 lines
- Complexity: 6 functions, statistical methods (ashr/mashr)
- Publication: Nature Genetics paper to cross-reference
- Estimated time: **6-10 hours** for complete analysis

**Breakdown:**
- Phase 0: ✅ DONE (1 hour)
- Phase 1: 3-4 hours (detailed reading + documentation)
- Phase 2: 1-2 hours (cross-references + examples)
- Phase 3: 2-3 hours (verification)
- Phase 4: 0-1 hour (gap filling if needed)
- Phase 5-6: 1 hour (audit report + git)

### Critical Success Factors

1. **Read code line-by-line** - Don't assume, verify every detail
2. **Understand statistical methods** - ashr and mashr are complex
3. **Cross-reference publication** - Nadig et al. is ground truth
4. **Verify vignette examples** - Real usage validates understanding
5. **Check Roxygen docs** - May be outdated or incomplete
6. **Report progress incrementally** - File-by-file updates
7. **Maintain audit trail** - Every verification documented
8. **Zero tolerance for assumptions** - If uncertain, read the code again

---

## Appendix: Repository Context

**Package Name:** TRADEtools

**Publication:**
- **Primary:** Nadig et al. "Transcriptome-wide characterization of genetic perturbations" *Nature Genetics* 57(5):1228-1237 (2025). doi: 10.1038/s41588-025-02169-3
- **Preprint:** bioRxiv 10.1101/2024.07.03.601903 (July 2024)

**Repository Structure:**
```
TRADEtools/
├── DESCRIPTION              - Package metadata
├── NAMESPACE                - Exports: exportPattern("^[[:alpha:]]+")
├── README.md                - Overview and installation
├── R/                       ← **DOCUMENTATION TARGET**
│   ├── TRADE.R              (116 lines) - Main API
│   ├── TRADE_univariate.R   (290 lines) - Univariate engine
│   └── TRADE_bivariate.R    (270 lines) - Bivariate engine
├── man/
│   └── TRADE.rd             - Existing documentation
├── vignettes/
│   ├── TRADEtools-intro.Rmd - Tutorial source
│   └── TRADEtools-intro.html - Compiled tutorial
└── example/                 - Example data
```

**Branch:** claude/perturbseq-library-documentation-01SUE1wkcvpPxXo6Ag1m5gVm

**Expected Outcome:** Complete, verified, accurate documentation matching Mixscale quality standards (100% coverage, 100% accuracy, full audit trail).

---

**END OF METHODOLOGY DOCUMENT**

**Next Action:** Execute Phase 1 - Systematic Code Reading & Initial Documentation

**Begin with:** Read R/TRADE.R and document the main TRADE() function
