# SV Pathogenicity Classifier with Short-Read vs Long-Read Feature Ablation
## AI Development Context Pack

**Project domain:** DNA Sequencing Technology / Structural Variation / Machine Learning  
**Project type:** Classification + controlled feature ablation  
**Primary task:** Predict whether a structural variant (SV) is pathogenic or benign  
**Novelty axis:** Quantify how much long-read-exclusive information improves pathogenicity prediction over information obtainable from short-read data alone  
**Default implementation stack:** Python, pandas, scikit-learn, XGBoost, SHAP, cyvcf2/pysam, bedtools/pybedtools  
**Compute target:** Laptop/CPU-friendly; no wet lab; no GPU required  
**Recommended project duration:** ~2 weeks part-time for a course implementation  
**Status of this document:** Master context file for AI-assisted development

---

# 0. Instructions to Any AI Agent Reading This File

You are helping develop a DNA Sequencing Technology course project titled:

> **SV Pathogenicity Classifier with Short-Read vs Long-Read Feature Ablation**

Treat this file as the authoritative project scope unless the student explicitly asks to redesign the project.

Your job is to help implement, test, document, explain, and present the project.

## 0.1 Core rule

Do **not** reduce this project to:

> "Train a Random Forest/XGBoost model to classify pathogenic vs benign SVs."

That is only the baseline.

The actual project is:

> Train the same pathogenicity classifier under two information regimes and measure the gain from long-read-specific features.

The core controlled comparison is:

- **Arm A:** features realistically obtainable from short-read SV data
- **Arm B:** Arm A + features enabled by long-read SV characterization

The headline result must be a measured difference such as:

- ΔAUROC
- ΔAUPRC
- changes in recall/precision
- improvement by SV type or size band

The project is only complete if this ablation is completed and interpreted.

## 0.2 Novelty statement

Use this statement consistently in the report, presentation, README, and viva:

> **We quantify how much long-read-exclusive features — exact breakpoints, resolved insertion sequence, microhomology length, and phasing — improve structural variant pathogenicity prediction over features obtainable from short-read data alone.**

Do not dilute this by adding unrelated secondary tasks.

## 0.3 Important scientific constraint

Do **not** assume that a sufficiently large dataset of **labelled pathogenic/benign SVs with both short-read and long-read characterizations** automatically exists.

Before building the final Arm A vs Arm B experiment, perform a **feasibility check**:

1. Identify labelled pathogenic/benign SVs.
2. Identify comparable short-read and long-read SV callsets.
3. Match them using genomic coordinates / reciprocal overlap.
4. Determine how many labelled SVs survive.
5. Report counts by:
   - label
   - SV type
   - chromosome
   - size band
6. Only then decide whether the full long-read ablation is statistically viable.

Never fabricate Arm B feature values.

Never set "long-read feature = 0" merely because the feature is unavailable in short reads unless that encoding is scientifically justified.

If the matched labelled cohort is too small, stop and report the problem to the student before changing project scope.

## 0.4 Do not silently redesign

If a dataset, accession, genome build, or feature is unavailable:

- explain the issue,
- provide the closest valid alternative,
- state exactly what changes,
- preserve the novelty claim if possible,
- ask for approval before materially changing the experimental question.

---

# 1. Project in One Sentence

Given a structural variant in the human genome, predict whether it is **pathogenic** or **benign**, and then determine whether adding information obtainable from long-read sequencing improves that prediction compared with short-read-derived information alone.

---

# 2. Simple Conceptual Explanation

A human genome contains large genomic changes called **structural variants (SVs)**.

Examples include:

- deletion (DEL)
- duplication (DUP)
- insertion (INS)
- inversion (INV)
- translocation (TRA)

A structural variant is generally considered a genomic change of approximately **50 bp or larger**.

Many SVs are harmless population variants.

Some are disease-causing because they may:

- delete an important gene,
- duplicate a dosage-sensitive region,
- interrupt a coding exon,
- disrupt regulatory DNA,
- create a fusion,
- alter gene dosage,
- change genomic structure.

The machine-learning task is therefore:

```text
Structural Variant
        ↓
Genomic + sequencing-derived features
        ↓
Machine-learning classifier
        ↓
Pathogenic or Benign
```

The research question is not only whether classification is possible.

The research question is:

```text
How much better can we classify pathogenicity
when long-read-specific SV information is available?
```

---

# 3. Why Short Reads vs Long Reads Matters

## 3.1 Short-read sequencing

Illumina-style short reads are typically around ~150 bp.

Structural variants can be thousands or millions of base pairs long.

Short-read pipelines often infer SVs indirectly using signals such as:

- abnormal read depth,
- discordant paired-end mappings,
- split reads,
- soft clipping.

This can result in:

- uncertain breakpoint locations,
- unresolved insertion sequences,
- ambiguity in repetitive regions,
- weak phasing,
- difficulty resolving complex SVs.

## 3.2 Long-read sequencing

Long-read technologies such as:

- Oxford Nanopore
- PacBio HiFi

produce reads that may span thousands to tens of thousands of bases or more.

A long read may directly span a structural variant and its breakpoints.

This can provide:

- more precise breakpoints,
- full insertion sequences,
- repeat structure,
- phasing,
- breakpoint sequence context,
- microhomology information.

## 3.3 Experimental idea

Instead of merely repeating the statement:

> "Long reads are better for structural variant analysis."

this project tries to **measure the value of that extra information**.

---

# 4. Main Research Question

> **How much do long-read-exclusive structural-variant features improve pathogenicity prediction compared with short-read-derivable features alone?**

Secondary questions:

1. Which long-read features contribute most to the improvement?
2. Does the gain differ by SV type?
3. Does the gain differ by SV size?
4. Does performance remain after removing obvious confounders such as size and allele frequency?
5. Are insertions or repeat-associated SVs helped more by long-read information?

---

# 5. Hypotheses

## H1 — Primary hypothesis

Adding long-read-enabled features will improve pathogenicity classification performance over short-read-derived features alone.

Expected direction:

```text
AUPRC_ArmB > AUPRC_ArmA
AUROC_ArmB > AUROC_ArmA
```

## H2 — Breakpoint hypothesis

Better breakpoint precision will improve the ability to correctly determine whether an SV:

- intersects a gene,
- intersects an exon,
- affects coding sequence,
- interrupts regulatory regions.

## H3 — Insertion hypothesis

Resolved insertion sequence features will be particularly valuable for INS variants because short-read callers often cannot reconstruct long inserted sequences.

## H4 — Mechanism hypothesis

Breakpoint microhomology and sequence-context information may contribute predictive information because they reflect mechanisms of SV formation.

## H5 — Confound-control hypothesis

Performance will decrease when SV size and allele frequency are removed, but biological/contextual features should retain non-trivial predictive ability.

---

# 6. Course / Syllabus Alignment

This project directly addresses:

- structural variation,
- applications of long-read sequencing,
- advantages of long reads for structural-variant analysis,
- genomic data analysis.

This makes it particularly suitable for a DNA Sequencing Technology course.

The project connects strongly to the following biological ideas:

- SV discovery,
- breakpoint resolution,
- long-read sequencing,
- population structural variation,
- disease-causing structural variation.

---

# 7. Base Paper and Literature Foundation

## 7.1 Primary ML base paper

**Sharo AG, Hu Z, Sunyaev SR, Brenner SE. (2022).**  
**"StrVCTVRE: A supervised learning method to predict the pathogenicity of human genome structural variants."**  
*The American Journal of Human Genetics*, 109(2), 195–209.  
DOI: `10.1016/j.ajhg.2021.12.007`

Role in this project:

- establishes that ML can be used for SV pathogenicity prediction,
- provides a strong methodological reference,
- acts as the closest base paper for the classification task.

Important distinction:

StrVCTVRE is the **foundation**, not the novelty.

Our novelty is the controlled short-read vs long-read feature comparison.

Official/public references verified on 2026-08-18:

- PubMed: `https://pubmed.ncbi.nlm.nih.gov/35032432/`
- Project site: `https://strvctvre.berkeley.edu/`

## 7.2 Structural-variation population reference

**Sudmant PH et al. (2015).**  
**"An integrated map of structural variation in 2,504 human genomes."**  
*Nature*, 526, 75–81.

Role:

- population SV reference,
- supports benign/common SV framing,
- directly aligned with the course references.

Useful source:

`https://www.internationalgenome.org/data-portal/data-collection/phase3/`

## 7.3 Disease framing

**Stankiewicz P, Lupski JR. (2010).**  
**"Structural variation in the human genome and its role in disease."**

Role:

- biological explanation of how SVs cause disease,
- useful for introduction and discussion,
- useful when explaining breakpoint mechanisms and disease relevance.

## 7.4 Benchmark / long-read evidence reference

Genome in a Bottle (GIAB) structural-variant benchmark work provides high-confidence SV resources derived from multiple technologies.

Important benchmark paper:

**Zook et al.**
"A robust benchmark for detection of germline large deletions and insertions."

GIAB/NIST resource:

`https://www.nist.gov/programs-projects/genome-bottle`

Note:

GIAB resources evolve. As of 2026, newer HG002 benchmark products exist. Do not hard-code an old filename without checking the official GIAB release information.

---

# 8. Classification Task

## 8.1 Unit of classification

One structural variant.

## 8.2 Input

A numerical/categorical feature vector describing:

- the SV itself,
- genomic context,
- population frequency,
- repeat context,
- sequencing-resolution information.

## 8.3 Output

Binary label:

```text
0 = benign
1 = pathogenic
```

## 8.4 Example row

```text
chrom                  = chr7
svtype                 = DEL
length                 = 120000
population_af          = 0.00002
n_genes_overlapped     = 3
n_exons_overlapped     = 7
overlaps_segdup        = False
breakpoint_uncertainty = 4
is_phased              = True
label                  = 1
```

---

# 9. Label Definition

Label construction is one of the most attackable parts of the project.

The AI agent must document the label policy exactly.

## 9.1 Benign class

Initial proposed definition:

> Common structural variants from population datasets such as the 1000 Genomes Project, using an allele-frequency threshold such as AF > 0.01.

Interpretation:

A variant present in more than ~1% of a population is unlikely to be a highly penetrant severe pathogenic variant.

Label:

```text
benign = 0
```

Potential source:

- 1000 Genomes Phase 3 integrated SV map
- optionally gnomAD-SV

## 9.2 Pathogenic class

Initial proposed definition:

> Structural variants / CNVs with clinical significance "Pathogenic" or "Likely pathogenic" from ClinVar/dbVar-derived resources.

Label:

```text
pathogenic = 1
```

## 9.3 Exclude ambiguous clinical classes

Do not include:

- VUS / uncertain significance,
- conflicting interpretation,
- likely benign,
- benign in the pathogenic set,
- records without clear significance.

## 9.4 Important caveat

"Common population variant = benign" is a pragmatic project label, not perfect clinical truth.

State this limitation explicitly.

---

# 10. Primary Data Sources

Dataset URLs and accessions can move. Verify before downloading.

## 10.1 1000 Genomes Phase 3 SV map

Purpose:

- common/benign population SVs,
- allele frequencies,
- SV coordinates/types.

Historical filename from the project catalogue:

`ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz`

Associated reference:

Sudmant et al. 2015.

Portal:

`https://www.internationalgenome.org/data-portal/data-collection/phase3/`

Typical build:

GRCh37/hg19 for the classic Phase 3 integrated SV map.

## 10.2 ClinVar / dbVar

Purpose:

- pathogenic / likely pathogenic CNVs and SVs.

Portal:

`https://www.ncbi.nlm.nih.gov/dbvar/`

The original project catalogue mentioned the ClinVar-related dbVar study commonly known as `nstd102`, but:

> **Do not assume this accession is still the best/current source. Verify before use.**

The AI agent must inspect:

- assembly/build,
- clinical significance,
- SV type,
- coordinate quality,
- accession version.

## 10.3 gnomAD-SV

Purpose:

- optional additional population SVs,
- population frequencies,
- additional benign/reference data.

Portal:

`https://gnomad.broadinstitute.org/downloads`

Use only if necessary.

Avoid adding data merely for size if it increases integration complexity.

## 10.4 GENCODE

Purpose:

- gene coordinates,
- transcript coordinates,
- exon coordinates,
- CDS coordinates,
- TSS-related annotations.

Official:

`https://www.gencodegenes.org/human/`

As of 2026-08-18:

- current GRCh38 human release: GENCODE Release 50
- GRCh37 mapping remains available via the GRCh37 release line

If the project is standardized to GRCh37, use the matching GENCODE GRCh37 annotation rather than mixing GRCh38 coordinates.

## 10.5 UCSC tracks

Purpose:

- RepeatMasker,
- segmental duplications,
- optionally mappability.

Portal:

`https://genome.ucsc.edu/cgi-bin/hgTables`

Export as BED when convenient.

## 10.6 GIAB / HG002

Purpose:

- long-read / multi-technology SV benchmark resources,
- possible source for the long-read arm,
- technology comparison.

Official GIAB/NIST:

`https://www.nist.gov/programs-projects/genome-bottle`

Important:

Do not automatically use "latest" benchmark without checking reference build and compatibility with the pathogenicity dataset.

## 10.7 HGSVC / long-read resources

Possible additional source of long-read SV calls.

Use only after verifying:

- sample identity,
- genome build,
- technology,
- SV call format,
- ability to match to the labelled cohort.

---

# 11. Genome Build Policy

This is a critical requirement.

Never intersect records from different genome builds directly.

Common builds:

- GRCh37 / hg19
- GRCh38 / hg38
- T2T-CHM13

The classic 1000 Genomes Phase 3 integrated SV map is GRCh37-oriented.

Therefore the simplest first implementation may standardize all data to:

> **GRCh37/hg19**

However, the final build choice must be made only after verifying all required sources.

## 11.1 Build decision checklist

Before feature engineering, record:

```text
Benign dataset build:
Pathogenic dataset build:
GENCODE build:
RepeatMasker build:
SegDup build:
Long-read callset build:
Reference FASTA build:
```

Every one must match.

If not, convert using:

- UCSC liftOver,
- CrossMap,
- another validated coordinate conversion workflow.

After liftOver:

- record failed mappings,
- record multi-mappings,
- drop ambiguous conversion when necessary,
- report counts before/after conversion.

---

# 12. Coordinate System Policy

Be careful:

- BED: 0-based, half-open
- VCF: 1-based positions
- GTF/GFF: generally 1-based inclusive

Do not manually mix these conventions.

Prefer:

- bedtools
- pybedtools
- cyvcf2
- pysam

Write unit tests for coordinate conversion.

---

# 13. Feature Engineering

Features are divided into two experimental arms.

---

# 14. Arm A — Short-Read-Derivable Features

These features should be available or derivable without depending on long-read-exclusive sequence resolution.

## 14.1 `svtype`

Type:

categorical

Possible values:

- DEL
- DUP
- INV
- INS
- CNV
- possibly TRA if supported consistently

Encoding:

one-hot encode for linear/tree models if necessary.

## 14.2 `log_length`

Definition:

```python
log_length = log10(abs(SVLEN))
```

Reason:

SV sizes span orders of magnitude.

## 14.3 `chrom`

Chromosome.

Potential concern:

chromosome may introduce unwanted positional shortcuts.

Evaluate performance with and without it.

## 14.4 `gc_content_flank`

GC percentage in a fixed window around the SV breakpoints.

Suggested:

±1 kb.

Requires:

reference FASTA.

## 14.5 `n_genes_overlapped`

Number of GENCODE genes intersected by the SV.

Tool:

`bedtools intersect`

## 14.6 `n_exons_overlapped`

Number of exons intersected.

Important because exon disruption is biologically relevant.

## 14.7 `frac_cds_overlapped`

Possible definition:

```text
overlap length with CDS
-----------------------
relevant total CDS length
```

Document the exact denominator used.

## 14.8 `dist_to_nearest_tss`

Distance to nearest transcription start site.

Tool:

`bedtools closest`

## 14.9 `population_af`

Population allele frequency.

Potential sources:

- 1000G
- gnomAD-SV

This is highly predictive but also a major potential confound.

## 14.10 `overlaps_segdup`

Boolean indicating overlap with segmental duplication.

## 14.11 `repeat_class_frac`

Fraction of the SV interval overlapping repeat classes.

Possible categories:

- Alu
- LINE
- SINE
- LTR
- simple repeats
- tandem repeats
- other

Use separate numeric columns when possible.

Example:

```text
frac_alu
frac_line
frac_sine
frac_ltr
frac_simple_repeat
```

---

# 15. Arm B — Long-Read-Enabled Features

Arm B must include all Arm A features plus long-read-specific or long-read-improved information.

Do not invent values where data are unavailable.

## 15.1 `breakpoint_precision`

Concept:

How uncertain are the start/end breakpoints?

Possible short-read representation:

- `CIPOS`
- `CIEND`

Possible measure:

```text
uncertainty_width_start + uncertainty_width_end
```

Long-read calls may have near-exact breakpoints.

Important:

Different callers encode breakpoint uncertainty differently.

Normalize carefully.

## 15.2 `has_resolved_insertion_seq`

Boolean.

```text
1 = actual inserted sequence available
0 = unresolved symbolic insertion
```

Only meaningful for insertion-like SVs.

## 15.3 `insertion_seq_length`

Length of actual resolved inserted sequence.

Use only where actual inserted bases are available.

## 15.4 `insertion_seq_gc`

GC content of resolved inserted sequence.

## 15.5 `insertion_repeat_class`

Repeat classification of inserted sequence.

## 15.6 `tandem_repeat_unit_length`

Repeat motif/unit size in repeat-associated SVs.

## 15.7 `copy_number_change`

Estimated repeat copy-number change if reliably available.

## 15.8 `is_phased`

Whether the SV is phased to a haplotype.

## 15.9 `in_phase_block`

Information describing membership in a phased block if available.

## 15.10 `breakpoint_microhomology_length`

Length of matching short sequence at the rearrangement junction.

Biological interpretation:

May provide information about the mechanism that created the SV.

This requires highly accurate breakpoint sequence context.

---

# 16. Feature Availability Matrix

Before modelling, create a table like:

| Feature | Arm A | Arm B | Source | Missingness | Notes |
|---|---:|---:|---|---:|---|
| svtype | yes | yes | VCF | ... | |
| log_length | yes | yes | VCF | ... | |
| population_af | yes | yes | 1000G/gnomAD | ... | |
| n_genes_overlapped | yes | yes | GENCODE | ... | |
| breakpoint_precision | limited | yes | VCF/callset | ... | |
| resolved insertion seq | no/limited | yes | long-read callset | ... | |
| microhomology | no/limited | yes | breakpoint sequence | ... | |
| phasing | limited | yes | long-read callset | ... | |

This table must be included in the final report or appendix.

---

# 17. Critical Feasibility Gate for the Ablation

This is the first major engineering milestone.

## 17.1 Goal

Determine whether enough **labelled SVs** can be represented under both feature regimes.

## 17.2 Matching procedure

Where comparable short-read and long-read calls exist for the same SV/sample:

Use reciprocal overlap.

Initial matching rule:

```bash
bedtools intersect -f 0.5 -r -wa -wb
```

Meaning:

at least 50% reciprocal overlap.

## 17.3 Sensitivity analysis

Try:

- 0.3
- 0.5
- 0.7

to determine whether conclusions are sensitive to matching threshold.

## 17.4 Matching is not enough

Two calls overlapping does not automatically mean they are the same event.

Also compare:

- SVTYPE,
- breakpoint distance,
- SV length similarity,
- inserted sequence where relevant.

## 17.5 Required feasibility report

Produce:

```text
Total labelled pathogenic SVs:
Total labelled benign SVs:

Matched to long-read characterization:
  pathogenic:
  benign:

By SV type:
  DEL:
  DUP:
  INS:
  INV:

By size:
  50 bp–1 kb:
  1–10 kb:
  10–100 kb:
  >100 kb:
```

If sample size is too small for reliable modelling, notify the student.

---

# 18. Dataset Assembly Pipeline

Recommended conceptual pipeline:

```text
1000 Genomes / population SVs
        ↓
Benign candidate table
        ↓
AF filtering
        ↓
label = 0
        │
        │
        ├──────────────┐
        │              │
        │              ↓
        │       unified SV table
        │              ↑
        │              │
ClinVar/dbVar          │
        ↓              │
Pathogenic SVs         │
        ↓              │
clinical filtering     │
        ↓              │
label = 1 ─────────────┘
                       ↓
              normalize genome build
                       ↓
               normalize SV fields
                       ↓
               annotation features
                       ↓
          match long-read characterization
                       ↓
             Arm A / Arm B feature tables
```

---

# 19. Unified Internal SV Schema

Every data source should be converted to a common internal schema before modelling.

Suggested columns:

```text
sv_id
source
chrom
start
end
svtype
svlen
abs_svlen
label
clinical_significance
population_af
genome_build
original_id
original_source
```

Add annotation columns later.

Never train directly from source-specific raw columns.

---

# 20. Source Leakage Risk

This is one of the most important pitfalls.

If:

```text
all benign SVs come from 1000G
all pathogenic SVs come from ClinVar
```

the model may learn database-specific artifacts instead of pathogenicity.

Examples:

- coordinate rounding,
- variant-size conventions,
- type distributions,
- caller-specific fields,
- missingness patterns,
- chromosome distributions.

## 20.1 Leakage test

Train an auxiliary classifier:

```text
features → source_database
```

If source database can be predicted nearly perfectly from supposedly biological features, investigate leakage.

Also ensure that source identifiers or source-specific INFO fields are never included as model features.

---

# 21. Major Confound: Size and Allele Frequency

Pathogenic SVs may systematically be:

- larger,
- rarer.

A classifier could learn:

```text
large + rare = pathogenic
```

and achieve impressive metrics without learning deeper biology.

Therefore run a mandatory confound ablation.

## Model C1 — Full baseline

Includes:

- length
- population AF
- context features

## Model C2 — Confound-reduced

Remove:

- `log_length`
- `population_af`

Evaluate again.

Report:

```text
Full model performance
vs
No-size/no-AF performance
```

Also stratify by size.

---

# 22. Class Imbalance

Benign variants may greatly outnumber pathogenic variants.

Do not rely on accuracy.

Possible strategy:

- subsample benign class to approximately 1:5 or 1:10,
- retain enough diversity,
- use `class_weight="balanced"` where appropriate,
- use XGBoost `scale_pos_weight`.

Do not oversample before splitting.

If SMOTE is ever used:

> It must occur inside the training fold only.

However, prefer class weighting over SMOTE for the initial implementation.

---

# 23. Train/Test Splitting

Do not use a naive random row split if nearby variants or related records can leak.

Preferred:

> **chromosome-held-out validation**

Example conceptual split:

```text
training:
chr1–chr16

validation/test:
chr17–chr22
```

A better implementation may use grouped cross-validation where chromosome is the group.

Keep sex chromosomes separate or exclude them initially unless explicitly justified.

Important:

Use the **same split** for Arm A and Arm B.

---

# 24. Models

Run all three.

## 24.1 Logistic Regression

Purpose:

- interpretable baseline,
- sanity check.

Needs:

- numeric scaling,
- one-hot encoding.

## 24.2 Random Forest

Purpose:

- nonlinear baseline,
- handles mixed feature relationships,
- useful feature importance.

## 24.3 XGBoost

Purpose:

- likely strongest tabular model,
- handles nonlinearity,
- compatible with SHAP.

Do not introduce deep learning unless the project is explicitly redesigned.

---

# 25. Main Ablation Protocol

This is the core experiment.

## Step 1

Construct a cohort where the relevant SVs have valid labels and sufficient features.

## Step 2

Create:

```text
X_A = Arm A features only
X_B = Arm A + Arm B features
```

## Step 3

Use:

- same samples,
- same labels,
- same train/test folds,
- same random seed,
- same model family,
- same hyperparameter policy.

## Step 4

Train Arm A.

## Step 5

Train Arm B.

## Step 6

Compare:

```text
ΔAUROC = AUROC_B - AUROC_A
ΔAUPRC = AUPRC_B - AUPRC_A
ΔF1    = F1_B - F1_A
```

## Step 7

Use repeated grouped validation if feasible.

Original project suggestion:

- repeated 5-fold style analysis,
- multiple repetitions,
- paired comparison of folds.

Because chromosomes are groups, adapt repetition carefully to avoid violating group separation.

## Step 8

Estimate uncertainty:

- mean,
- standard deviation,
- bootstrap confidence interval where appropriate.

## Step 9

Perform paired statistical comparison across matched folds/repetitions where defensible.

---

# 26. Metrics

Mandatory:

- AUROC
- AUPRC / Average Precision
- Precision
- Recall
- F1
- Confusion matrix

Recommended:

- calibration curve
- Brier score
- per-size-band performance
- per-SV-type performance

Primary metric:

> **AUPRC**

Reason:

The dataset is imbalanced.

---

# 27. Size-Stratified Analysis

Report metrics for:

```text
50 bp – 1 kb
1 kb – 10 kb
10 kb – 100 kb
>100 kb
```

This is important because pathogenic and benign sets may have very different size distributions.

---

# 28. SV-Type Analysis

Where sample sizes permit, evaluate separately for:

- DEL
- DUP
- INS
- INV

Do not report unstable per-class metrics for categories with extremely small sample counts.

For rare classes, state that sample size is insufficient.

---

# 29. SHAP / Interpretability

For the best-performing Arm B model:

Generate SHAP feature importance.

Questions:

1. Which long-read feature contributes most?
2. Is performance mostly driven by size/AF?
3. Are breakpoint features important?
4. Are insertion-specific features important only for INS?
5. Does phasing contribute?
6. Does repeat context matter?

Recommended figure:

SHAP beeswarm.

Also create a bar plot grouped by:

```text
basic SV properties
genomic context
population features
repeat context
long-read-only features
```

---

# 30. Headline Figure

The core presentation figure should compare Arm A and Arm B.

Preferred:

Two ROC curves and/or two PR curves.

Since AUPRC is especially important, strongly consider making the **precision-recall comparison the primary figure** and ROC secondary.

Example caption:

> Pathogenicity classification using short-read-derivable features compared with the same classifier augmented with long-read-enabled features.

Display:

```text
Arm A AUPRC = ...
Arm B AUPRC = ...
ΔAUPRC = ...
```

Do not manufacture the expected direction before running the experiment.

If Arm B does not improve performance, that is still a valid result.

---

# 31. Secondary Figures

Useful:

1. Class distribution
2. SV size distribution by label
3. SV type distribution by label
4. AUPRC Arm A vs Arm B
5. SHAP beeswarm
6. Performance by size band
7. Performance by SV type
8. Calibration curve
9. Source-leakage diagnostic
10. Missingness heatmap

---

# 32. Expected Final Feature Table

Preferred output:

`data/processed/features.parquet`

Possible columns:

```text
sv_id
chrom
start
end
svtype_DEL
svtype_DUP
svtype_INS
svtype_INV
log_length
gc_content_flank
n_genes_overlapped
n_exons_overlapped
frac_cds_overlapped
dist_to_nearest_tss
population_af
overlaps_segdup
frac_alu
frac_line
frac_sine
frac_ltr
frac_simple_repeat
breakpoint_precision
has_resolved_insertion_seq
insertion_seq_length
insertion_seq_gc
tandem_repeat_unit_length
copy_number_change
is_phased
breakpoint_microhomology_length
label
source
```

Do **not** train on:

```text
source
original_source
clinical_significance text
database accession patterns
```

unless explicitly part of a diagnostic analysis.

---

# 33. Project Repository Structure

Use:

```text
sv-pathogenicity-longread/
│
├── README.md
├── requirements.txt
├── environment.yml
├── context.md
│
├── data/
│   ├── raw/
│   ├── interim/
│   └── processed/
│
├── notebooks/
│   ├── 01_data_audit.ipynb
│   ├── 02_feature_engineering.ipynb
│   ├── 03_feasibility_ablation.ipynb
│   ├── 04_modelling.ipynb
│   └── 05_evaluation.ipynb
│
├── src/
│   ├── config.py
│   ├── download.py
│   ├── parse_sv.py
│   ├── harmonize.py
│   ├── annotate.py
│   ├── features.py
│   ├── match_callsets.py
│   ├── models.py
│   ├── evaluate.py
│   └── utils.py
│
├── tests/
│   ├── test_coordinates.py
│   ├── test_sv_length.py
│   ├── test_overlap.py
│   └── test_features.py
│
├── configs/
│   ├── data_sources.yaml
│   └── experiment.yaml
│
├── figures/
│
├── models/
│
└── report/
    ├── report.md
    ├── references.bib
    └── slides/
```

---

# 34. Configuration Philosophy

Avoid hard-coded URLs and thresholds in Python files.

Example `configs/data_sources.yaml`:

```yaml
genome_build: GRCh37

benign:
  source: 1000_genomes
  url: null
  af_threshold: 0.01

pathogenic:
  source: clinvar_dbvar
  url: null
  allowed_significance:
    - Pathogenic
    - Likely pathogenic

annotations:
  gencode_gtf: null
  repeatmasker_bed: null
  segdup_bed: null

long_read:
  source: null
  build: null
```

The AI agent should populate verified URLs only after checking them.

---

# 35. Python Environment

Recommended:

```text
python >= 3.10
pandas
numpy
scikit-learn
xgboost
shap
matplotlib
biopython
cyvcf2
pysam
pybedtools
pyfaidx
imbalanced-learn
jupyter
pyarrow
```

Command-line:

```text
bedtools
samtools
bcftools
CrossMap or liftOver
```

On Windows:

Use WSL2/Ubuntu for:

- bedtools
- samtools
- bcftools
- bioconda tools

Pure Python work can be done natively.

---

# 36. Reproducibility Rules

Always:

```python
random_state = 42
```

Record:

- package versions,
- dataset versions,
- access date,
- genome build,
- file checksums if practical,
- filtering counts,
- train/test chromosomes,
- feature names,
- hyperparameters.

Commit:

`features.parquet`

if licensing permits.

Do not commit huge raw genomic files.

---

# 37. Data Audit Requirements

Before modelling, generate an audit report.

Include:

```text
number of benign variants
number of pathogenic variants
number per chromosome
number per SV type
SV length distribution
AF distribution
missing values
duplicate coordinates
genome build
source distribution
matched long-read counts
```

Plot:

- size histogram/log distribution,
- label counts,
- SVTYPE by label,
- AF by label.

---

# 38. Duplicate / Overlapping Record Handling

Pathogenic databases may contain multiple records describing the same or similar event.

Define deduplication rules.

Possible strategy:

- same build,
- same chromosome,
- same SV type,
- highly overlapping interval,
- same clinical label.

Do not allow near-duplicate variants to appear in both train and test.

Chromosome-group splitting helps but does not solve duplicates within a chromosome.

Deduplicate before splitting.

---

# 39. Missing Data Policy

For each feature record:

- missing percentage,
- biological reason,
- source reason.

Do not blindly mean-impute biologically meaningful missingness.

Examples:

`insertion_seq_length` missing for a deletion is **not the same thing** as unknown insertion length.

Consider:

- type-specific features,
- NA + missingness indicator,
- separate per-SV-type analysis.

Tree models can sometimes handle missing values, but preprocessing must be consistent.

---

# 40. Hyperparameter Tuning

Keep tuning modest.

The scientific comparison is more important than squeezing out 0.5% performance.

For a fair Arm A vs Arm B comparison:

- same tuning strategy,
- tuning inside training folds only.

Never tune using the final test set.

Possible XGBoost grid:

```text
n_estimators: [200, 500]
max_depth: [3, 5, 7]
learning_rate: [0.03, 0.1]
subsample: [0.8, 1.0]
colsample_bytree: [0.8, 1.0]
```

Use a smaller grid if runtime is high.

---

# 41. Statistical Comparison

If Arm A and Arm B are evaluated on the same folds:

Store per-fold metrics.

Then calculate:

```text
metric_B - metric_A
```

for every paired fold.

Report:

- mean delta,
- std,
- confidence interval.

A paired significance test may be used if assumptions are justified.

Do not overstate statistical significance from a tiny number of folds.

---

# 42. Failure Modes and Pitfalls

## 42.1 Genome build mismatch

Symptom:

Gene/repeat overlaps become nonsensical.

Fix:

Resolve build before feature engineering.

## 42.2 Source leakage

Symptom:

Near-perfect performance with suspiciously simple features.

Fix:

Source classifier diagnostic + remove artifacts.

## 42.3 Size leakage/confounding

Symptom:

Length dominates SHAP.

Fix:

No-length/no-AF ablation + size-stratified metrics.

## 42.4 Coordinate bugs

Symptom:

Off-by-one or zero overlaps.

Fix:

Unit tests + bedtools.

## 42.5 Too few matched long-read SVs

Symptom:

Arm B contains very few labelled examples.

Fix allowed within current scope:

- inspect matching quality,
- try reciprocal overlap sensitivity at 0.3/0.5/0.7,
- analyze per SV type if sample counts allow.

Do not silently replace the experiment.

## 42.6 Missing insertion features

Do not fill unavailable insertion sequence with a fake sequence or 0.

Use correct missingness semantics.

## 42.7 Imbalanced evaluation

Do not report accuracy as the main metric.

Use AUPRC.

## 42.8 Feature leakage

Never include:

- clinical significance text,
- source database,
- disease name,
- submission IDs,
- labels encoded indirectly in column names.

## 42.9 Random train/test split

Avoid because genomic regions may be correlated.

Prefer grouped chromosome splitting.

## 42.10 Overclaiming clinical utility

This is a research/course model.

Do not claim:

- clinical diagnostic readiness,
- replacement of clinical geneticists,
- direct patient use.

---

# 43. Development Roadmap

## Phase 0 — Feasibility

Deliverable:

`data_feasibility.md`

Tasks:

1. verify source URLs,
2. verify genome builds,
3. download small subsets,
4. parse labels,
5. count classes,
6. inspect long-read matching feasibility.

Exit criterion:

Enough high-quality labelled examples exist for baseline and a credible Arm B experiment.

## Phase 1 — Baseline dataset

Create:

```text
benign_sv.parquet
pathogenic_sv.parquet
sv_master.parquet
```

## Phase 2 — Harmonization

Tasks:

- normalize chromosome names,
- normalize SVTYPE,
- normalize coordinates,
- convert builds,
- remove invalid/duplicate records.

## Phase 3 — Arm A feature engineering

Create short-read-derived feature table.

## Phase 4 — Baseline modelling

Train:

- Logistic Regression
- Random Forest
- XGBoost

Produce baseline metrics.

## Phase 5 — Long-read matching and Arm B

Match long-read characterization.

Build Arm B features.

## Phase 6 — Controlled ablation

Train same models / same folds.

Calculate Δ metrics.

## Phase 7 — Controls

Run:

- no length + no AF,
- size stratification,
- source leakage diagnostic,
- match-threshold sensitivity.

## Phase 8 — Interpretability

SHAP.

## Phase 9 — Report / presentation

Lead with:

> measured long-read gain

not generic classifier accuracy.

---

# 44. Suggested Two-Week Timeline

## Days 1–2

- verify data sources,
- choose genome build,
- create environment,
- inspect VCF schemas,
- perform long-read feasibility audit.

## Days 3–5

- download selected data,
- parse pathogenic/benign labels,
- harmonize coordinates,
- produce master SV table.

## Days 5–8

- GENCODE overlaps,
- repeat annotations,
- allele frequency,
- Arm A features.

## Days 7–9

- baseline models,
- grouped validation,
- first metrics.

## Days 9–11

- long-read matching,
- Arm B,
- main ablation.

## Days 11–12

- SHAP,
- confound tests,
- figures.

## Days 12–13

- report,
- README,
- slides.

## Day 14

- reproduce pipeline from clean environment,
- viva preparation.

---

# 45. Team Split for Four People

## Person 1 — Data engineering

Own:

- downloads,
- parsing,
- data provenance,
- genome-build harmonization.

## Person 2 — Feature engineering

Own:

- GENCODE,
- RepeatMasker,
- GC,
- overlap calculations,
- Arm B feature extraction.

## Person 3 — Modelling

Own:

- train/test strategy,
- models,
- hyperparameter tuning,
- metrics.

## Person 4 — Analysis / documentation

Own:

- SHAP,
- statistical comparison,
- plots,
- report,
- slides,
- reproducibility.

The main ablation is jointly owned by everyone.

---

# 46. README Structure

Recommended:

```text
# Project Title

## Problem
## Research Question
## Novelty
## Base Paper
## Dataset Sources
## Genome Build
## Label Definition
## Features
## Arm A vs Arm B
## Installation
## Data Download
## Running the Pipeline
## Evaluation
## Reproducing Results
## Limitations
## References
```

---

# 47. Report Structure

## 1. Introduction

Explain:

- structural variants,
- disease relevance,
- short-read limitations,
- long-read advantage.

End with research question.

## 2. Related Work

Discuss:

- StrVCTVRE,
- Sudmant 2015,
- structural variation and disease,
- GIAB benchmark literature.

## 3. Data

Include:

- source,
- version,
- date,
- genome build,
- filtering,
- label definition,
- final class counts.

## 4. Methods

Include:

- feature definitions,
- Arm A,
- Arm B,
- models,
- grouped splitting,
- ablation protocol.

## 5. Results

Order:

1. data characteristics,
2. baseline,
3. main Arm A vs Arm B result,
4. confound controls,
5. stratified analysis,
6. SHAP.

## 6. Discussion

Interpret:

- whether long-read features help,
- which features matter,
- where gains occur,
- limitations.

## 7. Conclusion

Lead with the actual measured result.

---

# 48. Presentation Structure

Suggested 12–14 slides:

1. Title
2. What is a structural variant?
3. Why pathogenic vs benign matters
4. Short reads vs long reads
5. Existing work: StrVCTVRE
6. Research gap + novelty
7. Data sources
8. Feature design: Arm A vs Arm B
9. Pipeline + validation
10. Baseline results
11. **Main ablation result**
12. SHAP / interpretation
13. Limitations
14. Conclusion

Spend most presentation time on Slide 11.

---

# 49. Viva Questions to Prepare

## Q: What is a structural variant?

A genomic change generally ≥50 bp, including deletions, duplications, insertions, inversions, and translocations.

## Q: Why can an SV be pathogenic?

It may disrupt genes/exons, alter dosage, affect regulatory regions, or change chromosome structure.

## Q: What is the project's novelty?

We do not simply classify pathogenicity. We measure the added predictive value of long-read-enabled SV information relative to short-read-derived features.

## Q: Why long reads?

Long reads can span complex SVs and repetitive regions and provide more precise breakpoints, inserted sequence, and phasing.

## Q: What is the base paper?

StrVCTVRE (Sharo et al., 2022), a supervised ML method for predicting SV pathogenicity.

## Q: How is this different from StrVCTVRE?

StrVCTVRE establishes pathogenicity classification. This project performs a controlled short-read vs long-read feature ablation to quantify the value of sequencing-resolution information.

## Q: Why not accuracy?

Class imbalance can make accuracy misleading. AUPRC is more informative for a rare positive class.

## Q: What is the biggest confound?

SV size and population allele frequency.

## Q: How do you control it?

Run models with those features removed and stratify evaluation by SV size.

## Q: Why not random split?

Nearby genomic events/features can be correlated. Chromosome-held-out evaluation better tests genomic generalization.

## Q: What is CIPOS/CIEND?

VCF fields that can represent confidence intervals around structural-variant breakpoints.

## Q: What is breakpoint microhomology?

A short matching sequence around a rearrangement junction that can reflect the DNA-repair mechanism responsible for the SV.

## Q: Is the model clinically deployable?

No. It is a research/course study designed to quantify information value, not a validated clinical diagnostic.

---

# 50. Definition of Done

The project is **not done** merely when a classifier runs.

It is done when all of the following exist:

- [ ] verified public data sources
- [ ] documented genome build
- [ ] reproducible pathogenic/benign label construction
- [ ] cleaned master SV table
- [ ] Arm A feature table
- [ ] validated Arm B feature table
- [ ] long-read feasibility counts
- [ ] Logistic Regression baseline
- [ ] Random Forest baseline
- [ ] XGBoost baseline
- [ ] chromosome/group-aware evaluation
- [ ] AUPRC/AUROC/F1 metrics
- [ ] Arm A vs Arm B paired comparison
- [ ] ΔAUPRC and ΔAUROC
- [ ] size/AF confound ablation
- [ ] size-stratified evaluation
- [ ] source-leakage diagnostic
- [ ] SHAP interpretation
- [ ] headline figure
- [ ] README
- [ ] report
- [ ] presentation
- [ ] limitations clearly stated
- [ ] results reproducible from saved feature table

---

# 51. Recommended First Prompt to an AI Coding Agent

Use:

> Read `context.md` completely before writing code. We are implementing the "SV Pathogenicity Classifier with Short-Read vs Long-Read Feature Ablation." Do not redesign the project. First perform Phase 0 only: verify each proposed data source, identify the current downloadable files/accessions, record genome builds and approximate file sizes, and determine whether a sufficiently large labelled cohort can be matched to long-read SV characterization. Produce a feasibility report before downloading large files or writing the modelling pipeline. Do not guess URLs or accessions.

---

# 52. Prompt for Data Engineering Agent

> Using `context.md` as the authoritative project specification, implement `src/download.py` and `src/parse_sv.py`. Verify URLs before adding them. Store raw files without modification under `data/raw/`. Convert each source into a documented common SV schema. Record source, version, genome build, original ID, filtering counts, and file checksum where practical. Do not mix genome builds.

---

# 53. Prompt for Feature Engineering Agent

> Read `context.md`. Implement Arm A features first in `src/features.py`. Every feature must have a docstring containing its biological definition, input source, coordinate convention, and expected missing-value behavior. Add unit tests for interval overlaps. Do not implement long-read feature placeholders with fake zeros. Only implement an Arm B feature if the required source field or sequence is actually available.

---

# 54. Prompt for Modelling Agent

> Read `context.md`. Use only the processed feature table. Implement Logistic Regression, Random Forest, and XGBoost. Use chromosome/group-aware splitting, class weighting, and identical folds for Arm A and Arm B. Report AUROC, AUPRC, precision, recall, F1, confusion matrices, and per-fold metrics. The primary experimental output is ΔAUPRC between Arm B and Arm A. No test-set tuning.

---

# 55. Prompt for Evaluation Agent

> Read `context.md`. Generate the main Arm A vs Arm B comparison, confidence intervals, size-stratified metrics, SV-type metrics where sample size permits, the no-size/no-AF confound ablation, source-leakage diagnostic, calibration analysis, and SHAP plots. Explicitly distinguish genuine biological findings from dataset artifacts.

---

# 56. Prompt for Report / Viva Agent

> Read `context.md` and the generated results. Explain the project in simple language first, then technical language. The presentation must lead with the measured value of long-read information, not generic classifier accuracy. Prepare answers to dataset-label validity, source leakage, genome-build harmonization, class imbalance, size confounding, chromosome-held-out validation, and clinical limitations.

---

# 57. Non-Negotiable Rules for Future AI Assistance

1. Do not fabricate genomic data.
2. Do not fabricate metrics.
3. Do not invent accessions.
4. Verify links before depending on them.
5. Never mix genome builds.
6. Never use source database as a predictive feature.
7. Never use clinical significance text as a feature.
8. Never run feature selection on the full dataset before splitting.
9. Never evaluate only with accuracy.
10. Never claim long-read benefit until the measured Arm A/Arm B result exists.
11. Never create fake long-read features for samples without long-read characterization.
12. Preserve identical cohorts/folds for the controlled ablation whenever possible.
13. Record every filtering step and sample count.
14. State limitations honestly.
15. Prefer one rigorous ablation over many weak extensions.

---

# 58. Optional Extensions — Only After the Core Project Works

## Extension A — Long-read-only detection

Estimate the fraction of SVs found by long reads that have no short-read counterpart.

If clinical labels exist, ask whether pathogenic variants are enriched in this set.

## Extension B — Per-SV-type gain

Measure:

```text
ΔAUPRC_DEL
ΔAUPRC_DUP
ΔAUPRC_INS
ΔAUPRC_INV
```

## Extension C — Precision degradation simulation

Artificially perturb breakpoint positions and measure pathogenicity-performance degradation.

This could isolate the value of breakpoint precision without requiring every long-read-specific feature.

Only do this as an explicitly labelled simulation, never present it as real long-read data.

## Extension D — External validation

Validate the final model on an independent labelled SV dataset if one is available and compatible.

---

# 59. Main Scientific Limitations

State these openly.

1. Population frequency is an imperfect benign label.
2. ClinVar/dbVar pathogenic annotations may include ascertainment bias.
3. Pathogenic and benign sources may differ systematically.
4. SV types and sizes are imbalanced.
5. Long-read feature availability may be incomplete.
6. Matching SVs across callers/technologies is imperfect.
7. Different variant callers may encode uncertainty differently.
8. Genome-build conversion may lose variants.
9. Some genomic-context annotations are correlated.
10. Results may not generalize to all SV classes or clinical populations.
11. A machine-learning prediction is not a clinical diagnosis.

---

# 60. Final Mental Model

The entire project can be remembered as:

```text
Known benign + pathogenic SVs
            ↓
      harmonize / annotate
            ↓
     same labelled cohort
            ↓
      ┌───────────────┐
      │               │
      ↓               ↓
   Arm A            Arm B
 short-read      short-read
  features          +
                long-read
                 features
      │               │
      ↓               ↓
   same ML         same ML
      │               │
      └───────┬───────┘
              ↓
          compare
              ↓
   ΔAUPRC / ΔAUROC
              ↓
            SHAP
              ↓
Which long-read information actually helped?
```

The classifier is the measuring instrument.

The **measured information gain from long-read SV characterization** is the project.

---

# 61. Project Identity Summary

**Title:**  
SV Pathogenicity Classifier with Short-Read vs Long-Read Feature Ablation

**Task:**  
Binary structural-variant pathogenicity classification.

**Baseline:**  
Short-read-derivable structural, population, and genomic-context features.

**Experimental arm:**  
Baseline + long-read-enabled breakpoint, insertion-sequence, phasing, repeat, and microhomology information.

**Primary metric:**  
AUPRC.

**Primary result:**  
ΔAUPRC between Arm B and Arm A.

**Base paper:**  
Sharo et al. 2022 — StrVCTVRE.

**Main population reference:**  
Sudmant et al. 2015.

**Main long-read / benchmark family:**  
GIAB / HG002 and compatible long-read SV resources.

**Core warning:**  
Verify that a sufficiently large labelled short-read/long-read matched cohort actually exists before committing to the final ablation.

**Core novelty:**  
Quantify the value of long-read-exclusive information rather than merely asserting that long reads are better.

---

# 62. Provenance Note

This context pack was created from the project's original DNA Sequencing Technology idea catalogue and then expanded into an AI-development specification.

External resource status was checked on **2026-08-18** for the major references/resources. Genomics portals, releases, accessions, and filenames can change; future AI agents must re-verify them before use.
