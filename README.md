 # Integrated PanGWAS–ML Analysis of Recurrent UTI in *Escherichia coli*

This repository contains the computational workflows, analysis scripts, processed datasets, intermediate results, and visualization outputs associated with the study:

**Integrated Machine Learning–PanGWAS Reveals Chromosome-Encoded Persistence Networks and Plasmid Plasticity in Recurrent Urinary Tract Infection in *Escherichia coli***

---

## Overview

Recurrent urinary tract infection (rUTI) represents a major clinical challenge, and the genomic features that contribute to recurrent infection remain incompletely understood.

This study integrates **pangenome-wide association analysis (PanGWAS)** and **machine learning (ML)** to investigate genomic features associated with recurrent and sporadic urinary tract infection in *Escherichia coli*.

The computational framework evaluates genomic signals across three genomic compartments:

* **Combined chromosome + plasmid dataset**
* **Chromosome-only dataset**
* **Plasmid-only dataset**

Multiple complementary approaches were integrated to identify genomic features associated with recurrent UTI, including pangenome analysis, population structure analysis, machine learning, Pan-GWAS, and association rule mining.

## Study Dataset

The primary study dataset consists of **491 *Escherichia coli* genomes** associated with recurrent UTI and non recurrent UTI.

An independent external validation dataset consisting of **63 *E. coli* genomes** was used to evaluate the robustness of the predictive genomic signatures.

The external validation dataset comprises:

* **31 recurrent UTI isolates**
* **32 sporadic UTI isolates**

Genome sequences were obtained from publicly available repositories, including:

* National Center for Biotechnology Information (NCBI) Assembly
* European Nucleotide Archive (ENA)

The accession numbers and associated phenotype annotations used in this study are provided in the repository.

---

## Computational Workflow

The overall computational workflow consists of the following major stages:

                    GENOME DATASET
                          │
                          ▼
                    MOB-suite
                          │
          ┌───────────────┼───────────────┐
          ▼               ▼               ▼
     CHROMOSOME        PLASMID         COMBINED
          │               │               │
          └───────┬───────┘               │
                  │                       │
                  ▼                       ▼
              Panaroo                 Panaroo
                  │                       │
                  │                       ▼
                  │               Core gene alignment
                  │                       │
                  │                       ▼
                  │                    Gubbins
                  │                       │
                  │                       ▼
                  │               Recombination-filtered
                  │                 core alignment
                  │                       │
                  │                       ▼
                  │                    IQ-TREE
                  │                       │
                  │                       ▼
                  │                   Phylogeny
                  │
                  │
                  └───────────────────────────────┐
                                                  │
                                                  ▼
                              MACHINE LEARNING + FEATURE SELECTION
                                                  │
                          ┌───────────────────────┼───────────────────────┐
                          ▼                       ▼                       ▼
                     Chromosome                Plasmid               Combined
                          │                       │                       │
                          └───────────────────────┼───────────────────────┘
                                                  │
                                                  ▼
                                      Feature selection
                                  ┌─────────┬─────────┬─────────┐
                                  ▼         ▼         ▼
                              Chi-square   MI      Boruta
                                  └─────────┬─────────┘
                                            ▼
                                     STRICT / RELAXED
                                            │
                                            ▼
                                     ML classification
                                            │
                                            ▼
                              Genomic compartment comparison
                                            │
                                            ▼
                         Identify recurrence-associated features
                       across chromosome / plasmid / combined data
                                            │
                                            ▼
                               COMBINED DATASET SELECTED
                                            │
                              ┌─────────────┴─────────────┐
                              ▼                           ▼
                           Scoary                       Pyseer
                              │                           │
                              └─────────────┬─────────────┘
                                            ▼
                                Integrate ML + Scoary +
                                      Pyseer features
                                            │
                                            ▼
                                  SEEDED GENE LIST
                                            │
                                            ▼
                                    Apriori analysis
                                            │
                                            ▼
                              Co-occurring gene networks
                                            │
                                            ▼
                              Biological interpretation
                                            │
                                            ▼
                                External validation
    
---

## Genome Quality Control and Processing

Genome datasets were subjected to quality assessment and filtering prior to downstream analysis. The processing strategy differed depending on whether an assembled genome or raw sequencing reads were available.

A subset of recurrent UTI (rUTI) genomes were obtained from NCBI/ENA as raw genome data. Assembled genomes were therefore used directly for genome-quality assessment and subsequent analyses without re-performing read-level quality control, trimming, or genome assembly.

For the remaining rUTI isolates for which raw sequencing reads were available, sequencing data were subjected to read-level quality assessment using FastQC v0.12.1, followed by quality trimming using Trimmomatic v0.39. The filtered reads were subsequently assembled using Unicycler v0.5.0, where applicable. The resulting genome assemblies were then assessed for genome quality using CheckM v1.2.2.

Genome assemblies were retained for downstream analysis according to the predefined quality criteria:

* N50 ≥ 150,000 bp
* Number of contigs < 100
* Completeness > 99%

Thus, FastQC, Trimmomatic, and Unicycler were applied only to datasets for which raw sequencing reads were available, whereas publicly available assembled genomes were subjected directly to assembly-level quality assessment and filtering. The accession sources, processing status, and quality-filtering results for individual genomes are documented in the corresponding metadata and analysis files.

## Chromosome and Plasmid Classification

Genome assemblies were compartmentalized into chromosome and plasmid components using **MOB-suite v3.1.9**.

The resulting genomic compartments were subsequently used to construct:

1. Chromosome-only datasets
2. Plasmid-only datasets
3. Combined chromosome + plasmid datasets

Plasmid-associated features including plasmid mobility, relaxase types, mating pair formation (MPF) types, and incompatibility groups were further investigated.

---

## Pangenome Reconstruction

Pangenome analysis was performed using **Panaroo v1.3.4**.

Pangenome reconstruction was performed for the relevant genomic compartments to generate gene presence–absence matrices for downstream population analysis, machine learning, and genome-wide association analysis.

Processed gene presence–absence matrices generated during the study are provided where appropriate.

---

## Phylogenetic Analysis

To reduce the influence of homologous recombination on phylogenetic inference:

* Core-genome alignments were generated from the pangenome analysis.
* Recombination regions were identified and filtered using **Gubbins v3.3.0**.
* Maximum-likelihood phylogenetic analysis was performed using **IQ-TREE v2.2.0**.
* Phylogenetic trees were subsequently visualized and annotated.

The final phylogenetic outputs and associated scripts are provided in the corresponding repository directory.

---

## Population Structure Analysis

Population structure was investigated using multiple complementary approaches.

### Principal Component Analysis

Principal component analysis (PCA) was performed using gene presence–absence matrices and implemented using Python and scikit-learn.

### Multidimensional Scaling

Multidimensional scaling (MDS) was performed to further evaluate relationships among isolates based on genomic feature profiles.

### Phylogroup Classification

*E. coli* phylogroups were determined using **EzClermont**.

### Multilocus Sequence Typing

Multilocus sequence typing (MLST) was performed according to the **Achtman *E. coli* MLST scheme**.

---

## Pangenome Structure Analysis

Pangenome characteristics were evaluated using the gene presence–absence information generated by Panaroo.

The analysis included:

* Pangenome openness
* Heap's Law analysis
* Core genome diversity
* Core genome shrinkage
* Gene accumulation patterns

Heap's Law parameters were calculated using custom Python scripts developed for this study.

---

# Machine Learning Analysis

Machine learning was used to evaluate whether recurrent and sporadic UTI isolates could be discriminated based on genomic features.

Three genomic datasets were independently evaluated:

* **Combined chromosome + plasmid**
* **Chromosome-only**
* **Plasmid-only**

## Feature Selection

Three complementary feature-selection approaches were applied:

* Chi-square test
* Mutual information
* Boruta

Feature selection was performed to identify genomic features associated with discrimination between recurrent and sporadic UTI isolates.

### STRICT Feature Set

The **STRICT** feature set represents high-confidence features identified consistently across the feature-selection approaches and cross-validation procedures.

### RELAXED Feature Set

The **RELAXED** feature set represents a broader collection of features identified by any of the feature-selection approaches.

STRICT and RELAXED feature sets were generated independently for chromosome, plasmid, and combined genomic datasets.

---

## Machine Learning Models

The following supervised machine learning classifiers were evaluated:

* Logistic Regression
* Random Forest
* ExtraTrees
* Support Vector Machine
* XGBoost

Model performance was evaluated using appropriate test-set and cross-validation procedures.

Performance metrics include:

* Accuracy
* Area under the receiver operating characteristic curve (AUC)
* Precision
* Recall
* F1-score
* Confusion matrices
* ROC curves

The scripts used for feature selection, model training, cross-validation, performance evaluation, and visualization are provided in the machine-learning directory.

---

# Pan-GWAS Analysis

Genome-wide association analysis was performed to identify genes associated with recurrent UTI.

Two complementary approaches were used:

### Scoary

**Scoary v1.6.16** was used for gene-level association analysis based on gene presence–absence patterns and phenotype classification.

### Pyseer

**Pyseer v1.3.1** was used for association analysis incorporating population structure through a kinship matrix.

Results from complementary analytical approaches were integrated to identify genes showing consistent association with the recurrent UTI phenotype.

---

# Association Rule Mining

Genes identified as associated with recurrent UTI through complementary analytical approaches were used to generate a seeded gene list.

Association rule mining was subsequently applied to the gene presence–absence matrices to identify co-occurring genomic features.

This analysis was used to investigate potential cooperative genomic patterns rather than focusing exclusively on individual gene effects.

The resulting association rules and derived gene networks are provided where appropriate.

---

# External Validation

An independent dataset of **63 *E. coli* genomes** was used for external validation.

The external dataset consisted of:

* 31 recurrent UTI isolates
* 32 sporadic UTI isolates

These genomes were not used during model training.

The same genomic feature-generation framework was applied to the external dataset, and predictive models trained on the primary dataset were evaluated using the independent genomes.

External validation metrics include:

* Accuracy
* Precision
* Recall
* F1-score
* AUC
* Confusion matrix
* ROC analysis

The external validation dataset, metadata, processed feature matrices, scripts, and results are organized separately from the primary training dataset.

---

# Repository Organization

The repository is organized according to the major stages of the computational workflow:

```text
Integrated_PanGWAS-ML_rUTI_Ecoli/
│
├── README.md
├── LICENSE
├── .gitignore
│
├── 01_data/
│   ├── accession_lists/
│   └── metadata/
│
├── 02_genome_QC/
│   └── scripts/
│
├── 03_MOB_suite/
│   └── scripts/
│
├── 04_Panaroo/
│   └── scripts/
│
├── 05_phylogeny/
│   └── scripts/
│
├── 06_population_structure/
│   └── scripts/
│
├── 07_pan_genome_analysis/
│   └── scripts/
│
├── 08_machine_learning/
│   └── scripts/
│
├── 09_PanGWAS/
│   └── scripts/
│
├── 10_association_rules/
│   └── scripts/
│
├── 11_figures/
│
├── 12_tables/
│
├── 13_external_validation/
│   └── scripts/
│
└── environment/
```

The repository structure may be updated as additional analysis scripts and processed files are organized.

---

# Software and Versions

The principal software used in the computational workflow includes:

|| Software / Tool |     Version |
| --------------- | ----------: |
| FastQC          |      0.12.1 |
| Trimmomatic     |        0.39 |
| Unicycler       |       0.5.0 |
| CheckM          |       1.2.2 |
| MOB-suite       |       3.1.9 |
| Panaroo         |       1.3.4 |
| Gubbins         |       3.3.0 |
| IQ-TREE         |       2.2.0 |
| Scoary          |      1.6.16 |
| Pyseer          |       1.3.1 |
| scikit-learn    |       1.7.2 |
| XGBoost         |       3.1.2 |
| Boruta          |       0.4.3 |
| Python          |     3.10.19 |


Exact software versions and computational dependencies will be documented in the `environment/` directory.

---

# Data Availability

All genome sequences analysed in this study were obtained from publicly available datasets deposited in the **NCBI Assembly database** and the **European Nucleotide Archive (ENA)**.

The accession numbers used for the primary and external datasets are provided in the corresponding accession files within this repository.

The original genome sequences are **not redistributed through this repository**. Users should retrieve the corresponding sequence data directly from NCBI, ENA, or the original public data record using the accession numbers provided and comply with the applicable terms and conditions associated with those records.

---

# Study-Derived Data

Where appropriate, this repository provides study-derived and processed data, including:

* Curated accession lists
* Phenotype annotations
* Processed gene presence–absence matrices
* Machine-learning feature lists
* STRICT feature sets
* RELAXED feature sets
* Pan-GWAS results
* Seeded gene lists
* Association-rule results
* External validation matrices
* Summary tables
* Analysis outputs

Files containing information that cannot appropriately be redistributed will not be included.

---

# Reproducibility

The repository is intended to provide the computational workflow required to reproduce the analyses described in the associated manuscript.

Scripts are organized according to the corresponding stages of the Materials and Methods section.

Where applicable, each analysis directory will contain:

* Input files
* Analysis scripts
* Command-line execution instructions
* Processed intermediate files
* Final analysis outputs
* Figure-generation scripts

Large raw sequencing datasets and publicly available genome assemblies are not redistributed and should instead be retrieved using the accession numbers provided.

---

# Licensing

## Original scripts and code

Original computational scripts developed for this study are made available under the **MIT License**.

## Study-derived datasets

Study-derived and curated datasets made available through this repository are intended to be released under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license, unless otherwise specified for an individual file or dataset.

## Public genome sequences

Genome sequences obtained from NCBI, ENA, or other public repositories are **not relicensed by this repository**. The original terms and conditions associated with the respective database records and source datasets apply.

## Third-party software

Third-party software and databases used in this workflow remain subject to their respective licenses and terms of use.

---

# Citation

If you use the computational scripts, processed datasets, or analytical workflows provided in this repository, please cite the associated publication:

**[Full publication citation will be added upon publication.]**

Repository DOI:

**[DOI will be added after repository archival.]**

---

# Contact

For questions regarding the computational workflow, data processing, or reproducibility of the analyses, please contact the corresponding author listed in the associated publication.
