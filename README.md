# 🧬 MDT-TM_CRC_Microbiome-Framework

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18760840.svg)](https://doi.org/10.5281/zenodo.18760840)

Genome-resolved workflow implementing Minimal Dominance Threshold–Transition Modeling (MDT-TM) to quantify tipping-point dynamics in colorectal cancer microbiome progression. This workflow builds on and extends the GRUMB framework.

This repository contains the full computational workflow used in:

**Ascandari A, Aminu S, Benhida R, Daoud R.**  
*Modeling compositional perturbations and predicted state transitions in the colorectal cancer microbiome.*

The study introduces **Minimal Dominance Threshold – Transition Modeling (MDT-TM)**, a genome-resolved ecological modeling framework that quantifies tipping-point dynamics along the colorectal cancer (CRC) continuum (Healthy → Adenoma → Carcinoma).


## Analytical Workflow

```mermaid
flowchart TB
    INPUT["Public SRA metagenomes"]

    subgraph BIO["01 · Genome-resolved bioinformatics"]
        direction TB
        B1["Assembly and binning"]
        B2["MAG dereplication and taxonomy"]
        B3["Abundance estimation and species mapping"]
        B1 --> B2 --> B3
    end

    subgraph PRE["02 · Quality control and preprocessing"]
        direction TB
        P1["MAG and assembly quality assessment"]
        P2["Contamination screening and coverage assessment"]
        P3["Profile cleanup and metadata harmonization"]
        P1 --> P2 --> P3
    end

    DATA["Species abundance matrix and harmonized metadata"]

    subgraph ECO["03 · Species ecology"]
        direction TB
        E1["Species matrix processing"]
        E2["CLR transformation and batch correction"]
        E3["Alpha diversity"]
        E4["Beta diversity and community structure"]
        E5["Indicator species and prevalence"]
        E1 --> E2
        E1 --> E3
        E2 --> E4
        E1 --> E5
    end

    subgraph MDT["04 · MDT-TM transition modeling"]
        direction TB
        M1["Model benchmarking and SVM observer"]
        M2["Synthetic disease-signal blending"]
        M3["MDT estimation and instability profiling"]
        M4["Directional perturbation analysis"]
        M5["MDT-based stratification and MaAsLin3"]
        M1 --> M2 --> M3
        M3 --> M4
        M3 --> M5
    end

    INPUT --> B1
    B3 --> P1
    P3 --> DATA
    DATA --> E1
    E2 --> M1

    classDef input fill:#eff6ff,stroke:#2563eb,color:#172554;
    classDef ecological fill:#ecfdf5,stroke:#059669,color:#064e3b;
    classDef modeling fill:#faf5ff,stroke:#9333ea,color:#581c87;

    class INPUT,DATA input;
    class E1,E2,E3,E4,E5 ecological;
    class M1,M2,M3,M4,M5 modeling;
```

*This diagram summarizes the analytical stages. Consult individual scripts and the manuscript Methods for input requirements and the data representation used in each analysis.*

## Repository Structure

```text
MDT-TM_CRC_Microbiome-Framework/
├── 01_Bioninformatics/
├── 02_Preprocessing/
├── 03_species_ecology/
└── 04_Synthetic_blending/
```

Directories and scripts are numbered to indicate the intended execution order. Individual scripts may require outputs from earlier stages.


---

# 01_Bioinformatics

Genome-resolved metagenomic reconstruction.

Scripts include:

- 01_run.sh  
- 02_Assembly_Binning.sh  
- 03_collect_files.sh  
- 04_global_drep_taxonomy_abundance.sh  
- 05_global_ARG.sh  
- 06_global_VFDB.sh  
- 07_clean_abundance_files.sh  
- 08_merge_abundances.sh  
- 09_mapping_to_species.sh  

Implements:

- De novo assembly (MEGAHIT)  
- Binning (MetaBAT2)  
- Dereplication (dRep)  
- Taxonomic assignment (GTDB-Tk)  
- Species-level abundance matrix construction  

---

#  02_Preprocessing

Quality control and harmonization.

Scripts include:

- 00_gather_script.sh  
- 01_mag_quality_crc_new.py  
- 02_quast_crc.py  
- 03_fastqscreen_crc_mags.py  
- 04_Nonpareil_CRC.R  
- 05_mags_vf_arg_cleanup.py  

Implements:

- MAG quality filtering  
- Assembly quality metrics  
- Host contamination screening  
- Coverage estimation (Nonpareil)  
- Metadata harmonization  

---

# 03_species_ecology

Baseline ecological characterization.

Scripts include:

- 01_species_matrix_processing.py  
- 02_CLR_batch_correction_species.R  
- 03_Species_alpha_diversity.R  
- 04_Beta_diversity_species.R  
- 05_Indicator_Specie_Analysis.R  
- 06_Species_Prevalence.R  

Implements:

- CLR transformation and batch correction  
- Alpha diversity (Shannon, Simpson, Richness, Evenness)  
- Beta diversity (Aitchison distance, NMDS, PERMANOVA)  
- Indicator species analysis  
- Core–periphery prevalence architecture  

---

# 04_Synthetic_blending

Minimal Dominance Threshold – Transition Modeling (MDT-TM)

Scripts include:

- 01_species_benchmark.py  
- 02_feature_importance.py  
- 03_blending_species.py  
- 04_metadata_threshold_mdt.py  
- 05_maaslin3_analysis.R
- LOCO_SVM.ipynb
- Sensitivity_analysis_MDT.ipynb

Implements:

- SVM observer model training
- Leave-one-cohort-out validation of the SVM model using CLR batch-corrected and non-batch-corrected matrix
- Synthetic disease-signal blending  
- Minimal Dominance Threshold (MDT) estimation
- Sensitivity analysis of the MDT threshold stability
- Shannon entropy instability profiling
- Sensitivity analysis of the Shannon entropy instability at different set thresholds
- Directional perturbation analysis  
- MDT-based sample stratification  
- MaAsLin3 differential abundance testing  

---

# Software Requirements

Core tools:

- MEGAHIT  
- MetaBAT2  
- dRep  
- GTDB-Tk  
- CheckM2  
- QUAST  
- BBTools  
- FastQ Screen  
- Nonpareil  
- CoverM  
- R (vegan, limma, indicspecies, MaAsLin3)  
- Python (scikit-learn, pandas, numpy, scipy)  

Exact versions are specified in the manuscript Methods section.

---

# Data Availability

All metagenomic datasets are publicly available from the NCBI Sequence Read Archive (SRA).  
Accession numbers are provided in Supplementary Table 1 of the manuscript.

---

# Citation

If you use this workflow, please cite:
Ascandari A et al.
Modeling compositional perturbations and predicted state transitions in the colorectal cancer microbiome. (under review).

Software archive: Zenodo: https ://doi.org/10.5281/zenodo.18760840

---

# License

Released for academic research use. Please cite appropriately when reusing the framework or scripts.
