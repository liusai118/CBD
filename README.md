# An AI-Integrated Pharmacophore and Transcriptomic Framework for Rapid Discovery of Neurotherapeutic Leads

This repository contains code, data, and analyses for an AI-integrated framework that combines pharmacophore-inspired modeling and transcriptomic prediction to rapidly identify neurotherapeutic leads.

---

## Table of Contents
- [Results 1. Construction of a Machine-Learning Model Integrating Stroke-Associated Biological Pathways](#results-1-construction-of-a-machine-learning-model-integrating-stroke-associated-biological-pathways)
- [Results 2. Development and Application of the LGE-GNN Model for Predicting Gene-Expression Profiles](#results-2-development-and-application-of-the-lge-gnn-model-for-predicting-gene-expression-profiles)
- [Results 3. Establishment of the DDN Model for Predicting Bioactive Substructures](#results-3-establishment-of-the-ddn-model-for-predicting-bioactive-substructures)
- [Results 4. AI-Aided Drug Repurposing Identifies CBD as a Promising Therapeutic Candidate for Ischemic Stroke (IS)](#results-4-ai-aided-drug-repurposing-identifies-cbd-as-a-promising-therapeutic-candidate-for-ischemic-stroke-is)
- [Results 5. LGE-GNN-Guided Identification of the NRF2 Pathway as a Key Mediator of CBD’s Neuroprotective Effects](#results-5-lge-gnn-guided-identification-of-the-nrf2-pathway-as-a-key-mediator-of-cbds-neuroprotective-effects)
- [Results 6. DeepD2V-Enhanced Identification of BMAL1 as a Centrally NRF2-Regulated Target Driving CBD-Induced Angiogenesis](#results-6-deepd2v-enhanced-identification-of-bmal1-as-a-centrally-nrf2-regulated-target-driving-cbd-induced-angiogenesis)

---

## Results 1. Construction of a Machine-Learning Model Integrating Stroke-Associated Biological Pathways

- **Population: GBD analysis**
  - `pre.R` — Script for GBD analysis
- **MR** — Mendelian randomization of stroke and anxiety
- **compound_pathway** — Machine-learning model integrating stroke-associated biological pathways
  - `new_test.csv` — List of anxiolytics
  - `training2.csv` — List of compounds used in the training set of LGE-GNN and DDN
  - `Tanimoto.ipynb` — Script for computing Tanimoto similarities
  - `other` — Scripts and datasets for compound–pathway models

---

## Results 2. Development and Application of the LGE-GNN Model for Predicting Gene-Expression Profiles

- `LGE-GNN.py` — Implementation of the LGE-GNN model

---

## Results 3. Establishment of the DDN Model for Predicting Bioactive Substructures

- `DDNmodel.py` — Implementation of the DDN model

---

## Results 4. AI-Aided Drug Repurposing Identifies CBD as a Promising Therapeutic Candidate for Ischemic Stroke (IS)

- **All results in this section are derived from biological experiments.**

---

## Results 5. LGE-GNN-Guided Identification of the NRF2 Pathway as a Key Mediator of CBD’s Neuroprotective Effects

- **Step 1: Whole-transcriptome prediction model**
  - Required file: `bgedv2_QNORM.gctx` (download: <https://cbcl.ics.uci.edu/public_data/D-GEX/bgedv2_QNORM.gctx>)
- **Step 2: Downstream analysis of Step 1**

---

## Results 6. DeepD2V-Enhanced Identification of BMAL1 as a Centrally NRF2-Regulated Target Driving CBD-Induced Angiogenesis

- **step1_datapre** — Data preparation for training DeepD2V
  - `all2.bed` — Consolidated BED file from ENCODE accessions ENCSR000ESK, ENCSR584GHV, and ENCSR707IUN
  - `datapre.R` — Script for preparing the training set of DeepD2V
  - `testpre.R` — Script for preparing the test set of DeepD2V
- **step2_DEEPD2V** — Training DeepD2V
  - `main.py` — Script for training DeepD2V
  - `predict.py` — Script for predicting NRF2 downstream targets
- **step3_RNA-seq** — RNA-seq analysis
  - `RNAseq.sh` — Script for upstream analysis of RNA-seq
  - `RNAseq.R` — Script for downstream analysis of RNA-seq
- **step4_CUT&Tag** — CUT&Tag analysis
  - `CUT&Tag.sh` — Script for upstream analysis of CUT&Tag
  - `CUT&Tag.R` — Script for downstream analysis of CUT&Tag
