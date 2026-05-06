# GOR Cross-Validation Pipeline

Java-based utilities for generating cross-validation folds and running automated evaluation workflows for GOR secondary structure prediction models.

The repository contains two main tools:

- `TrainSetCreator.jar`
- `CrossValidationRunner.jar`

Together, they provide a complete workflow for:
- generating k-fold train/test splits
- training GOR models
- predicting secondary structures
- validating predictions
- computing mean and standard deviation for Q3 and SOV scores across folds

---

# Workflow Overview

The cross-validation workflow consists of two steps:

1. Generate fold datasets using `TrainSetCreator`
2. Run the full cross-validation pipeline using `CrossValidationRunner`

---

# Repository Contents

```text
.
├── TrainSetCreator.jar
├── CrossValidationRunner.jar
├── src/
├── examples/
├── folds/
└── README.md
```

---

# TrainSetCreator

`TrainSetCreator` creates reproducible k-fold train/test splits from a Seclib dataset.

---

## Usage

```bash
java -jar TrainSetCreator.jar \
  --db <seclib-file> \
  --out <output-prefix> \
  --cv <k> \
  [--seed <long>] \
  [--mafSource <alignment-folder>]
```

---

## Example

```bash
java -jar TrainSetCreator.jar \
  --db data/pc25_train_set.db \
  --out folds/dataset \
  --cv 5 \
  --seed 42
```

---

## Arguments

| Argument | Description |
|---|---|
| `--db` | Input Seclib dataset |
| `--out` | Output prefix for generated folds |
| `--cv` | Number of folds |
| `--seed` | Random seed for reproducible shuffling |
| `--mafSource` | Optional alignment folder for GOR5 workflows |

---

## Generated Files

For `k = 5`, the tool creates:

```text
dataset_fold1_train.db
dataset_fold1_test.db
dataset_fold2_train.db
dataset_fold2_test.db
...
dataset_fold5_train.db
dataset_fold5_test.db
```

---

## Alignment Filtering

If `--mafSource` is provided:

- proteins without matching alignment files are excluded
- supported alignment extensions:
  - `.aln`
  - `.ma`
  - `.maf`

Alignment files are matched by:
1. exact filename
2. unique prefix match

This is primarily used for alignment-based GOR5 workflows.

---

# CrossValidationRunner

`CrossValidationRunner` runs the complete cross-validation workflow using previously generated folds.

For each fold, it:

- trains a model
- predicts secondary structures
- validates predictions
- computes Q3 and SOV metrics
- stores fold-specific outputs
- aggregates overall statistics

---

## Usage

```bash
java -jar CrossValidationRunner.jar \
  --foldPrefix <prefix> \
  --cv <k> \
  --trainJar <train.jar> \
  --predictJar <predict.jar> \
  --validateJar <validateGor.jar> \
  --out <output-folder> \
  --method <method> \
  [--mafSource <alignment-folder>] \
  [--format <txt|html>]
```

---

## Standard GOR Example

```bash
java -jar CrossValidationRunner.jar \
  --foldPrefix folds/dataset \
  --cv 5 \
  --trainJar trainBruno.jar \
  --predictJar predictBruno.jar \
  --validateJar validateGor.jar \
  --out results \
  --method gor4 \
  --format txt
```

---

## GOR5 Example

```bash
java -jar CrossValidationRunner.jar \
  --foldPrefix folds/dataset \
  --cv 5 \
  --trainJar trainBruno.jar \
  --predictJar predictBruno.jar \
  --validateJar validateGor.jar \
  --out results \
  --method gor5gor4 \
  --mafSource data/maf \
  --format txt
```

---

## Supported Methods

| Method | Description |
|---|---|
| `gor1` | Standard GOR I prediction |
| `gor3` | Standard GOR III prediction |
| `gor4` | Standard GOR IV prediction |
| `gor5gor1` | Alignment-based prediction using GOR I |
| `gor5gor3` | Alignment-based prediction using GOR III |
| `gor5gor4` | Alignment-based prediction using GOR IV |

---

## Arguments

| Argument | Description |
|---|---|
| `--foldPrefix` | Prefix used by `TrainSetCreator` |
| `--cv` | Number of folds |
| `--trainJar` | Path to training jar |
| `--predictJar` | Path to prediction jar |
| `--validateJar` | Path to validation jar |
| `--out` | Output directory |
| `--method` | Prediction method |
| `--mafSource` | Alignment folder (required for GOR5 methods) |
| `--format` | Validation output format (`txt` or `html`) |

---

# Generated Outputs

For each fold, the pipeline creates:

```text
fold1_gor4_model.txt
fold1_gor4_pred.txt
fold1_gor4_summary.txt
fold1_gor4_detail.txt
```

Additionally, the final summary file:

```text
cv_summary.txt
```

contains:

```text
fold_metric    mean    stddev
Q3             ...
SOV            ...
```

---

# Validation Metrics

## Q3 Accuracy

Q3 measures the percentage of correctly predicted residues across:

- helix (`H`)
- sheet (`E`)
- coil (`C`)

---

## SOV Score

The Segment Overlap (SOV) score evaluates overlap quality between predicted and reference secondary structure segments.

The workflow computes:
- overall SOV
- helix SOV
- sheet SOV
- coil SOV

---

# Notes

- Cross-validation folds must be generated before running `CrossValidationRunner`
- GOR5 methods require alignment files
- Fold shuffling is reproducible using a fixed random seed
---
