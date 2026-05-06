# GOR IV Protein Secondary Structure Prediction Pipeline

Java-based implementation of a GOR IV protein secondary structure prediction pipeline.

The project provides an end-to-end workflow for:
- model training
- secondary structure prediction
- alignment-based prediction
- optional probability output
- post-processing of predictions

The program automatically decides whether to run in training or prediction mode based on the provided command-line arguments.

---

## Features

- GOR IV-based secondary structure prediction
- Model training from labeled protein sequence data
- Prediction from FASTA sequences
- Alignment-based prediction from multiple sequence alignments
- Optional probability/confidence output
- Post-processing to remove unlikely short secondary structure segments

---

## Tech Stack

- Java
- SQL / MySQL
- HTML
- Git

---

## Input Formats

### Training Data (Seclib Format)

Training data must follow the format:

```text
>protein_id
AS AMINOACIDSEQUENCE
SS SECONDARYSTRUCTURE
```

Example:

```text
>1ABC
AS MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMF
SS CCCCCCHHHHHHHHHHHCCCCCCEEEECCCCCCC
```

---

### FASTA Input

```text
>protein_id
AMINOACIDSEQUENCE
```

Example:

```text
>ExampleProtein
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMF
```

---

### Multiple Alignment Input

Alignment files must contain:
- one target sequence
- homologous aligned sequences

Supported file extensions:
- `.aln`
- `.ma`
- `.maf`

Example:

```text
>TargetProtein
AS MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMF
1 MVLSPADKTNIKAAWGKVGAHAGEYGAEALERMF
2 MVLSPADKTNVKAAWGKVGPHAGEYGAEALERMF
```

---

## Usage

The program switches between training and prediction mode automatically:

- If `--db` is provided → training mode
- Otherwise → prediction mode

Prediction mode requires exactly one of:
- `--seq`
- `--maf`

---

## Training

### Train a Model

```bash
java -jar gor4.jar --db <seclib-file> --method <gor1|gor3|gor4> --model <model-file>
```

Example:

```bash
java -jar gor4.jar --db data/train.db --method gor4 --model models/gor4.model
```

### Required Arguments

| Argument | Description |
|---|---|
| `--db` | Path to training database |
| `--model` | Output model file |

---

## Prediction

### Predict from FASTA

```bash
java -jar gor4.jar --model <model-file> --format txt --seq <fasta-file>
```

Example:

```bash
java -jar gor4.jar --model models/gor4.model --format txt --seq data/example.fasta
```

---

### Predict from Multiple Alignments

```bash
java -jar gor4.jar --model <model-file> --format txt --maf <alignment-folder>
```

Example:

```bash
java -jar gor4.jar --model models/gor4.model --format txt --maf data/maf/
```

---

### Enable Probability Output

```bash
java -jar gor4.jar --probabilities --model <model-file> --format txt --seq data/example.fasta
```

---

### Prediction Arguments

| Argument | Description |
|---|---|
| `--model` | Path to trained model |
| `--format` | Output format (`txt` or `html`) |
| `--seq` | FASTA input file |
| `--maf` | Folder containing multiple alignment files |
| `--probabilities` | Enables probability output |

---

## Output Format

### Standard Output

```text
>protein_id
AS AMINOACIDSEQUENCE
PS PREDICTEDSTRUCTURE
```

Example:

```text
>ExampleProtein
AS MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMF
PS CCCCCCHHHHHHHHHHHCCCCCCEEEECCCCCCC
```

---

### Probability Output

When `--probabilities` is enabled:

```text
>protein_id
AS AMINOACIDSEQUENCE
PS PREDICTEDSTRUCTURE
PH HELIX_PROBABILITIES
PE SHEET_PROBABILITIES
PC COIL_PROBABILITIES
```

Example:

```text
>ExampleProtein
AS MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMF
PS CCCCCCHHHHHHHHHHHCCCCCCEEEECCCCCCC
PH 0001236789987654321000000000000000
PE 0000000000000000012345678999999999
PC 9998763210000000009876543211111111
```

Probability values are scaled from `0` to `9`:
- `0` = low confidence
- `9` = high confidence

---

## How It Works

### Training

- Reads labeled protein sequences
- Computes frequency-based matrices (4D and 6D)
- Stores model parameters

### Prediction

- Uses a sliding window of 17 amino acids
- Computes log-odds scores for:
  - coil (`C`)
  - sheet (`E`)
  - helix (`H`)
- Assigns the most probable structure class to each residue

### Alignment-Based Prediction

- Aggregates predictions across homologous sequences
- Averages probabilities per alignment column
- Uses evolutionary information to improve prediction robustness

### Post-Processing

- Enforces minimum segment lengths
- Removes isolated predictions
- Smooths local inconsistencies

---

## Notes

- Window size: 17 amino acids
- Sequence edges are marked with `-` because the full prediction window is not available
- Only `txt` output is fully supported
- The `html` output flag exists but may not be fully implemented
- Alignment-based prediction ignores gaps and unknown amino acids
- Developed as part of the *Programming Praktikum Bioinformatics* course

---

