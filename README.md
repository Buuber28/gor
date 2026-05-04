# GOR IV 

Java-based implementation of a GOR IV protein secondary structure prediction pipeline.  
The project provides an end-to-end workflow for training, prediction, and validation, including support for FASTA input, multiple sequence alignments (GOR V), and optional probability outputs.

## Features

- Implementation of the GOR IV algorithm
- Model training from labeled protein sequence data
- Secondary structure prediction from FASTA sequences
- Prediction from multiple sequence alignments (alignment-based extension)
- k-fold cross-validation support (project-level)
- Optional probability output for helix (H), sheet (E), and coil (C)
- Post-processing to remove short, unlikely structural segments

## Input Formats

### Training Data (Seclib format)

```text
>protein_id
AS AMINOACIDSEQUENCE
SS SECONDARYSTRUCTURE
```

### Prediction input (FASTA)

```text
>protein_id
AMINOACIDSEQUENCE
```

### Multiple Alignment Input

- supported file extensions:
  .aln
  .ma
  .maf

## Usage
The program uses the provided command-line arguments to decide whether to run in training or prediction mode.

- If `--db` is provided, the program runs in training mode.
- Otherwise, the program runs in prediction mode.
- Prediction requires exactly one of `--seq` or `--maf`.

### Train a Model
```bash
java -jar gor4.jar --db <seclib-file> --method gor4 --model <model-file>
```
### Predict from FASTA
```bash
java -jar gor4.jar --model <model-file> --format txt --seq <fasta-file>
```

### Predict from Multiple Alignments
```bash
java -jar gor4.jar --model <model-file> --format txt --maf <alignment-folder>
```

### Enable Probability Output
```bash
java -jar gor4.jar --probabilities --model <model-file> --format txt --seq data/example.fasta
```

## Output Format

### Standard Ouput
```text
>protein_id
AS AMINOACIDSEQUENCE
PS PREDICTEDSTRUCTURE
```

### With probability output enabled:
```text
>protein_id
AS AMINOACIDSEQUENCE
PS PREDICTEDSTRUCTURE
PH HELIX_PROBABILITIES
PE SHEET_PROBABILITIES
PC COIL_PROBABILITIES
```
Probability values are scaled from 0 to 9
Higher values indicate higher confidence

## How It Works

### Training

- Reads labeled protein sequences
- Computes frequency-based matrices (4D and 6D)
- Stores model parameters

### Prediction

- Computes log-odds scores for each structure class: coil (C), sheet (E), and helix (H)
- Assigns the most probable structure class to each residue

### Alignment-Based Prediction

- Aggregates predictions across homologous sequences
- Averages probabilities per alignment column
- Uses evolutionary information to improve prediction robustness

### Post-Processing

- Enforces minimum segment lengths
- Removes isolated predictions
- Smooths local inconsistencies

## Notes

- Window size: 17 amino acids
- Sequence edges are marked with `-` because the full prediction window is not available
- Only `txt` output is fully supported
