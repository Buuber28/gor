# GOR Validation Tool

Java-based validation tool for evaluating predicted protein secondary structures against reference secondary structures.

The tool compares prediction files with reference Seclib files and computes:

- Q3 accuracy
- Class-specific Q3 values for helix, sheet, and coil
- SOV scores
- Class-specific SOV values for helix, sheet, and coil

---

## Usage

```bash
java -jar validateGor.jar -p <predictions> -r <seclib-file> -s <summary-file> -d <detailed-file> -f <txt|html>
```

Example:

```bash
java -jar validateGor.jar \
  -p examples/predictions.txt \
  -r data/reference.db \
  -s results/summary.txt \
  -d results/details.txt \
  -f txt
```

---

## Arguments

| Argument | Description |
|---|---|
| `-p` | Prediction file |
| `-r` | Reference file in Seclib format |
| `-s` | Output path for summary results |
| `-d` | Output path for detailed per-protein results |
| `-f` | Output format: `txt` or `html` |

---

## Prediction Input Format

The prediction file is expected to contain entries in the following format:

```text
>protein_id
AS AMINOACIDSEQUENCE
PS PREDICTEDSTRUCTURE
```

Probability lines are ignored if present:

```text
PH ...
PE ...
PC ...
```

---

## Reference Input Format

The reference file is expected to follow Seclib-style formatting:

```text
>protein_id
AS AMINOACIDSEQUENCE
SS REFERENCESTRUCTURE
```

---

## Output

### Summary Output

The summary file contains five-number summaries for each metric:

```text
q3: min q25 median q75 max
qObs_H: min q25 median q75 max
qObs_E: min q25 median q75 max
qObs_C: min q25 median q75 max
SOV: min q25 median q75 max
SOV_H: min q25 median q75 max
SOV_E: min q25 median q75 max
SOV_C: min q25 median q75 max
```

### Detailed Output

The detailed output contains one entry per matched protein:

```text
>protein_id q3 SOV qH qE qC SOV_H SOV_E SOV_C
AS AMINOACIDSEQUENCE
PS PREDICTEDSTRUCTURE
SS REFERENCESTRUCTURE
```

---

## Metrics

### Q3 Accuracy

Q3 measures the percentage of correctly predicted residues across the three secondary structure classes:

- `H` = helix
- `E` = sheet
- `C` = coil

Positions marked with `-` in the prediction are ignored.

### SOV Score

The Segment Overlap (SOV) score evaluates how well predicted secondary structure segments overlap with reference segments.

The implementation computes:

- overall SOV
- SOV for helix segments
- SOV for sheet segments
- SOV for coil segments

---

## Notes

- Prediction and reference entries are matched by protein ID.
- Probability lines are ignored during validation.
- Positions marked with `-` are excluded from Q3 and SOV calculation.
- This tool validates existing predictions; cross-validation workflow documentation is provided in ../CrossValidation directoy.
