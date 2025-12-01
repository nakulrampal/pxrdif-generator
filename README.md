# PXRDIF Generator

Converts CSV files to PXRDIF format.

## Features

- Parses experimental conditions (precursors, solvents, catalysts, temperature, duration)
- Converts PXRD data (theta, intensity)
- Handles inline and column-based CSV formats
- Supports "no catalyst" experiments
- Handles dash types: -, –, —
- Batch processing

## Files

- `excel_to_pxrdif_converter.py` - Main converter
- `batch_convert_pxrdif.py` - Batch processor

## Requirements

```
pandas
numpy
```

## Usage

### Single Conversion

```python
from excel_to_pxrdif_converter import excel_to_pxrdif_multiple

excel_to_pxrdif_multiple(
    'conditions.csv',
    'pxrd_data.csv',
    output_dir='output',
    verbose=False
)
```

### Batch Conversion

```python
python batch_convert_pxrdif.py
```

## CSV Format

### Conditions CSV
- Experiment ID
- Precursor 1 amount
- Precursor 2 amount
- Solvent 1 (name and volume)
- Solvent 2 (name and volume)
- Solvent 3 (optional)
- Catalyst (name and volume/mass) or "-"
- Total volume
- Rationale

### PXRD CSV
- Column 1: theta values
- Other columns: intensity (labeled by experiment ID)

## Output

PXRDIF files with:
- Experimental conditions
- Synthesis parameters
- PXRD diffraction data

---

*Developed with assistance from GitHub Copilot (Claude Sonnet 4.5)*
