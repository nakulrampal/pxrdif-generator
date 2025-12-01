# PXRDIF Generator

Convert Excel CSV files containing experimental conditions and PXRD data into PXRDIF format files.

## Features

- Parses experimental conditions from CSV files (precursors, solvents, catalysts, temperature, duration)
- Converts PXRD data (theta, intensity) to PXRDIF format
- Handles inline and column-based CSV formats
- Supports "no catalyst" experiments
- Handles multiple dash types (-, –, —)
- Batch processing for multiple datasets

## Files

- `excel_to_pxrdif_converter.py` - Main converter with parsing logic
- `batch_convert_pxrdif.py` - Batch processing script for multiple datasets

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
Can be inline format (entire row as quoted string) or column-based with headers:
- Experiment ID
- Precursor 1 amount
- Precursor 2 amount
- Solvent 1 (name and volume)
- Solvent 2 (name and volume)
- Solvent 3 (optional)
- Catalyst (name and volume/mass) or "-" for no catalyst
- Total volume
- Rationale

### PXRD CSV
- First column: theta values
- Subsequent columns: intensity values for each experiment (labeled by experiment ID)

## Output Format

PXRDIF files following crystallographic information file standards with:
- Experimental conditions (precursors, solvents, catalysts)
- Synthesis parameters (temperature, duration)
- PXRD diffraction data (theta, intensity)

## Special Handling

- "no catalyst" experiments: No catalyst fields in output
- Multiple dash types: Regular dash (-), en-dash (–), em-dash (—) all treated as empty
- Force catalyst mass: Converts volume to mass for specific catalyst types
