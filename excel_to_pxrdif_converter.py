import pandas as pd
import numpy as np
from datetime import datetime
import re
import os
import csv
import unicodedata
import builtins
import json


class DebugFilter:
    """
    Context manager to suppress lines that start with 'DEBUG:' when enable=False.
    """
    def __init__(self, enable: bool = True):
        self.enable = enable
        self._orig_print = None

    def __enter__(self):
        self._orig_print = builtins.print

        def filtered_print(*args, **kwargs):
            if not self.enable and args and isinstance(args[0], str) and args[0].startswith('DEBUG:'):
                return
            return self._orig_print(*args, **kwargs)

        builtins.print = filtered_print
        return self

    def __exit__(self, exc_type, exc, tb):
        if self._orig_print is not None:
            builtins.print = self._orig_print
        self._orig_print = None

def find_condition_rows(df):
    """
    Find rows that contain experiment condition data.
    Supports two formats:
      1) Inline strings in a cell starting with A-Z number (e.g., "A1,3.15,4.23,...").
      2) Column-based format with headers like 'Condition ID', 'Tp mass (mg)', etc.
    Returns list of tuples (row_idx, col_idx, value_or_dict). For column-based rows,
    value_or_dict is a dict of the row keyed by column name.
    """
    condition_rows = []

    # Heuristic: detect column-based format
    lower_cols = [str(c).strip().lower() for c in df.columns]
    col_based_keys = ['condition id', 'condition_id']
    has_condition_id_col = any(k in lower_cols for k in col_based_keys)

    if has_condition_id_col:
        # Map headers to normalized keys
        col_map = {}
        for i, c in enumerate(df.columns):
            cl = str(c).strip().lower()
            col_map[i] = cl
        for idx, row in df.iterrows():
            row_dict = {col_map[i]: row[i] for i in range(len(df.columns))}
            cond_id_raw = row_dict.get('condition id') or row_dict.get('condition_id')
            if pd.isna(cond_id_raw):
                continue
            cond_id = str(cond_id_raw).strip()
            # Accept IDs like 'A1' or '6A1' (leading batch digit)
            if re.search(r'[A-Z]\d+', cond_id):
                condition_rows.append((idx, 0, row_dict))
                print(f"DEBUG: Found column-based condition at row {idx}: {cond_id}")
        return condition_rows

    # Fallback: inline-string format
    for idx, row in df.iterrows():
        for col_idx, cell in enumerate(row):
            if pd.isna(cell):
                continue
            cell_str = str(cell).strip()
            # Accept inline strings containing codes like 'A1' or '6A1' anywhere in the text
            if re.search(r'[A-Z]\d+', cell_str):
                condition_rows.append((idx, col_idx, cell_str))
                print(f"DEBUG: Found inline condition at row {idx}, col {col_idx}: {cell_str}")
                break

    return condition_rows

def normalize_str(s):
    """
    Normalize unicode string and replace common problematic characters.
    """
    if not isinstance(s, str):
        return s
    # Replace problematic characters BEFORE NFKC normalization to avoid unwanted conversions
    s = s.replace('�', 'μ')
    s = s.replace('ʵ', 'μ')  # U+02B5 modifier letter small turned r (appears in some CSV files)
    s = s.replace('ɻ', 'μ')  # U+027B latin small letter turned r (result of NFKC on U+02B5)
    # Handle 2-byte UTF-8 sequences read as latin1
    s = s.replace('\xca\xb5', 'μ')  # 0xCA 0xB5 = U+02B5 read as latin1
    s = s.replace('Êµ', 'μ')  # Same as above
    s = s.replace('Ê', 'μ')  # 0xCA alone (first byte of UTF-8 μ sequences)
    s = unicodedata.normalize('NFKC', s)
    s = s.replace('uL', 'μL')
    s = s.replace('μl', 'μL')
    s = s.replace('μL', 'μL')
    s = s.replace('·', '·')
    # Handle non-breaking spaces and additional Unicode variants
    s = s.replace('\u00A0', ' ')  # Replace non-breaking space
    s = s.replace('\u202F', ' ')  # Replace narrow no-break space
    s = re.sub(r'\s+', ' ', s)    # Normalize all whitespace to single space
    return s.strip()

def normalize_catalyst_name(name: str) -> str:
    if not isinstance(name, str):
        return name
    n = name.strip()
    low = n.lower()
    if 'bf3' in low:
        return 'BF3·OEt2'
    # p-TsOH variants
    if low.replace('-', '') in {'ptsoh', 'p-tsoh'} or 'p-toluenesulfonic' in low:
        return 'p-TsOH'
    return n

def normalize_solvent_name(name: str) -> str:
    if not isinstance(name, str):
        return name
    low = name.strip().lower()
    if 'mesitylene' in low:
        return 'Mesitylene'
    if 'butanol' in low:
        return 'n-Butanol'
    if 'dioxane' in low:
        return '1_4-Dioxane'
    if 'dmf' in low or 'dimethylformamide' in low:
        return 'DMF'
    return name.strip()

def parse_catalyst_info(catalyst_part: str) -> dict:
    """
    Parse catalyst info into type, concentration (M), solvent, and volume (µL) with robust handling.
    Supports patterns like:
      - "p-TsOH (0.15 M in dioxane, 30 µL)"
      - "0.5M pTsOH in DMF (30 µL)"
      - "BF3·OEt2 (0.1 M in o-DCB, 20 uL)"
      - Solid catalysts like "ZnCl2 (30 µL)" or just "Zeolite"
    """
    out = {}
    if not catalyst_part:
        return out
    s = normalize_str(str(catalyst_part))
    # Remove leading/trailing quotes that can break regexes and leak into names
    s = re.sub(r'^[\'\"]+|[\'\"]+$', '', s)

    # Volume (µL)
    vm = re.search(r'(\d+\.?\d*)\s*[μu]?L', s, re.IGNORECASE)
    if vm:
        try:
            out['catalyst_volume'] = float(vm.group(1))
        except ValueError:
            pass

    # Mass (mg/g) for solids; convert g -> mg
    mm = re.search(r'(\d+\.?\d*)\s*(mg|g)\b', s, re.IGNORECASE)
    if mm:
        try:
            mass = float(mm.group(1))
            unit = mm.group(2).lower()
            if unit == 'g':
                mass *= 1000.0
            out['catalyst_mass'] = mass
        except ValueError:
            pass

    # Concentration (M)
    # Handle both "3 M" and "2μM" (micromolar) formats
    cm = re.search(r'(\d+\.?\d*)\s*[μu]?\s*M\b', s, re.IGNORECASE)
    if cm:
        try:
            out['catalyst_concentration'] = float(cm.group(1))
        except ValueError:
            pass

    # Solvent: try "M in solvent" first, then generic "in solvent"
    sm = re.search(r'M\s+in\s+([^,\)]+)', s, re.IGNORECASE)
    if not sm:
        sm = re.search(r'\bin\s+([^,\)]+)', s, re.IGNORECASE)
    if sm:
        out['catalyst_solvent'] = normalize_solvent_name(sm.group(1))

    # Type:
    # Case 1: Leading concentration then type: "0.5 M pTsOH in DMF ..."
    tm = re.search(r'^\s*(\d+\.?\d*)\s*M\s*([^\(,]+?)(?:\s+in\b|\(|$)', s, re.IGNORECASE)
    if tm:
        t = tm.group(2).strip()
        out['catalyst_type'] = normalize_catalyst_name(t)
    else:
        # Case 2: Type before parentheses: "p-TsOH (0.15 M in dioxane, 30 µL)"
        tm2 = re.search(r'^\s*([^\(,]+?)\s*(?:\(|,|$)', s)
        if tm2:
            out['catalyst_type'] = normalize_catalyst_name(tm2.group(1).strip())

    # Final normalization for BF3 variants and common names
    if 'catalyst_type' in out:
        # Clean up any stray quotes around the name
        ct = str(out['catalyst_type']).strip().strip('\"').strip("'")
        out['catalyst_type'] = normalize_catalyst_name(ct)

    return out

def parse_condition_string(condition_string, force_catalyst_mass: bool = False):
    """
    Parse the comma-delimited condition string with improved parsing.
    Expected format: A8,3.15,4.23,Mesitylene (100 μL),1,4-Dioxane (200 μL),-,p-Toluenesulfonic acid (0.5 M, 30 μL),330,Dioxane-rich binary
    """
    conditions = {}
    
    if pd.isna(condition_string):
        return conditions
    
    condition_str = str(condition_string).strip()
    condition_str = normalize_str(condition_str)
    print(f"DEBUG: Parsing condition string: {condition_str}")
    
    # Split by commas but preserve content within parentheses and quotes more robustly
    parts = []
    current_part = ""
    paren_count = 0
    in_quotes = False
    
    for i, char in enumerate(condition_str):
        if char == '"':
            in_quotes = not in_quotes
        elif char == '(' and not in_quotes:
            paren_count += 1
        elif char == ')' and not in_quotes:
            paren_count -= 1
        elif char == ',' and paren_count == 0 and not in_quotes and not (current_part.strip() and current_part[-1] == '('):
            parts.append(current_part.strip())
            current_part = ""
            continue
        current_part += char
    
    if current_part:
        current_part = normalize_str(current_part)
        parts.append(current_part.strip())
    
    # Normalize all parts and strip surrounding quotes
    parts = [normalize_str(p).strip().strip('"').strip("'").strip() for p in parts]
    print(f"DEBUG: Parsed condition parts: {parts}")
    
    if len(parts) == 0:
        return conditions
    
    # [0] Extract condition ID (A1, E3, etc.) possibly with leading batch digit like '6A1'.
    # Allow the letter-digit code anywhere in the token so '6E7' -> 'E7'.
    exp_match = re.search(r'([A-Z]\d+)', parts[0])
    if exp_match:
        conditions['condition_id'] = exp_match.group(1)
    
    # [1] Extract TP mass
    if len(parts) >= 2:
        try:
            mass_value = float(parts[1])
            conditions['tp_mass'] = mass_value
            print(f"DEBUG: Found TP mass: {mass_value}")
        except ValueError:
            print(f"DEBUG: Could not parse TP mass from: {parts[1]}")
    
    # [2] Extract second precursor mass (Pa-SO₃H, TABTA, or DABA)
    # Try to detect which one from the original condition_string
    if len(parts) >= 3:
        try:
            mass_value = float(parts[2])
            # Infer precursor type from condition_string header if present
            lower_str = condition_str.lower()
            if 'tabta' in lower_str:
                conditions['tabta_mass'] = mass_value
                print(f"DEBUG: Found TABTA mass: {mass_value}")
            elif 'daba' in lower_str:
                conditions['daba_mass'] = mass_value
                print(f"DEBUG: Found DABA mass: {mass_value}")
            else:
                conditions['pa_so3h_mass'] = mass_value
                print(f"DEBUG: Found Pa-SO₃H mass: {mass_value}")
        except ValueError:
            print(f"DEBUG: Could not parse second precursor mass from: {parts[2]}")
    
    # Optional Temperature column (e.g., after masses). If present, capture and shift positions.
    start_pos = 3
    if len(parts) > start_pos:
        try:
            temp_val = float(parts[start_pos])
            # temperatures are generally small integers like 90, 120, etc.
            # If numeric, treat as temperature in °C
            conditions['reaction_temperature'] = temp_val
            conditions['temperature_c'] = temp_val  # backward-compatible alias
            start_pos += 1
        except ValueError:
            pass

    # [start_pos], [start_pos+1], [start_pos+2] Extract solvents - be more flexible with parsing
    solvent_positions = [start_pos, start_pos + 1, start_pos + 2]
    solvent_index = 1
    solvent_volumes = []
    last_solvent_pos = start_pos - 1  # Track last position where solvent was found
    
    for pos in solvent_positions:
        if len(parts) > pos and parts[pos] != '-':
            part = normalize_str(parts[pos])
            print(f"DEBUG: Processing solvent at position {pos}: {part}")
            
            # Check if this is a combined solvent string (contains semicolon)
            if ';' in part:
                print(f"DEBUG: Detected combined solvent string: {part}")
                # Parse as combined string: "Mesitylene 200; 1,4-Dioxane 100"
                # Pattern: Name (vol µL) OR Name vol µL, allow commas in name
                pattern = re.compile(r"\s*([^();]+?)\s*\((\d+\.?\d*)\s*[μu]?l\)|\s*([^();]+?)\s+(\d+\.?\d*)\s*[μu]?l",
                                      re.IGNORECASE)
                matches = list(pattern.findall(part))
                for match in matches:
                    if match[0]:  # Pattern with parentheses
                        name, vol = match[0].strip(), match[1]
                    else:  # Pattern without parentheses
                        name, vol = match[2].strip(), match[3]
                    # Normalize name
                    if 'mesitylene' in name.lower():
                        name = 'Mesitylene'
                    elif 'butanol' in name.lower():
                        name = 'n-Butanol'
                    elif 'dioxane' in name.lower():
                        name = '1_4-Dioxane'
                    elif 'dmf' in name.lower():
                        name = 'DMF'
                    conditions[f'solvent_{solvent_index}_type'] = name
                    conditions[f'solvent_{solvent_index}_volume'] = float(vol)
                    solvent_volumes.append(float(vol))
                    print(f"DEBUG: Found combined solvent {solvent_index}: {name} - {vol} μL")
                    solvent_index += 1
                last_solvent_pos = pos
                continue
            
            # Single solvent parsing
            solvent_patterns = [
                r'^([^(]+)\s*\((\d+\.?\d*)\s*[μu]?L?\)',  # Standard pattern
                r'^([^(]+)\s*\((\d+\.?\d*)\s*μL\)',       # Explicit μL
                r'^([^(]+)\s*\((\d+\.?\d*)\)',            # No unit specified
                r'^([^(]+?)\s+(\d+\.?\d*)\s*(μL|uL|mL|L)'
            ]
            
            solvent_found = False
            for pattern in solvent_patterns:
                solvent_match = re.search(pattern, part, re.IGNORECASE)
                if solvent_match:
                    solvent_name = solvent_match.group(1).strip()
                    solvent_volume = float(solvent_match.group(2))
                    if 'mesitylene' in solvent_name.lower():
                        solvent_name = 'Mesitylene'
                    elif 'butanol' in solvent_name.lower():
                        solvent_name = 'n-Butanol'
                    elif 'dioxane' in solvent_name.lower():
                        solvent_name = '1_4-Dioxane'
                    elif 'dmf' in solvent_name.lower() or 'dimethylformamide' in solvent_name.lower():
                        solvent_name = 'DMF'
                    
                    conditions[f'solvent_{solvent_index}_type'] = solvent_name
                    conditions[f'solvent_{solvent_index}_volume'] = solvent_volume
                    solvent_volumes.append(solvent_volume)
                    print(f"DEBUG: Found solvent {solvent_index}: {solvent_name} - {solvent_volume} μL")
                    solvent_index += 1
                    last_solvent_pos = pos
                    solvent_found = True
                    break
            
            if not solvent_found:
                print(f"DEBUG: Could not parse solvent from: {part}")
    
    # Find catalyst (next position after last solvent, or skip dashes/en-dashes/em-dashes)
    catalyst_pos = last_solvent_pos + 1
    while catalyst_pos < len(parts) and (parts[catalyst_pos].strip() in ['-', '–', '—', ''] or not parts[catalyst_pos].strip()):
        catalyst_pos += 1
    if catalyst_pos < len(parts):
        catalyst_part = parts[catalyst_pos].strip()
        print(f"DEBUG: Checking catalyst at position {catalyst_pos}: '{catalyst_part}'")
        # Skip if the part says "no catalyst" or if the rationale field says "no catalyst"
        if catalyst_part and catalyst_part not in ['-', '–', '—']:
            # Check if this catalyst field or the rationale explicitly says "no catalyst"
            # Rationale is typically at catalyst_pos + 2, but check both +1 and +2 positions
            has_no_catalyst = 'no catalyst' in catalyst_part.lower()
            for offset in [1, 2]:
                rationale_pos = catalyst_pos + offset
                if not has_no_catalyst and rationale_pos < len(parts):
                    has_no_catalyst = 'no catalyst' in str(parts[rationale_pos]).lower()
            
            if not has_no_catalyst:
                catalyst_part = normalize_str(catalyst_part)
                print(f"DEBUG: Processing catalyst: {catalyst_part}")
                cat = parse_catalyst_info(catalyst_part)
                # If header dictates mass mg, override any volume to be mass
                if force_catalyst_mass:
                    # If the catalyst text looks like a solution (has µL or concentration or 'in solvent'), do not force mass
                    has_ul = re.search(r"[μu]\s*L", catalyst_part, re.IGNORECASE)
                    has_conc = re.search(r"\b\d+\.?\d*\s*M\b", catalyst_part, re.IGNORECASE)
                    has_in = re.search(r"\bin\s+[^,\)]", catalyst_part, re.IGNORECASE)
                    if not (has_ul or has_conc or has_in):
                        # Try explicit mg/g first
                        mm = re.search(r"(\d+\.?\d*)\s*(mg|g)\b", catalyst_part, re.IGNORECASE)
                        forced_mass = None
                        if mm:
                            try:
                                forced_mass = float(mm.group(1))
                                if mm.group(2).lower() == 'g':
                                    forced_mass *= 1000.0
                            except ValueError:
                                forced_mass = None
                        if forced_mass is None:
                            # Fall back to the first bare number as mg
                            mn = re.search(r"\b(\d+\.?\d*)\b", catalyst_part)
                            if mn:
                                try:
                                    forced_mass = float(mn.group(1))
                                except ValueError:
                                    forced_mass = None
                        if forced_mass is not None:
                            cat['catalyst_mass'] = forced_mass
                        # Remove volume if present
                        if 'catalyst_volume' in cat:
                            cat.pop('catalyst_volume', None)
                conditions.update(cat)
    
    # [7 or next] Extract total volume with fallback
    total_vol_pos = catalyst_pos + 1 if catalyst_pos < len(parts) else 7
    if len(parts) > total_vol_pos:
        try:
            total_vol = float(parts[total_vol_pos])
            conditions['total_volume'] = total_vol
            print(f"DEBUG: Found total volume: {total_vol}")
        except ValueError:
            print(f"DEBUG: Could not parse total volume from: {parts[total_vol_pos]}")
            # Fallback: Sum solvent volumes
            if solvent_volumes:
                total_vol_fallback = sum(solvent_volumes)
                if 'catalyst_volume' in conditions:
                    total_vol_fallback += conditions['catalyst_volume']
                conditions['total_volume'] = total_vol_fallback
                print(f"DEBUG: Fallback total volume (sum of solvents + catalyst): {total_vol_fallback}")
    
    # [8 or next] Extract ratio from the last part
    ratio_pos = total_vol_pos + 1 if total_vol_pos < len(parts) else 8
    if len(parts) > ratio_pos:
        ratio_part = parts[ratio_pos]
        ratio_match = re.search(r'((?:Mes|Bu|Diox):[^,\s]*\s+\d+:\d+)', ratio_part)
        if ratio_match:
            conditions['ratio'] = ratio_match.group(1).strip()
            print(f"DEBUG: Found ratio: {conditions['ratio']}")
    
    return conditions

def extract_temperature_from_name(path_or_name: str):
    """Extract temperature in °C from filenames like 'GPT-3_120C.csv' or '..._90C_PXRD.csv'.
    Returns a float or None if not found.
    """
    try:
        name = os.path.basename(str(path_or_name))
    except Exception:
        name = str(path_or_name)
    name = normalize_str(name)
    # Match patterns like '120C', '120 C', '120°C' even if followed by '_' or other non-letters
    m = re.search(r'(\d+)\s*°?\s*[cC](?![A-Za-z])', name)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            return None
    return None

def load_fixed_params(config_path: str = None) -> dict:
    """Load fixed instrument/reaction parameters from JSON file.
    Keys expected:
      - wavelength_A (float)
      - scanning_speed_deg_per_min (float)
      - xray_power_W (float)
      - mode_of_operation (str)
      - reaction_time_hours (float)
    Returns a dict mapped to internal field names with units normalized (kW, hours).
    """
    defaults = {}
    try:
        if config_path is None:
            here = os.path.dirname(__file__)
            config_path = os.path.join(here, 'fixed_params.json')
        if os.path.exists(config_path):
            with open(config_path, 'r', encoding='utf-8') as f:
                cfg = json.load(f)
            # Map to internal fields
            if 'wavelength_A' in cfg:
                defaults['wavelength'] = float(cfg['wavelength_A'])
            if 'scanning_speed_deg_per_min' in cfg:
                defaults['scanning_speed'] = float(cfg['scanning_speed_deg_per_min'])
            if 'xray_power_W' in cfg:
                # Convert W -> kW for PXRDIF
                defaults['xray_power'] = float(cfg['xray_power_W']) / 1000.0
            if 'mode_of_operation' in cfg and cfg['mode_of_operation']:
                defaults['mode_of_operation'] = str(cfg['mode_of_operation'])
            if 'reaction_time_hours' in cfg:
                defaults['reaction_time'] = float(cfg['reaction_time_hours'])
    except Exception:
        # Silently ignore config load errors to avoid breaking conversions
        pass
    return defaults

def parse_condition_row_dict(row_dict):
    """
    Parse a column-based row dict into the conditions structure.
    Expects keys like:
      - 'condition id'
      - 'tp mass (mg)' / 'tp mass'
      - 'pa-so₃h mass (mg)' / 'tabta mass (mg)' / 'pa-so3h mass (mg)'
      - solvent columns such as 'solvents & volumes (µl)' or per-solvent columns
        'solvent 1 (type, volume (µl))', 'solvent 2 ...', 'solvent 3 ...'
      - catalyst info like 'catalyst info' or 'catalyst (type, conc. m, vol. µl)'
      - 'total volume (µl)' / 'total vol. (µl)'
    """
    conditions = {}

    def get_key(*candidates):
        for k in candidates:
            if k in row_dict and not pd.isna(row_dict[k]):
                return row_dict[k]
        return None

    # Normalize row_dict keys to lowercase and strip
    norm = {str(k).strip().lower(): v for k, v in row_dict.items()}

    cond_id = get_key('condition id', 'condition_id')
    if cond_id is not None:
        cond_id = str(cond_id).strip()
        # Accept IDs like 'A1' or '6A1' (leading batch digit). Store normalized letter-digit code.
        m = re.search(r'([A-Z]\d+)', cond_id)
        if m:
            conditions['condition_id'] = m.group(1)

    # Masses
    tp_mass = get_key('tp mass (mg)', 'tp mass')
    if tp_mass is not None:
        try:
            conditions['tp_mass'] = float(str(tp_mass).replace(',', ''))
            conditions['precursor_1_mass'] = conditions['tp_mass']
            conditions['precursor_1_name'] = '1,3,5-Triformylphloroglucinol (Tp)'
        except ValueError:
            pass

    other_mass = get_key('pa-so₃h mass (mg)', 'pa-so3h mass (mg)', 'pa-so₃h mass',
                         'tabta mass (mg)', 'tabta mass',
                         'daba mass (mg)', 'daba mass')
    if other_mass is not None:
        try:
            val = float(str(other_mass).replace(',', ''))
            # Heuristic: choose name by which column exists
            if 'tabta mass' in norm or 'tabta mass (mg)' in norm:
                conditions['precursor_2_name'] = 'Tris(4-aminophenyl)benzene-1,3,5-triscarboxamide (TABTA)'
                conditions['tabta_mass'] = val
            elif 'daba mass' in norm or 'daba mass (mg)' in norm:
                conditions['precursor_2_name'] = '2,5-Diaminobenzenesulfonic acid (Pa-SO₃H)'
                conditions['pa_so3h_mass'] = val
                conditions['daba_mass'] = val  # Store both for compatibility
            else:
                conditions['precursor_2_name'] = '2,5-Diaminobenzenesulfonic acid (Pa-SO₃H)'
                conditions['pa_so3h_mass'] = val
            conditions['precursor_2_mass'] = val
        except ValueError:
            pass

    # Solvents: either a combined column or separate solvent 1/2/3 columns
    solvents = []
    combined = get_key(
        'solvents & volumes (µl)', 'solvents & volumes (�l)',
        'solvents', 'solvent list', 'solvent(s)',
        'solvent list (µl)', 'solvent list (�l)',
        'solvent volumes', 'solvents (type and volume)'
    )
    if combined:
        combined_str = normalize_str(str(combined))
        # Prefer pattern-based global extraction to tolerate commas/semicolons/plus signs
        # Matches: Name (123 uL) OR Name 123 uL
        pattern = re.compile(r"\s*([^(),;+]+?)\s*\((\d+\.?\d*)\s*[μu]?l\)|\s*([^(),;+]+?)\s+(\d+\.?\d*)\s*[μu]?l",
                              re.IGNORECASE)
        matches = list(pattern.findall(combined_str))
        for mm in matches:
            name = mm[0] or mm[2]
            vol = mm[1] or mm[3]
            if name and vol:
                try:
                    solvents.append((name.strip(), float(vol)))
                except ValueError:
                    pass
        # Fallback: split if regex found nothing
        if not solvents:
            for part in re.split(r'[;,+]', combined_str):
                part = part.strip()
                if not part or part == '-':
                    continue
                m2 = re.search(r'^([^\(]+)\((\d+\.?\d*)\s*[μu]?l\)', part, re.IGNORECASE)
                if not m2:
                    m2 = re.search(r'^([^\d]+?)\s+(\d+\.?\d*)\s*[μu]?l', part, re.IGNORECASE)
                if m2:
                    name = m2.group(1).strip()
                    try:
                        vol = float(m2.group(2))
                        solvents.append((name, vol))
                    except ValueError:
                        continue
    else:
        # Separate solvent columns - be flexible about header wording like 'Solvent 1 (type vol µL)'
        for i in [1, 2, 3]:
            # Try known keys first
            val = get_key(
                f'solvent {i} (type, volume (µl))', f'solvent {i} (type, volume (�l))',
                f'solvent {i} (vol)', f'solvent {i} (type vol µl)', f'solvent {i} (type vol �l)'
            )
            # If not found, search any column key that contains 'solvent {i}'
            if val is None:
                for k in norm.keys():
                    if re.search(rf'\bsolvent\s*{i}\b', str(k), re.IGNORECASE):
                        candidate = norm.get(k)
                        if candidate and str(candidate).strip() not in ['-']:
                            val = candidate
                            break
            if not val or str(val).strip() == '-':
                continue
            s = normalize_str(str(val))
            # Accept patterns: Name (123 µL), Name 123 µL, Name 123
            m = re.match(r'^([^\(]+)\((\d+\.?\d*)\s*([μu]?l|ml|l)?\)', s, re.IGNORECASE)
            if not m:
                m = re.match(r'^([^,]+?)\s+(\d+\.?\d*)\s*([μu]?l|ml|l)?\b', s, re.IGNORECASE)
            if m:
                name = m.group(1).strip()
                vol = float(m.group(2))
                unit = (m.group(3) or '').lower()
                # Normalize units to µL
                if unit == 'ml':
                    vol *= 1000.0
                elif unit == 'l':
                    vol *= 1_000_000.0
                # normalize common names
                lname = name.lower()
                if 'mesitylene' in lname:
                    name = 'Mesitylene'
                elif 'butanol' in lname:
                    name = 'n-Butanol'
                elif 'dioxane' in lname:
                    name = '1_4-Dioxane'
                elif 'dmf' in lname or 'n_n-dimethylformamide' in lname or 'dimethylformamide' in lname:
                    name = 'DMF'
                solvents.append((name, vol))

    for idx, (name, vol) in enumerate(solvents, start=1):
        # normalize common names
        lname = name.lower()
        if 'mesitylene' in lname:
            name = 'Mesitylene'
        elif 'butanol' in lname:
            name = 'n-Butanol'
        elif 'dioxane' in lname:
            name = '1_4-Dioxane'
        elif 'dmf' in lname or 'n_n-dimethylformamide' in lname or 'dimethylformamide' in lname:
            name = 'DMF'
        conditions[f'solvent_{idx}_type'] = name
        conditions[f'solvent_{idx}_volume'] = vol

    # Catalyst (header-aware parsing to distinguish mass mg vs volume µL)
    catalyst_value = None
    catalyst_key_used = None
    # Prefer explicit known keys
    catalyst_keys_priority = [
        'catalyst info',
        'catalyst (type, conc. m, vol. µl)', 'catalyst (type, conc. m, vol. �l)', 'catalyst (conc vol)',
        'catalyst (type conc. m mass mg)', 'catalyst (type conc. m mass (mg))', 'catalyst (type conc m mass mg)'
    ]
    for ck in catalyst_keys_priority:
        if ck in norm and norm[ck] and str(norm[ck]).strip().lower() not in ['-', 'none']:
            catalyst_value = norm[ck]
            catalyst_key_used = ck
            break
    # If still not found, pick any column containing 'catalyst' with a non-empty value
    if catalyst_value is None:
        for k, v in norm.items():
            if 'catalyst' in str(k) and v and str(v).strip().lower() not in ['-', 'none']:
                catalyst_value = v
                catalyst_key_used = str(k)
                break
    if catalyst_value is not None and str(catalyst_value).strip().lower() not in ['none', '-']:
        cap = normalize_str(str(catalyst_value))
        # Check if the value explicitly indicates no catalyst
        if 'no catalyst' not in cap.lower():
            # First, parse generally to extract type/conc/solvent plus any quantities present
            cat = parse_catalyst_info(cap)
            # Header-driven override: if header indicates mass mg, force quantity to mass (mg)
            used_key_lc = (catalyst_key_used or '').lower()
            mass_header = bool(re.search(r"mass\s*\(?\s*mg\)?", used_key_lc))
            if mass_header:
                # If the text clearly indicates a solution catalyst (has concentration/solvent or µL), keep volume semantics
                has_ul = re.search(r"[μu]\s*L", cap, re.IGNORECASE)
                has_conc = re.search(r"\b\d+\.?\d*\s*M\b", cap, re.IGNORECASE)
                has_in = re.search(r"\bin\s+[^,\)]", cap, re.IGNORECASE)
                if not (has_ul or has_conc or has_in):
                    # Try to extract mg/g explicitly from the cell first
                    mm = re.search(r"(\d+\.?\d*)\s*(mg|g)\b", cap, re.IGNORECASE)
                    forced_mass = None
                    if mm:
                        try:
                            forced_mass = float(mm.group(1))
                            if mm.group(2).lower() == 'g':
                                forced_mass *= 1000.0
                        except ValueError:
                            forced_mass = None
                    # If no mg/g unit present, take the first standalone number as mg
                    if forced_mass is None:
                        # Prefer the last number (often at the end) but ensure it's standalone
                        nums = re.findall(r"(?<![A-Za-z])([0-9]+\.?[0-9]*)(?![A-Za-z])", cap)
                        if nums:
                            try:
                                forced_mass = float(nums[-1])
                            except ValueError:
                                forced_mass = None
                    if forced_mass is not None:
                        cat['catalyst_mass'] = forced_mass
                    # Suppress volume when interpreting as mass
                    if 'catalyst_volume' in cat:
                        cat.pop('catalyst_volume', None)
            conditions.update(cat)

    # Total volume
    tot = get_key('total volume (µl)', 'total vol. (µl)', 'total volume (�l)', 'total vol. (�l)')
    if tot is not None:
        try:
            conditions['total_volume'] = float(str(tot).replace(',', ''))
        except ValueError:
            pass

    # Ratio/rationale column sometimes contains "Mes:Bu 2:1" etc.
    rationale = get_key('notes', 'rationale')
    if rationale:
        txt = normalize_str(str(rationale))
        ratio_match = re.search(r'(Mes|Bu|Diox|Hex|o-DCB|Benz)[^\d]*\d+\s*:\s*\d+', txt, re.IGNORECASE)
        if ratio_match:
            conditions['ratio'] = ratio_match.group(0)

    return conditions

def find_parameter_columns(df, header_row_idx=None):
    """
    Search for additional parameter columns in the dataframe.
    """
    column_mappings = {}
    
    # Define keyword mappings for parameters, hardcoded for specific columns
    keyword_mappings = {
        'wavelength': ['wavelength', 'wave length', 'λ', 'lambda', 'Wavelength', 'Wave Length', 'WAVELENGTH'],
        'scanning_speed': ['scanning speed', 'scan speed', 'speed', 'deg/min', 'Scanning Speed', 'Scan Speed', 'SCANNING SPEED', 'scanning_speed', 'scan_speed'],
        'xray_power': ['x-ray power', 'xray power', 'power', 'x-ray intensity', 'X-ray Power', 'XRay Power', 'X-RAY POWER', 'xray_power', 'x_ray_power'],
        'mode_of_operation': ['mode', 'operation mode', 'reflection', 'transmission', 'Mode', 'Operation Mode', 'MODE OF OPERATION', 'mode_of_operation'],
        'pxrd_model': ['model', 'instrument', 'pxrd model', 'diffractometer', 'Model', 'Instrument', 'PXRD Model', 'Diffractometer', 'pxrd_model'],
        'reaction_time': ['time', 'reaction time', 'hours', 'duration', 'Time', 'Reaction Time', 'REACTION TIME', 'reaction_time', 'rxn_time', 'rxn time'],
        'reaction_temperature': ['temperature', 'temp', '°c', 'reaction temp', 'Temperature', 'Reaction Temperature', 'REACTION TEMPERATURE', 'reaction_temperature', 'rxn_temp', 'rxn temperature'],
        'catalyst_solvent': ['catalyst solvent', 'cat solvent', 'catalyst medium', 'Catalyst Solvent', 'Cat Solvent', 'CATALYST SOLVENT', 'catalyst_solvent'],
        # Hardcoded precursor keywords for specific columns
        'precursor_1_name': ['tp mass', 'Tp mass', 'TP MASS', 'tp_mass'],
        'precursor_2_name': ['pa-so₃h mass', 'Pa-SO₃H mass', 'PA-SO3H MASS', 'pa_so3h_mass', 'pa-so3h_mass'],
        'precursor_1_mass': ['tp mass', 'Tp mass', 'TP MASS', 'tp_mass'],
        'precursor_2_mass': ['pa-so₃h mass', 'Pa-SO₃H mass', 'PA-SO3H MASS', 'pa_so3h_mass', 'pa-so3h_mass']
    }
    
    # Search through all columns
    for col_idx, col in enumerate(df.columns):
        col_str = str(col).lower().strip()
        for param_type, keywords in keyword_mappings.items():
            if any(keyword in col_str for keyword in keywords):
                column_mappings[param_type] = col_idx
                print(f"DEBUG: Found parameter column {col_idx} ('{col}') for {param_type}")
                break
    
    return column_mappings

def parse_experiment_conditions(df):
    """
    Parse experiment conditions using the improved approach.
    Returns a list of conditions dicts in order of appearance.
    """
    all_conditions = []
    condition_rows = find_condition_rows(df)
    if not condition_rows:
        print("DEBUG: No condition rows found")
        return all_conditions

    header_row_idx = None
    if condition_rows:
        first_condition_row = condition_rows[0][0]
        if first_condition_row > 0:
            header_row_idx = first_condition_row - 1

    column_mappings = find_parameter_columns(df, header_row_idx)

    # Determine if this CSV encodes conditions as one quoted CSV string per row.
    # Also detect which precursor type from header
    header_force_mass = False
    is_tabta_dataset = False
    is_daba_dataset = False
    try:
        header_row = df.columns[0]
        hdr = normalize_str(str(header_row)).lower()
        if 'catalyst' in hdr and re.search(r'mass\s*\(?\s*mg\)?', hdr):
            header_force_mass = True
        if 'tabta' in hdr:
            is_tabta_dataset = True
        elif 'daba' in hdr:
            is_daba_dataset = True
    except Exception:
        header_force_mass = False

    for row_idx, col_idx, payload in condition_rows:
        # payload can be inline string or a row dict depending on detection
        if isinstance(payload, dict):
            conditions = parse_condition_row_dict(payload)
        else:
            conditions = parse_condition_string(payload, force_catalyst_mass=header_force_mass)
            # For inline parsing, set the correct precursor mass key based on dataset type
            if is_tabta_dataset and 'pa_so3h_mass' in conditions:
                conditions['tabta_mass'] = conditions.pop('pa_so3h_mass')
            elif is_daba_dataset and 'pa_so3h_mass' in conditions:
                conditions['daba_mass'] = conditions.pop('pa_so3h_mass')
        if not conditions:
            continue
        # Ensure defaults if not provided
        if 'precursor_1_name' not in conditions:
            conditions['precursor_1_name'] = '1,3,5-Triformylphloroglucinol (Tp)'
        if 'precursor_1_mass' not in conditions and 'tp_mass' in conditions:
            conditions['precursor_1_mass'] = conditions.get('tp_mass')
        if 'precursor_2_name' not in conditions:
            # Infer precursor name from which mass field is present
            if 'tabta_mass' in conditions:
                conditions['precursor_2_name'] = 'Tris(4-aminophenyl)benzene-1,3,5-triscarboxamide (TABTA)'
            elif 'daba_mass' in conditions or 'pa_so3h_mass' in conditions:
                conditions['precursor_2_name'] = '2,5-Diaminobenzenesulfonic acid (Pa-SO₃H)'
            else:
                conditions['precursor_2_name'] = '2,5-Diaminobenzenesulfonic acid (Pa-SO₃H)'  # ultimate fallback
        if 'precursor_2_mass' not in conditions:
            conditions['precursor_2_mass'] = conditions.get('tabta_mass') or conditions.get('pa_so3h_mass') or conditions.get('daba_mass')
        all_conditions.append(conditions)
    return all_conditions

def parse_xrd_data_by_order(df, num_conditions):
    """
    Parse PXRD data columns by order, mapping last condition to first two columns, first condition to last two columns.
    Returns a list of dicts with theta/intensity arrays in the order of conditions.
    """
    xrd_data_list = []
    num_cols = len(df.columns)
    # Each experiment has two columns: 2θ and Intensity
    num_experiments = num_cols // 2  # Assuming even number of columns, all data

    # Limit to available experiments
    effective_num = min(num_conditions, num_experiments)

    for i in range(effective_num):
        # Map: condition 0 (first) -> last pair (columns num_cols-2, num_cols-1)
        # condition num_conditions-1 (last) -> first pair (columns 0,1)
        angle_col = num_cols - 2 - (i * 2)
        intensity_col = num_cols - 1 - (i * 2)
        if angle_col < 0 or intensity_col < 0:
            break
        theta_values = []
        intensity_values = []
        for row_idx in range(1, len(df)):  # Skip header
            try:
                theta_val = df.iloc[row_idx, angle_col]
                intensity_val = df.iloc[row_idx, intensity_col]
                if pd.isna(theta_val) or pd.isna(intensity_val):
                    continue
                theta_float = float(theta_val)
                intensity_float = float(intensity_val)
                if not (np.isnan(theta_float) or np.isnan(intensity_float)):
                    theta_values.append(theta_float)
                    intensity_values.append(intensity_float)
            except (ValueError, TypeError):
                continue
        if theta_values and intensity_values:
            min_length = min(len(theta_values), len(intensity_values))
            xrd_data_list.append({
                'theta': theta_values[:min_length],
                'intensity': intensity_values[:min_length]
            })
            print(f"DEBUG: Found {min_length} XRD data points for condition index {i}")
    return xrd_data_list

def create_pxrdif_content_single_experiment(conditions, xrd_data, operator="Yaghi group", suppress_nulls: bool = False):
    """
    Create PXRDIF format content for a single experiment.
    """
    pxrdif_lines = []
    # Header
    pxrdif_lines.append("data_pxrdif")
    pxrdif_lines.append("")
    
    # Experimental information
    pxrdif_lines.append(f"_exptl_operator '{operator}'")
    
    # Default to 'Rigaku Smartlab' if not specified in conditions
    pxrd_model = conditions.get('pxrd_model', 'Rigaku Smartlab')
    if pxrd_model:
        pxrdif_lines.append(f"_exptl_pxrd_model '{pxrd_model}'")
    
    # Compute experiment number string from normalized condition id (e.g., 'A1').
    condition_id_raw = conditions.get('condition_id', 'Unknown')
    m = re.search(r'([A-Z]\d+)', str(condition_id_raw))
    exp_num = m.group(1) if m else str(condition_id_raw)
    pxrdif_lines.append(f"_exptl_experiment_number '{exp_num}'")
    
    if conditions.get('condition_id'):
        pxrdif_lines.append(f"_exptl_condition_id '{exp_num}'")
    
    pxrdif_lines.append("")
    
    # XRD measurement parameters
    if conditions.get('wavelength') is not None:
        pxrdif_lines.append(f"_diffrn_radiation_wavelength {conditions['wavelength']}")
    else:
        if not suppress_nulls:
            pxrdif_lines.append("_diffrn_radiation_wavelength null")
    
    if conditions.get('scanning_speed') is not None:
        pxrdif_lines.append(f"_diffrn_measurement_scanning_speed {conditions['scanning_speed']}")
    else:
        if not suppress_nulls:
            pxrdif_lines.append("_diffrn_measurement_scanning_speed null")
    
    if conditions.get('xray_power') is not None:
        pxrdif_lines.append(f"_diffrn_radiation_power {conditions['xray_power']}")
    else:
        if not suppress_nulls:
            pxrdif_lines.append("_diffrn_radiation_power null")
    
    if conditions.get('mode_of_operation'):
        pxrdif_lines.append(f"_diffrn_measurement_mode '{conditions['mode_of_operation']}'")
    else:
        if not suppress_nulls:
            pxrdif_lines.append("_diffrn_measurement_mode null")
    
    pxrdif_lines.append("")
    
    # Reaction conditions
    if conditions.get('reaction_time') is not None:
        pxrdif_lines.append(f"_reaction_time {conditions['reaction_time']}")
    else:
        if not suppress_nulls:
            pxrdif_lines.append("_reaction_time null")
    
    if conditions.get('reaction_temperature') is not None:
        pxrdif_lines.append(f"_reaction_temperature {conditions['reaction_temperature']}")
    else:
        if not suppress_nulls:
            pxrdif_lines.append("_reaction_temperature null")
    
    pxrdif_lines.append("")

    # Precursor information data (used for units and optional loop)
    # Determine precursor_2 name, CAS, and mass key
    # If not explicitly set, infer from which mass field is present
    precursor_2_name = conditions.get('precursor_2_name')
    if not precursor_2_name:
        if conditions.get('tabta_mass') is not None:
            precursor_2_name = 'Tris(4-aminophenyl)benzene-1,3,5-triscarboxamide (TABTA)'
        elif conditions.get('daba_mass') is not None or conditions.get('pa_so3h_mass') is not None:
            precursor_2_name = '2,5-Diaminobenzenesulfonic acid (Pa-SO₃H)'
        else:
            precursor_2_name = '2,5-Diaminobenzenesulfonic acid (Pa-SO₃H)'
    
    # CAS number lookup
    cas_numbers = {
        '1,3,5-Triformylphloroglucinol (Tp)': '34374-88-4',
        '2,5-Diaminobenzenesulfonic acid (Pa-SO₃H)': '88-45-9',
        'Tris(4-aminophenyl)benzene-1,3,5-triscarboxamide (TABTA)': '205653-12-9'
    }
    
    # Mass key lookup based on precursor name
    mass_keys = {
        '2,5-Diaminobenzenesulfonic acid (Pa-SO₃H)': 'pa_so3h_mass',
        'Tris(4-aminophenyl)benzene-1,3,5-triscarboxamide (TABTA)': 'tabta_mass'
    }
    
    # DABA is Pa-SO₃H, so fallback to checking daba_mass if pa_so3h_mass not found
    if precursor_2_name == '2,5-Diaminobenzenesulfonic acid (Pa-SO₃H)':
        precursor_2_mass_key = 'pa_so3h_mass' if 'pa_so3h_mass' in conditions else 'daba_mass'
    
    precursor_2_cas = cas_numbers.get(precursor_2_name, 'N/A')
    if precursor_2_name not in mass_keys:
        precursor_2_mass_key = 'tabta_mass'
    
    precursor_data = {
        'precursor_1': {
            'name': conditions.get('precursor_1_name', '1,3,5-Triformylphloroglucinol (Tp)'),
            'cas': cas_numbers.get(conditions.get('precursor_1_name', '1,3,5-Triformylphloroglucinol (Tp)'), '34374-88-4'),
            'target_mass': conditions.get('tp_mass'),
            'actual_mass': conditions.get('tp_mass')
        },
        'precursor_2': {
            'name': precursor_2_name,
            'cas': precursor_2_cas,
            'target_mass': conditions.get('tabta_mass') or conditions.get('daba_mass') or conditions.get('pa_so3h_mass'),
            'actual_mass': conditions.get('tabta_mass') or conditions.get('daba_mass') or conditions.get('pa_so3h_mass')
        }
    }

    # Units
    if not suppress_nulls:
        pxrdif_lines.append("_units_angle 'deg'")
        pxrdif_lines.append("_units_intensity 'cps'")
        pxrdif_lines.append("_units_mass 'mg'")
        pxrdif_lines.append("_units_volume 'μL'")
        pxrdif_lines.append("_units_concentration 'M'")
        pxrdif_lines.append("_units_wavelength 'Å'")
        pxrdif_lines.append("_units_scanning_speed 'deg/min'")
        pxrdif_lines.append("_units_power 'kW'")
        pxrdif_lines.append("_units_time 'hours'")
        pxrdif_lines.append("_units_temperature '°C'")
        pxrdif_lines.append("")
    else:
        used_units = set()
        # Always include angle/intensity if data present
        if xrd_data and xrd_data.get('theta') and xrd_data.get('intensity'):
            used_units.update({'angle', 'intensity'})
        # Instrument/reaction
        if conditions.get('wavelength') is not None:
            used_units.add('wavelength')
        if conditions.get('scanning_speed') is not None:
            used_units.add('scanning_speed')
        if conditions.get('xray_power') is not None:
            used_units.add('power')
        if conditions.get('reaction_time') is not None:
            used_units.add('time')
        if conditions.get('reaction_temperature') is not None:
            used_units.add('temperature')

        # Precursors
        any_precursor_mass = any(
            (v.get('target_mass') is not None) or (v.get('actual_mass') is not None)
            for v in precursor_data.values()
        )
        if any_precursor_mass:
            used_units.add('mass')

        # Solvents
        for i in range(1, 6):
            if conditions.get(f'solvent_{i}_type') and conditions.get(f'solvent_{i}_volume') is not None:
                used_units.add('volume')

        # Catalyst
        if conditions.get('catalyst_volume') is not None:
            used_units.add('volume')
        if conditions.get('catalyst_mass') is not None:
            used_units.add('mass')
        if conditions.get('catalyst_concentration') is not None:
            used_units.add('concentration')

        # Emit only used units in canonical order
        order = [
            ('angle', "_units_angle 'deg'"),
            ('intensity', "_units_intensity 'cps'"),
            ('mass', "_units_mass 'mg'"),
            ('volume', "_units_volume 'μL'"),
            ('concentration', "_units_concentration 'M'"),
            ('wavelength', "_units_wavelength 'Å'"),
            ('scanning_speed', "_units_scanning_speed 'deg/min'"),
            ('power', "_units_power 'kW'"),
            ('time', "_units_time 'hours'"),
            ('temperature', "_units_temperature '°C'"),
        ]
        for key, line in order:
            if key in used_units:
                pxrdif_lines.append(line)
        pxrdif_lines.append("")
    
    # Precursor information loop - using parsed masses from condition string
    # Decide whether to include precursor loop
    include_precursors = True
    if suppress_nulls:
        any_mass_present = any(
            (v.get('target_mass') is not None) or (v.get('actual_mass') is not None)
            for v in precursor_data.values()
        )
        if not any_mass_present:
            include_precursors = False

    if include_precursors:
        pxrdif_lines.append("loop_")
        pxrdif_lines.append("_precursor_name")
        pxrdif_lines.append("_precursor_cas_number")
        pxrdif_lines.append("_precursor_target_mass")
        pxrdif_lines.append("_precursor_actual_mass")

        for precursor_id, data in precursor_data.items():
            name = data.get('name', 'Unknown')
            cas = data.get('cas', 'N/A')
            target_mass = data.get('target_mass')
            actual_mass = data.get('actual_mass')

            # Format masses - use actual values or 'null' if not available
            target_str = str(target_mass) if target_mass is not None else 'null'
            actual_str = str(actual_mass) if actual_mass is not None else 'null'

            pxrdif_lines.append(f"'{name}' '{cas}' {target_str} {actual_str}")

        pxrdif_lines.append("")
    
    # Solvent information loop
    solvents_found = []
    for i in range(1, 6):
        if conditions.get(f'solvent_{i}_type'):
            solvents_found.append(i)
    
    if solvents_found:
        pxrdif_lines.append("loop_")
        pxrdif_lines.append("_solvent_name")
        pxrdif_lines.append("_solvent_volume")
        
        for i in solvents_found:
            solvent_type = conditions[f'solvent_{i}_type']
            solvent_volume = conditions.get(f'solvent_{i}_volume')
            volume_str = str(solvent_volume) if solvent_volume is not None else 'null'
            pxrdif_lines.append(f"'{solvent_type}' {volume_str}")
        
        pxrdif_lines.append("")
    
    # Catalyst information
    if conditions.get('catalyst_type'):
        pxrdif_lines.append(f"_catalyst_type '{conditions['catalyst_type']}'")
        
        if conditions.get('catalyst_concentration') is not None:
            pxrdif_lines.append(f"_catalyst_concentration {conditions['catalyst_concentration']}")
        
        if conditions.get('catalyst_solvent'):
            pxrdif_lines.append(f"_catalyst_solvent '{conditions['catalyst_solvent']}'")
        
        if conditions.get('catalyst_volume') is not None:
            pxrdif_lines.append(f"_catalyst_volume {conditions['catalyst_volume']}")
        
        # For solid catalysts, include mass when available (mg)
        if conditions.get('catalyst_mass') is not None:
            pxrdif_lines.append(f"_catalyst_mass {conditions['catalyst_mass']}")
        
        pxrdif_lines.append("")
    
    # Ratio information
    if conditions.get('ratio'):
        pxrdif_lines.append(f"_sample_ratio '{conditions['ratio']}'")
        pxrdif_lines.append("")
    
    # XRD Data
    if xrd_data:
        pxrdif_lines.append("loop_")
        pxrdif_lines.append("_diffrn_reflns_theta")
        pxrdif_lines.append("_diffrn_reflns_intensity")
        for i in range(len(xrd_data['theta'])):
            pxrdif_lines.append(f"{xrd_data['theta'][i]} {xrd_data['intensity'][i]}")
        pxrdif_lines.append("")
    return "\n".join(pxrdif_lines)

def excel_to_pxrdif_multiple(condition_file="51-55Conditions.csv", pxrd_file="51-55Conditions_PXRD.csv", operator="Unknown", verbose: bool = True):
    """
    Main function to convert two CSV files (conditions and PXRD) to multiple PXRDIF files.
    Uses order-based mapping.
    """
    try:
        print(f"Reading condition CSV file: {condition_file}")
        encodings = ['utf-8', 'latin1', 'iso-8859-1', 'cp1252']
        df_cond = None
        for encoding in encodings:
            try:
                df_cond = pd.read_csv(condition_file, encoding=encoding)
                print(f"Successfully read conditions file with {encoding} encoding")
                break
            except UnicodeDecodeError:
                continue
        if df_cond is None:
            raise ValueError("Could not read condition CSV file with any encoding")

        print(f"Reading PXRD CSV file: {pxrd_file}")
        df_pxrd = None
        for encoding in encodings:
            try:
                df_pxrd = pd.read_csv(pxrd_file, encoding=encoding)
                print(f"Successfully read PXRD file with {encoding} encoding")
                break
            except UnicodeDecodeError:
                continue
        if df_pxrd is None:
            raise ValueError("Could not read PXRD CSV file with any encoding")

        print(f"Condition file shape: {df_cond.shape}, PXRD file shape: {df_pxrd.shape}")

    # Parse experiment conditions from condition file (ordered list)
        print("Parsing experiment conditions...")
        with DebugFilter(enable=verbose):
            all_conditions = parse_experiment_conditions(df_cond)
        print(f"Parsed {len(all_conditions)} conditions")

        # Ensure reaction_temperature is set: use parsed value if present, else infer from filename
        file_temp = extract_temperature_from_name(condition_file) or extract_temperature_from_name(pxrd_file)
        if file_temp is not None:
            for cond in all_conditions:
                # Map alias if used earlier
                if 'reaction_temperature' not in cond and 'temperature_c' in cond:
                    cond['reaction_temperature'] = cond.get('temperature_c')
                if cond.get('reaction_temperature') is None:
                    cond['reaction_temperature'] = file_temp

        # Apply fixed defaults if provided via config
        fixed_defaults = load_fixed_params()
        if fixed_defaults:
            for cond in all_conditions:
                for k, v in fixed_defaults.items():
                    if cond.get(k) is None:
                        cond[k] = v

        # Parse PXRD data from PXRD file (ordered list)
        print("Parsing PXRD data...")
        with DebugFilter(enable=verbose):
            xrd_data_list = parse_xrd_data_by_order(df_pxrd, len(all_conditions))
        print(f"Found {len(xrd_data_list)} PXRD data sets")

        # Create folder named after condition file
        condition_folder = os.path.splitext(os.path.basename(condition_file))[0]
        output_dir = condition_folder
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")

        # Create PXRDIF files in the folder
        created_files = []
        for i, conditions in enumerate(all_conditions):
            # Compute exp_num for filename from condition_id.
            # If condition_id like '6A1', normalize to 'A1'. If missing, use sequential index.
            raw_id = conditions.get('condition_id')
            exp_num = None
            if isinstance(raw_id, str):
                m = re.search(r'([A-Z]\d+)', raw_id)
                if m:
                    exp_num = m.group(1)
            if not exp_num:
                exp_num = f"Unknown_{i+1}"
            if i < len(xrd_data_list):
                xrd_data = xrd_data_list[i]
            else:
                xrd_data = None
            print(f"Creating PXRDIF file for experiment {exp_num}...")
            pxrdif_content = create_pxrdif_content_single_experiment(conditions, xrd_data, operator)
            output_file = os.path.join(output_dir, f"experiment_{exp_num}.pxrdif")
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(pxrdif_content)
            created_files.append(output_file)
            print(f"Created {output_file}")

        print(f"\nSuccessfully created {len(created_files)} PXRDIF files:")
        for file in created_files:
            print(f"  - {file}")

        return True, created_files

    except FileNotFoundError as e:
        print(f"Error: {e}")
        return False, []
    except Exception as e:
        print(f"Error during conversion: {str(e)}")
        import traceback
        traceback.print_exc()
        return False, []

def preview_pxrdif_multiple(condition_file="51-55Conditions.csv", pxrd_file="51-55Conditions_PXRD.csv", num_lines=50, verbose: bool = True):
    """
    Preview PXRDIF output for each experiment using two separate files.
    Uses order-based mapping.
    """
    try:
        encodings = ['utf-8', 'latin1', 'iso-8859-1', 'cp1252']
        df_cond = None
        for encoding in encodings:
            try:
                df_cond = pd.read_csv(condition_file, encoding=encoding)
                print(f"Successfully read conditions file with {encoding} encoding")
                break
            except UnicodeDecodeError:
                continue
        if df_cond is None:
            raise ValueError("Could not read condition CSV file with any encoding")

        df_pxrd = None
        for encoding in encodings:
            try:
                df_pxrd = pd.read_csv(pxrd_file, encoding=encoding)
                print(f"Successfully read PXRD file with {encoding} encoding")
                break
            except UnicodeDecodeError:
                continue
        if df_pxrd is None:
            raise ValueError("Could not read PXRD CSV file with any encoding")

        with DebugFilter(enable=verbose):
            all_conditions = parse_experiment_conditions(df_cond)
            # Populate reaction temperature from filenames if missing
            file_temp = extract_temperature_from_name(condition_file) or extract_temperature_from_name(pxrd_file)
            if file_temp is not None:
                for cond in all_conditions:
                    if 'reaction_temperature' not in cond and 'temperature_c' in cond:
                        cond['reaction_temperature'] = cond.get('temperature_c')
                    if cond.get('reaction_temperature') is None:
                        cond['reaction_temperature'] = file_temp
            # Apply fixed defaults
            fixed_defaults = load_fixed_params()
            if fixed_defaults:
                for cond in all_conditions:
                    for k, v in fixed_defaults.items():
                        if cond.get(k) is None:
                            cond[k] = v
            xrd_data_list = parse_xrd_data_by_order(df_pxrd, len(all_conditions))

        for i, conditions in enumerate(all_conditions):
            # Compute exp_num
            condition_id = conditions.get('condition_id', f"Unknown_{i+1}")
            if condition_id.startswith('E'):
                try:
                    num = int(condition_id[1:])
                    exp_num = str(num + 48)
                except ValueError:
                    exp_num = condition_id
            else:
                exp_num = condition_id
            if i < len(xrd_data_list):
                xrd_data = xrd_data_list[i]
            else:
                xrd_data = None
            pxrdif_content = create_pxrdif_content_single_experiment(conditions, xrd_data)
            lines = pxrdif_content.split('\n')
            preview_lines = lines[:num_lines]
            print("=" * 60)
            print(f"PXRDIF PREVIEW for Experiment {exp_num} (first {num_lines} lines):")
            print("=" * 60)
            for line in preview_lines:
                print(line)
            print("=" * 60)
            print(f"Total lines in this PXRDIF: {len(lines)}")
            print()

    except Exception as e:
        print(f"Error during preview: {str(e)}")
        import traceback
        traceback.print_exc()

def pxrd_splitter(pxrd_filepath, output_dir):
    import os
    import csv
    with open(pxrd_filepath, newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        rows = list(reader)

    # First row: sample names, second row: column headers
    sample_names = rows[0]
    col_headers = rows[1]
    data_rows = rows[2:]

    # Each experiment has two columns: 2θ and Intensity
    num_experiments = len(col_headers) // 2

    base_name = os.path.splitext(os.path.basename(pxrd_filepath))[0]

    for i in range(num_experiments):
        # Get last two digits from sample name (e.g., "96" from "250422_001_96/96")
        sample = sample_names[2*i]
        exp_num = sample.split('/')[-1][-2:]
        filename = f"{base_name}exp_{exp_num}.csv"
        filepath = os.path.join(output_dir, filename)

        # Prepare data for this experiment
        exp_data = []
        for row in data_rows:
            theta = row[2*i]
            intensity = row[2*i+1]
            exp_data.append([theta, intensity])

        # Write to file
        with open(filepath, 'w', newline='', encoding='utf-8') as out_f:
            writer = csv.writer(out_f)
            writer.writerow(['2theta', 'Intensity'])
            writer.writerows(exp_data)

if __name__ == "__main__":
    condition_filename = "SO3H-COF\\GPT-1_120C.csv"
    pxrd_filename = "SO3H-COF\\GPT-1_120C_PXRD.csv"
    operator_name = "Researcher"

    if not os.path.exists(condition_filename):
        print(f"Condition file '{condition_filename}' not found in current directory.")
    elif not os.path.exists(pxrd_filename):
        print(f"PXRD file '{pxrd_filename}' not found in current directory.")
    else:
        print("Generating preview for all experiments...")
        preview_pxrdif_multiple(condition_filename, pxrd_filename)
        proceed = input("\nDo you want to save these as separate PXRDIF files? (y/n): ")
        if proceed.lower().startswith('y'):
            success, created_files = excel_to_pxrdif_multiple(condition_filename, pxrd_filename, operator_name)
            if success:
                print(f"\nConversion complete! Created {len(created_files)} PXRDIF files:")
                for file in created_files:
                    print(f"  - {file}")
        else:
            print("Conversion cancelled.")
