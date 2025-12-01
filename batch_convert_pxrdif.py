import os
import sys
from typing import List, Tuple


def find_pairs(dataset_dir: str) -> List[Tuple[str, str]]:
    pairs = []
    for name in sorted(os.listdir(dataset_dir)):
        if not name.lower().endswith('.csv'):
            continue
        if name.lower().endswith('_pxrd.csv'):
            continue
        cond = name
        pxrd = os.path.splitext(name)[0] + '_PXRD.csv'
        if os.path.exists(os.path.join(dataset_dir, pxrd)):
            pairs.append((cond, pxrd))
    return pairs


def main() -> int:
    # Ensure repo root is on sys.path to import Research.excel_to_pxrdif_converter
    this_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(this_dir, '..'))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)

    try:
        from Research.excel_to_pxrdif_converter import excel_to_pxrdif_multiple
    except Exception as e:
        print(f"Failed to import converter: {e}")
        return 1

    datasets = [
        os.path.join(repo_root, 'Research', 'Dataset', 'SO3H-COF'),
        os.path.join(repo_root, 'Research', 'Dataset', 'TABTA-COF'),
    ]

    overall_created: List[str] = []
    for ds in datasets:
        if not os.path.isdir(ds):
            print(f"Dataset directory missing: {ds}")
            continue

        print(f"\n=== Processing dataset: {os.path.relpath(ds, repo_root)} ===")
        pairs = find_pairs(ds)
        if not pairs:
            print("No condition/PXRD CSV pairs found.")
            continue

        for cond, pxrd in pairs:
            print(f"\nConverting pair: {cond}  +  {pxrd}")
            cwd = os.getcwd()
            try:
                os.chdir(ds)  # ensure outputs go inside the dataset folder
                ok, created = excel_to_pxrdif_multiple(cond, pxrd, operator='Yaghi group', verbose=False)
                if ok:
                    # Prepend absolute path for reporting
                    created_abs = [os.path.join(ds, os.path.basename(p)) for p in created]
                    overall_created.extend(created_abs)
                    print(f"Created {len(created)} files under {ds}")
                else:
                    print("Conversion returned failure for this pair.")
            finally:
                os.chdir(cwd)

    print("\n==== Summary ====")
    if overall_created:
        print(f"Total PXRDIF files created: {len(overall_created)}")
        # Print a few sample paths
        for p in overall_created[:10]:
            print(f" - {p}")
        if len(overall_created) > 10:
            print(" ...")
        return 0
    else:
        print("No PXRDIF files were created.")
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
