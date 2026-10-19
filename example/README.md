# Worked example: E. coli DPS (EMD-17992)

## Why this map

- Deposited with a voxel size of **0.646 Å** (nominal).
- The same data were calibrated on the gold foil lattice: **0.6484 Å** (Dickerson et al. 2024, Table 1).
- Eight X-ray structures of E. coli DPS are available (the ones used in the paper).
- The scan should therefore find ~0.648 Å rather than 0.646 Å.

## Run

```bash
cd example
bash download_dps.sh
chimerax --nogui --exit --script "../chimerax/px_scan.py data/emd_17992.map models/*.cif --contour_level 0.0275"
```

- Default range: 0.646 Å −5% / +5% → 0.614–0.678 Å, step 0.001 Å.
- `--contour_level 0.0275` (optional) is the author-recommended contour of EMD-17992; without it ChimeraX sets the level automatically.

Quick test first (one model, narrow range):

```bash
chimerax --nogui --exit --script "../chimerax/px_scan.py data/emd_17992.map models/1DPS.cif --min 0.640 --max 0.656 --contour_level 0.0275"
```

Results are written to `data/`:

- `emd_17992_px_summary.txt`: highest-correlation voxel size per model, mean ± SD, fit correlation per model
- `emd_17992_px_scan.png`: correlation curves, nominal voxel size, best voxel size (mean ± SD)
- `emd_17992_px_scan.csv`: all correlations
- `px_scan_log.txt` (in `example/`): full log of the run
- `emd_17992_px_fits/`: fitted models; open them with the map to check the fits

## What the script does

- Fits the X-ray models into the map at the nominal voxel size: global search for each, the best-fitting one becomes the reference, and the others are superposed on it with `matchmaker`, then refined locally.
  - DPS has tetrahedral symmetry, so a global search alone can land on partially overlapping, wrong placements. Superposing on one good fit avoids this.
- Changes the voxel size step by step, from the nominal value down to `--min`, then from the nominal value up to `--max`.
- At each step: re-fits the model, simulates a map from it (`molmap`, 4 Å) on the grid of the EM map, and fits the EM map into the simulated map. The correlation reported by `fitmap` is recorded.
- Takes the voxel size with the highest correlation for each model.
- Averages over models: each X-ray structure has its own small errors (unit cell, refinement), so the spread across models is the uncertainty.
