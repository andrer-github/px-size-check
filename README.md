# px-size-check

Estimate the true pixel size of a cryo-EM map by comparing it with X-ray structures, in ChimeraX.

Based on the "comparison to a known structure" approach in Dickerson et al. 2024 ([doi:10.1016/j.ultramic.2023.113883](https://doi.org/10.1016/j.ultramic.2023.113883)).

## Usage

```bash
chimerax --nogui --exit --script "chimerax/px_scan.py map.mrc models/*.cif"
chimerax --nogui --exit --script "chimerax/px_scan.py maps/*.mrc models/*.cif"   # several maps
```

Note: `chimerax` should be replaced with the full path of the chimerax executable. E.g. on a Mac it could be: `/Applications/ChimeraX-1.11.1.app/Contents/bin/ChimeraX`

- Maps (`.mrc`, `.map`, `.ccp4`) and models (`.cif`, `.pdb`) are recognised by extension.

- Prints and saves the estimated pixel size (mean ± SD across models).
- No preparation needed: models are fitted into the map by the script.
- Check the fits (`<map>_px_fits/`) before trusting the result.

## Options

- `--min`, `--max`: voxel size range to test (Å), same for all maps; default: nominal voxel size of each map −5% / +5%, rounded to the step
- `--step 0.001`: step (Å)
- `--res 4.0`: resolution (Å) of the map simulated from the models
- `--contour_level`: contour level of the EM maps (default: ChimeraX automatic)
- `--prefitted`: models are already fitted in the map; skip the automatic fitting
- `--reference NAME`: model used as reference for the automatic fitting (file name without extension); default: the best-fitting one

## Method

Fitting (at the nominal voxel size):

- Global search for each model; the best-fitting one is the reference.
- All other models are superposed on the reference (`matchmaker`), then refined locally.
- Requires models of the same protein. Otherwise fit them by hand and use `--prefitted`.

Scan, for each model, from the nominal voxel size down to `--min`, then up to `--max`:

- `volume MAP voxelSize S`
- `fitmap MODEL inMap MAP`
- `molmap MODEL RES onGrid MAP`
- `fitmap MAP inMap MOLMAP` to obtain the  correlation

The voxel size with the highest correlation is taken for each model.

## Output

Next to each map:

- `<map>_px_summary.txt`: best voxel size per model, mean ± SD, fit correlation per model (a clearly lower value means a bad fit)
- `<map>_px_scan.png`: correlation curves, with the nominal voxel size (grey) and the best voxel size, mean ± SD across models (red)
- `<map>_px_scan.csv`: all correlations
- `<map>_px_fits/`: fitted models, to check the fits

In the folder where ChimeraX was started:

- `px_scan_log.txt`: full log (arguments, reference model, fit correlations, progress, summaries)
- With several maps: `px_scan_combined_summary.txt` (estimate per map, change vs nominal) and `px_scan_combined.png` (one panel per map)

## Several maps

- Use the same models and settings for all maps.
- Use the same kind of map throughout (all unsharpened, or all post-processed).
- Maps from the same microscope, detector and nominal magnification should agree within the model-to-model spread; their mean is the best estimate for that setup.
- A map that clearly deviates points to that map (processing pixel size, magnification anisotropy, energy filter, different session): check it rather than averaging it in.

## Input requirements

- Map processed at its **nominal** pixel size.
- **Several X-ray models** of the same complex.
  - Not cryo-EM models: they inherit the pixel size of the map they were built in.
- If the map contains more than the models (e.g. 70S map, 50S models), crop it first.

## Independent check: Cs from RELION

```bash
python relion/relion_cs_px.py CtfRefine/jobXXX/particles_ctf_refine.star
```

- Estimates the pixel size from the apparent Cs after CtfRefine with 4th-order aberrations (formula from the RELION documentation, with the factor 12 missing in eq. 5 of Dickerson et al.).
- Reports the result per optics group. Only the Python standard library is needed.
- Caveats: needs successful 4th-order aberration refinement (less reliable at lower resolution); invalid if the pixel size was rescaled during processing, or if voltage or Cs are wrong. Use it as a check, not as a calibration.

## Design choices

Details and test results: [`method_notes.md`](method_notes.md).

- **Resample at every step.** The simulated map is regenerated on the grid of the rescaled map. On EMD-17992, skipping this gave the same mean but ~3× broader, noisier peaks, and two of eight models stuck at the nominal voxel size.
- **Re-fit at every step**, scanning outward from the nominal voxel size, so the fit follows the map.
- **Fit via a reference model.** A global search alone can misplace models in symmetric particles; superposing all models on the best-fitting one avoids this.
- **Several X-ray models**, reported as mean ± SD.
- **Fixed contour level** (`--contour_level`) when comparing runs or maps.

## Example

- E. coli DPS, EMD-17992: deposited at 0.646 Å, gold calibration 0.6484 Å.
- See [`example/README.md`](example/README.md).

## Requirements

- [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/) (includes matplotlib, used for the plot). Tested on version 1.11.1.

## Limitations

- Only as accurate as the X-ray models: use several.
- The correlation depends on the contour level of the EM map: set `--contour_level` for reproducible results.
- Measures the pixel size of the map as processed. For absolute calibration on gold-foil grids, use [magCalEM](https://pypi.org/project/magCalEM/).
- Alternatives: magCalEM (gold lattice), EMDA (map-to-model magnification refinement), and the Cs check above. See [`method_notes.md`](method_notes.md).

## Citation

If you use this tool, please consider citing this repository and the original paper that inspired it:
- Dickerson JL, Leahy E, Peet MJ, Naydenova K, Russo CJ (2024). Accurate magnification determination for cryoEM using gold. *Ultramicroscopy* 256, 113883. [doi:10.1016/j.ultramic.2023.113883](https://doi.org/10.1016/j.ultramic.2023.113883), which describes the approach.

## Acknowledgements

- I adapted my original scripts (written in R and ChimeraX script) using Claude Opus 5.5 (Anthropic), who wrote also most of the documentation.

