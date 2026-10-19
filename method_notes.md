# Method notes

Tests behind the choices made in `px_scan.py`. Most were run during development on unpublished ribosome maps; the resampling comparison was repeated on the DPS example (EMD-17992).

## Resample the X-ray map at every step

- The map simulated from the X-ray model (`molmap`) is regenerated `onGrid` of the tested map **after every voxel size change**.
- Without resampling, the rescaled map and the simulated map are on different grids, so the comparison is not meaningful.
- Resampling only once, before the scan, produced a false local maximum at the original voxel size of the map.
- Resampling on the grid of the tested map follows the recommendation of the EMDA tutorial (Murshudov group, MRC-LMB).

### Test on EMD-17992

Same 8 fitted models, contour level 0.025, range 0.614–0.678 Å, step 0.001 Å, molmap 4 Å.

- **With resampling** (`px_scan.py`): molmap regenerated `onGrid` and model re-fitted at every step.
- **Without resampling**: molmap generated once at the nominal voxel size, on its own grid; no re-fit.

| | With resampling | Without resampling |
|---|---|---|
| Estimate | 0.6485 ± 0.0016 Å | 0.6484 ± 0.0022 Å |
| Correlation drop within ±0.01 Å of the peak | 1.8–2.0% | 0.5–1.1% |
| Scatter around the peak | 0.00005–0.00014 | 0.00011–0.00033 |
| Models with the best value exactly at the nominal size | 0 | 2 |

- The mean is the same here (gold calibration: 0.6484 Å; Dickerson et al. 2024).
- Without resampling, peaks are ~3× broader, noisier, and the spread across models is larger.
- Two models stayed at the nominal voxel size without resampling: the same bias towards the original size seen in the ribosome tests. With a lower-resolution map and flatter peaks, this bias can decide the result.
- Absolute correlation values are not comparable between the two (different grids); only the curve shapes are.
- With resampling, curves get noisier far from the peak (beyond ~±3%), where re-fitting on a strongly rescaled grid is harder. The region around the peak is not affected.

## Display step (`step 1` vs `step 2`)

- `fitmap` and `molmap onGrid` use the grid at the map's current display step.
- `step 1` and `step 2` gave very similar results; `step 1` was much slower.
- The script uses the ChimeraX default.

## Order of the voxel sizes

- Ascending and random order gave similar maxima, once the X-ray map was resampled at every step.
- The script scans from the nominal voxel size outwards (down, then up), re-fitting the model at each step, so the fit follows the map gradually.

## Fitting the models

- A global search alone can place a model wrongly in symmetric particles (e.g. tetrahedral DPS), giving flat or jumping curves.
- Superposing all models on the best-fitting one (`matchmaker`) solved this.
- The fit correlation in the summary flags bad fits: one model clearly lower than the others should be checked.

## Contour level

- `fitmap MAP inMap MOLMAP` uses the EM voxels above the contour level, so the level affects the correlation.
- Set `--contour_level` explicitly when comparing runs or maps.

## Several X-ray models

- Single models can disagree by a few thousandths of an Å (unit cell and refinement errors).
- Use several and report mean ± SD, as in Dickerson et al. 2024.

## Independent check: Cs from RELION CtfRefine (`relion/relion_cs_px.py`)

- After CtfRefine with 4th-order aberrations, a pixel size error appears as an error in the apparent Cs.
- Formula from the RELION documentation (Reference > Pixel size issues):
  - Cs_apparent = Cs + 12 · Z₄⁰ · 10⁻⁷ / (π λ³), with Z₄⁰ the 7th value of `rlnEvenZernike`
  - true pixel size = nominal pixel size · (Cs / Cs_apparent)^(1/4)
- Eq. 5 of Dickerson et al. 2024 omits the factor 12. Without it, the correction is ~12× too small.
- The script computes λ from the voltage and takes Cs from the star file, per optics group.

Caveats:

- Valid only if 4th-order aberrations were refined successfully. This becomes less reliable as resolution decreases.
- Invalid if the pixel size was rescaled during processing, or if the voltage or the true Cs are wrong.
- Cs and defocus partly compensate each other: treat the result as a check, not as a calibration.

## RELION-rescaled maps

- During development, maps rescaled in RELION (PostProcess with a new pixel size) were also scanned. Their maximum correlations were lower than for the same maps rescaled in ChimeraX, and the best values were not identical.
- For the scan, rescale in ChimeraX (as the script does).
- If a new pixel size is adopted, re-run CTF refinement at that pixel size rather than only rescaling the final map.

## Alternatives

- **Gold lattice (magCalEM)**: absolute calibration from gold-foil grids (Dickerson et al. 2024). Preferred when gold is in the images.
- **EMDA** (Murshudov group): includes map-to-model magnification refinement. During development it could not be used because of a hard-to-install dependency (proshade); this may have changed in newer versions.
