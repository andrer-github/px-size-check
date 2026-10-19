"""
px_scan.py - estimate the true pixel size of cryo-EM maps from X-ray models.

Usage (one command):
  chimerax --nogui --exit --script "px_scan.py MAP [MAP ...] MODEL [MODEL ...] [options]"
  Maps (.mrc, .map, .ccp4) and models (.cif, .pdb) are recognised by extension.
  Wildcards are allowed, e.g. "px_scan.py maps/*.mrc models/*.cif"

Options:
  --min, --max      voxel size range to test (A), same for all maps; default:
                    nominal voxel size of each map -5% / +5%, rounded to the step
  --step 0.001      step (A)
  --res 4.0         resolution (A) of the map simulated from the model (molmap)
  --contour_level X contour level of the EM maps (default: ChimeraX automatic)
  --prefitted       models are already fitted in the map: skip the automatic fitting
  --reference NAME  model used as reference for the automatic fitting (file name
                    without extension); default: the one that fits best

Steps, for each map:
  1. fit the models into the map: global search for each model, the best-fitting
     one is the reference, all others are superposed on it (matchmaker) and
     refined locally (skipped with --prefitted)
  2. for each model, starting from the nominal voxel size, step down to --min,
     then from the nominal size step up to --max. At each voxel size:
        volume MAP voxelSize S
        fitmap MODEL inMap MAP
        molmap MODEL RES onGrid MAP
        fitmap MAP inMap MOLMAP        -> correlation
  3. take the voxel size with the highest correlation, per model; mean +/- SD

Outputs, next to each map:
  <map>_px_summary.txt   best voxel size per model, mean +/- SD, fit correlations
  <map>_px_scan.csv      all correlations
  <map>_px_scan.png      correlation curves, nominal voxel size, best (mean +/- SD)
  <map>_px_fits/         fitted models (.cif), to check the fits in ChimeraX
In the folder where ChimeraX was started:
  px_scan_log.txt        full log of the run
  px_scan_combined_summary.txt, px_scan_combined.png   (only with several maps)
"""

import argparse
import csv
import glob
import os
import sys
from pathlib import Path

import numpy as np

MAP_EXT = (".mrc", ".map", ".ccp4", ".mrc.gz", ".map.gz", ".ccp4.gz")
MODEL_EXT = (".cif", ".pdb", ".mmcif", ".ent", ".cif.gz", ".pdb.gz")
LOG_FILE = None


# ------------------------------------------------------------------ helpers

def log(msg):
    """Print, show in the ChimeraX log, and append to the log file."""
    session.logger.info(msg)
    print(msg, flush=True)
    if LOG_FILE is not None:
        with open(LOG_FILE, "a") as fh:
            fh.write(msg + "\n")


def stem_of(path):
    return Path(path).name.split(".")[0]


def run_new(command, cls):
    """Run a command and return the first new model of the requested class."""
    from chimerax.core.commands import run
    before = set(session.models.list())
    run(session, command, log=False)
    new = [m for m in session.models.list() if m not in before and isinstance(m, cls)]
    if not new:
        raise RuntimeError(f"No {cls.__name__} created by: {command}")
    return new[0]


def cmd(command):
    from chimerax.core.commands import run
    return run(session, command, log=False)


def fit_correlation(command):
    """Run a fitmap command and return the correlation it reports."""
    fits = cmd(command)
    fit = fits[0] if isinstance(fits, (list, tuple)) else fits
    corr = fit.correlation() if hasattr(fit, "correlation") else fit.stats["correlation"]
    if corr is None:
        raise RuntimeError(f"fitmap reported no correlation: {command}")
    return float(corr)


def mean_sd(values):
    values = np.asarray(values, dtype=float)
    sd = values.std(ddof=1) if len(values) > 1 else float("nan")
    return values.mean(), sd


# ------------------------------------------------------------------ plotting

def model_colours(model_names):
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap("tab10" if len(model_names) <= 10 else "tab20")
    return {n: cmap(i % cmap.N) for i, n in enumerate(model_names)}


def plot_panel(ax, res, colours):
    """Curves normalised to each model's maximum; nominal size, mean best size and SD band."""
    for name, _ in res["peaks"]:
        pts = sorted((r["voxel_size"], r["correlation"]) for r in res["rows"] if r["model"] == name)
        x, y = zip(*pts)
        top = max(y)
        ax.plot(x, [v / top for v in y], "o", markersize=3, color=colours[name], label=name)
    mean, sd = mean_sd([b for _, b in res["peaks"]])
    ax.axvline(res["nominal"], color="grey", linestyle=":", linewidth=1.5,
               label=f"nominal  {res['nominal']:.4f}")
    if not np.isnan(sd):
        ax.axvspan(mean - sd, mean + sd, color="red", alpha=0.15, linewidth=0)
        best_label = f"best  {mean:.4f} \u00b1 {sd:.4f}"
    else:
        best_label = f"best  {mean:.4f}"
    ax.axvline(mean, color="red", linewidth=1.5, label=best_label)
    ax.set_title(res["map"], fontsize=10)
    ax.set_xlabel("Voxel size (\u00c5)")
    ax.set_ylabel("Correlation (normalised)")
    ax.spines[["top", "right"]].set_visible(False)


def save_map_plot(res, colours, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    plot_panel(ax, res, colours)
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def save_combined_plot(results, colours, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(results)
    ncol = min(n, 3)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.5 * ncol, 3.6 * nrow), squeeze=False)
    for ax, res in zip(axes.flat, results):
        plot_panel(ax, res, colours)
        mean, sd = mean_sd([b for _, b in res["peaks"]])
        ax.set_title(f"{res['map']}\nnominal {res['nominal']:.4f}, best {mean:.4f} \u00b1 {sd:.4f} \u00c5",
                     fontsize=9)
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    keep = [(h, l) for h, l in zip(handles, labels) if not l.startswith(("nominal", "best"))]
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    keep += [(Line2D([], [], color="grey", linestyle=":"), "nominal"),
             (Line2D([], [], color="red"), "best (mean of models)"),
             (Patch(color="red", alpha=0.15), "\u00b1 SD")]
    fig.legend(*zip(*keep), frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ text output

def map_summary_lines(res, args):
    mean, sd = mean_sd([b for _, b in res["peaks"]])
    lines = [f"Map: {res['map']}",
             f"Nominal voxel size: {res['nominal']:.4f} A",
             f"Scan: {res['min']}-{res['max']} A, step {args.step} A, molmap {args.res} A",
             "",
             "Best voxel size per model (fit correlation at nominal size in brackets):"]
    lines += [f"  {n:<20} {b:.4f} A   ({res['fitq'][n]:.3f})" for n, b in res["peaks"]]
    lines += ["",
              f"Estimated pixel size: {mean:.4f} +/- {sd:.4f} A (SD, {len(res['peaks'])} models)",
              f"Change vs nominal: {100 * (mean / res['nominal'] - 1):+.2f}%"]
    return lines


def combined_summary_lines(results):
    lines = ["Combined summary", "",
             f"  {'map':<24}{'nominal':>9}{'best (mean)':>13}{'SD':>9}{'change':>9}"]
    changes = []
    for res in results:
        mean, sd = mean_sd([b for _, b in res["peaks"]])
        ch = 100 * (mean / res["nominal"] - 1)
        changes.append(ch)
        lines.append(f"  {res['map']:<24}{res['nominal']:>9.4f}{mean:>13.4f}{sd:>9.4f}{ch:>+8.2f}%")
    lines.append("")
    nominals = {round(r["nominal"], 6) for r in results}
    if len(nominals) == 1:
        m, s = mean_sd([mean_sd([b for _, b in r["peaks"]])[0] for r in results])
        lines.append(f"Same nominal voxel size: mean across maps {m:.4f} +/- {s:.4f} A (SD, {len(results)} maps)")
    m, s = mean_sd(changes)
    lines.append(f"Change vs nominal, mean across maps: {m:+.2f} +/- {s:.2f}% (SD, {len(results)} maps)")
    return lines


# ------------------------------------------------------------------ one map

def process_map(map_path, model_files, args):
    from chimerax.atomic import AtomicStructure
    from chimerax.map import Volume

    stem = stem_of(map_path)
    out_dir = map_path.parent
    fit_dir = out_dir / f"{stem}_px_fits"
    fit_dir.mkdir(exist_ok=True)

    vol = run_new(f'open "{map_path}"', Volume)
    if args.contour_level is not None:
        cmd(f"volume {vol.atomspec} level {args.contour_level}")
    nominal = round(float(vol.data.step[0]), 6)
    pos0_map = vol.position

    on_step = lambda x: round(round(x / args.step) * args.step, 6)
    vmin = args.min if args.min is not None else on_step(nominal * 0.95)
    vmax = args.max if args.max is not None else on_step(nominal * 1.05)
    n_left = int(round((nominal - vmin) / args.step))
    n_right = int(round((vmax - nominal) / args.step))
    left = [round(nominal - i * args.step, 6) for i in range(0, n_left + 1)]
    right = [round(nominal + i * args.step, 6) for i in range(1, n_right + 1)]

    log(f"\n=== Map {map_path.name}: nominal voxel size {nominal:.4f} A")
    log(f"Scanning {vmin}-{vmax} A, step {args.step} A, {len(left) + len(right)} sizes, "
        f"{len(model_files)} model(s), molmap {args.res} A")

    # ---- Fit the models at the nominal voxel size
    structs = {stem_of(f): run_new(f'open "{f}"', AtomicStructure) for f in model_files}

    def fit_quality(s):
        return fit_correlation(f"fitmap {s.atomspec} inMap {vol.atomspec} "
                               f"resolution {args.res} shift false rotate false")

    if not args.prefitted:
        if args.reference:
            if args.reference not in structs:
                raise RuntimeError(f"--reference {args.reference} is not among the models")
            candidates = [args.reference]
        else:
            candidates = list(structs)
        quality = {}
        for name in candidates:
            s = structs[name]
            cmd(f"fitmap {s.atomspec} inMap {vol.atomspec} resolution {args.res} search 50")
            cmd(f"fitmap {s.atomspec} inMap {vol.atomspec}")
            quality[name] = fit_quality(s)
        ref = max(quality, key=quality.get)
        log(f"Reference model: {ref} (fit correlation {quality[ref]:.4f})")
        for name, s in structs.items():
            if name != ref:
                cmd(f"matchmaker {s.atomspec} to {structs[ref].atomspec}")

    pos0_model, fitq = {}, {}
    for name, s in structs.items():
        cmd(f"fitmap {s.atomspec} inMap {vol.atomspec}")
        fitq[name] = fit_quality(s)
        pos0_model[name] = s.position
        cmd(f'save "{fit_dir / (name + "_fit.cif")}" models {s.atomspec} relModel {vol.atomspec}')
        log(f"{name}: fit correlation {fitq[name]:.4f}")

    # ---- Scan
    rows, peaks = [], []
    for name, struct in structs.items():
        log(f"{name}: scanning")
        scores = {}
        for side, sizes in (("left", left), ("right", right)):
            # each side starts from the fitted state at the nominal size
            cmd(f"volume {vol.atomspec} voxelSize {nominal}")
            vol.position = pos0_map
            struct.position = pos0_model[name]
            for size in sizes:
                cmd(f"volume {vol.atomspec} voxelSize {size}")
                cmd(f"fitmap {struct.atomspec} inMap {vol.atomspec}")
                sim = run_new(f"molmap {struct.atomspec} {args.res} onGrid {vol.atomspec}", Volume)
                corr = fit_correlation(f"fitmap {vol.atomspec} inMap {sim.atomspec}")
                session.models.close([sim])
                scores[size] = corr
                rows.append({"map": stem, "model": name, "voxel_size": size,
                             "correlation": corr, "side": side})
        best = max(scores, key=scores.get)
        peaks.append((name, best))
        log(f"  best {best:.4f} A (correlation {scores[best]:.4f})")

    session.models.close(list(structs.values()) + [vol])
    return {"map": stem, "out_dir": out_dir, "nominal": nominal, "min": vmin, "max": vmax,
            "rows": rows, "peaks": peaks, "fitq": fitq}


def write_map_outputs(res, args, colours):
    stem, out_dir = res["map"], res["out_dir"]
    with open(out_dir / f"{stem}_px_scan.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res["rows"][0].keys()))
        w.writeheader()
        w.writerows(res["rows"])
    lines = map_summary_lines(res, args)
    (out_dir / f"{stem}_px_summary.txt").write_text("\n".join(lines) + "\n")
    log("\n" + "\n".join(lines))
    try:
        save_map_plot(res, colours, out_dir / f"{stem}_px_scan.png")
        log(f"Plot saved: {stem}_px_scan.png")
    except Exception as e:   # the results above are already saved
        log(f"Plot not saved: {e}")


# ------------------------------------------------------------------ main

def main():
    global LOG_FILE
    p = argparse.ArgumentParser(prog="px_scan.py")
    p.add_argument("files", nargs="+")
    p.add_argument("--min", type=float, default=None)
    p.add_argument("--max", type=float, default=None)
    p.add_argument("--step", type=float, default=0.001)
    p.add_argument("--res", type=float, default=4.0)
    p.add_argument("--contour_level", type=float, default=None)
    p.add_argument("--prefitted", action="store_true")
    p.add_argument("--reference", default=None)
    args = p.parse_args(sys.argv[1:])

    # Relative paths and wildcards are resolved from the folder ChimeraX was started in
    base = Path(os.environ.get("PWD", os.getcwd()))
    resolve = lambda f: Path(f) if Path(f).is_absolute() else base / f
    files = sorted({Path(m).resolve() for pat in args.files for m in glob.glob(str(resolve(pat)))})
    maps = [f for f in files if f.name.lower().endswith(MAP_EXT)]
    models = [f for f in files if f.name.lower().endswith(MODEL_EXT)]
    other = [f for f in files if f not in maps and f not in models]
    if not maps or not models:
        raise RuntimeError(f"Need at least one map and one model. Got maps: {maps}, models: {models}")

    LOG_FILE = base / "px_scan_log.txt"
    LOG_FILE.write_text("")
    from datetime import datetime
    try:
        from chimerax.core import version as cx_version
    except ImportError:
        cx_version = "unknown"
    log(f"px_scan.py  {datetime.now():%Y-%m-%d %H:%M}  ChimeraX {cx_version}")
    log("Arguments: " + " ".join(sys.argv[1:]))
    log("Maps: " + ", ".join(m.name for m in maps))
    log("Models: " + ", ".join(m.name for m in models))
    if other:
        log("Ignored (unknown extension): " + ", ".join(o.name for o in other))

    try:
        colours = model_colours([stem_of(m) for m in models])
    except Exception:
        colours = {stem_of(m): None for m in models}

    results = []
    for map_path in maps:
        res = process_map(map_path, models, args)
        write_map_outputs(res, args, colours)
        results.append(res)

    if len(results) > 1:
        lines = combined_summary_lines(results)
        (base / "px_scan_combined_summary.txt").write_text("\n".join(lines) + "\n")
        log("\n" + "\n".join(lines))
        try:
            save_combined_plot(results, colours, base / "px_scan_combined.png")
            log("Plot saved: px_scan_combined.png")
        except Exception as e:
            log(f"Combined plot not saved: {e}")


main()
