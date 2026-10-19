#!/usr/bin/env python3
"""
relion_cs_px.py - pixel size from the apparent Cs refined by RELION CtfRefine.

After CtfRefine with 4th-order aberrations, an error in the pixel size shows up
as an error in the spherical aberration (Cs). From the RELION documentation
(Reference > Pixel size issues):

    Cs_apparent = Cs + 12 * Z40 * 1e-7 / (pi * lambda^3)      (Cs in mm, lambda in A)
    true pixel size = nominal pixel size * (Cs / Cs_apparent)^(1/4)

Z40 is the 7th value of rlnEvenZernike. Note: eq. 5 of Dickerson et al. 2024
omits the factor 12.

Usage:
  python relion_cs_px.py CtfRefine/job123/particles_ctf_refine.star
  python relion_cs_px.py CtfRefine/job123/            # a job folder also works
  python relion_cs_px.py job1/ job2/ --cs 2.7         # override the true Cs (mm)

Only the Python standard library is needed.
"""

import argparse
import math
import sys
from pathlib import Path


def electron_wavelength(kv):
    """Relativistic electron wavelength (A) for an accelerating voltage in kV."""
    v = kv * 1e3
    return 12.2643247 / math.sqrt(v * (1 + 0.978466e-6 * v))


def read_optics(star):
    """Return the data_optics table of a RELION star file as a list of dicts."""
    lines = Path(star).read_text().splitlines()
    try:
        i = next(k for k, l in enumerate(lines) if l.strip() == "data_optics")
    except StopIteration:
        raise ValueError(f"No data_optics block in {star}")
    cols, rows = [], []
    for line in lines[i + 1:]:
        s = line.strip()
        if s.startswith("data_"):
            break
        if not s or s.startswith("#") or s == "loop_":
            continue
        if s.startswith("_rln"):
            cols.append(s.split()[0][4:])          # e.g. EvenZernike
        elif cols:
            rows.append(dict(zip(cols, s.split())))
    return rows


def parse_list(value):
    return [float(x) for x in value.strip("[]").split(",") if x]


def true_pixel_size(group, cs_true=None):
    if "EvenZernike" not in group:
        raise ValueError("no rlnEvenZernike: run CtfRefine with 4th-order aberrations")
    z = parse_list(group["EvenZernike"])
    if len(z) < 7:
        raise ValueError("rlnEvenZernike has fewer than 7 values: 4th-order aberrations not refined")
    z40 = z[6]
    kv = float(group["Voltage"])
    cs = cs_true if cs_true is not None else float(group["SphericalAberration"])
    lam = electron_wavelength(kv)
    cs_app = cs + 12 * z40 * 1e-7 / (math.pi * lam ** 3)
    factor = (cs / cs_app) ** 0.25
    return {"z40": z40, "kv": kv, "cs": cs, "cs_app": cs_app, "factor": factor}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("star", nargs="+", help="particles_ctf_refine.star, or the CtfRefine job folder")
    p.add_argument("--cs", type=float, default=None,
                   help="true Cs (mm); default: rlnSphericalAberration from the star file")
    args = p.parse_args()

    for item in args.star:
        star = Path(item)
        if star.is_dir():
            star = star / "particles_ctf_refine.star"
        print(f"\n{star}")
        try:
            groups = read_optics(star)
        except (OSError, ValueError) as e:
            print(f"  error: {e}")
            continue
        print(f"  {'group':<14}{'kV':>6}{'Cs':>7}{'Cs app.':>9}{'factor':>9}"
              f"{'mic. px nom':>13}{'-> true':>9}{'part. px nom':>14}{'-> true':>9}")
        for g in groups:
            name = g.get("OpticsGroupName", g.get("OpticsGroup", "?"))
            try:
                r = true_pixel_size(g, args.cs)
            except ValueError as e:
                print(f"  {name:<14}error: {e}")
                continue
            mic = float(g["MicrographOriginalPixelSize"]) if "MicrographOriginalPixelSize" in g else float("nan")
            img = float(g["ImagePixelSize"]) if "ImagePixelSize" in g else float("nan")
            print(f"  {name:<14}{r['kv']:>6.0f}{r['cs']:>7.2f}{r['cs_app']:>9.3f}{r['factor']:>9.5f}"
                  f"{mic:>13.4f}{mic * r['factor']:>9.4f}{img:>14.4f}{img * r['factor']:>9.4f}")


if __name__ == "__main__":
    sys.exit(main())
