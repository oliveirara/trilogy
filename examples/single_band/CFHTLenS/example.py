#!/usr/bin/env python
"""Example: Single band grayscale image from CFHTLenS."""

from pathlib import Path
from trilogy import Trilogy

trilogy = Trilogy(
    images="J021408.1-053532.4_W1m1p1_g.V2.2A.swarp.cut.10032_11323_19034_20325.fits",
    outname="J021408.1-053532.4_W1m1p1_g.V2.2A.swarp.cut.10032_11323_19034_20325",
    indir=Path(__file__).parent,
    # Optimized parameters from manual tuning
    noiselum=0.5,
    satpercent=0.005,
    noisesig=50,
    noisesig0=0,
    correctbias=True,
    samplesize=20000,
    stampsize=20000,
    maxstampsize=20000,
    combine="sum",
)

print("Creating CFHTLenS single band image...")
output = trilogy.save()
print(f"✓ Saved to {output}")
