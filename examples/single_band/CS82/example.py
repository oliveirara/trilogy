#!/usr/bin/env python
"""Example: Single band grayscale image from CS82."""

from trilogy import Trilogy

trilogy = Trilogy(
    images="J001424.3+004145.3_S82p4p_y.V2.7A.swarp.cut.fits",
    outname="J001424.3+004145.3_S82p4p_y.V2.7A.swarp",
    # Optimized parameters from manual tuning
    noiselum=3.0,
    satpercent=0.01,
    noisesig=1,
    noisesig0=0.009,
    correctbias=True,
    samplesize=20000,
    stampsize=20000,
    maxstampsize=20000,
    combine="sum",
)

print("Creating CS82 single band image...")
output = trilogy.save()
print(f"✓ Saved to {output}")
