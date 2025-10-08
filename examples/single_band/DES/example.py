#!/usr/bin/env python
"""Example: Single band grayscale image from DES."""

from trilogy import Trilogy

trilogy = Trilogy(
    images="J000316.4-334804.3_DES0002-3332_r4907p01_r.fits.fz",
    outname="J000316.4-334804.3_DES0002-3332_r4907p01_r",
    # Optimized parameters from manual tuning
    noiselum=1.01,
    satpercent=0.001,
    noisesig=0.05,
    noisesig0=1,
    correctbias=True,
    samplesize=20000,
    stampsize=20000,
    maxstampsize=20000,
    combine="sum",
)

print("Creating DES single band image...")
output = trilogy.save()
print(f"✓ Saved to {output}")
