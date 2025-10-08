#!/usr/bin/env python
"""Example: RGB image from RCSLenS survey."""

from trilogy import Trilogy

images = {
    "R": ["J105750.9+573026.0_CDE1040A2_z.V2.7A.swarp.cut.release.12157_13448_14303_15594.fits"],
    "G": ["J105750.9+573026.0_CDE1040A2_r.V2.7A.swarp.cut.release.12157_13448_14303_15594.fits"],
    "B": ["J105750.9+573026.0_CDE1040A2_g.V2.7A.swarp.cut.release.12157_13448_14303_15594.fits"],
}

trilogy = Trilogy(
    images=images,
    outname="J105750.9+573026.0_CDE1040A2_g.V2.7A.swarp.cut.release.12157_13448_14303_15594",
    # Optimized parameters from manual tuning
    noiselum=0.5,
    satpercent=0.0005,
    noisesig=35,
    noisesig0=20,
    correctbias=False,
    samplesize=20000,
    stampsize=20000,
    maxstampsize=20000,
    combine="sum",
)

print("Creating RCSLenS RGB image...")
output = trilogy.save()
print(f"✓ Saved to {output}")
