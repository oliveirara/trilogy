#!/usr/bin/env python
"""Example: RGB image from KIDS survey."""

from trilogy import Trilogy

images = {
    "R": ["J120447.2-024312.3_KiDS_DR4.0_181.5_-2.5_i_sci.fits"],
    "G": ["J120447.2-024312.3_KiDS_DR4.0_181.5_-2.5_r_sci.fits"],
    "B": ["J120447.2-024312.3_KiDS_DR4.0_181.5_-2.5_g_sci.fits"],
}

trilogy = Trilogy(
    images=images,
    outname="J120447.2-024312.3_KiDS_DR4.0_181.5_-2.5_sci",
    # Optimized parameters from manual tuning
    noiselum=0.5,
    satpercent=0.0009,
    noisesig=35,
    noisesig0=35,
    correctbias=False,
    samplesize=20000,
    stampsize=20000,
    maxstampsize=20000,
    combine="sum",
)

print("Creating KIDS RGB image...")
output = trilogy.save()
print(f"✓ Saved to {output}")
