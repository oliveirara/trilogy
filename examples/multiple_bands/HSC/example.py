#!/usr/bin/env python
"""Example: RGB image from HSC survey."""

from trilogy import Trilogy

images = {
    "R": ["J114444.8+001347.0-HSC-I-pdr3_wide.fits[1]"],
    "G": ["J114444.8+001347.0-HSC-R-pdr3_wide.fits[1]"],
    "B": ["J114444.8+001347.0-HSC-G-pdr3_wide.fits[1]"],
}

trilogy = Trilogy(
    images=images,
    outname="J114444.8+001347.0-HSC-pdr3_wide",
    # Optimized parameters from manual tuning
    noiselum=0.5,
    satpercent=0.0009,
    noisesig=50,
    noisesig0=10,
    correctbias=False,
    samplesize=20000,
    stampsize=20000,
    maxstampsize=20000,
    combine="sum",
)

print("Creating HSC RGB image...")
output = trilogy.save()
print(f"✓ Saved to {output}")
