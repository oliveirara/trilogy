#!/usr/bin/env python
"""Example: RGB image from Legacy Survey."""

from trilogy import Trilogy

images = {
    "R": ["J054358.0-303449.5_legacysurvey-0860m305-image-z.fits.fz"],
    "G": ["J054358.0-303449.5_legacysurvey-0860m305-image-r.fits.fz"],
    "B": ["J054358.0-303449.5_legacysurvey-0860m305-image-g.fits.fz"],
}

trilogy = Trilogy(
    images=images,
    outname="J054358.0-303449.5_legacysurvey-0860m305-image",
    # Optimized parameters from manual tuning
    noiselum=0.5,
    satpercent=0.005,
    noisesig=45,
    noisesig0=0,
    correctbias=True,
    samplesize=20000,
    stampsize=20000,
    maxstampsize=20000,
    combine="sum",
)

print("Creating Legacy Survey RGB image...")
output = trilogy.save()
print(f"✓ Saved to {output}")
