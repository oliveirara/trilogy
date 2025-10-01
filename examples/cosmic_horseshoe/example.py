#!/usr/bin/env python
"""Example: Cosmic Horseshoe - RGB image from Hubble WFC3."""

from pathlib import Path
from trilogy import Trilogy

# RGB composite from three HST filters
images = {
    "R": ["hst_mos_1020228_wfc3_uvis_f814w_drz_cutout_20asec.fits"],
    "G": ["hst_mos_1020228_wfc3_uvis_f606w_drz_cutout_20asec.fits"],
    "B": ["hst_mos_1020228_wfc3_uvis_f475w_drz_cutout_20asec.fits"],
}

trilogy = Trilogy(
    images=images,
    outname="cosmic_horseshoe",
    # Optimized parameters from manual tuning
    noiselum=0.0009,
    satpercent=0.0009,
    noisesig=0.0009,
    noisesig0=0.009,
    correctbias=True,
    bscale=0.0009,
    samplesize=10000,
    stampsize=10000,
    maxstampsize=10000,
    combine="average",
)

print("Creating Cosmic Horseshoe RGB image...")
output = trilogy.save()
print(f"✓ Saved to {output}")
