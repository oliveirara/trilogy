#!/usr/bin/env python
"""Example: Cosmic Horseshoe - RGB image from Hubble WFC3."""

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
    # Corrected parameters - original had bscale too small
    noiselum=0.10,
    satpercent=0.0001,
    noisesig=1,
    noisesig0=2,
    correctbias=False,
    bscale=1.0,  # Fixed: was 0.0009, making image too dark
    colorsatfac=1.3,  # Boost color saturation
    samplesize=10000,
    stampsize=10000,
    maxstampsize=10000,
    combine="average",
)

print("Creating Cosmic Horseshoe RGB image...")
output = trilogy.save()
print(f"✓ Saved to {output}")
