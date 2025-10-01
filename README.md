# trilogy

Python library for converting astronomical FITS images into beautiful color or grayscale images. Originally written by [Dan Coe](https://www.stsci.edu/~dcoe), now refactored and optimized for Python 3.11+.

Modified by [Renan Alves de Oliveira](https://github.com/oliveirara)

<a href="https://ascl.net/1508.009"><img src="https://img.shields.io/badge/ascl-1508.009-blue.svg?colorB=262255" alt="ascl:1508.009" /></a>

<p align="center">
  <img width="300" src="https://raw.githubusercontent.com/oliveirara/trilogy/main/examples/cosmic_horseshoe/cosmic_horseshoe.png" alt="Cosmic Horseshoe" title="Cosmic Horseshoe">
  <br>
  <em><strong>Cosmic Horseshoe</strong> image taken by Hubble WFC3, using filters F475W, F606W, and F814W.
              <br>We use <b>trilogy</b> to combine all FITS files into this beautiful RGB image!
  </em>
</p>

Try it! [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/oliveirara/trilogy/HEAD)

## Installation

### From PyPI:

```bash
pip install trilogy
```

### From GitHub:

```bash
git clone https://github.com/oliveirara/trilogy.git
cd trilogy
pip install -e .
```

## Quick Start

### Single Band (Grayscale) Image

```python
from trilogy import Trilogy

# Simple single image
trilogy = Trilogy(
    images="path/to/image.fits",
    outname="output",
    noiselum=0.15,
    satpercent=0.001
)
img = trilogy.run()
img.show()  # Display image
```

### Multi-Band (RGB) Image

```python
from trilogy import Trilogy

# RGB composite from three filters
images = {
    "R": ["path/to/red_filter.fits"],
    "G": ["path/to/green_filter.fits"], 
    "B": ["path/to/blue_filter.fits"]
}

trilogy = Trilogy(
    images=images,
    outname="rgb_output",
    noiselum=0.15,
    satpercent=0.001,
    colorsatfac=1.5  # Boost color saturation
)
img = trilogy.save("output.png")  # Save directly
```

### Advanced Configuration

```python
from trilogy import Trilogy, TrilogyConfig
from pathlib import Path

# Use configuration dataclass for complex settings
config = TrilogyConfig(
    indir=Path("data/fits"),
    outdir=Path("output"),
    outname="advanced_image",
    
    # Scaling parameters
    satpercent=0.01,
    noiselum=0.2,
    noiselums={"R": 0.2, "G": 0.15, "B": 0.1},  # Per-channel
    colorsatfac=1.3,
    
    # Processing
    samplesize=2000,
    stampsize=2000,
    combine="average",  # or "sum"
    
    # Visual
    invert=False,
    legend=True  # Add filter legend to image
)

trilogy = Trilogy(
    images={
        "R": ["i_band.fits"],
        "G": ["r_band.fits"],
        "B": ["g_band.fits"]
    },
    config=config
)

img = trilogy.make_image()
img.save("output.png")
```

## Usage in Notebooks

Trilogy is designed to work seamlessly in Jupyter notebooks:

```python
from trilogy import Trilogy
from pathlib import Path

# Process and display inline
t = Trilogy(
    images={"R": ["i.fits"], "G": ["r.fits"], "B": ["g.fits"]},
    outname="notebook_output",
    noiselum=0.15
)

# Returns PIL Image that displays automatically in notebooks
img = t.run()
```

See `examples/with_notebook/notebook.ipynb` for complete examples.

## Configuration Parameters

### Input/Output
- `indir`: Input directory (default: current directory)
- `outdir`: Output directory (default: current directory)  
- `outname`: Output filename without extension

### Image Scaling
- `satpercent`: Percentage of pixels to saturate (default: 0.001)
- `noiselum`: Noise luminosity level 0-1 (default: 0.15)
- `noiselums`: Per-channel noise luminosity dict
- `noisesig`: Noise sigma for output (default: 1)
- `noisesig0`: Noise sigma for measurement (default: 2)
- `correctbias`: Measure noise mean vs assume 0 (default: False)
- `colorsatfac`: Color saturation factor, >1 boosts (default: 1)

### Processing
- `samplesize`: Sample region size for scaling (default: 1000)
- `sampledx`, `sampledy`: Sample offsets (default: 0)
- `stampsize`: Processing stamp size (default: 1000)
- `maxstampsize`: Maximum stamp size (default: 6000)
- `combine`: "average" or "sum" for multiple images (default: "average")
- `bscale`: Multiply all image values (default: 1)
- `bzero`: Add to all image values (default: 0)

### Advanced
- `noise`: Manual noise level (overrides automatic)
- `saturate`: Manual saturation level (overrides automatic)
- `invert`: Invert luminosity (default: False)
- `legend`: Add filter legend to RGB images (default: False)

## What's New in v1.0

- ✅ **Removed CLI** - Focus on Python/notebook usage
- ✅ **Modern Python 3.11+** - Type hints, dataclasses, pathlib
- ✅ **Performance improvements** - Caching, better numpy usage
- ✅ **Cleaner API** - Simplified interface, better error handling
- ✅ **No interactive prompts** - Fully scriptable
- ✅ **Better documentation** - Type-safe and IDE-friendly

## Requirements

- Python >= 3.11
- astropy >= 6.1.7
- numpy >= 2.2.4
- pillow >= 11.2.1
- scipy >= 1.15.2

## Algorithm

Trilogy uses **Lupton's method** for creating color images from astronomical data:
- Logarithmic scaling for wide dynamic range
- Robust statistics with sigma clipping
- Automatic or manual level determination
- Color saturation adjustment
- Per-channel noise luminosity control

## Resources

- [Lupton's method](http://www.astro.princeton.edu/~rhl/PrettyPictures/)
- [STIFF](https://github.com/astromatic/stiff)
- [Original trilogy](https://www.stsci.edu/~dcoe/trilogy/)

## Examples

Check the `examples/` directory for sample FITS files and configurations:
- `cosmic_horseshoe/` - HST WFC3 multi-band example
- `single_band/` - Single filter examples
- `multiple_bands/` - Various survey examples (HSC, KIDS, Legacy, etc.)
- `with_notebook/` - Jupyter notebook examples

## License

See LICENSE file.

## Citation

If you use trilogy in your research, please cite:
```
Coe, D. 2012, "Trilogy", Astrophysics Source Code Library, record ascl:1508.009
```
