# Trilogy

Python script that converts astronomical FITS images into color/grayscale images. [Trilogy](https://www.stsci.edu/~dcoe/trilogy/) was originally written in Python 2.

<a  href="https://ascl.net/1508.009"><img  src="https://img.shields.io/badge/ascl-1508.009-blue.svg?colorB=262255"  alt="ascl:1508.009" /></a>

Author: [Dan Coe](https://www.stsci.edu/~dcoe)

Modified by [Renan Alves de Oliveira](https://github.com/oliveirara)

## Installation

1. From **PyPi**:

```bash
pip install trilogy
```

2. From **Github**:

```bash
git clone https://github.com/oliveirara/trilogy.git
cd trilogy
python setup.py build
python setup.py install
```

## Usage

1. Command line:

```bash
trilogy-cl -params
```

2. With input file (e.g. see *.in in ~/examples/):

```bash
trilogy-cl single.in
```

3. With input file and command line parameters:

```bash
trilogy-cl single.in -deletefilters 0 -showwith PIL -sampledx 300
```

4. Check notebooks in `~/notebooks` for examples.

## Requirements

* `Pillow`
* `astropy`
* `numpy`
* `scipy`

## Parameters and default values

```python
"indir" = ''          # Input directory.
"outname" = None      # Output filename.
"outdir" = ''         # Output directory.
"saturate" = None     # Determined automatically if None: image data value allowed to saturate.
"satpercent" = 0.001  # Percentage of pixels which will be saturated.
"colorsatfac" = 1     # \> 1 to boost color saturation.
"noise" = None        # Noise luminosity is determined automatically if None.
"noiselum" = 0.15     # Noise luminosity for single channel (between 0 - 1).
"noiselums" = {}      # Noise luminosity for each channel (between 0 - 1).
"noisesig" = 1        # Data noise level output to noiselum: measured sigma above the measured mean.
"noisesig0": 2        # Data noise level: measured sigma above the measured mean.
"correctbias"= 0      # Measure data noise mean (otherwise assume = 0).
"combine" = 'average' # "average" or "sum" combine images.
"invert" =  0         # Invert luminosity (black on white).
"scaling" = None      # Determined automatically if None: image scaling.
"bscale" = 1          # Multiply all images by this value.
"bzero" = 0           # Add this value to all images.
"samplesize" = 1000   # Determine number of levels.
"stampsize" = 1000    # Making final color image (memory issue).
"sampledx" = 0        # Offset in x direction.
"sampledy" = 0        # Offset in y direction.
"show" = 0            # Show final image at the end.
"showstamps" = 0      # Show image config stamps.
"showwith" = 'open'   # Command to display images.
"thumbnail" = None    # Show thumbnail.
"legend" = 0          # Adds legend to top-left corner indicating which filters were used (only for RGB).
"maxstampsize" = 6000 # Memory fix.
"testfirst" = 1       # Test some options before making the final image.
"deletetests" = 0     # Delete testing files.
"deletefilters" = 1   # Delete filter files.
```

## Resources

* [Lupton's method](http://www.astro.princeton.edu/~rhl/PrettyPictures/)
* [STIFF](https://github.com/astromatic/stiff)
