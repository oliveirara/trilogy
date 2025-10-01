"""
Trilogy: Modern Python library for converting astronomical FITS images to color/grayscale images.

Based on the original trilogy by Dan Coe, modernized and optimized.
Uses Lupton's method for image scaling.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Literal

import astropy.io.fits as fits
import numpy as np
from numpy.typing import NDArray
from PIL import Image, ImageDraw
from scipy.optimize import golden

# Luminance weights (D65: red boosted, blue muted)
LUMINANCE_WEIGHTS = {"R": 0.212671, "G": 0.715160, "B": 0.072169}

ChannelType = Literal["R", "G", "B", "L"]
ModeType = Literal["RGB", "L"]
CombineMethod = Literal["average", "sum"]


@dataclass
class TrilogyConfig:
    """Configuration for Trilogy image processing."""

    # Input/Output
    indir: Path = field(default_factory=lambda: Path("."))
    outdir: Path = field(default_factory=lambda: Path("."))
    outname: str | None = None

    # Image scaling parameters
    satpercent: float = 0.001  # Percentage of pixels to saturate
    noiselum: float = 0.15  # Noise luminosity (0-1)
    noiselums: dict[str, float] = field(default_factory=dict)  # Per-channel noise
    noisesig: float = 1  # Noise sigma for output
    noisesig0: float = 2  # Noise sigma for measurement
    correctbias: bool = False  # Measure noise mean vs assume 0
    colorsatfac: float = 1  # Color saturation factor (>1 boosts)

    # Processing parameters
    samplesize: int = 1000  # Sample size for determining scaling
    sampledx: int = 0  # X offset for sampling
    sampledy: int = 0  # Y offset for sampling
    stampsize: int = 1000  # Stamp size for final image
    maxstampsize: int = 6000  # Maximum stamp size

    # Image manipulation
    invert: bool = False  # Invert luminosity
    combine: CombineMethod = "average"  # How to combine multiple images
    bscale: float = 1  # Multiply all images
    bzero: float = 0  # Add to all images

    # Advanced
    noise: float | None = None  # Manual noise level
    saturate: float | None = None  # Manual saturation level
    scaling: Path | None = None  # Load pre-computed scaling

    # Output options
    legend: bool = False  # Add legend to image
    
    # Auto-adjustment
    auto_adjust: bool = True  # Automatically adjust problematic parameters


@dataclass
class ImageStats:
    """Robust statistics for image data using sigma clipping."""

    mean: float
    std: float
    median: float
    min_val: float
    max_val: float
    percentile_99: float


def compute_robust_stats(
    data: NDArray[np.floating],
    n_sigma: float = 3.0,
    n_iterations: int = 5,
    presorted: bool = False,
) -> ImageStats:
    """
    Compute robust mean and std using sigma clipping.

    Args:
        data: Input data array
        n_sigma: Number of sigmas for clipping
        n_iterations: Number of clipping iterations
        presorted: Whether data is already sorted

    Returns:
        ImageStats with robust mean and standard deviation
    """
    data_sorted = data if presorted else np.sort(data.ravel())
    data_sorted = data_sorted[~np.isnan(data_sorted)]

    if len(data_sorted) == 0 or data_sorted[0] == data_sorted[-1]:
        return ImageStats(
            mean=0.0, 
            std=1.0, 
            median=0.0, 
            min_val=0.0, 
            max_val=0.0,
            percentile_99=0.0
        )

    ilo, ihi = 0, len(data_sorted)

    for _ in range(n_iterations):
        subset = data_sorted[ilo:ihi]
        imed = (ilo + ihi) // 2
        median_val = data_sorted[imed]

        # Use RMS instead of std for better robustness
        rms = np.sqrt(np.mean((subset - median_val) ** 2))

        lo_val = median_val - n_sigma * rms
        hi_val = median_val + n_sigma * rms

        new_ilo = np.searchsorted(data_sorted, lo_val)
        new_ihi = np.searchsorted(data_sorted, hi_val, side="right")

        if new_ihi - new_ilo == ihi - ilo:
            break

        ilo, ihi = new_ilo, new_ihi

    clipped = data_sorted[ilo:ihi]
    mean = float(np.mean(clipped))
    std = float(np.sqrt(np.mean((clipped - mean) ** 2)))
    median = float(np.median(clipped))
    
    # Additional stats for parameter estimation
    min_val = float(data_sorted[0])
    max_val = float(data_sorted[-1])
    percentile_99 = float(data_sorted[int(0.99 * len(data_sorted))])

    return ImageStats(
        mean=mean, 
        std=std, 
        median=median, 
        min_val=min_val, 
        max_val=max_val,
        percentile_99=percentile_99
    )


def determine_scaling(
    data: NDArray[np.floating],
    unsatpercent: float,
    noisesig: float = 1.0,
    correctbias: bool = False,
    noisesig0: float = 2.0,
) -> tuple[float, float, float]:
    """
    Determine data values (x0, x1, x2) to scale to (0, noiselum, 1).

    Args:
        data: Input image data
        unsatpercent: Percentile for saturation (1 - satpercent/100)
        noisesig: Sigma level for noise brightness
        correctbias: Whether to correct for bias
        noisesig0: Sigma level for bias measurement

    Returns:
        Tuple of (x0, x1, x2) scaling levels
    """
    data_flat = data.ravel()
    data_sorted = np.sort(data_flat)
    data_sorted[np.isnan(data_sorted)] = 0

    if data_sorted[0] == data_sorted[-1]:
        return (0.0, 1.0, 100.0)

    stats = compute_robust_stats(data_sorted, presorted=True)

    x0 = stats.mean - noisesig0 * stats.std if correctbias else 0.0
    x1 = stats.mean + noisesig * stats.std

    # Find saturation level
    idx = int(unsatpercent * len(data_sorted))
    idx = np.clip(idx, 0, len(data_sorted) - 1)
    x2 = float(data_sorted[idx])

    return (x0, x1, x2)


def auto_adjust_parameters(
    config: TrilogyConfig,
    data_stats: ImageStats,
    channel: str = "L",
) -> TrilogyConfig:
    """
    Automatically adjust parameters based on data statistics.
    
    This function detects problematic parameter combinations and suggests
    better values while preserving user intent.
    
    Args:
        config: Current configuration
        data_stats: Statistics from the sample data
        channel: Channel being processed
    
    Returns:
        Adjusted configuration
    """
    adjusted = False
    
    # Get current noiselum for this channel
    noiselum = config.noiselums.get(channel, config.noiselum)
    
    # Detect extreme noiselum values that cause numerical issues
    if noiselum > 1.5:
        print(f"  ⚠️  noiselum={noiselum:.2f} is very high, may cause numerical issues")
        print(f"  ℹ️  Auto-adjusting to noiselum=0.5 for better results")
        if config.noiselums:
            config.noiselums[channel] = 0.5
        else:
            config.noiselum = 0.5
        adjusted = True
    elif 0.95 < noiselum <= 1.5:
        # Values slightly above 1.0 often cause issues too
        print(f"  ⚠️  noiselum={noiselum:.2f} is close to 1.0, may cause numerical issues")
        print(f"  ℹ️  Auto-adjusting to noiselum=0.15 for better results")
        if config.noiselums:
            config.noiselums[channel] = 0.15
        else:
            config.noiselum = 0.15
        adjusted = True
    elif noiselum < 0.01:
        print(f"  ⚠️  noiselum={noiselum:.4f} is very low, image may be too bright")
        print(f"  ℹ️  Auto-adjusting to noiselum=0.15 for better results")
        if config.noiselums:
            config.noiselums[channel] = 0.15
        else:
            config.noiselum = 0.15
        adjusted = True
    
    # Detect extreme noisesig values
    if config.noisesig < 0.1:
        print(f"  ⚠️  noisesig={config.noisesig:.3f} is very low")
        print(f"  ℹ️  Auto-adjusting to noisesig=1.0 for better results")
        config.noisesig = 1.0
        adjusted = True
    elif config.noisesig > 100:
        print(f"  ⚠️  noisesig={config.noisesig:.1f} is very high")
        print(f"  ℹ️  Auto-adjusting to noisesig=50 for better results")
        config.noisesig = 50.0
        adjusted = True
    
    # Detect if data range is very small (all zeros or very uniform)
    data_range = data_stats.max_val - data_stats.min_val
    if data_range < 1e-10:
        print(f"  ⚠️  Data range is extremely small ({data_range:.2e})")
        print(f"  ℹ️  Image may appear blank - check input data")
    elif abs(data_stats.mean) < 1e-10 and abs(data_stats.std) < 1e-10:
        print(f"  ⚠️  Data appears to be all zeros or near-zero")
        print(f"  ℹ️  Auto-adjusting satpercent for better contrast")
        config.satpercent = 0.1  # More aggressive clipping for low-contrast data
        adjusted = True
    
    # Detect very high bscale values (often used incorrectly)
    if config.bscale < 0.001 and config.bscale != 1.0:
        print(f"  ⚠️  bscale={config.bscale:.6f} is very small")
        print(f"  ℹ️  This will make the image very dark")
        print(f"  💡 Consider using bscale=1.0 unless intentional")
    
    if adjusted:
        print(f"  ✓ Parameters auto-adjusted for optimal results")
    
    return config


def _scaling_objective(k: float, x0: float, x1: float, x2: float, n: float) -> float:
    """Objective function for golden section search in image scaling."""
    if abs(k) < 1e-10:
        return _scaling_objective(1e-10, x0, x1, x2, n)

    a1 = k * (x1 - x0) + 1
    a2 = k * (x2 - x0) + 1
    a1n = abs(a1**n)
    da1 = abs(a1n - a2)

    return da1 / abs(k)


def scale_image(
    data: NDArray[np.floating],
    levels: tuple[float, float, float],
    noiselum: float,
) -> NDArray[np.uint8]:
    """
    Scale image data to 0-255 using logarithmic scaling.

    Args:
        data: Input image data
        levels: Tuple of (x0, x1, x2) scaling levels
        noiselum: Target luminosity for noise level

    Returns:
        Scaled image as uint8 array
    """
    x0, x1, x2 = levels

    # Calculate scaling factor k
    if abs(noiselum - 0.5) < 1e-6:
        k = (x2 - 2 * x1 + x0) / max((x1 - x0) ** 2, 1e-30)
    else:
        n = 1.0 / noiselum
        # Use closure to pass parameters to objective
        try:
            k = abs(golden(lambda k_val: _scaling_objective(k_val, x0, x1, x2, n)))
        except (ValueError, RuntimeError):
            # Fallback to simple estimation if golden search fails
            # This can happen with extreme noiselum values
            k = 1.0 / max(x2 - x0, 1e-30)

    # Apply logarithmic scaling
    r1 = np.log10(k * (x2 - x0) + 1)
    if r1 <= 0:
        r1 = 1.0  # Avoid division by zero
    
    data_clipped = np.clip(data, 0, None)
    d = k * (data_clipped - x0) + 1
    d = np.clip(d, 1e-30, None)

    scaled = np.log10(d) / r1
    scaled = np.clip(scaled, 0, 1)
    scaled = (scaled * 255).astype(np.uint8)

    return scaled


def adjust_saturation(rgb: NDArray[np.floating], saturation: float) -> NDArray[np.floating]:
    """
    Adjust color saturation of RGB image.

    Args:
        rgb: RGB image array of shape (3, ny, nx)
        saturation: Saturation factor (>1 boosts, <1 reduces)

    Returns:
        Saturation-adjusted RGB array
    """
    if abs(saturation - 1.0) < 1e-6:
        return rgb

    # Build saturation adjustment matrix
    k = saturation
    rw, gw, bw = LUMINANCE_WEIGHTS["R"], LUMINANCE_WEIGHTS["G"], LUMINANCE_WEIGHTS["B"]

    sat_matrix = np.array(
        [
            [rw * (1 - k) + k, gw * (1 - k), bw * (1 - k)],
            [rw * (1 - k), gw * (1 - k) + k, bw * (1 - k)],
            [rw * (1 - k), gw * (1 - k), bw * (1 - k) + k],
        ]
    )

    # Reshape and apply transformation
    original_shape = rgb.shape
    rgb_flat = rgb.reshape(3, -1)
    rgb_adjusted = sat_matrix @ rgb_flat
    return rgb_adjusted.reshape(original_shape)


def rgb_to_image(rgb: NDArray[np.uint8], invert: bool = False) -> Image.Image:
    """
    Convert RGB array to PIL Image.

    Args:
        rgb: RGB array of shape (3, ny, nx) or (1, ny, nx) for grayscale
        invert: Whether to invert luminosity

    Returns:
        PIL Image
    """
    data = np.transpose(rgb, (1, 2, 0))  # (3, ny, nx) -> (ny, nx, 3)
    data = np.clip(data, 0, 255).astype(np.uint8)

    if invert:
        data = 255 - data

    n_channels = data.shape[-1]
    if n_channels == 3:
        img = Image.fromarray(data, mode="RGB")
    elif n_channels == 1:
        img = Image.fromarray(data[:, :, 0], mode="L")
    else:
        raise ValueError(f"Unexpected number of channels: {n_channels}")

    return img.transpose(Image.FLIP_TOP_BOTTOM)


class Trilogy:
    """Main class for creating color/grayscale images from FITS files."""

    def __init__(
        self,
        images: str | Path | list[str | Path] | dict[str, list[str | Path]] | None = None,
        config: TrilogyConfig | None = None,
        **kwargs,
    ):
        """
        Initialize Trilogy image processor.

        Args:
            images: Input images. Can be:
                - Single file path for grayscale
                - List of paths for grayscale (will be combined)
                - Dict with 'R', 'G', 'B' keys for color image
            config: Configuration object
            **kwargs: Additional config parameters as keyword arguments
        """
        # Initialize configuration
        if config is None:
            config = TrilogyConfig(**kwargs)
        else:
            # Merge kwargs into config
            for key, value in kwargs.items():
                if hasattr(config, key):
                    setattr(config, key, value)

        self.config = config
        self.images: dict[ChannelType, list[Path]] = {"R": [], "G": [], "B": [], "L": []}
        self.filters: dict[ChannelType, list[str]] = {"R": [], "G": [], "B": [], "L": []}
        self.mode: ModeType = "L"
        self.shape: tuple[int, int] | None = None  # (ny, nx)
        self.levels: dict[str, tuple[float, float, float]] = {}
        self._data_cache: dict[tuple[str, int, int, int, int], NDArray[np.floating]] = {}

        # Process input images
        if images is not None:
            self._process_images(images)

    def _process_images(self, images: str | Path | list[str | Path] | dict[str, list[str | Path]]) -> None:
        """Process and categorize input images."""
        if isinstance(images, (str, Path)):
            # Single image -> grayscale
            self.images["L"] = [self._normalize_path(images)]
            self.mode = "L"

        elif isinstance(images, list):
            # List of images -> grayscale
            self.images["L"] = [self._normalize_path(img) for img in images]
            self.mode = "L"

        elif isinstance(images, dict):
            # Dict with channels -> RGB
            self.mode = "RGB"
            for channel in ["R", "G", "B"]:
                if channel in images:
                    self.images[channel] = [self._normalize_path(img) for img in images[channel]]

    def _normalize_path(self, path: str | Path) -> Path:
        """Normalize image path and extract FITS extension if present."""
        path_str = str(path)

        # Handle FITS extensions like "image.fits[1]"
        if path_str.endswith("]"):
            # Keep the bracket notation for later
            pass

        # Add .fits extension if missing
        if not any(path_str.endswith(ext) for ext in [".fits", ".fits.gz", ".fits.fz", "]"]):
            path_str += ".fits"

        return Path(path_str)

    def _load_fits_data(self, filepath: Path) -> tuple[NDArray[np.floating], str]:
        """
        Load FITS image data and extract filter information.

        Args:
            filepath: Path to FITS file (may include [ext] notation)

        Returns:
            Tuple of (data array, filter name)
        """
        path_str = str(filepath)
        extension = None

        # Handle FITS extension notation
        if path_str.endswith("]"):
            ext_start = path_str.rfind("[")
            extension = int(path_str[ext_start + 1 : -1])
            path_str = path_str[:ext_start]

        full_path = self.config.indir / path_str

        with fits.open(full_path, memmap=True) as hdul:
            # Auto-detect extension with data if not specified
            if extension is None:
                extension = 0
                # If primary HDU has no data, try to find first extension with data
                if hdul[0].data is None:
                    for i in range(1, len(hdul)):
                        if hdul[i].data is not None:
                            extension = i
                            break
            
            hdu = hdul[extension]
            
            if hdu.data is None:
                raise ValueError(f"No image data found in {full_path}")
            
            data = hdu.data.astype(np.float64)

            # Extract filter information
            filter_name = hdu.header.get("FILTER", "")
            if not filter_name:
                filter1 = hdu.header.get("FILTER1", "")
                if filter1 and not filter1.startswith("CLEAR"):
                    filter_name = filter1
                else:
                    filter_name = hdu.header.get("FILTER2", "")

        return data, filter_name

    def _get_image_shape(self) -> None:
        """Determine output image shape from first image."""
        for channel in self.mode:
            if self.images[channel]:
                data, _ = self._load_fits_data(self.images[channel][0])
                self.shape = data.shape
                return

        raise ValueError("No images provided")

    def _load_channel_data(
        self,
        channel: ChannelType,
        yslice: slice,
        xslice: slice,
    ) -> NDArray[np.floating]:
        """
        Load and combine all images for a channel in the specified region.

        Args:
            channel: Channel to load ('R', 'G', 'B', or 'L')
            yslice: Y slice of image to load
            xslice: X slice of image to load

        Returns:
            Combined channel data
        """
        if not self.images[channel]:
            raise ValueError(f"No images for channel {channel}")

        # Check cache
        cache_key = (channel, yslice.start, yslice.stop, xslice.start, xslice.stop)
        if cache_key in self._data_cache:
            return self._data_cache[cache_key]

        combined = None
        filter_added = False

        for img_path in self.images[channel]:
            data, filt = self._load_fits_data(img_path)

            # Apply scaling
            data = self.config.bscale * data + self.config.bzero

            # Extract region
            region = data[yslice, xslice]

            if combined is None:
                combined = region
                if not filter_added:
                    self.filters[channel].append(filt)
                    filter_added = True
            else:
                combined = combined + region
                if not filter_added:
                    self.filters[channel].append(filt)

        # Average if requested
        if self.config.combine == "average":
            combined = combined / len(self.images[channel])

        # Cache result
        self._data_cache[cache_key] = combined

        return combined

    def _determine_levels(self) -> None:
        """Determine scaling levels for each channel based on sample region."""
        if self.shape is None:
            self._get_image_shape()

        ny, nx = self.shape
        yc, xc = ny // 2, nx // 2

        # Determine sample region
        sample_size = self.config.samplesize if self.config.samplesize > 0 else min(nx, ny, self.config.maxstampsize)

        y0 = max(0, yc - sample_size // 2 + self.config.sampledy)
        y1 = min(ny, yc + sample_size // 2 + self.config.sampledy)
        x0 = max(0, xc - sample_size // 2 + self.config.sampledx)
        x1 = min(nx, xc + sample_size // 2 + self.config.sampledx)

        print(f"Determining scaling from {x1-x0}x{y1-y0} sample region")
        if self.config.sampledx or self.config.sampledy:
            print(f"  offset by ({self.config.sampledx}, {self.config.sampledy})")

        unsatpercent = 1.0 - 0.01 * self.config.satpercent

        # Determine levels for each channel
        for channel in self.mode:
            data = self._load_channel_data(channel, slice(y0, y1), slice(x0, x1))

            # Get statistics for auto-adjustment
            data_stats = compute_robust_stats(data)
            
            # Auto-adjust parameters if needed (unless disabled by user)
            if self.config.auto_adjust:
                self.config = auto_adjust_parameters(self.config, data_stats, channel)

            if self.config.noise is not None and self.config.saturate is not None:
                # Manual levels
                self.levels[channel] = (0.0, self.config.noise, self.config.saturate)
            else:
                # Automatic levels
                self.levels[channel] = determine_scaling(
                    data,
                    unsatpercent,
                    noisesig=self.config.noisesig,
                    correctbias=self.config.correctbias,
                    noisesig0=self.config.noisesig0,
                )

            x0_lev, x1_lev, x2_lev = self.levels[channel]
            print(f"  {channel}: {x0_lev:.6f}  {x1_lev:.6f}  {x2_lev:.6f}")

    def _process_region(
        self,
        yslice: slice,
        xslice: slice,
    ) -> Image.Image:
        """
        Process a region of the image.

        Args:
            yslice: Y slice to process
            xslice: X slice to process

        Returns:
            PIL Image of the processed region
        """
        n_channels = len(self.mode)
        y0, y1 = yslice.start, yslice.stop
        x0, x1 = xslice.start, xslice.stop
        ny, nx = y1 - y0, x1 - x0

        # Load and scale each channel
        scaled = np.zeros((n_channels, ny, nx), dtype=np.uint8)

        noiselums = self.config.noiselums or {ch: self.config.noiselum for ch in self.mode}

        for i, channel in enumerate(self.mode):
            data = self._load_channel_data(channel, yslice, xslice)
            noiselum = noiselums.get(channel, self.config.noiselum)
            scaled[i] = scale_image(data, self.levels[channel], noiselum)

        # Adjust saturation for RGB
        if self.mode == "RGB" and abs(self.config.colorsatfac - 1.0) > 1e-6:
            scaled = adjust_saturation(scaled.astype(np.float64), self.config.colorsatfac)
            scaled = np.clip(scaled, 0, 255).astype(np.uint8)

        return rgb_to_image(scaled, invert=self.config.invert)

    def make_image(self) -> Image.Image:
        """
        Create the final image.

        Returns:
            PIL Image of the final result
        """
        if self.shape is None:
            self._get_image_shape()

        # Determine scaling levels
        if not self.levels:
            self._determine_levels()

        ny, nx = self.shape
        stamp_size = self.config.stampsize if self.config.stampsize > 0 else min(nx, ny, self.config.maxstampsize)

        # Create full image
        full_image = Image.new(self.mode, (nx, ny))

        print(f"\nCreating {self.mode} image...")

        # Process in stamps
        for y in range(0, ny, stamp_size):
            y1 = min(y + stamp_size, ny)
            for x in range(0, nx, stamp_size):
                x1 = min(x + stamp_size, nx)

                print(f"  Processing region ({x}, {y}) -> ({x1}, {y1})")

                region = self._process_region(slice(y, y1), slice(x, x1))

                # Paste into full image (flip vertically for astronomical convention)
                full_image.paste(region, (x, ny - y1, x1, ny - y))

        return full_image

    def save(self, output_path: str | Path | None = None) -> Path:
        """
        Create and save the image.

        Args:
            output_path: Output file path. If None, uses config.outname

        Returns:
            Path to saved file
        """
        if output_path is None:
            if self.config.outname is None:
                raise ValueError("No output path specified")
            output_path = self.config.outdir / f"{self.config.outname}.png"
        else:
            output_path = Path(output_path)

        img = self.make_image()

        # Add legend if requested
        if self.config.legend and self.mode == "RGB":
            img = self._add_legend(img)

        print(f"\nSaving to {output_path}")
        img.save(output_path)

        return output_path

    def _add_legend(self, img: Image.Image) -> Image.Image:
        """Add filter legend to image."""
        draw = ImageDraw.Draw(img)

        x, y = 20, 20
        dy = 15

        colors = {
            "R": (255, 0, 0),
            "G": (0, 255, 0),
            "B": (0, 128, 255),
        }

        for i, channel in enumerate(["B", "G", "R"]):
            filters = " + ".join(self.filters[channel])
            text = f"{channel} = {filters}"
            draw.text((x, y + i * dy), text, fill=colors[channel])

        return img

    def run(self) -> Image.Image:
        """
        Convenience method that creates and returns the image.

        Returns:
            PIL Image
        """
        return self.make_image()

