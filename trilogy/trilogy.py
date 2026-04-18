"""
Trilogy: Python library for converting astronomical FITS images to color/grayscale images.
"""

from __future__ import annotations

from dataclasses import dataclass, field
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

    When `presorted=True`, `data` must be 1-D, finite, and sorted ascending.
    """
    if presorted:
        data_sorted = data
    else:
        data_flat = data.ravel()
        data_flat = data_flat[np.isfinite(data_flat)]
        data_sorted = np.sort(data_flat) if len(data_flat) else data_flat

    if len(data_sorted) == 0 or data_sorted[0] == data_sorted[-1]:
        return ImageStats(mean=0.0, std=1.0, median=0.0, min_val=0.0, max_val=0.0, percentile_99=0.0)

    ilo, ihi = 0, len(data_sorted)

    for _ in range(n_iterations):
        subset = data_sorted[ilo:ihi]
        imed = (ilo + ihi) // 2
        median_val = data_sorted[imed]

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

    min_val = float(data_sorted[0])
    max_val = float(data_sorted[-1])
    percentile_99 = float(data_sorted[int(0.99 * len(data_sorted))])

    return ImageStats(mean=mean, std=std, median=median, min_val=min_val, max_val=max_val, percentile_99=percentile_99)


def _levels_from_stats(
    data_sorted: NDArray[np.floating],
    stats: ImageStats,
    unsatpercent: float,
    noisesig: float,
    correctbias: bool,
    noisesig0: float,
) -> tuple[float, float, float]:
    """Compute (x0, x1, x2) levels from pre-sorted data and already-computed stats."""
    if len(data_sorted) == 0 or data_sorted[0] == data_sorted[-1]:
        return (0.0, 1.0, 100.0)

    x0 = stats.mean - noisesig0 * stats.std if correctbias else 0.0
    x1 = stats.mean + noisesig * stats.std

    idx = int(unsatpercent * len(data_sorted))
    idx = int(np.clip(idx, 0, len(data_sorted) - 1))
    x2 = float(data_sorted[idx])

    return (x0, x1, x2)


def determine_scaling(
    data: NDArray[np.floating],
    unsatpercent: float,
    noisesig: float = 1.0,
    correctbias: bool = False,
    noisesig0: float = 2.0,
) -> tuple[float, float, float]:
    """Determine data values (x0, x1, x2) to scale to (0, noiselum, 1)."""
    data_flat = data.ravel()
    data_valid = data_flat[np.isfinite(data_flat)]

    total_pixels = len(data_flat)
    valid_pixels = len(data_valid)
    nan_percentage = ((total_pixels - valid_pixels) / max(total_pixels, 1)) * 100.0

    if nan_percentage > 10.0:
        print(f"  ⚠️  Warning: {nan_percentage:.1f}% of pixels are NaN/Inf")

    if valid_pixels == 0:
        print("  ❌ Error: All pixels are NaN/Inf - cannot determine scaling")
        return (0.0, 1.0, 100.0)

    data_sorted = np.sort(data_valid)
    stats = compute_robust_stats(data_sorted, presorted=True)
    return _levels_from_stats(data_sorted, stats, unsatpercent, noisesig, correctbias, noisesig0)


def auto_adjust_parameters(
    config: TrilogyConfig,
    data_stats: ImageStats,
    channel: str = "L",
) -> TrilogyConfig:
    """Automatically adjust problematic parameters while preserving user intent."""
    adjusted = False

    noiselum = config.noiselums.get(channel, config.noiselum)

    if noiselum > 1.5:
        print(f"  ⚠️  noiselum={noiselum:.2f} is very high, may cause numerical issues")
        print("  ℹ️  Auto-adjusting to noiselum=0.5 for better results")  # noqa: RUF001
        if config.noiselums:
            config.noiselums[channel] = 0.5
        else:
            config.noiselum = 0.5
        adjusted = True
    elif 0.95 < noiselum <= 1.5:
        print(f"  ⚠️  noiselum={noiselum:.2f} is close to 1.0, may cause numerical issues")
        print("  ℹ️  Auto-adjusting to noiselum=0.15 for better results")  # noqa: RUF001
        if config.noiselums:
            config.noiselums[channel] = 0.15
        else:
            config.noiselum = 0.15
        adjusted = True
    elif noiselum < 0.01:
        print(f"  ⚠️  noiselum={noiselum:.4f} is very low, image may be too bright")
        print("  ℹ️  Auto-adjusting to noiselum=0.15 for better results")  # noqa: RUF001
        if config.noiselums:
            config.noiselums[channel] = 0.15
        else:
            config.noiselum = 0.15
        adjusted = True

    if config.noisesig < 0.1:
        print(f"  ⚠️  noisesig={config.noisesig:.3f} is very low")
        print("  ℹ️  Auto-adjusting to noisesig=1.0 for better results")  # noqa: RUF001
        config.noisesig = 1.0
        adjusted = True
    elif config.noisesig > 100:
        print(f"  ⚠️  noisesig={config.noisesig:.1f} is very high")
        print("  ℹ️  Auto-adjusting to noisesig=50 for better results")  # noqa: RUF001
        config.noisesig = 50.0
        adjusted = True

    data_range = data_stats.max_val - data_stats.min_val
    if data_range < 1e-10:
        print(f"  ⚠️  Data range is extremely small ({data_range:.2e})")
        print("  ℹ️  Image may appear blank - check input data")  # noqa: RUF001
    elif abs(data_stats.mean) < 1e-10 and abs(data_stats.std) < 1e-10:
        print("  ⚠️  Data appears to be all zeros or near-zero")
        print("  ℹ️  Auto-adjusting satpercent for better contrast")  # noqa: RUF001
        config.satpercent = 0.1
        adjusted = True

    if config.bscale < 0.001 and config.bscale != 1.0:
        print(f"  ⚠️  bscale={config.bscale:.6f} is very small")
        print("  ℹ️  This will make the image very dark")  # noqa: RUF001
        print("  💡 Consider using bscale=1.0 unless intentional")

    if adjusted:
        print("  ✓ Parameters auto-adjusted for optimal results")

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


def compute_scaling_k(levels: tuple[float, float, float], noiselum: float) -> float:
    """Compute the scaling factor `k` for given levels and target noise luminosity."""
    x0, x1, x2 = levels

    if abs(noiselum - 0.5) < 1e-6:
        return (x2 - 2 * x1 + x0) / max((x1 - x0) ** 2, 1e-30)

    n = 1.0 / noiselum
    try:
        return abs(golden(lambda k_val: _scaling_objective(k_val, x0, x1, x2, n)))
    except (ValueError, RuntimeError):
        return 1.0 / max(x2 - x0, 1e-30)


def apply_scaling(
    data: NDArray[np.floating],
    levels: tuple[float, float, float],
    k: float,
) -> NDArray[np.uint8]:
    """Apply logarithmic scaling with a precomputed k. NaN/Inf pixels become 0."""
    x0, _x1, x2 = levels

    valid_mask = np.isfinite(data)
    scaled = np.zeros(data.shape, dtype=np.uint8)

    if not np.any(valid_mask):
        return scaled

    data_valid = data[valid_mask]

    r1 = np.log10(k * (x2 - x0) + 1)
    if r1 <= 0:
        r1 = 1.0

    data_clipped = np.clip(data_valid, 0, None)
    d = k * (data_clipped - x0) + 1
    np.clip(d, 1e-30, None, out=d)

    scaled_valid = np.log10(d) / r1
    np.clip(scaled_valid, 0, 1, out=scaled_valid)
    scaled[valid_mask] = (scaled_valid * 255).astype(np.uint8)

    return scaled


def scale_image(
    data: NDArray[np.floating],
    levels: tuple[float, float, float],
    noiselum: float,
) -> NDArray[np.uint8]:
    """Scale image data to 0-255 using logarithmic scaling."""
    k = compute_scaling_k(levels, noiselum)
    return apply_scaling(data, levels, k)


def adjust_saturation(rgb: NDArray[np.floating], saturation: float) -> NDArray[np.floating]:
    """Adjust color saturation of an RGB image of shape (3, ny, nx)."""
    if abs(saturation - 1.0) < 1e-6:
        return rgb

    k = saturation
    rw, gw, bw = LUMINANCE_WEIGHTS["R"], LUMINANCE_WEIGHTS["G"], LUMINANCE_WEIGHTS["B"]

    sat_matrix = np.array(
        [
            [rw * (1 - k) + k, gw * (1 - k), bw * (1 - k)],
            [rw * (1 - k), gw * (1 - k) + k, bw * (1 - k)],
            [rw * (1 - k), gw * (1 - k), bw * (1 - k) + k],
        ],
        dtype=np.float32,
    )

    original_shape = rgb.shape
    rgb_flat = rgb.reshape(3, -1)
    rgb_adjusted = sat_matrix @ rgb_flat
    return rgb_adjusted.reshape(original_shape)


def rgb_to_image(rgb: NDArray[np.uint8], invert: bool = False) -> Image.Image:
    """Convert (C, ny, nx) array to PIL Image, flipping vertically for astro convention."""
    data = np.transpose(rgb, (1, 2, 0))  # (C, ny, nx) -> (ny, nx, C)
    data = np.ascontiguousarray(np.clip(data, 0, 255).astype(np.uint8))

    if invert:
        data = 255 - data

    n_channels = data.shape[-1]
    if n_channels == 3:
        img = Image.fromarray(data, mode="RGB")
    elif n_channels == 1:
        img = Image.fromarray(data[:, :, 0], mode="L")
    else:
        raise ValueError(f"Unexpected number of channels: {n_channels}")

    return img.transpose(Image.Transpose.FLIP_TOP_BOTTOM)


class Trilogy:
    """Main class for creating color/grayscale images from FITS files."""

    def __init__(
        self,
        images: str | Path | list[str | Path] | dict[str, list[str | Path]] | None = None,
        config: TrilogyConfig | None = None,
        **kwargs,
    ):
        if config is None:
            config = TrilogyConfig(**kwargs)
        else:
            for key, value in kwargs.items():
                if hasattr(config, key):
                    setattr(config, key, value)

        self.config = config
        self.images: dict[ChannelType, list[Path]] = {"R": [], "G": [], "B": [], "L": []}
        self.filters: dict[ChannelType, list[str]] = {"R": [], "G": [], "B": [], "L": []}
        self.mode: ModeType = "L"
        self.shape: tuple[int, int] | None = None  # (ny, nx)
        self.levels: dict[str, tuple[float, float, float]] = {}
        self.ks: dict[str, float] = {}

        if images is not None:
            self._process_images(images)

    def _process_images(self, images: str | Path | list[str | Path] | dict[str, list[str | Path]]) -> None:
        if isinstance(images, str | Path):
            self.images["L"] = [self._normalize_path(images)]
            self.mode = "L"
        elif isinstance(images, list):
            self.images["L"] = [self._normalize_path(img) for img in images]
            self.mode = "L"
        elif isinstance(images, dict):
            self.mode = "RGB"
            for channel in ["R", "G", "B"]:
                if channel in images:
                    self.images[channel] = [self._normalize_path(img) for img in images[channel]]

    def _normalize_path(self, path: str | Path) -> Path:
        path_str = str(path)
        if not any(path_str.endswith(ext) for ext in [".fits", ".fits.gz", ".fits.fz", "]"]):
            path_str += ".fits"
        return Path(path_str)

    @staticmethod
    def _split_extension(path_str: str) -> tuple[str, int | None]:
        """Split a 'path.fits[ext]' string into (path, ext_or_None)."""
        if path_str.endswith("]"):
            ext_start = path_str.rfind("[")
            return path_str[:ext_start], int(path_str[ext_start + 1 : -1])
        return path_str, None

    def _load_fits_data(
        self,
        filepath: Path,
        yslice: slice | None = None,
        xslice: slice | None = None,
    ) -> tuple[NDArray[np.float32], str]:
        """
        Load a FITS region as float32 and extract the filter name.

        If `yslice`/`xslice` are given, only that region is materialized — the rest
        stays on disk via memmap. Without slices, the full image is loaded.
        """
        path_str, extension = self._split_extension(str(filepath))
        full_path = self.config.indir / path_str

        with fits.open(full_path, memmap=True) as hdul:
            if extension is None:
                extension = 0
                if hdul[0].data is None:
                    for i in range(1, len(hdul)):
                        if hdul[i].data is not None:
                            extension = i
                            break

            hdu = hdul[extension]

            if hdu.data is None:
                raise ValueError(f"No image data found in {full_path}")

            if yslice is None and xslice is None:
                region = np.asarray(hdu.data, dtype=np.float32)
            else:
                ys = yslice if yslice is not None else slice(None)
                xs = xslice if xslice is not None else slice(None)
                # Slicing a memmap materializes only the requested region.
                region = np.ascontiguousarray(hdu.data[ys, xs], dtype=np.float32)

            filter_name = hdu.header.get("FILTER", "")
            if not filter_name:
                filter1 = hdu.header.get("FILTER1", "")
                filter_name = filter1 if filter1 and not filter1.startswith("CLEAR") else hdu.header.get("FILTER2", "")

        return region, filter_name

    def _get_image_shape(self) -> None:
        """Determine output image shape from the first image's header."""
        for channel in self.mode:
            if self.images[channel]:
                path_str, extension = self._split_extension(str(self.images[channel][0]))
                full_path = self.config.indir / path_str
                with fits.open(full_path, memmap=True) as hdul:
                    if extension is None:
                        extension = 0
                        if hdul[0].data is None:
                            for i in range(1, len(hdul)):
                                if hdul[i].data is not None:
                                    extension = i
                                    break
                    hdr = hdul[extension].header
                    self.shape = (int(hdr["NAXIS2"]), int(hdr["NAXIS1"]))
                return

        raise ValueError("No images provided")

    def _load_channel_data(
        self,
        channel: ChannelType,
        yslice: slice,
        xslice: slice,
    ) -> NDArray[np.float32]:
        """Load and combine all images for a channel over the given region."""
        if not self.images[channel]:
            raise ValueError(f"No images for channel {channel}")

        populate_filters = not self.filters[channel]
        arrays: list[NDArray[np.float32]] = []

        for img_path in self.images[channel]:
            data, filt = self._load_fits_data(img_path, yslice, xslice)
            if self.config.bscale != 1.0 or self.config.bzero != 0.0:
                data = self.config.bscale * data + self.config.bzero
            arrays.append(data)
            if populate_filters:
                self.filters[channel].append(filt)

        if len(arrays) == 1:
            return arrays[0]

        stack = np.stack(arrays, axis=0)
        # Use nan-aware reduction so a single bad pixel in one frame doesn't
        # poison the combined result.
        reduce_fn = np.nanmean if self.config.combine == "average" else np.nansum
        combined = reduce_fn(stack, axis=0)
        return combined.astype(np.float32, copy=False)

    def _determine_levels(self) -> None:
        """Determine scaling levels and the per-channel `k` from a sample region."""
        if self.shape is None:
            self._get_image_shape()

        ny, nx = self.shape
        yc, xc = ny // 2, nx // 2

        sample_size = self.config.samplesize if self.config.samplesize > 0 else min(nx, ny, self.config.maxstampsize)

        y0 = max(0, yc - sample_size // 2 + self.config.sampledy)
        y1 = min(ny, yc + sample_size // 2 + self.config.sampledy)
        x0 = max(0, xc - sample_size // 2 + self.config.sampledx)
        x1 = min(nx, xc + sample_size // 2 + self.config.sampledx)

        print(f"Determining scaling from {x1 - x0}x{y1 - y0} sample region")
        if self.config.sampledx or self.config.sampledy:
            print(f"  offset by ({self.config.sampledx}, {self.config.sampledy})")

        unsatpercent = 1.0 - 0.01 * self.config.satpercent

        for channel in self.mode:
            data = self._load_channel_data(channel, slice(y0, y1), slice(x0, x1))

            # Sort once and share between robust-stats and level determination.
            data_flat = data.ravel()
            data_valid = data_flat[np.isfinite(data_flat)]

            total = len(data_flat)
            valid = len(data_valid)
            nan_pct = ((total - valid) / max(total, 1)) * 100.0
            if nan_pct > 10.0:
                print(f"  ⚠️  Warning: {nan_pct:.1f}% of pixels are NaN/Inf")

            if valid == 0:
                print("  ❌ Error: All pixels are NaN/Inf - cannot determine scaling")
                data_sorted = np.zeros(1, dtype=np.float32)
                stats = ImageStats(mean=0.0, std=1.0, median=0.0, min_val=0.0, max_val=0.0, percentile_99=0.0)
            else:
                data_sorted = np.sort(data_valid)
                stats = compute_robust_stats(data_sorted, presorted=True)

            if self.config.auto_adjust:
                self.config = auto_adjust_parameters(self.config, stats, channel)

            if self.config.noise is not None and self.config.saturate is not None:
                self.levels[channel] = (0.0, self.config.noise, self.config.saturate)
            else:
                self.levels[channel] = _levels_from_stats(
                    data_sorted,
                    stats,
                    unsatpercent,
                    noisesig=self.config.noisesig,
                    correctbias=self.config.correctbias,
                    noisesig0=self.config.noisesig0,
                )

            noiselum_ch = self.config.noiselums.get(channel, self.config.noiselum)
            self.ks[channel] = compute_scaling_k(self.levels[channel], noiselum_ch)

            x0_lev, x1_lev, x2_lev = self.levels[channel]
            print(f"  {channel}: {x0_lev:.6f}  {x1_lev:.6f}  {x2_lev:.6f}")

    def _process_region(
        self,
        yslice: slice,
        xslice: slice,
    ) -> NDArray[np.uint8]:
        """Process a region and return a (C, ny, nx) uint8 array."""
        n_channels = len(self.mode)
        ny = yslice.stop - yslice.start
        nx = xslice.stop - xslice.start

        scaled = np.zeros((n_channels, ny, nx), dtype=np.uint8)

        for i, channel in enumerate(self.mode):
            data = self._load_channel_data(channel, yslice, xslice)
            scaled[i] = apply_scaling(data, self.levels[channel], self.ks[channel])

        if self.mode == "RGB" and abs(self.config.colorsatfac - 1.0) > 1e-6:
            adjusted = adjust_saturation(scaled.astype(np.float32), self.config.colorsatfac)
            np.clip(adjusted, 0, 255, out=adjusted)
            scaled = adjusted.astype(np.uint8)

        return scaled

    def make_image(self) -> Image.Image:
        """Create the final image by processing tiles into a single array."""
        if self.shape is None:
            self._get_image_shape()

        if not self.levels:
            self._determine_levels()

        ny, nx = self.shape
        stamp_size = self.config.stampsize if self.config.stampsize > 0 else min(nx, ny, self.config.maxstampsize)
        n_channels = len(self.mode)

        # Assemble directly into a NumPy array; wrap as PIL once at the end.
        full_array = np.zeros((n_channels, ny, nx), dtype=np.uint8)

        print(f"\nCreating {self.mode} image...")

        for y in range(0, ny, stamp_size):
            y1 = min(y + stamp_size, ny)
            for x in range(0, nx, stamp_size):
                x1 = min(x + stamp_size, nx)
                print(f"  Processing region ({x}, {y}) -> ({x1}, {y1})")
                region = self._process_region(slice(y, y1), slice(x, x1))
                full_array[:, y:y1, x:x1] = region

        return rgb_to_image(full_array, invert=self.config.invert)

    def save(self, output_path: str | Path | None = None) -> Path:
        """Create and save the image."""
        if output_path is None:
            if self.config.outname is None:
                raise ValueError("No output path specified")
            output_path = self.config.outdir / f"{self.config.outname}.png"
        else:
            output_path = Path(output_path)

        img = self.make_image()

        if self.config.legend and self.mode == "RGB":
            img = self._add_legend(img)

        print(f"\nSaving to {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
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
        """Convenience method that creates and returns the image."""
        return self.make_image()
