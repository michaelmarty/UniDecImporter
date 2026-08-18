"""Small numerical helpers used by the standalone importers.

These functions intentionally cover only the operations required by this package. Keeping
them here avoids importing the much larger :mod:`unidec.tools` module.
"""

from __future__ import annotations

import numpy as np


def isempty(value) -> bool:
    """Return ``True`` when *value* is ``None`` or has no elements."""
    return value is None or np.asarray(value, dtype=object).size == 0


def safedivide(numerator, denominator):
    """Divide arrays while returning zero where the denominator is zero."""
    numerator = np.asarray(numerator)
    denominator = np.asarray(denominator)
    shape = np.broadcast_shapes(numerator.shape, denominator.shape)
    result = np.zeros(shape, dtype=np.result_type(numerator, denominator, float))
    return np.divide(numerator, denominator, out=result, where=denominator != 0)


def nearestunsorted(values, target) -> int:
    """Return the index of the value nearest to *target*."""
    return int(np.argmin(np.abs(np.asarray(values) - target)))


def nonlinear_axis(start: float, end: float, resolution: float) -> np.ndarray:
    """Create an m/z axis with approximately constant resolving power."""
    start, end, resolution = float(start), float(end), float(resolution)
    if start <= 0 or end < start or resolution <= 0:
        raise ValueError("start and resolution must be positive and end must be >= start")
    values = [start]
    current = start + start / resolution
    while current < end:
        values.append(current)
        current += current / resolution
    return np.asarray(values)


def mergedata(template: np.ndarray, data: np.ndarray) -> np.ndarray:
    """Interpolate two-column *data* onto a two-column template axis."""
    template = np.asarray(template)
    data = np.asarray(data)
    if template.ndim != 2 or template.shape[1] < 1:
        raise ValueError("template must be an N x 2 array")
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("data must be an N x 2 array")
    order = np.argsort(data[:, 0])
    intensities = np.interp(template[:, 0], data[order, 0], data[order, 1], left=0, right=0)
    return np.column_stack((template[:, 0], intensities))


def lintegrate(data: np.ndarray, axis: np.ndarray, fastmode: bool = False) -> np.ndarray:
    """Conservatively bin intensities from two-column *data* onto *axis*."""
    del fastmode
    data = np.asarray(data)
    axis = np.asarray(axis)
    if axis.ndim != 1 or axis.size == 0:
        raise ValueError("axis must be a non-empty one-dimensional array")
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("data must be an N x 2 array")
    if axis.size == 1:
        return np.column_stack((axis, [np.sum(data[:, 1])]))
    midpoints = axis[:-1] + np.diff(axis) / 2
    first = axis[0] - (midpoints[0] - axis[0])
    last = axis[-1] + (axis[-1] - midpoints[-1])
    edges = np.concatenate(([first], midpoints, [last]))
    intensity, _ = np.histogram(data[:, 0], bins=edges, weights=data[:, 1])
    return np.column_stack((axis, intensity))


def get_autocorr_ratio(data: np.ndarray) -> float:
    """Return the lag-one/lag-zero intensity autocorrelation ratio."""
    intensities = np.asarray(data, dtype=float)[:, 1]
    if intensities.size < 2:
        return 0.0
    denominator = float(np.dot(intensities, intensities))
    if denominator == 0:
        return 0.0
    return float(np.dot(intensities[:-1], intensities[1:]) / denominator)


def data_extract(data: np.ndarray, mz: float, method: int, window: float | None = None):
    """Extract a peak height (method 1) or its m/z position (method 4)."""
    data = np.asarray(data)
    if data.ndim != 2 or data.shape[1] < 2 or len(data) == 0:
        return 0
    candidates = (np.arange(len(data)) if window is None else
                  np.flatnonzero(np.abs(data[:, 0] - mz) <= window))
    if candidates.size == 0:
        return 0
    peak_index = candidates[np.argmax(data[candidates, 1])]
    if method == 1:
        return data[peak_index, 1]
    if method == 4:
        return data[peak_index, 0]
    raise ValueError(f"unsupported extraction method: {method}")

