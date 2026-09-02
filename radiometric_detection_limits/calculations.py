"""Core radiometric and gamma-spectrum ROI calculations.

The equations use a normal approximation to Poisson counting statistics.
They are intended for teaching and reproducible laboratory calculations, not
as a replacement for a validated regulatory measurement procedure.
"""

from dataclasses import dataclass
from math import sqrt
from statistics import NormalDist
from typing import Iterable, Mapping, Optional, Tuple

ChannelRange = Tuple[int, int]


def normal_quantile(probability: float) -> float:
    """Return the standard-normal quantile for ``0 < probability < 1``."""
    if not 0.0 < probability < 1.0:
        raise ValueError("probability must be between 0 and 1")
    return NormalDist().inv_cdf(probability)


@dataclass(frozen=True)
class RadiometryResult:
    """Calculated count-rate, decision-threshold, and activity quantities."""

    background_rate: float
    gross_rate: float
    net_rate: float
    k_alpha: float
    k_beta: float
    k_confidence: float
    decision_threshold_rate: float
    detection_limit_rate: float
    net_rate_lower: float
    net_rate_upper: float
    activity_bq: float
    activity_detection_limit_bq: float
    activity_lower_bq: float
    activity_upper_bq: float
    specific_activity_bq_per_kg: Optional[float]
    specific_detection_limit_bq_per_kg: Optional[float]
    signal_detected: bool
    method_suitable: Optional[bool]


@dataclass(frozen=True)
class PeakROIResult:
    """Calculated quantities for a peak and its side-background regions."""

    peak_range: ChannelRange
    background_ranges: Tuple[ChannelRange, ...]
    peak_width_channels: int
    background_width_channels: int
    gross_peak_counts: float
    side_background_counts: float
    background_under_peak_counts: float
    net_peak_counts: float
    net_peak_rate: float
    net_counts_standard_uncertainty: float
    net_rate_standard_uncertainty: float
    k_alpha: float
    k_beta: float
    k_confidence: float
    decision_threshold_counts: float
    detection_limit_counts: float
    decision_threshold_rate: float
    detection_limit_rate: float
    activity_bq: float
    activity_detection_limit_bq: float
    signal_detected: bool


def _validate_probability(value: float, name: str) -> None:
    if not 0.0 < value < 1.0:
        raise ValueError(f"{name} must be between 0 and 1")


def calc_radiometry(
    background_counts: float,
    background_time_s: float,
    gross_counts: float,
    gross_time_s: float,
    efficiency: float,
    *,
    gamma_yield: float = 1.0,
    mass_kg: Optional[float] = None,
    alpha: float = 0.05,
    beta: float = 0.05,
    confidence: float = 0.95,
    regulatory_level: Optional[float] = None,
) -> RadiometryResult:
    """Calculate radiometric detection limits and activity.

    ``regulatory_level`` is interpreted as Bq when ``mass_kg`` is omitted and
    as Bq/kg when a mass is supplied.
    """
    if background_time_s <= 0 or gross_time_s <= 0:
        raise ValueError("measurement times must be positive")
    if background_counts < 0 or gross_counts < 0:
        raise ValueError("counts cannot be negative")
    if efficiency <= 0 or gamma_yield <= 0:
        raise ValueError("efficiency and gamma_yield must be positive")
    if mass_kg is not None and mass_kg <= 0:
        raise ValueError("mass_kg must be positive")
    if regulatory_level is not None and regulatory_level < 0:
        raise ValueError("regulatory_level cannot be negative")
    _validate_probability(alpha, "alpha")
    _validate_probability(beta, "beta")
    _validate_probability(confidence, "confidence")

    background_rate = background_counts / background_time_s
    gross_rate = gross_counts / gross_time_s
    net_rate = gross_rate - background_rate

    k_alpha = normal_quantile(1.0 - alpha)
    k_beta = normal_quantile(1.0 - beta)
    k_confidence = normal_quantile((1.0 + confidence) / 2.0)

    time_factor = 1.0 / gross_time_s + 1.0 / background_time_s
    decision_threshold_rate = (k_alpha**2 / (2.0 * background_time_s)) * (
        1.0
        + sqrt(
            1.0
            + background_rate
            * (2.0 * background_time_s / k_alpha) ** 2
            * time_factor
        )
    )
    detection_limit_rate = (
        0.25 * (k_alpha + k_beta) ** 2 * time_factor
        + (k_alpha + k_beta) * sqrt(background_rate * time_factor)
    )

    net_rate_uncertainty = sqrt(
        gross_rate / gross_time_s + background_rate / background_time_s
    )
    net_rate_lower = net_rate - k_confidence * net_rate_uncertainty
    net_rate_upper = net_rate + k_confidence * net_rate_uncertainty

    response = efficiency * gamma_yield
    activity_bq = net_rate / response
    activity_detection_limit_bq = detection_limit_rate / response
    activity_lower_bq = net_rate_lower / response
    activity_upper_bq = net_rate_upper / response

    specific_activity = None
    specific_detection_limit = None
    if mass_kg is not None:
        specific_activity = activity_bq / mass_kg
        specific_detection_limit = activity_detection_limit_bq / mass_kg

    method_suitable = None
    if regulatory_level is not None:
        quantity = (
            activity_detection_limit_bq
            if mass_kg is None
            else specific_detection_limit
        )
        method_suitable = quantity is not None and quantity < regulatory_level

    return RadiometryResult(
        background_rate=background_rate,
        gross_rate=gross_rate,
        net_rate=net_rate,
        k_alpha=k_alpha,
        k_beta=k_beta,
        k_confidence=k_confidence,
        decision_threshold_rate=decision_threshold_rate,
        detection_limit_rate=detection_limit_rate,
        net_rate_lower=net_rate_lower,
        net_rate_upper=net_rate_upper,
        activity_bq=activity_bq,
        activity_detection_limit_bq=activity_detection_limit_bq,
        activity_lower_bq=activity_lower_bq,
        activity_upper_bq=activity_upper_bq,
        specific_activity_bq_per_kg=specific_activity,
        specific_detection_limit_bq_per_kg=specific_detection_limit,
        signal_detected=net_rate > decision_threshold_rate,
        method_suitable=method_suitable,
    )


def _range_width(channel_range: ChannelRange) -> int:
    start, end = channel_range
    if start > end:
        raise ValueError("channel range start cannot exceed its end")
    return end - start + 1


def _sum_range(counts: Mapping[int, float], channel_range: ChannelRange) -> float:
    start, end = channel_range
    return sum(counts.get(channel, 0.0) for channel in range(start, end + 1))


def _contains(channel_range: ChannelRange, channel: int) -> bool:
    return channel_range[0] <= channel <= channel_range[1]


def calc_peak_roi(
    counts_by_channel: Mapping[int, float],
    peak_range: ChannelRange,
    background_ranges: Iterable[ChannelRange],
    measurement_time_s: float,
    efficiency: float,
    *,
    gamma_yield: float = 1.0,
    alpha: float = 0.025,
    beta: float = 0.025,
    confidence: float = 0.95,
) -> PeakROIResult:
    """Estimate a net gamma-peak area using scaled side-background regions."""
    if measurement_time_s <= 0:
        raise ValueError("measurement_time_s must be positive")
    if efficiency <= 0 or gamma_yield <= 0:
        raise ValueError("efficiency and gamma_yield must be positive")
    if any(value < 0 for value in counts_by_channel.values()):
        raise ValueError("channel counts cannot be negative")
    _validate_probability(alpha, "alpha")
    _validate_probability(beta, "beta")
    _validate_probability(confidence, "confidence")

    backgrounds = tuple(background_ranges)
    if not backgrounds:
        raise ValueError("at least one background range is required")

    peak_width = _range_width(peak_range)
    background_width = sum(_range_width(item) for item in backgrounds)
    for channel in range(peak_range[0], peak_range[1] + 1):
        if any(_contains(item, channel) for item in backgrounds):
            raise ValueError("peak and background ranges cannot overlap")

    gross_peak = _sum_range(counts_by_channel, peak_range)
    side_background = sum(_sum_range(counts_by_channel, item) for item in backgrounds)
    scale = peak_width / background_width
    background_under_peak = scale * side_background
    net_peak = gross_peak - background_under_peak

    uncertainty_counts = sqrt(gross_peak + scale**2 * side_background)
    net_peak_rate = net_peak / measurement_time_s
    uncertainty_rate = uncertainty_counts / measurement_time_s

    k_alpha = normal_quantile(1.0 - alpha)
    k_beta = normal_quantile(1.0 - beta)
    k_confidence = normal_quantile((1.0 + confidence) / 2.0)
    null_uncertainty = sqrt(background_under_peak + scale**2 * side_background)
    decision_threshold_counts = k_alpha * null_uncertainty
    detection_limit_counts = (k_alpha + k_beta) * null_uncertainty
    decision_threshold_rate = decision_threshold_counts / measurement_time_s
    detection_limit_rate = detection_limit_counts / measurement_time_s

    response = efficiency * gamma_yield
    return PeakROIResult(
        peak_range=peak_range,
        background_ranges=backgrounds,
        peak_width_channels=peak_width,
        background_width_channels=background_width,
        gross_peak_counts=gross_peak,
        side_background_counts=side_background,
        background_under_peak_counts=background_under_peak,
        net_peak_counts=net_peak,
        net_peak_rate=net_peak_rate,
        net_counts_standard_uncertainty=uncertainty_counts,
        net_rate_standard_uncertainty=uncertainty_rate,
        k_alpha=k_alpha,
        k_beta=k_beta,
        k_confidence=k_confidence,
        decision_threshold_counts=decision_threshold_counts,
        detection_limit_counts=detection_limit_counts,
        decision_threshold_rate=decision_threshold_rate,
        detection_limit_rate=detection_limit_rate,
        activity_bq=net_peak_rate / response,
        activity_detection_limit_bq=detection_limit_rate / response,
        signal_detected=net_peak > decision_threshold_counts,
    )


