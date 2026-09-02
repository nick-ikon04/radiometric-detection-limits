"""Reproduce the three laboratory examples and optionally export a CSV."""

import argparse
import csv
from dataclasses import asdict
from pathlib import Path
from typing import Any

from .calculations import calc_peak_roi, calc_radiometry


def calculate_examples() -> list[dict[str, Any]]:
    """Return normalized result rows for all bundled examples."""
    rows: list[dict[str, Any]] = []

    for measurement_time_s, gross_counts in (
        (60, 2861),
        (600, 28612),
        (1800, 85863),
    ):
        result = calc_radiometry(
            background_counts=1652,
            background_time_s=1800,
            gross_counts=gross_counts,
            gross_time_s=measurement_time_s,
            efficiency=0.47,
            regulatory_level=200,
        )
        rows.append(
            {
                "example": f"surface_smear_{measurement_time_s}s",
                "activity_bq": result.activity_bq,
                "detection_limit_bq": result.activity_detection_limit_bq,
                "specific_activity_bq_per_kg": "",
                "specific_detection_limit_bq_per_kg": "",
                "signal_detected": result.signal_detected,
                "method_suitable": result.method_suitable,
            }
        )

    waste = calc_radiometry(
        background_counts=256,
        background_time_s=7200,
        gross_counts=779,
        gross_time_s=3600,
        efficiency=0.018,
        mass_kg=0.750,
        beta=0.001,
        confidence=0.998,
        regulatory_level=50,
    )
    rows.append(
        {
            "example": "radioactive_waste",
            "activity_bq": waste.activity_bq,
            "detection_limit_bq": waste.activity_detection_limit_bq,
            "specific_activity_bq_per_kg": waste.specific_activity_bq_per_kg,
            "specific_detection_limit_bq_per_kg": waste.specific_detection_limit_bq_per_kg,
            "signal_detected": waste.signal_detected,
            "method_suitable": waste.method_suitable,
        }
    )

    counts = dict(
        zip(
            range(378, 398),
            [22, 19, 21, 12, 15, 11, 13, 17, 20, 25, 26, 36, 20, 20, 14, 11, 13, 18, 18, 19],
        )
    )
    peak = calc_peak_roi(
        counts_by_channel=counts,
        peak_range=(384, 390),
        background_ranges=((380, 383), (391, 394)),
        measurement_time_s=50_000,
        efficiency=0.01315,
    )
    rows.append(
        {
            "example": "be7_roi",
            "activity_bq": peak.activity_bq,
            "detection_limit_bq": peak.activity_detection_limit_bq,
            "specific_activity_bq_per_kg": "",
            "specific_detection_limit_bq_per_kg": "",
            "signal_detected": peak.signal_detected,
            "method_suitable": "",
        }
    )
    return rows


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    """Write example result rows to a UTF-8 CSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, help="optional output CSV path")
    args = parser.parse_args()

    rows = calculate_examples()
    for row in rows:
        print(
            f"{row['example']:<24} "
            f"activity={float(row['activity_bq']):.6g} Bq  "
            f"detection_limit={float(row['detection_limit_bq']):.6g} Bq  "
            f"detected={row['signal_detected']}"
        )
    if args.csv:
        write_csv(rows, args.csv)
        print(f"Results written to {args.csv}")


if __name__ == "__main__":
    main()


