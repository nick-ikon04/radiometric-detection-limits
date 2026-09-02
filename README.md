# Radiometric Detection Limits

[Українська версія](README_UA.md)

An educational Python project for calculating decision thresholds, detection
limits, confidence intervals, activity, and specific activity from radiometric
measurements. It also contains a side-background region-of-interest (ROI)
calculation for a gamma-spectrum peak.

## What the project demonstrates

- background, gross, and net count-rate calculation;
- decision threshold and detection limit under a normal approximation;
- confidence interval for the true net count rate;
- conversion from count rate to activity and specific activity;
- comparison of a detection limit with a regulatory level;
- gamma-peak net-area estimation using two side-background regions;
- explicit input validation and automated tests.

The implementation uses only the Python standard library.

## Main results

| Example | Activity | Detection limit | Conclusion |
|---|---:|---:|---|
| Surface smear, 60 s | 99.5012 Bq | 0.9791 Bq | Signal detected; method suitable |
| Surface smear, 600 s | 99.5083 Bq | 0.3289 Bq | Signal detected; method suitable |
| Surface smear, 1800 s | 99.5402 Bq | 0.2299 Bq | Signal detected; method suitable |
| Radioactive waste | 10.0463 Bq | 1.1423 Bq | 13.3951 Bq/kg; signal detected |
| Be-7 peak ROI | 0.08308 Bq | 0.08260 Bq | Peak detected by the approximate ROI criterion |

Machine-readable values are stored in
[`results/example_results.csv`](results/example_results.csv).

## Method overview

For background counts `N0` measured during `t0` and gross counts `Nb` measured
during `tb`, the net count rate is

```text
r_net = Nb / tb - N0 / t0
```

Activity is calculated from the detector response:

```text
A = r_net / (efficiency * gamma_yield)
```

For the peak ROI, side-background counts are scaled by the ratio of the peak
window width to the total background-window width. Poisson count variances are
then propagated through the subtraction.

## Run the examples

Python 3.10 or newer is required.

```bash
python -m radiometric_detection_limits.examples
python -m radiometric_detection_limits.examples --csv results/example_results.csv
```

No third-party runtime dependencies are required.

## Run the tests

```bash
python -m unittest discover -s tests -v
```

## Repository structure

```text
radiometric_detection_limits/  calculation library and examples
tests/                          automated numerical and validation tests
results/                        reproduced example results
README.md                       English documentation
README_UA.md                    Ukrainian documentation
```

## Limitations

The equations use a normal approximation to Poisson counting statistics. The
ROI decision and detection thresholds are simplified educational estimates.
This software has not been validated for regulatory, safety-critical, medical,
or dosimetric decisions. A real measurement procedure must follow the relevant
standard, calibration method, uncertainty budget, and detector model.

## License

MIT License. See [LICENSE](LICENSE).

