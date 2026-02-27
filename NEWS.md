# experdyn (development version)

# experdyn 0.1.0

## New features

* `resample_timeseries()` — converts irregularly sampled ESM/EMA time series to
  a regular grid using natural cubic spline interpolation. Automatically detects
  and splits on large temporal gaps (e.g., overnight breaks) before interpolating.
