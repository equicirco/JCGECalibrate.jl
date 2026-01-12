# Usage

`JCGECalibrate` loads SAM data and derives calibrated parameters and starting values.

## Load SAM

```julia
using JCGECalibrate

sam = load_sam_table("path/to/sam.csv"; goods=["A","B"], factors=["K","L"])
```

## Calibrate

```julia
start = compute_starting_values(sam)
params = compute_calibration_params(sam, start)
```

## Notes

Calibration assumes the canonical SAM structure and consistent labels.

