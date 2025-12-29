# HLSP Linelist Cross-Check Report

## Summary

- HLSP E140H linelist total lines: 1155
- Wavelength range: 1164.96 - 1647.23 Å

## Species Breakdown in HLSP

- Other: 1155

## Cross-Match with q-Bearing Atomic Table

- Atomic table lines: 351
- HLSP lines matched (within 0.02 Å): 29
- HLSP matches with q values: 29
- HLSP lines NOT in atomic table: 1126

## Key Finding: Pipeline Detection vs HLSP

The HLSP linelist represents lines that were independently identified in this spectrum.
When our fitter fails to detect a line that HLSP identifies, it suggests:
- Window positioning error (z_guess incorrect)
- Fit initialization finding wrong minimum
- Line is blended/weak for our fitting approach

## Example Lines HLSP Identifies but Fitter Misses

- line_0018: Ni V
- line_0019: Ni V
- line_0025: Ni V
- line_0041: Ni V
- line_0047: Ni V
- line_0048: Ni V
- line_0049: Ni V
- line_0054: Ni V
- line_0057: Ni V
- line_0061: Ni V
- line_0062: Fe V
- line_0065: Ni V
- line_0069: Fe V
- line_0070: Fe V
- line_0074: Ni V
- line_0078: Ni V
- line_0081: Fe V
- line_0082: Fe V
- line_0083: Ni V
- line_0090: Fe V
