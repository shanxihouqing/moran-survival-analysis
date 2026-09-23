# Moran's I habitat analysis

## Overview

This project uses local spatial autocorrelation to characterize intratumoural heterogeneity in contrast-enhanced CT. Each voxel inside the tumour region of interest is classified according to its intensity relative to the tumour mean and the similarity between that voxel and its spatial neighbours.

The resulting habitat map contains four spatial-association patterns plus a non-significant category:

- High–High (HH): a high-intensity voxel surrounded by high-intensity neighbours.
- Low–Low (LL): a low-intensity voxel surrounded by low-intensity neighbours.
- High–Low (HL): a high-intensity voxel surrounded by low-intensity neighbours.
- Low–High (LH): a low-intensity voxel surrounded by high-intensity neighbours.
- Non-significant: a voxel that does not meet the prespecified spatial-association criterion.

## Local Moran's I

For a tumour region containing `n` voxels, let `x_i` denote the intensity of voxel `i`, and let `x_bar` denote the mean intensity within the tumour. The variance term is

```text
M2 = (1 / n) * sum((x_i - x_bar)^2).
```

The local Moran statistic for voxel `i` is

```text
I_i = ((x_i - x_bar) / M2) * sum_j(w_ij * (x_j - x_bar)),
```

where `w_ij` is the spatial weight between voxels `i` and `j`. In the referenced implementation, voxels within a three-dimensional neighbourhood distance of `sqrt(3)` are treated as neighbours (`w_ij = 1`); otherwise, `w_ij = 0`.

A positive `I_i` indicates local similarity, whereas a negative `I_i` indicates local dissimilarity. Classification combines the sign of `I_i` with whether `x_i` is above or below `x_bar`.

| Pattern | Local association | Voxel intensity | Interpretation |
|---|---|---|---|
| High–High | `I_i > 0` | `x_i > x_bar` | Spatial cluster of high CT intensity |
| Low–Low | `I_i > 0` | `x_i < x_bar` | Spatial cluster of low CT intensity |
| High–Low | `I_i < 0` | `x_i > x_bar` | High-intensity spatial outlier |
| Low–High | `I_i < 0` | `x_i < x_bar` | Low-intensity spatial outlier |

Voxels failing the prespecified association or significance criterion are assigned to the non-significant class. The exact threshold and inferential procedure must be reported with the study. If statistical significance is claimed, a permutation-derived p-value or another explicitly defined inferential method should be used rather than interpreting the magnitude of `I_i` alone as a p-value.

## Potential interpretation in contrast-enhanced CT

These patterns describe image intensity and spatial organization; they are not direct histopathological measurements.

- High–High regions may correspond to spatially coherent enhancing tissue, potentially reflecting relatively high perfusion, vascularity, inflammation, or contrast-agent concentration.
- Low–Low regions may correspond to spatially coherent low-enhancement tissue, potentially reflecting necrosis, low perfusion, ischaemia, cystic change, or fluid-containing regions.
- High–Low regions represent focal high-intensity outliers within a lower-intensity neighbourhood and may occur near enhancing boundaries, small vessels, or focal contrast accumulation.
- Low–High regions represent focal low-intensity outliers within a higher-intensity neighbourhood and may occur in non-enhancing foci surrounded by enhancing tissue, including possible necrotic or fluid-rich centres.

These interpretations should be presented as biological hypotheses. Confirmation requires correlation with acquisition phase, pathology, perfusion measurements, treatment information, and other clinical data.

## Outputs

The workflow produces a NIfTI habitat mask that preserves the source image geometry and assigns each tumour voxel to HH, LL, HL, LH, or the non-significant category. Global Moran's I may additionally be exported as a patient-level measure of overall spatial autocorrelation.

For reproducibility, a public analysis should report:

1. CT acquisition and enhancement phase.
2. Image resampling and intensity preprocessing.
3. Tumour segmentation procedure.
4. Neighbourhood definition and spatial-weight normalization.
5. Local-association threshold or permutation-testing procedure.
6. Habitat label encoding in the exported NIfTI files.
7. Software package and version used to calculate Moran's I.


