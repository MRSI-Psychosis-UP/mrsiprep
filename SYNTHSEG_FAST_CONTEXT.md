# SynthSeg + FAST Tissue Backend Context

SynthSeg + FAST replaced the experimental CAT12-like Python port as the
preferred tissue segmentation path.

## Rationale

The CAT12 replication work did not reach acceptable CSF agreement. Even when
seeded with CAT12 p0/p1/p2/p3, the Python AMAP/cleanup path topped out below
the desired CSF Dice range. SynthSeg handled interior and exterior CSF,
including ventricular CSF, more reliably. FAST then provided partial-volume
GM/WM/CSF maps from the original T1 intensities.

## Current Backend

Backend name:

```text
synthseg-fast
```

High-level steps:

1. Run `mri_synthseg --robust` on the raw T1.
2. Resample SynthSeg labels to the raw T1 grid with nearest-neighbor
   interpolation using `nibabel.processing.resample_from_to`.
3. Build a CSF/ventricle mask from SynthSeg labels:

```text
4   left lateral ventricle
5   left inferior lateral ventricle
14  third ventricle
15  fourth ventricle
24  CSF
43  right lateral ventricle
44  right inferior lateral ventricle
```

4. Run HD-BET on the raw T1.
5. Build the FAST working mask in Python:

```text
HD-BET brain mask OR SynthSeg CSF/ventricle mask
```

6. Mask the raw T1 with that combined mask using NumPy.
7. Run FSL FAST on the masked T1.
8. Zero all GM/WM/CSF probabilities where the native-space SynthSeg label is
   `0`.
9. Correct small GM/WM leakage into the SynthSeg exterior CSF layer: where
   SynthSeg label `24` overlaps FAST GM or WM probability `> 0.01`, set GM
   and WM to `0` and CSF to `1`.
10. Export FAST partial-volume maps as:

```text
p1 = GM = FAST pve_1
p2 = WM = FAST pve_2
p3 = CSF = FAST pve_0
```

## Manual Validation Result

The manual HD-BET + SynthSeg CSF/ventricle + FAST run produced better CSF maps
than the CAT12-like implementation. The Atropos equivalent was also tested and
had improved CSF relative to the CAT12-like path, but FAST remained the better
choice for partial-volume maps.

## Retired Code

The experimental CAT12-like Python port and vendored CAT12 AMAP/template assets
were removed. The `existing` backend still supports consuming precomputed
CAT12-style p1/p2/p3 maps from `derivatives/cat12` for compatibility.
