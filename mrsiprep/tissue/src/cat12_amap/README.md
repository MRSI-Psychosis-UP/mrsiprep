# CAT12 AMAP C sources

This folder vendors the CAT12 C sources needed by the `cat_amap` AMAP/PVE
segmentation step.

Source checkout:

`/home/flucchetti/Connectome/Dev/cat12`

Copied files:

- `cat_amap.c`: original MATLAB MEX entry point kept for reference
- `Amap.c`
- `Amap.h`
- `Kmeans.c`
- `MrfPrior.c`
- `Pve.c`
- `vollib.c`
- `COPYING.CAT12`

The original CAT12 compile hint is:

```matlab
mex -O cat_amap.c Kmeans.c Amap.c MrfPrior.c Pve.c vollib.c
```

`cat_amap_core.c` and `cat_amap_core.h` are local non-MEX wrapper files for a
future Python extension. They preserve the same core AMAP call sequence without
depending on MATLAB headers.

CAT12 is GPL-2-or-later. Any distributed binary or Python package that links or
vendors this code needs to respect that license boundary.
