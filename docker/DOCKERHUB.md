# MRSIPrep

BIDS App for preprocessing already-quantified whole-brain MRSI maps:
biharmonic spike repair, MRSI↔T1w registration, tissue segmentation
(SynthSeg/FAST), partial-volume correction, T1w↔MNI normalization, and
Chimera/MNI-atlas regional profile extraction.

Full documentation, CLI options, and output layout:
https://github.com/MRSI-Psychosis-UP/mrsiprep

## Pull

```bash
docker pull fedlucchetti/mrsiprep:cpu
```

Bundles ANTs, FSL (FAST only), FreeSurfer (`recon-all`, `mri_synthseg`,
`mri_vol2vol`), PETPVC, and Chimera. No FreeSurfer license is included —
mount your own and set `FS_LICENSE`.

## Run

```bash
docker run --rm \
  -v /path/to/bids:/data:ro \
  -v /path/to/bids/derivatives:/out \
  -v /path/to/freesurfer/license.txt:/opt/freesurfer/license.txt:ro \
  -e FS_LICENSE=/opt/freesurfer/license.txt \
  fedlucchetti/mrsiprep:cpu \
  /data /out participant \
  --participant-label S001 \
  --session-label V1 \
  --mode light \
  --nthreads 8
```

`mrsiprep` does not perform spectral fitting; it expects quantified MRSI
maps, quality maps (CRLB/SNR/linewidth), and T1w images already in BIDS
layout. Use `--validate-only` to check inputs before an expensive batch run.

## License

CHUV academic non-commercial research license. See the LICENSE file in the
source repository.
