#!/usr/bin/env bash
set -euo pipefail

: "${FS_LICENSE:?Set FS_LICENSE to your FreeSurfer license file before running this script.}"

docker build -t mrsiprep:latest /home/flucchetti/Connectome/Dev/mrsiprep

docker run --rm \
  -v /home/flucchetti/Connectome/BIDS/Dummy-Project:/data:ro \
  -v /home/flucchetti/Connectome/BIDS/Dummy-Project/derivatives:/out \
  -v "${FS_LICENSE}:/opt/freesurfer/license.txt:ro" \
  -e FS_LICENSE=/opt/freesurfer/license.txt \
  mrsiprep:latest \
  /data /out participant \
  --participant-label CHUVUP013 \
  --session-label V1 \
  --tissue-backend freesurfer \
  --registration-t1-target brain-csf \
  --parcellation-mode mni \
  --nthreads 4 \
  --verbose
