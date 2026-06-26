#!/usr/bin/env bash
set -euo pipefail

FSL_INSTALLER_URL="${FSL_INSTALLER_URL:-https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/releases/getfsl.sh}"
FREESURFER_DEB_URL="${FREESURFER_DEB_URL:-https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/8.2.0/freesurfer_ubuntu22-8.2.0_amd64.deb}"

printf 'Installing FSL into /opt/fsl...\n'
rm -rf /opt/fsl
curl -fL --retry 3 "${FSL_INSTALLER_URL}" -o /tmp/getfsl.sh
sh /tmp/getfsl.sh /opt/fsl
rm -f /tmp/getfsl.sh

printf 'Installing FreeSurfer Ubuntu package...\n'
curl -fL --retry 3 "${FREESURFER_DEB_URL}" -o /tmp/freesurfer.deb
apt-get update
apt-get install -y /tmp/freesurfer.deb
rm -f /tmp/freesurfer.deb
rm -rf /var/lib/apt/lists/*

fs_root=""
for candidate in /usr/local/freesurfer/8.2.0 /usr/local/freesurfer; do
  if [[ -x "${candidate}/bin/mri_synthseg" ]]; then
    fs_root="${candidate}"
    break
  fi
done

if [[ -z "${fs_root}" ]]; then
  printf 'FreeSurfer installed, but mri_synthseg was not found.\n' >&2
  exit 1
fi

rm -rf /opt/freesurfer
ln -s "${fs_root}" /opt/freesurfer

printf '\nFSL and FreeSurfer installed. Running verification...\n'
REQUIRE_FSL=1 REQUIRE_FREESURFER=1 REQUIRE_PETPVC=1 REQUIRE_CHIMERA=1 \
  /usr/local/bin/mrsiprep-check-neurodeps
