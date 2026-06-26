#!/usr/bin/env bash
set -euo pipefail

BOOTSTRAP_IMAGE="${BOOTSTRAP_IMAGE:-mrsiprep-bootstrap:ubuntu22.04-cpu}"
CONTAINER_NAME="${DEPS_CONTAINER:-mrsiprep-deps-install}"

if docker container inspect "${CONTAINER_NAME}" >/dev/null 2>&1; then
  exec docker start -ai "${CONTAINER_NAME}"
fi

printf 'Inside the container, run:\n'
printf '  bash /root/manual_install_fsl_freesurfer.sh\n'
printf 'Then exit and run docker/finalize_manual_deps.sh\n\n'

exec docker run -it \
  --name "${CONTAINER_NAME}" \
  --entrypoint /bin/bash \
  "${BOOTSTRAP_IMAGE}"
