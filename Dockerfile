FROM freesurfer/freesurfer:7.4.1

LABEL org.opencontainers.image.title="MRSIPrep"
LABEL org.opencontainers.image.description="BIDS App for preprocessing quantified whole-brain MRSI derivatives"
LABEL org.opencontainers.image.licenses="CHUV academic non-commercial research license"

ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    MPLBACKEND=Agg

RUN if command -v apt-get >/dev/null 2>&1; then \
        export DEBIAN_FRONTEND=noninteractive; \
        apt-get update; \
        apt-get install -y --no-install-recommends \
            build-essential ca-certificates git libgl1 libglib2.0-0 libgomp1 python3-pip; \
        rm -rf /var/lib/apt/lists/*; \
    elif command -v dnf >/dev/null 2>&1; then \
        dnf install -y \
            gcc gcc-c++ make ca-certificates git mesa-libGL glib2 libgomp python3-pip; \
        dnf clean all; \
    elif command -v yum >/dev/null 2>&1; then \
        yum install -y \
            gcc gcc-c++ make ca-certificates git mesa-libGL glib2 libgomp python3-pip; \
        yum clean all; \
    elif command -v microdnf >/dev/null 2>&1; then \
        microdnf install -y \
            gcc gcc-c++ make ca-certificates git mesa-libGL glib2 libgomp python3-pip; \
        microdnf clean all; \
    else \
        command -v python3 >/dev/null 2>&1; \
        command -v git >/dev/null 2>&1 || true; \
    fi

WORKDIR /opt/mrsiprep
COPY pyproject.toml README.md LICENSE ./
COPY mrsiprep ./mrsiprep

RUN python3 -m ensurepip --upgrade || true \
    && python3 -m pip install --upgrade pip \
    && python3 -m pip install ".[ants]"

ENTRYPOINT ["mrsiprep"]
