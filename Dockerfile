FROM freesurfer/freesurfer:7.4.1

LABEL org.opencontainers.image.title="MRSIPrep"
LABEL org.opencontainers.image.description="BIDS App for preprocessing quantified whole-brain MRSI derivatives"
LABEL org.opencontainers.image.licenses="CHUV academic non-commercial research license"

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    MPLBACKEND=Agg

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        git \
        libgl1 \
        libglib2.0-0 \
        libgomp1 \
        python3-pip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/mrsiprep
COPY pyproject.toml README.md LICENSE ./
COPY mrsiprep ./mrsiprep

RUN python3 -m pip install --upgrade pip \
    && python3 -m pip install ".[ants]"

ENTRYPOINT ["mrsiprep"]
