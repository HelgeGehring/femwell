FROM ghcr.io/mamba-org/micromamba:latest

COPY  --chmod=777 . ${HOME}

RUN micromamba install -y -n base -f environment.yml julia && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-cache notebook jupyterlab
RUN pip install .

RUN /bin/sh -c JULIA_PROJECT="" julia -e "using Pkg; Pkg.add(\"IJulia\"); using IJulia; installkernel(\"Julia\", \"--project=.\");" && julia --project=${REPO_DIR} -e 'using Pkg; Pkg.instantiate(); Pkg.resolve(); pkg"precompile"'