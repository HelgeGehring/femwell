FROM jupyter/scipy-notebook:latest

COPY environment.yml ${HOME}
COPY Project.toml ${HOME}
COPY Manifest.toml ${HOME}

RUN conda env update --file environment.yml --prune --quiet
RUN conda install julia --quiet

RUN /bin/sh -c JULIA_PROJECT="" julia -e "using Pkg; Pkg.add(\"IJulia\"); using IJulia; installkernel(\"Julia\", \"--project=.\");" && julia --project=${REPO_DIR} -e 'using Pkg; Pkg.instantiate(); Pkg.resolve(); pkg"precompile"'