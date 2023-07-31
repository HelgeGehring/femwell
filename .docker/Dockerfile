FROM jupyter/scipy-notebook:latest

COPY . ${HOME}

RUN conda env update --file environment.yml --prune
RUN conda install julia

RUN /bin/sh -c JULIA_PROJECT="" julia -e "using Pkg; Pkg.add(\"IJulia\"); using IJulia; installkernel(\"Julia\", \"--project=${REPO_DIR}\");" && julia --project=${REPO_DIR} -e 'using Pkg; Pkg.instantiate(); Pkg.resolve(); pkg"precompile"'