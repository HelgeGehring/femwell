FROM python:3.11

RUN apt-get install -y libblas-dev liblapack-dev mpich

RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache notebook jupyterlab

ENV PETSC_CONFIGURE_OPTIONS="--with-scalar-type=complex"
RUN pip install --no-cache petsc &&\
    pip install --no-cache petsc4py &&\
    pip install --no-cache slepc &&\
    pip install --no-cache slepc4py

RUN pip install --no-cache git+https://github.com/kinnala/scikit-fem.git
RUN pip install --no-cache git+https://github.com/HelgeGehring/femwell.git

# create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}
USER ${USER}