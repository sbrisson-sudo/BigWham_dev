FROM ubuntu:jammy
LABEL maintainer="Brice Lecampion <brice.lecampion@epfl.ch>"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y wget cpio sudo
RUN wget -O - https://registrationcenter-download.intel.com/akdlm/IRC_NAS/992857b9-624c-45de-9701-f6445d845359/l_BaseKit_p_2023.2.0.49397.sh > script.sh

# RUN wget -O- https://registrationcenter-download.intel.com/akdlm/IRC_NAS/992858b9-624c-45de-9701-f6445d845359/l_BaseKit_p_2023.2.0.49397.sh > script.sh
RUN chmod +x script.sh && ./script.sh -a -s --eula accept
ENV ONEAPI_ROOT=/opt/intel/oneapi
ENV MKLROOT=/opt/intel/oneapi
ENV PATH=$ONEAPI_ROOT/compiler/latest/linux/bin/intel64:$PATH
ENV LD_LIBRARY_PATH=$ONEAPI_ROOT/compiler/latest/linux/lib/intel64_lin:$LD_LIBRARY_PATH


# RUN apt-get -qq update && apt-get -qq -y install wget
# RUN wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/992857b9-624c-45de-9701-f6445d845359/l_BaseKit_p_2023.2.0.49397.sh
# RUN sh ./l_BaseKit_p_2023.2.0.49397.sh
# RUN wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \ | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
# RUN echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list
# RUN apt update


# Make sure the package repository is up to date.
RUN apt-get -qq update && apt-get -qq -y install \
    g++ \
    gfortran \
    libtbb-dev \
    # intel-oneapi-mkl-devel \
    bc \
    libhugetlbfs-bin \
    libboost-all-dev \
    libmetis-dev \
    libhdf5-serial-dev \
    libfftw3-dev \
    libopenblas-dev \
    libeigen3-dev \
    libgtest-dev\
    googletest \
    valgrind \
    parallel \
    gdb \
    cmake-curses-gui \
    libgsl-dev \
    libbz2-dev \
    python3-pip \
    vim \
    git \
    bash-completion \
    build-essential \
    libssl-dev \
    libffi-dev \
    python3-dev \
    wget \
    less \
    net-tools \
    openssh-server \
    tmux \
    ranger \
    libglu1-mesa \
    libgl1 \
    libxrender1 \
    libxcursor1 \
    libxfixes3 \
    libxext6 \
    libxft2 \
    libfontconfig1 \
    libxinerama1 \
    libx11-6 \
    emacs \
    libpthread-stubs0-dev \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install pip cmake ipython setuptools wheel -U --force
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install --upgrade setuptools
RUN python3 -m pip install --upgrade wheel
RUN python3 -m pip install numpy matplotlib pandas cython -U --force
RUN python3 -m pip install jupyter scipy matplotlib numpy jupyterlab
RUN python3 -m pip install mpmath
RUN python3 -m pip install cmake
RUN python3 -m pip install python-lsp-server
RUN python3 -m pip install py-ubjson
RUN python3 -m pip install dataclasses-json
RUN python3 -m pip install numba
RUN python3 -m pip install meshio
RUN python3 -m pip install numpy
RUN python3 -m pip install dill
RUN python3 -m pip install pygmsh
RUN python3 -m pip install multimethod

# install latest version of gmsh
WORKDIR /opt
RUN wget https://gmsh.info/bin/Linux/gmsh-4.5.4-Linux64.tgz
RUN tar -xf /opt/gmsh-4.5.4-Linux64.tgz
RUN ln -s /opt/gmsh-4.5.4-Linux64/bin/gmsh /bin/gmsh
