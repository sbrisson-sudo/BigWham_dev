FROM ubuntu:focal
MAINTAINER Brice Lecampion <brice.lecampion@epfl.ch>

ENV DEBIAN_FRONTEND=noninteractive

# Make sure the package repository is up to date.
RUN apt-get -qq update && apt-get -qq -y install \
    g++ \
    gfortran \
    libtbb-dev \
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
RUN python3 -m pip install -Iv pygmsh==6.1.1
RUN python3 -m pip install -Iv meshio==4.0.11
RUN python3 -m pip install -Iv atomman==1.3.0
RUN python3 -m pip uninstall -y scipy
RUN python3 -m pip uninstall -y cython
RUN python3 -m pip install git+https://github.com/congma/scipy.git@ce61e02e912d15ea28b37ea036334b9ba266ebb5#egg=scipy

# install latest version of gmsh
WORKDIR /opt
RUN wget https://gmsh.info/bin/Linux/gmsh-4.5.4-Linux64.tgz
RUN tar -xf /opt/gmsh-4.5.4-Linux64.tgz
RUN ln -s /opt/gmsh-4.5.4-Linux64/bin/gmsh /bin/gmsh
