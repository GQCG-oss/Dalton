FROM ubuntu:18.04

RUN apt-get update && apt-get install -y wget \
    build-essential \
    git \ 
    gcc \
    g++ \
    gfortran \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && apt-get autoremove -y

ENV PATH="/usr/local/miniconda3/bin:${PATH}"
ARG PATH="/usr/local/miniconda3/bin:${PATH}"
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh \
    && bash Miniconda3-py37_4.8.2-Linux-x86_64.sh -p /usr/local/miniconda3 -b \
    && rm -f Miniconda3-py37_4.8.2-Linux-x86_64.sh
RUN conda --version

RUN conda install -c conda-forge openblas=0.3.10 lapack=3.6.1 openmpi-mpicc=4.0.5 openmpi-mpifort=4.0.5 hdf5=1.12.0 cmake=3.18.2 gcc_linux-64=7.5 gxx_linux-64=7.5 gfortran_linux-64=7.5

ENV DALTON_LAUNCHER="mpirun -n 4 --allow-run-as-root"
ARG DALTON_LAUNCHER="mpirun -n 4 --allow-run-as-root"
# Debugging
RUN ls 
RUN python setup --mpi --explicit-libs="-lopenblas -llapack" --prefix="/usr/local/miniconda3"
RUN cd build && make -j3 VERBOSE=1 && ctest --output-on-failure -L essential && make install

RUN ldconfig

ENTRYPOINT bash
