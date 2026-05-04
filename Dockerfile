# syntax=docker/dockerfile:1
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ \
    cmake \
    make \
    libtbb-dev \
    python3 \
    python3-pip \
    ca-certificates \
 && rm -rf /var/lib/apt/lists/*
RUN python3 -m pip install --break-system-packages --no-cache-dir ViennaRNA==2.7.2
COPY . /opt/mRNAfold
WORKDIR /opt/mRNAfold
RUN mkdir build \
 && cmake -S cpp_src/ -B build/ -DCMAKE_BUILD_TYPE=Release \
 && cmake --build build -j"$(nproc)"
ENV PATH="/opt/mRNAfold/build/exe:${PATH}"
CMD ["fold_codon_graph"]