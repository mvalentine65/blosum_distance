# pull from Alma8 container for working dependencies
FROM quay.io/pypa/manylinux_2_28_x86_64

# install libclang and librocks
# cmake builds libz-ng-sys, pulled in by the zlib-ng feature on flate2
RUN dnf install -y clang rocksdb-devel cmake

# set global home variables
# github actions run under a separate user profile which can't access root's rust bin
ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH
RUN curl --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

ENV PATH=/opt/python/cp310-cp310/bin:/opt/python/cp311-cp311/bin:/opt/python/cp312-cp312/bin:/opt/python/cp313-cp313/bin:/opt/python/cp314-cp314/bin:/root/.cargo/bin:$PATH

RUN python3.10 -m pip install maturin && \
  python3.11 -m pip install maturin && \
  python3.12 -m pip install maturin && \
  python3.13 -m pip install maturin && \
  python3.14 -m pip install maturin

#WORKDIR /io

