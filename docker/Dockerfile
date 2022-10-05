# use manylinux2014 as base
FROM quay.io/pypa/manylinux2014_x86_64:latest

# install llvm to get libclang to satisfy librocks dependencies
RUN yum -y install centos-release-scl && yum install -y llvm-toolset-7
ENV LIBCLANG_PATH=/opt/rh/llvm-toolset-7/root/usr/lib64/ \
  LIBCLANG_STATIC_PATH=/opt/rh/llvm-toolset-7/root/usr/lib64/ \
  CLANG_PATH=/opt/rh/llvm-toolset-7/root/usr/bin/clang

# install default rust toolchain
RUN curl --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

# add cargo and python interpreters to PATH
ENV PATH /opt/python/cp38-cp38/bin:/opt/python/cp39-cp39/bin:/opt/python/cp310-cp310/bin:/opt/python/cp311-cp311/bin:/root/.cargo/bin:$PATH

# install maturin in each supported version of python
RUN python3.8 -m pip install maturin && \
  python3.9 -m pip install maturin && \
  python3.10 -m pip install maturin && \
  python3.11 -m pip install maturin


# uncomment next line for local development, not for github actions
# WORKDIR /io
