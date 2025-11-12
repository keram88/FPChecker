# Documentation for Docker Container Usage
#
# Build with: docker build -t fpchecker .
#
# This container can be used to run commands on files from your local directory
# by mounting the directory inside the container.
#
# To mount your current directory and run commands:
#
# 1. Run interactively to execute multiple commands (preferred):
#    docker run -it -v $(pwd):/workspace -w /workspace fpchecker /bin/bash
#
# 2. Mount current directory to /workspace:
#    docker run -v $(pwd):/workspace -w /workspace fpchecker <command>
#
# 3. Mount specific directory:
#    docker run -v /path/to/your/files:/workspace -w /workspace fpchecker <command>
#
# 4. Example - run a specific command on mounted files:
#    docker run -v $(pwd):/workspace -w /workspace fpchecker ls -la
#    docker run -v $(pwd):/workspace -w /workspace fpchecker cat file.txt
#
# 5. For Windows PowerShell, use ${PWD} instead of $(pwd):
#    docker run -v ${PWD}:/workspace -w /workspace fpchecker <command>
#

# Use miniconda3 as base image for conda package management
FROM continuumio/miniconda3:latest

# Set working directory
WORKDIR /tmp/tutorial

# Update conda and install mamba for faster package resolution
RUN conda update -n base -c defaults conda && \
    conda install -n base -c conda-forge mamba

# Create and activate conda environment, then install all dependencies
RUN mamba create --name tutorial_env -y && \
    echo "source activate tutorial_env" > ~/.bashrc

# Install system dependencies that might be needed
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install conda packages in the tutorial_env environment
RUN mamba install -n tutorial_env -c conda-forge -y \
    cmake \
    llvmdev=19.1.7 \
    clangxx=19.1.7 \
    python=3.12.9 \
    openmpi=5.0.7

# Activate environment and install pip packages
RUN /bin/bash -c "source activate tutorial_env && pip3 install matplotlib"

# Copy the entire FPChecker repository into the container
COPY . .

# Set up environment variables and build FPChecker
RUN /bin/bash -c "source activate tutorial_env && \
    mkdir -p build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=../install .. && \
    make && make install"

# Add the installation path to PATH permanently
ENV PATH="/tmp/tutorial/install/bin:${PATH}"

# Set the conda environment to be activated by default
ENV CONDA_DEFAULT_ENV=tutorial_env
ENV CONDA_PREFIX=/opt/conda/envs/tutorial_env
ENV PATH="/opt/conda/envs/tutorial_env/bin:${PATH}"

# Set the working directory to the build directory for convenience
WORKDIR /tmp/tutorial

# Default command to activate the environment and start a bash shell
CMD ["/bin/bash", "-c", "source activate tutorial_env && /bin/bash"]