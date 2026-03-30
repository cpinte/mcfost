# Use a lightweight Debian base image
FROM debian:bookworm-slim

# Set environment variables for MCFOST installation and runtime
ENV MCFOST_INSTALL="/opt/mcfost" \
    MCFOST_UTILS="/opt/mcfost_utils" \
    MCFOST_GIT=1 \
    MCFOST_AUTO_UPDATE=0 \
    OMP_STACKSIZE="512M" \
    PATH="/opt/mcfost/bin:${PATH}"

# Create necessary directories
RUN mkdir -p ${MCFOST_INSTALL}/bin \
           ${MCFOST_INSTALL}/lib \
           ${MCFOST_INSTALL}/include \
           ${MCFOST_UTILS}

# Install build dependencies:
# - git: To clone repositories
# - gfortran: The Fortran compiler for MCFOST
# - make, build-essential: For compiling code
# - wget: Potentially used by some scripts for downloads
# - libhdf5-dev: HDF5 library development files, essential for h5py (pymcfost dependency)
# - python3, python3-pip: For Python and its package manager
# - python3-venv: Useful for managing Python environments, though not strictly required for this direct install
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git \
    gfortran \
    make \
    wget \
    build-essential \
    libhdf5-dev \
    python3 \
    python3-pip \
    python3-venv && \
    rm -rf /var/lib/apt/lists/*

# Clone the MCFOST repository from GitHub
RUN git clone https://github.com/cpinte/mcfost.git /usr/local/src/mcfost

# Navigate to the MCFOST source directory's 'lib' folder
# Run the install script which compiles and installs MCFOST into MCFOST_INSTALL
WORKDIR /usr/local/src/mcfost/lib
RUN ./install.sh

# Navigate back to the root directory
WORKDIR /

# Perform MCFOST initial setup to download required data files
# These files are stored in the directory specified by MCFOST_UTILS
RUN mcfost -setup

# --- pymcfost installation ---
# Clone the pymcfost repository
RUN git clone https://github.com/cpinte/pymcfost.git /usr/local/src/pymcfost

# Install pymcfost using pip.
# This will also install its Python dependencies (numpy, scipy, matplotlib, astropy, h5py).
# Using --no-cache-dir to prevent caching pip packages, which helps keep the image size down.
WORKDIR /usr/local/src/pymcfost
RUN pip3 install . --no-cache-dir

# --- End of pymcfost installation ---

# Set the default command to execute when the container runs
# This will display the MCFOST help message.
# You could change this to start a Python interpreter if you primarily plan to use pymcfost,
# e.g., CMD ["python3"]
CMD ["mcfost", "--help"]
