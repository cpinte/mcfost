# Use a lightweight Debian base image
FROM debian:bookworm-slim

# Set environment variables for MCFOST installation and runtime
# MCFOST_INSTALL: Directory where MCFOST will be installed
# MCFOST_UTILS: Directory for MCFOST's data files (stellar spectra, etc.)
# MCFOST_GIT, MCFOST_AUTO_UPDATE: Settings for source installation as per MCFOST docs
# OMP_STACKSIZE: Recommended OpenMP stack size for performance
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
# - git: To clone the MCFOST repository
# - gfortran: The Fortran compiler for MCFOST
# - make, build-essential: For compiling the code
# - wget: Potentially used by some scripts for downloads
# - libhdf5-dev: HDF5 library development files, a common dependency for scientific software
# - python3, python3-pip: For Python utilities like pre-commit or if xgboost is used
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git \
    gfortran \
    make \
    wget \
    build-essential \
    libhdf5-dev \
    python3 \
    python3-pip && \
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

# Set the default command to execute when the container runs
# This will display the MCFOST help message
CMD ["mcfost", "--help"]
