# Use an official Python runtime as the base image
FROM python:3.8.15-bullseye

COPY ./python_requirements.txt ./python_requirements.txt

# Install required packages
RUN apt-get update && apt-get install -y \
    build-essential \
    libpcre3-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    gcc \ 
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Pandas and NumPy
RUN pip install -r ./python_requirements.txt

# set a working directory to install cmdstan
WORKDIR /cmdstan
# Install CmdStan
RUN python -c "import cmdstanpy; cmdstanpy.install_cmdstan(dir='.')"

# Set CMDSTAN environment variable
ENV CMDSTAN=/cmdstan
