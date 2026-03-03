# Use the official Dataproc Serverless base image (includes Spark, Java, Hadoop)
# Explicit platform required: Dataproc workers are linux/amd64; building on Apple Silicon without
# this flag produces an arm64 image that cannot be pulled on Dataproc.
FROM --platform=linux/amd64 europe-west1-docker.pkg.dev/cloud-dataproc/spark/dataproc_2.2:latest

# Switch to root to install system dependencies
USER root

# Install compilers for CmdStan
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY Docker_container_VIDRA_image/python_requirements.txt /tmp/
RUN pip install -r /tmp/python_requirements.txt
RUN pip install cmdstanpy==1.2.4

# Install CmdStan
# We install it to a global path and set permissions
ENV CMDSTAN_INSTALL_PATH=/opt/cmdstan
RUN python3 -c "import cmdstanpy; cmdstanpy.install_cmdstan(dir='/opt/cmdstan', version='2.34.1')"
ENV CMDSTAN=/opt/cmdstan/cmdstan-2.34.1

# Copy and COMPILE the Stan models securely into the image
# This ensures workers don't need to compile them at runtime
COPY stan_models/ /opt/vidra/stan_models/
WORKDIR /opt/vidra/stan_models/
RUN python3 -c "from cmdstanpy import CmdStanModel; \
    CmdStanModel(stan_file='VIDRA.stan').compile(); \
    CmdStanModel(stan_file='VIDRA_single_variant.stan').compile()"

# Switch back to the 'spark' user (Required for Dataproc Serverless)
USER spark