FROM continuumio/miniconda:4.6.14
LABEL authors="Susanne Jodoin" \
      description="Docker image containing all software requirements for batch-correction and pca analysis"

# Create directories
RUN mkdir -p /data
RUN mkdir -p /results
RUN mkdir -p /app

# Define working directory
WORKDIR .

# Copy files to container
COPY data/* /data/
COPY app/* /app/

# Install the conda environment
COPY environment.yml .
RUN conda env create -f environment.yml && conda clean -a

# Activate the environment
RUN echo "conda activate susijo-batch-correction" >> ~/.bashrc
ENV PATH=/opt/conda/envs/susijo-batch-correction/bin:$PATH

# Set shell to bash
# SHELL ["/bin/bash", "-c"]

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name susijo-batch-correction > susijo-batch-correction_environment.yml

RUN touch .Rprofile
RUN touch .Renviron
