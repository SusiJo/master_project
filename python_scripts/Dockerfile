FROM nfcore/base:1.14
LABEL authors="Susanne Jodoin" \
      description="Docker image containing SusiJo/master_project python command-line tools and machine-learning"

# Create directories
RUN mkdir -p /data
RUN mkdir -p /results
RUN mkdir -p /app

# Define working directory
# WORKDIR .

# Copy files to container 
COPY data/* /data/
COPY app/* /app/

# Install the conda environment
COPY environment.yml .
RUN conda env create -f environment.yml && conda clean -a

# Activate environment
RUN echo "conda activate python_scripts" >> ~/.bashrc 

# Set environment variables
ENV PATH=/opt/conda/envs/python_scripts/bin:$PATH
ENV PYTHONPATH "${PYTHONPATH}:/app/"

# Set shell to bash
SHELL ["/bin/bash", "-c"]

