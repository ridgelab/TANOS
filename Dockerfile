# I just added this dockerfile to have it and start working on it, it is not finalized, I copied it from chat
# INSTRUCTION args
# Use a lightweight base image
#FROM ubuntu:22.04

# Avoid interactive prompts
#ENV DEBIAN_FRONTEND=noninteractive

# Install IQ-TREE and basic utilities
FROM mambaorg/micromamba:latest

RUN micromamba install -y -n base -c conda-forge -c bioconda iqtree && \
    micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH

# Set working directory
WORKDIR /app

# Copy your program into the container
COPY . /app

# Make sure scripts are executable if needed
# RUN chmod +x your_script.py

# Default command - this will change and will probably depend on how we connect the steps in the pathway
CMD ["iqtree", "-h"]
