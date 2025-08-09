FROM continuumio/miniconda3

# Copy your environment definition
COPY brava_hits_common_condition_check_conda_env.yaml /tmp/environment.yaml

# Create the environment and clean up
RUN conda env create -f /tmp/environment.yaml && conda clean -afy

# Use this env for all subsequent commands, including Snakemake
SHELL ["conda", "run", "--no-capture-output", "-n", "brava_hits_common_condition_check", "/bin/bash", "-c"]

# install dependencies, cleanup apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
 wget \
 unzip \
 ca-certificates \
 perl \
 bzip2 \
 autoconf \
 automake \
 make \
 gcc \
 zlib1g-dev \
 libbz2-dev \
 liblzma-dev \
 libcurl4-gnutls-dev \
 libssl-dev \
 libperl-dev \
 libgsl0-dev && \
 rm -rf /var/lib/apt/lists/* && apt-get autoclean

# for easy upgrade later. ARG variables only persist during build time
ARG bcftoolsVer="1.12"

RUN wget https://github.com/samtools/bcftools/releases/download/${bcftoolsVer}/bcftools-${bcftoolsVer}.tar.bz2 && \
 tar -vxjf bcftools-${bcftoolsVer}.tar.bz2 && \
 rm bcftools-${bcftoolsVer}.tar.bz2 && \
 cd bcftools-${bcftoolsVer} && \
 make && \
 make install && \
 mkdir /data

ARG plink2Ver="v2.0.0-a.6.21"
RUN wget -O /tmp/plink2.zip https://github.com/chrchang/plink-ng/releases/download/${plink2Ver}/plink2_linux_x86_64.zip && unzip /tmp/plink2.zip

# set $PATH (honestly unnecessary here, lol) and locale settings for singularity compatibility
ENV PATH="$PATH" \
 LC_ALL=C

# set working directory
WORKDIR /data