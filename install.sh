#!/usr/bin/env bash
set -e
unset PYTHONPATH

chimeras="$( cd .; pwd -P )"
venv="${chimeras}/venv"
conda create --prefix="${venv}" python=3.8

conda install -c bioconda ruffus cutadapt umi_tools fastq-tools --prefix="${venv}"
conda install -c bioconda seqtk bowtie samtools bedtools star=2.4.0j --prefix="${venv}"
conda install -c bioconda perl perl-app-cpanminus --prefix="${venv}"

cd "${venv}"
wget https://github.com/chaolinzhanglab/czplib/archive/v1.0.8.zip
unzip v1.0.8.zip
rm v1.0.8.zip

wget https://github.com/chaolinzhanglab/ctk/archive/v1.1.4.zip
unzip v1.1.4.zip
rm v1.1.4.zip

cd "${chimeras}"

