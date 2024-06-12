FROM ubuntu:20.04
MAINTAINER archeney@ucsc.edu

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    gcc \
    git \
    make \
    bzip2 \
    tabix \
    python3 \
    python3-pip \
    libncurses5-dev \
    libncursesw5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    autoconf \
    build-essential \
    pkg-config \
    apt-transport-https software-properties-common dirmngr gpg-agent \
    && rm -rf /var/lib/apt/lists/*

RUN yes yes | sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)" \
    && echo "export PATH=/home/jupyter/edirect:\${PATH}" >> ${HOME}/.bashrc \
    && echo "export PATH=/home/jupyter/packages/bin/edirect:\${PATH}" >> ${HOME}/.bashrc 


# python packages
RUN pip3 install pandas

# add script
#WORKDIR /opt/scripts
#ADD cv_merge_script.py /opt/scripts

WORKDIR /home
ADD cv_merge_script.py /home
