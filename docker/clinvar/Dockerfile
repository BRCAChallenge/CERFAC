FROM ubuntu:22.04

ARG version
ARG username=docker
ARG grpname=docker

LABEL Description="NCBI EDirect"
LABEL Vendor="NCBI/NLM/NIH"
LABEL Version=${version}
LABEL Maintainer=archeney@ucsc.edu


# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    curl libxml-simple-perl libwww-perl libnet-perl build-essential \
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
    pkg-config \
    apt-transport-https software-properties-common dirmngr gpg-agent \
    && rm -rf /var/lib/apt/lists/*

RUN groupadd docker && \
    useradd --no-log-init -r -g docker docker && \
    mkdir -vp /home/docker && \
    chown -Rc docker:docker /home/docker

USER docker
WORKDIR /home/docker/

RUN sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

ENV PATH="/home/docker/edirect:${PATH}"
ENV EMAIL=docker@ncbi.nlm.nih.gov

CMD ["/bin/bash"]
USER root
# python packages
RUN pip3 install pandas

# add script
#WORKDIR /opt/scripts
#ADD cv_merge_script.py /opt/scripts

WORKDIR /home
ADD cv_merge_script.py /home

