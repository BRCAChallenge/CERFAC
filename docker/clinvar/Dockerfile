FROM ncbi/edirect:latest
MAINTAINER archeney@ucsc.edu

USER root
# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    curl \
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


# python packages
RUN pip3 install pandas

# add script
#WORKDIR /opt/scripts
#ADD cv_merge_script.py /opt/scripts

WORKDIR /home
ADD cv_merge_script.py /home

