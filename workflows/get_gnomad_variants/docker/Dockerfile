FROM hailgenetics/hail:0.2.131-py3.11

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update --fix-missing && \
  apt-get install -y wget bzip2 build-essential postgresql libssl-dev libssl3 postgresql-doc-14 libssl-doc libpq-dev \
  less vim git &&  apt-get clean -y

USER root
ENV HOME=/root
ENV PATH="${PATH}:${HOME}/.local/bin:${HOME}/packages/bin"

RUN pip install psycopg2-binary
RUN pip3 install pandas
RUN pip3 --no-cache-dir install gnomad==0.8.0

ADD workflows/get_gnomad_variants/get_gnomad_variants.py /root/

WORKDIR $HOME
CMD ["/bin/bash"]
