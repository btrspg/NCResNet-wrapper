FROM ubuntu:18.04
MAINTAINER CHEN Yuelong <yuelong.chen.btr@gmail.com>

ENV DEBIAN_FRONTEND noninteractive

ARG depends="build-essential zip wget gnupg ca-certificates python3 python3-dev python3-pip git r-base"

# update



RUN apt update && \
    apt install -y  $depends && \
    Rscript -e 'install.packages("LncFinder",repos="https://cloud.r-project.org/");\
    install.packages("seqinr",repos="https://cloud.r-project.org/")'

WORKDIR /opt/tmp
ADD . ./
RUN pip3 install -r requirements.txt && \
    pip3 install .

WORKDIR /opt
RUN rm -rf /opt/tmp

CMD bash

