FROM ubuntu:18.04

RUN apt-get update \
    && apt-get install -y build-essential wget zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev unzip python openjdk-8-jre \
    && rm -rf /var/lib/apt/lists/* \
    && cd /tmp \
    && wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar -vxjf samtools-1.9.tar.bz2 \
    && cd samtools-1.9 \
    && make \
    && make install \
    && cd /tmp \
    && wget https://github.com/broadinstitute/gatk/releases/download/4.1.2.0/gatk-4.1.2.0.zip \
    && unzip gatk-4.1.2.0.zip \
    && mv gatk-4.1.2.0/* /usr/local/bin/ \
    && cd /tmp \
    && wget https://github.com/alexdobin/STAR/archive/2.7.1a.tar.gz \
    && tar -xvzf 2.7.1a.tar.gz \
    && mv STAR-2.7.1a/bin/Linux_x86_64/* /usr/local/bin/ \
    && rm -rf /tmp/*
