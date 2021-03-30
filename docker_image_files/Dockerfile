FROM debian:stretch-slim@sha256:d27bbe1259aeb6bed459440649ce5bda9083ab9782274c7bc469f02f283a9e18

#Copy local plotReads python program
COPY bamMetrics_modified.cpp /tmp
COPY setup.sh /opt

#Install necessary libraries, matplotlib for python, and samtools
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
    wget \
    make \
    gcc \
    g++ \
    ca-certificates \
    lbzip2 \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    perl \
    python2.7 \                                
    && apt-get clean -y \
    && cd /tmp \
    && g++ -O2 -o /usr/bin/bamMetrics /tmp/bamMetrics_modified.cpp \
    && wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 \
    && tar -xf samtools-1.11.tar.bz2 \
    && cd /tmp/samtools-1.11 && ./configure --prefix=/opt && make all all-htslib && make install install-htslib \
    && ln -s /opt/bin/samtools /usr/bin/samtools \
    && rm /tmp/samtools-1.11.tar.bz2 \
    && rm -rf /tmp/samtools-1.11

    ENTRYPOINT ["/opt/setup.sh"]
