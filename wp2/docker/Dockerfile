FROM rocker/r-ver:3.6.0
USER root


RUN apt-get update && apt-get install -y \
        git\
	tree\
        ssh\
        zlib1g-dev\
        sudo

RUN sudo apt-get install -y make build-essential libssl-dev zlib1g-dev libbz2-dev \
        libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev \
        xz-utils tk-dev libffi-dev liblzma-dev


RUN apt-get install -y python-pip
RUN sudo apt-get install -y libbz2-dev
RUN sudo apt-get install -y liblzma-dev
#RUN sudo apt-get install -y python3
#RUN sudo apt-get install -y python3-pip
RUN sudo apt-get install bwa

#RUN pip3 install snaptools
#RUN apt-get update && apt-get install -y software-properties-common gcc && \
#    add-apt-repository -y ppa:deadsnakes/ppa

#RUN apt-get update && apt-get install -y python3.6 python3-distutils python3-pip python3-apt
#RUN wget https://www.python.org/ftp/python/3.6.9/Python-3.6.9.tgz
#RUN tar xvf Python-3.6.9.tgz
#RUN Python-3.6.9/configure --enable-optimizations --enable-shared --with-ensurepip=install

#RUN sudo make altinstall

#RUN python3 -V

ARG BUILDDIR="/tmp/build"
ARG PYTHON_VER="3.6.8"
WORKDIR ${BUILDDIR}

RUN apt-get update -qq && \
apt-get upgrade -y  > /dev/null 2>&1 && \
apt-get install wget gcc make zlib1g-dev -y -qq > /dev/null 2>&1 && \
wget --quiet https://www.python.org/ftp/python/${PYTHON_VER}/Python-${PYTHON_VER}.tgz > /dev/null 2>&1 && \
tar zxf Python-${PYTHON_VER}.tgz && \
cd Python-${PYTHON_VER} && \
./configure  > /dev/null 2>&1 && \
make > /dev/null 2>&1 && \
make install > /dev/null 2>&1 && \
rm -rf ${BUILDDIR} 
#RUN curl https://bootstrap.pypa.io/pip/3.6/get-pip.py | sudo -H python3.6
RUN python3 -m pip install numpy
#RUN pip3 install numpy
RUN pip3 install macs2
RUN install2.r devtools
RUN install2.r usethis 
RUN install2.r doSNOW 
RUN install2.r plot3D
RUN install2.r --error BiocManager
RUN Rscript -e 'BiocManager::install("GenomicRanges")'
RUN Rscript -e 'BiocManager::install("rtracklayer")'

RUN sudo apt-get install -y make cmake gcc git libcurl4-openssl-dev
RUN sudo apt-get install -y libxml2-dev


RUN installGithub.r r3fang/SnapATAC


RUN Rscript -e 'BiocManager::install("rtracklayer")'
RUN pip3 install snaptools
RUN pip install pysam

RUN echo $FOO
COPY ./base.sh /
RUN chmod +x /base.sh
COPY ./fasta-to-snap.sh /
COPY ./start-fasta-to-snap.sh /
COPY ./preprocessing_WP2_2.R /
COPY ./parsebam.py /
COPY ./narrowpeaks_to_bed.py /

ENTRYPOINT ["/base.sh"]
