################## BASE IMAGE ######################

FROM biocontainers/biocontainers:latest

################## METADATA ######################
LABEL base.image="biocontainers:latest"
LABEL version="1"
LABEL software="dreem"
LABEL software.version="0.1.0"
LABEL about.summary="analysis of DMS reactivity of RNA by the Rouskin lab (https://www.rouskinlab.com/)"
LABEL about.tags="RNA bioinformatics,RNA structure"

################## MAINTAINER ######################
MAINTAINER Joe Yesselman <jyesselm@unl.edu>

################## INSTALLATION ######################

ENV ZIP=bowtie2-2.2.9-linux-x86_64.zip
ENV URL=https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/
ENV FOLDER=bowtie2-2.2.9
ENV DST=/home/biodocker/bin

USER root

## install bowtie2
RUN wget $URL/$ZIP -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    mv $DST/$FOLDER/* $DST && \
    rmdir $DST/$FOLDER

## install fastqc
RUN apt-get update && \
    apt-get -y install fastqc

## setup correct python and packages
RUN conda create --name py3 python=3.7 --yes
ENV PATH /opt/conda/envs/py3/bin:$PATH
RUN /bin/bash -c "source activate py3"
RUN echo "source activate py3" > ~/.bashrc
RUN pip install cython colorlog click future biopython plotly matplotlib pandas

## install cutadapt
RUN wget https://github.com/marcelm/cutadapt/archive/refs/tags/v1.18.zip -O $DST/v1.18.zip && \
    unzip $DST/v1.18.zip -d $DST && \
    rm $DST/v1.18.zip
RUN /bin/bash -c "cd $DST/cutadapt-1.18 && python setup.py install"

## install trim_galore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz && \
    tar xvzf trim_galore.tar.gz && \
    cp TrimGalore-0.6.6/trim_galore $DST && \
    rm -rf TrimGalore-0.6.6 trim_galore.tar.gz

RUN git clone https://github.com/jyesselm/dreem && \
   /bin/bash -c "cd dreem && python setup.py install"

WORKDIR /data

# Change user
USER biodocker


CMD ["dreem"]
