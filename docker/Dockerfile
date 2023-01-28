################## BASE IMAGE ######################

FROM continuumio/miniconda3:latest

################## METADATA ######################
LABEL base.image="continuumio/miniconda3:latest"
LABEL version="1"
LABEL software="dreem"
LABEL software.version="0.1.0"
LABEL about.summary="analysis of mutational profiling (MaP) of RNA"
LABEL about.tags="RNA bioinformatics,RNA structure"

################## MAINTAINER ######################
MAINTAINER Joe Yesselman <jyesselm@unl.edu>

################## INSTALLATION ######################

USER root

RUN /bin/bash -c "mkdir /data && mkdir /data2 && mkdir /data3 && mkdir /data4 "

## install fastqc and unzip
RUN apt-get update && \
    apt-get -y -q install fastqc gcc unzip curl

ENV ZIP=bowtie2-2.2.9-linux-x86_64.zip
ENV URL=https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/
ENV FOLDER=bowtie2-2.2.9
ENV DST=/bin

## install bowtie2
RUN wget -q $URL/$ZIP -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    mv $DST/$FOLDER/* $DST && \
    rmdir $DST/$FOLDER

 # check that bowtie2 is accessible
RUN if [$(which bowtie2) = ""]; then exit 1; fi

## setup correct python and packages
RUN pip install cython colorlog click future plotly matplotlib pandas kaleido

## install cutadapt
RUN wget -q https://github.com/marcelm/cutadapt/archive/refs/tags/v1.18.zip -O $DST/v1.18.zip && \
    unzip $DST/v1.18.zip -d $DST && \
    rm $DST/v1.18.zip
RUN /bin/bash -c "cd $DST/cutadapt-1.18 && python setup.py install"

## install trim_galore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz && \
    tar xvzf trim_galore.tar.gz && \
    cp TrimGalore-0.6.6/trim_galore $DST && \
    rm -rf TrimGalore-0.6.6 trim_galore.tar.gz

 # check that trim_galore is accessible
RUN if [$(which trim_galore) = ""]; then exit 1; fi

RUN /bin/bash -c "mkdir -p /bin/rna_map"
COPY . /bin/rna_map
RUN /bin/bash -c "cd /bin/rna_map && python setup.py install"
RUN /bin/bash -c "cd /bin/rna_map && pytest -m quick"


WORKDIR /data

CMD ["rna-map"]
