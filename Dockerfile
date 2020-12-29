FROM continuumio/miniconda

RUN conda install python=3.8
RUN apt-get update -y
RUN apt-get -y install gcc make libuv1-dev

##Grab BLAST
WORKDIR /blast
RUN wget -O blast_source ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.16.0-x64-linux.tar.gz
RUN tar -xvf blast_source && cp ncbi-igblast-1.16.0/bin/igblastn /usr/local/bin/igblastn && rm -rf ncbi-igblast-1.16.0

#COPY sadie from local
COPY sadie /source/sadie/
COPY setup.py /source/.
COPY MANIFEST.in /source/.
COPY requirements.txt /source/
WORKDIR /source/
RUN pip3 install --upgrade pip
RUN pip3 install .
