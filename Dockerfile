FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y git wget cpanminus build-essential python r-base

RUN apt-get install default-jdk -y
RUN git clone --single-branch --branch OPERA-MS-0.9.0 https://github.com/CSB5/OPERA-MS.git operams

WORKDIR /operams
RUN make
RUN perl OPERA-MS.pl check-dependency
RUN perl OPERA-MS.pl\
    --contig-file test_files/contigs.fasta\
    --short-read1 test_files/R1.fastq.gz\
    --short-read2 test_files/R2.fastq.gz\
    --long-read test_files/long_read.fastq\
    --no-ref-clustering \
    --out-dir test_files/RESULTS 2> test_files/log.err
ENTRYPOINT ["perl", "OPERA-MS.pl"]
