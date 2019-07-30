FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y git wget cpanminus build-essential python r-base

RUN apt-get install default-jdk -y
RUN git clone https://github.com/CSB5/OPERA-MS.git operams

WORKDIR /operams
RUN make
RUN perl OPERA-MS.pl CHECK_DEPENDENCY
RUN perl OPERA-MS.pl\
    --contig-file test_files/contigs.fasta\
    --short-read1 test_files/R1.fastq.gz\
    --short-read2 test_files/R2.fastq.gz\
    --long-read test_files/long_read.fastq\
    --out-dir test_files/RESULTS 2> test_files/log.err
ENTRYPOINT ["perl", "OPERA-MS.pl"]
