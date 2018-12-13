FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y git wget cpanminus build-essential python r-base

ARG UNAME=opera
ARG UID=1000
ARG GID=1000
RUN groupadd -g $GID -o $UNAME
RUN useradd -m -u $UID -g $GID -o -s /bin/bash $UNAME
USER $UNAME

WORKDIR /home/$UNAME/
RUN git clone https://github.com/CSB5/OPERA-MS.git operams

WORKDIR /home/$UNAME/operams
USER root
RUN make
USER $UNAME
ENTRYPOINT ["perl", "OPERA-MS.pl"]

