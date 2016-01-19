FROM debian:jessie
MAINTAINER Thomas Bigot <tomed.bigot@gmail.com>
RUN apt-get update && apt-get install -y git wget
RUN echo "deb http://biopp.univ-montp2.fr/repos/apt/ Trusty main" >> /etc/apt/sources.list
RUN wget http://biopp.univ-montp2.fr/repos/apt/conf/biopp.gpg.key
RUN apt-key add biopp.gpg.key
RUN apt-get update
RUN apt-get install -y libbpp-phyl-dev