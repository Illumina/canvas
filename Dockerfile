FROM ubuntu:14.04

RUN apt-get update
RUN apt-get upgrade -y

# install python and PIP
RUN apt-get install openssl
RUN apt-get install python2.7 python2.7-dev -y
RUN apt-get install cython -y
RUN apt-get install python-setuptools -y
RUN apt-get install python-pip -y
RUN apt-get install python-numpy -y
# ideally, we don't want to use the ubuntu version
# here because this bloats the Docker image
# RUN apt-get install python-pandas -y
RUN easy_install -U distribute

# gcc, numpy and such
RUN apt-get install build-essential -y
RUN apt-get install gfortran -y
RUN apt-get install -y libatlas-base-dev
RUN apt-get install pkg-config -y
RUN apt-get install software-properties-common python-software-properties -y
RUN apt-get install cmake -y
RUN apt-get install ncurses-dev -y
RUN apt-get install zlib1g-dev -y
RUN apt-get install bedtools -y

# upgrade to newest versions / install abd clean-up
RUN pip install --upgrade cython
RUN pip install --upgrade numpy
RUN apt-get install zlib1g-dev
RUN pip install pysam
RUN pip install bx-python
RUN apt-get clean -y

# copy git repository into the image
RUN mkdir -p /opt/thapmix
COPY . /opt/thapmix/

# run HapMix installer in the image
WORKDIR /opt/thapmix/
