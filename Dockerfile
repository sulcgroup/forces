FROM debian:bullseye

# update 
#RUN apt-get update -y 

# install required packages from repositories
#RUN sed -i -e "s/ main[[:space:]]*\$/ main contrib non-free/" /etc/apt/sources.list
RUN DEBIAN_FRONTEND=noninteractive apt-get update -y && apt-get upgrade -y && apt-get install -y --no-install-recommends python3 python3-biopython python3-numpy python3-scipy gawk sed python3-dev python3-numpy-dev gcc g++

# environment variables
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
ENV TERM xterm

WORKDIR /
RUN mkdir /build
ADD forces.tgz /build
WORKDIR /build/src/python
RUN python3 setup.py build && python3 setup.py install
WORKDIR /build/src/xds_calculation
RUN g++ -O3 -o fasta_xds calculate_xds.cpp && g++ -O3 -o window_scan scan_window.cpp && install -m 755 fasta_xds /usr/local/bin && install -m 755 window_scan /usr/local/bin

# cleanup
WORKDIR /
RUN rm -rf /build

# default interactive shell
CMD /bin/bash
