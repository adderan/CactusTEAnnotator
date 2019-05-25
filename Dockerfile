FROM ubuntu:18.04
RUN apt-get update && apt-get install -y git build-essential zlib1g-dev libboost-dev wget libtext-soundex-perl hmmer

RUN mkdir /home/CactusTEAnnotator
COPY . /home/CactusTEAnnotator/

RUN cd /home/CactusTEAnnotator/ && make

ADD wrapper.sh /opt/wrapper.sh

ENV PATH /home/CactusTEAnnotator/RepeatMasker/:/home/CactusTEAnnotator/RepeatScout/:/home/CactusTEAnnotator/bin/:/home/CactusTEAnnotator/cactus/bin/:$PATH
RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["bash", "/opt/wrapper.sh"]
