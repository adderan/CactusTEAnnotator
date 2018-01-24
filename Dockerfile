FROM ubuntu:16.04

RUN apt-get update && apt-get install -y git gcc g++ build-essential

RUN mkdir -p /home/CactusTEAnnotator
COPY . /home/CactusTEAnnotator

RUN cd /home/CactusTEAnnotator && make

RUN mkdir /data
WORKDIR /data

RUN mkdir -p /opt/CactusTEAnnotator/wrapper.sh
COPY ./wrapper.sh /opt/CactusTEAnnotator/
RUN chmod +x /opt/CactusTEAnnotator/wrapper.sh

ENTRYPOINT ["bash", "/opt/CactusTEAnnotator/wrapper.sh"]
