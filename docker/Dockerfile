FROM ubuntu:18.04 AS builder
RUN apt-get update && apt-get install -y git build-essential zlib1g-dev libboost-dev

RUN git clone https://github.com/adderan/CactusTEAnnotator /home/CactusTEAnnotator

RUN cd /home/CactusTEAnnotator/ && make hal cactusTEAnnotator RepeatScout


FROM ubuntu:18.04
RUN apt-get update && apt-get install -y hmmer wget

ADD wrapper.sh /opt/wrapper.sh
COPY --from=builder /home/CactusTEAnnotator/bin/* /usr/local/bin

COPY --from=builder /home/RepeatScout/RepeatScout /usr/local/bin

WORKDIR /home/
RUN wget http://www.repeatmasker.org/RepeatMasker-open-4-0-7.tar.gz && tar -xvf RepeatMasker-open-4-0-7.tar.gz && rm RepeatMasker-open-4-0-7.tar.gz

RUN wget http://tandem.bu.edu/trf/downloads/trf407b.linux64
RUN mv trf407b.linux64 /usr/local/bin/trf
RUN wget http://www.dfam.org/releases/Dfam_3.0/families/Dfam.embl.gz && gunzip Dfam.embl.gz && mv Dfam.embl /home/RepeatMasker/Libraries
RUN wget http://www.dfam.org/releases/Dfam_3.0/families/Dfam.hmm.gz && gunzip Dfam.hmm.gz && mv Dfam.hmm /home/RepeatMasker/Libraries
RUN /home/RepeatMasker/configure --hmmerbin=/usr/bin --trfbin=/usr/local/bin/trf

ENTRYPOINT ["bash", "/opt/wrapper.sh"]
