FROM ubuntu:18.04

RUN apt-get update && apt-get install -y wget unzip perl build-essential

RUN mkdir /home/bin


RUN wget https://tandem.bu.edu/trf/trf409.linux64.download.html /home/bin/trf && chmod +x /home/bin/trf

RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/2.2.28 /home/rmblast && cp /home/blast+/bin/* /home/bin/

RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/ /home/blast+ && cp /home/blast+/bin/* /home/bin


RUN RUN wget http://www.repeatmasker.org/RepeatMasker-open-4-0-8.tar.gz /home/RepeatMasker.tar.gz && tar -xvf /home/RepeatMasker.tar.gz

WORKDIR /home/RepeatMasker

RUN perl -0p -e 's/\/usr\/local\/rmblast/\/home/bin/g;' \
	-e 's/DEFAULT_SEARCH_ENGINE = "crossmatch"/DEFAULT_SEARCH_ENGINE = "ncbi"/g;' \
	-e '/s/TRF_PRGM = ""/TRF_PRGRM = \/home/bin/trf/g;' RepeatMaskerConfig.tmpl > RepeatMaskerConfig.pm

ENTRYPOINT /home/RepeatMasker/RepeatMasker

