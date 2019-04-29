rootPath = .

libHeaders = $(wildcard impl/*.h)

murmurHashSources=${PWD}/smhasher/src/MurmurHash3.cpp

cpp=g++

objs=Minhash.o
sources=impl/Minhash.cpp

all: hal cactusTEAnnotator RepeatScout RepeatMasker

cactusTEAnnotator: bin/poa ./bin/minhash ${PWD}/hal/lib/halLib.a ./bin/neighborJoining ./bin/denseBundles ./bin/clusterByAlignmentDistances ./bin/getThreadPartitions ./bin/tests ./bin/getHeaviestBundles

hal: ./hal/lib/halLib.a
	cp hal/bin/* bin

bin/poa: poa


RepeatScout: bin/RepeatScout

${objs}: ${sources} sonLib/lib/sonLib.a
	g++ -I sonLib/lib -I smhasher/src -c ${sources}

./bin/tests: impl/tests.cpp ${PWD}/sonLib/lib/sonLib.a ${objs}
	g++ -g -o bin/tests -I sonLib/lib impl/tests.cpp ${murmurHashSources} ${objs} ${PWD}/sonLib/lib/sonLib.a -lm

./bin/minhash: impl/MinhashMain.cpp ${PWD}/sonLib/lib/sonLib.a ${objs}
	PATH=${PWD}/hdf5/bin:${PATH} ${cpp} ${cppflags} -I smhasher/src -I hal/lib -I sonLib/lib -o ./bin/minhash impl/MinhashMain.cpp ${murmurHashSources} ${objs} sonLib/lib/sonLib.a

./bin/poToGraphViz: impl/poToGraphViz.c poaV2/liblpo.a ${PWD}/sonLib/lib/sonLib.a
	gcc -g -o bin/poToGraphViz -I poaV2/ -I sonLib/lib impl/poToGraphViz.c ${PWD}/poaV2/liblpo.a ${PWD}/sonLib/lib/sonLib.a -lm


./bin/neighborJoining: impl/neighborJoining.c ${PWD}/sonLib/lib/sonLib.a
	gcc -g -o bin/neighborJoining -I sonLib/lib impl/neighborJoining.c ${PWD}/sonLib/lib/sonLib.a -lm

./bin/denseBundles: impl/denseBundles.c ${PWD}/sonLib/lib/sonLib.a
	gcc -g -o bin/denseBundles -I sonLib/lib -I poaV2/ impl/denseBundles.c ${PWD}/sonLib/lib/sonLib.a ${PWD}/poaV2/liblpo.a -lm

./bin/clusterByAlignmentDistances: impl/clusterByAlignmentDistances.cpp ${PWD}/sonLib/lib/sonLib.a poaV2/liblpo.a ${objs}
	g++ -g -o bin/clusterByAlignmentDistances -I sonLib/lib -I poaV2/ impl/clusterByAlignmentDistances.cpp ${objs} ${PWD}/sonLib/lib/sonLib.a ${murmurHashSources} ${PWD}/poaV2/liblpo.a -lm

./bin/getThreadPartitions: impl/getThreadPartitions.c poaV2/liblpo.a
	gcc -g -o bin/getThreadPartitions -I poaV2/ impl/getThreadPartitions.c ${PWD}/poaV2/liblpo.a -lm

bin/getHeaviestBundles: impl/getHeaviestBundles.cpp poaV2/liblpo.a
	g++ -g -o bin/getHeaviestBundles -I poaV2/ impl/getHeaviestBundles.cpp ${PWD}/poaV2/liblpo.a -lm

sonLib/lib/sonLib.a:
	cd ${PWD}/sonLib/ && make

hdf5/bin/h5c++:
	cd ${PWD}/hdf5 && ./configure --enable-shared --enable-cxx --prefix=${PWD}/hdf5 && CFLAGS=-std=c99 make -j4 -e && make install

./hal/lib/halLib.a: ${PWD}/hdf5/bin/h5c++ ${PWD}/sonLib/lib/sonLib.a
	cd ${PWD}/hal && PATH=${PWD}/hdf5/bin:${PATH} make

poa:
	wget https://downloads.sourceforge.net/project/poamsa/poamsa/2.0/poaV2.tar.gz
	tar -xvf poaV2.tar.gz poaV2/
	rm poaV2.tar.gz
	patch poaV2/black_flag.h poa.patch
	cd ${PWD}/poaV2 && make poa
	cp ${PWD}/poaV2/poa ./bin

bin/RepeatScout:
	cd RepeatScout && make
	mv RepeatScout/RepeatScout bin/RepeatScout

./RepeatMasker/Libraries/Dfam.hmm:
	wget http://www.dfam.org/releases/Dfam_3.0/families/Dfam.hmm.gz
	gunzip Dfam.hmm.gz
	mv Dfam.hmm ./RepeatMasker/Libraries

./RepeatMasker/Libraries/Dfam.embl:
	wget http://www.dfam.org/releases/Dfam_3.0/families/Dfam.embl.gz
	gunzip Dfam.embl.gz
	mv Dfam.embl ./RepeatMasker/Libraries


RepeatMasker: ./RepeatMasker/Libraries/Dfam.hmm ./RepeatMasker/Libraries/Dfam.embl
	wget http://tandem.bu.edu/trf/downloads/trf407b.linux64
	mv trf407b.linux64 ./bin/trf
	chmod +x ./bin/trf
	cd RepeatMasker && ./configure -trfbin=${PWD}/bin/trf --hmmerbin=/usr/bin/ -defaultengine=hmmer
	
clean:
	rm -f repeats.a cactusRepeatAnnotator repeatAnnotatorTests *.o
