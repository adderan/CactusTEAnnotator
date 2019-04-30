rootPath = .

libHeaders = $(wildcard impl/*.h)

murmurHashSources=smhasher/src/MurmurHash3.cpp

cpp=g++

objs=Minhash.o
sources=impl/Minhash.cpp

all: sonLib poa RepeatScout RepeatMasker cte

docker: sonLib poa RepeatScout cte

cte: bin/neighborJoining bin/denseBundles bin/clusterByAlignmentDistances bin/getThreadPartitions bin/tests bin/getHeaviestBundles

poa:
	wget https://downloads.sourceforge.net/project/poamsa/poamsa/2.0/poaV2.tar.gz
	tar -xvf poaV2.tar.gz poaV2/
	rm poaV2.tar.gz
	patch poaV2/black_flag.h poa.patch
	cd poaV2 && make poa
	cp poaV2/poa ./bin

poaV2/liblpo.a: poa

sonLib:
	cd sonLib/ && make

sonLib/lib/sonLib.a: sonLib

hdf5:
	cd hdf5 && ./configure --enable-shared --enable-cxx --prefix=${PWD}/hdf5 && CFLAGS=-std=c99 make -j4 -e && make install

hal:
	cd hal && PATH=${PWD}/hdf5/bin:${PATH} make
	cp hal/bin/* bin

${objs}: ${sources} sonLib/lib/sonLib.a
	g++ -I sonLib/lib -I smhasher/src -c ${sources}

bin/tests: impl/tests.cpp sonLib/lib/sonLib.a ${objs}
	g++ -g -o bin/tests -I sonLib/lib impl/tests.cpp ${murmurHashSources} ${objs} sonLib/lib/sonLib.a -lm

bin/minhash: impl/MinhashMain.cpp sonLib/lib/sonLib.a ${objs}
	PATH=${PWD}/hdf5/bin:${PATH} ${cpp} ${cppflags} -I smhasher/src -I hal/lib -I sonLib/lib -o ./bin/minhash impl/MinhashMain.cpp ${murmurHashSources} ${objs} sonLib/lib/sonLib.a

bin/poToGraphViz: impl/poToGraphViz.c poaV2/liblpo.a sonLib/lib/sonLib.a
	gcc -g -o bin/poToGraphViz -I poaV2/ -I sonLib/lib impl/poToGraphViz.c poaV2/liblpo.a sonLib/lib/sonLib.a -lm


bin/neighborJoining: impl/neighborJoining.c sonLib/lib/sonLib.a
	gcc -g -o bin/neighborJoining -I sonLib/lib impl/neighborJoining.c sonLib/lib/sonLib.a -lm

bin/denseBundles: impl/denseBundles.c sonLib/lib/sonLib.a
	gcc -g -o bin/denseBundles -I sonLib/lib -I poaV2/ impl/denseBundles.c sonLib/lib/sonLib.a poaV2/liblpo.a -lm

bin/clusterByAlignmentDistances: impl/clusterByAlignmentDistances.cpp sonLib/lib/sonLib.a poaV2/liblpo.a ${objs}
	g++ -g -o bin/clusterByAlignmentDistances -I sonLib/lib -I poaV2/ impl/clusterByAlignmentDistances.cpp ${objs} sonLib/lib/sonLib.a ${murmurHashSources} poaV2/liblpo.a -lm

bin/getThreadPartitions: impl/getThreadPartitions.c poaV2/liblpo.a
	gcc -g -o bin/getThreadPartitions -I poaV2/ impl/getThreadPartitions.c poaV2/liblpo.a -lm

bin/getHeaviestBundles: impl/getHeaviestBundles.cpp poaV2/liblpo.a
	g++ -g -o bin/getHeaviestBundles -I poaV2/ impl/getHeaviestBundles.cpp poaV2/liblpo.a -lm

RepeatScout:
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
	cd hdf5 && make clean
	cd hal && make clean
	cd sonLib && make clean
