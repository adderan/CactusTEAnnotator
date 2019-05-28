murmurHashSources=smhasher/src/MurmurHash3.cpp

cpp=g++
cflags=-g -O0 -Wall -Werror

objs=Minhash.o
sources=impl/Minhash.cpp

libSonLib = cactus/submodules/sonLib/lib/sonLib.a
sonLibInc = cactus/submodules/sonLib/lib/

liblpo = poaV2/liblpo.a

all: cactus poa bin/RepeatScout RepeatMaskerRule halBinaries cte


cte: bin/neighborJoining bin/denseBundles bin/clusterByAlignmentDistances bin/getThreadPartitions bin/tests bin/getHeaviestBundles bin/build_families bin/minhash

halBinaries:
	cd cactus && make
	cp cactus/submodules/hal/bin/* bin/

${libSonLib}:
	cd cactus && make

${liblpo}:
	wget https://downloads.sourceforge.net/project/poamsa/poamsa/2.0/poaV2.tar.gz
	tar -xvf poaV2.tar.gz poaV2/
	rm poaV2.tar.gz
	patch poaV2/black_flag.h poa.patch
	cd poaV2 && make poa
	cp poaV2/poa ./bin

poa: ${liblpo}
	cp poaV2/poa bin/


${objs}: ${sources} ${libSonLib}
	g++ -I ${sonLibInc} -I smhasher/src -c ${sources}

bin/tests: impl/tests.cpp ${libSonLib} ${objs}
	g++ -g -o bin/tests -I ${sonLibInc} impl/tests.cpp ${murmurHashSources} ${objs} ${libSonLib} -lm

bin/minhash: impl/MinhashMain.cpp ${libSonLib} ${objs}
	${cpp} ${cppflags} -I smhasher/src -I ${sonLibInc} -o ./bin/minhash impl/MinhashMain.cpp ${murmurHashSources} ${objs} ${libSonLib}

bin/build_families: impl/build_families.c ${libSonLib}
	gcc ${cflags} -I ${sonLibInc} -o bin/build_families impl/build_families.c ${libSonLib} -lm

bin/poToGraphViz: impl/poToGraphViz.c ${liblpo} ${libSonLib}
	gcc -g -o bin/poToGraphViz -I poaV2/ -I cactus/sonLib/lib impl/poToGraphViz.c poaV2/liblpo.a cactus/sonLib/lib/sonLib.a -lm


bin/neighborJoining: impl/neighborJoining.c ${libSonLib}
	gcc -g -o bin/neighborJoining -I ${sonLibInc} impl/neighborJoining.c ${libSonLib} -lm

bin/denseBundles: impl/denseBundles.c ${libSonLib}
	gcc ${cflags} -o bin/denseBundles -I ${sonLibInc} -I poaV2/ impl/denseBundles.c ${libSonLib} ${liblpo} -lm

bin/clusterByAlignmentDistances: impl/clusterByAlignmentDistances.cpp ${libSonLib} ${liblpo} ${objs}
	g++ -g -o bin/clusterByAlignmentDistances -I ${sonLibInc} -I poaV2/ impl/clusterByAlignmentDistances.cpp ${objs} ${libSonLib} ${murmurHashSources} ${liblpo} -lm

bin/getThreadPartitions: impl/getThreadPartitions.c ${liblpo}
	gcc -g -o bin/getThreadPartitions -I poaV2/ impl/getThreadPartitions.c ${liblpo} -lm

bin/getHeaviestBundles: impl/getHeaviestBundles.cpp poaV2/liblpo.a
	g++ -g -o bin/getHeaviestBundles -I poaV2/ impl/getHeaviestBundles.cpp ${liblpo} -lm

bin/RepeatScout:
	cd RepeatScout && make
	mv RepeatScout/RepeatScout bin/RepeatScout

./RepeatMasker/Libraries/Dfam.hmm:
	wget http://www.dfam.org/releases/Dfam_3.0/families/Dfam.hmm.gz
	gunzip Dfam.hmm.gz
	mv Dfam.hmm ./RepeatMasker/Libraries

RepeatMasker/Libraries/Dfam.embl:
	wget http://www.dfam.org/releases/Dfam_3.0/families/Dfam.embl.gz
	gunzip Dfam.embl.gz
	mv Dfam.embl ./RepeatMasker/Libraries

bin/trf:
	wget http://tandem.bu.edu/trf/downloads/trf407b.linux64
	mv trf407b.linux64 ./bin/trf
	chmod +x ./bin/trf

RepeatMasker:
	wget http://repeatmasker.org/RepeatMasker-open-4-0-9-p2.tar.gz 
	tar -xvf RepeatMasker-open-4-0-9-p2.tar.gz 
	rm RepeatMasker-open-4-0-9-p2.tar.gz 

rmblast-2.9.0:
	wget http://www.repeatmasker.org/rmblast-2.9.0+-x64-linux.tar.gz

RepeatMaskerRule: RepeatMasker rmblast-2.9.0 bin/trf
	cd RepeatMasker && sed \
    	-e "s/DEFAULT_SEARCH_ENGINE = \"crossmatch\"/DEFAULT_SEARCH_ENGINE = 'ncbi'/g;" \
		-e "s|RMBLAST_DIR\\s\+=\\s\+\"\/usr\/local\/rmblast\"|RMBLAST_DIR = \"${PWD}/rmblast-2.9.0\"|g;" \
		-e "s|TRF_PRGM = \"\"|TRF_PRGM = \"${PWD}/bin/trf\"|g;" \
		RepeatMaskerConfig.tmpl > RepeatMaskerConfig.pm

clean:
	rm -f repeats.a cactusRepeatAnnotator repeatAnnotatorTests *.o
	cd cactus && make clean
