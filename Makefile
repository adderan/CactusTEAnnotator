murmurHashSources=smhasher/src/MurmurHash3.cpp

cpp=g++
cc=gcc

cflags_opt=-Wall -Werror -DNDEBUG
cflags_debug=-g -O0 -Wall -Werror

cppflags_opt=-Wall -DNDEBUG
cppflags_debug=-g -O0 -Wall

ifndef TE_ANNOTATOR_DEBUG
	cflags=${cflags_opt}
	cppflags=${cppflags_opt}
else
	cflags=${cflags_debug}
	cppflags=${cppflags_debug}
endif

objs=Minhash.o
sources=impl/Minhash.cpp

cObjs=repeatGraph.o
cSources=impl/repeatGraph.c

libSonLib = cactus/submodules/sonLib/lib/sonLib.a
libCu = cactus/submodules/sonLib/lib/cuTest.a
sonLibInc = cactus/submodules/sonLib/lib/
sonLibTestInc = cactus/submodules/sonLib/C/tests/

libHal = cactus/submodules/hal/lib/halLib.a
halInc = cactus/submodules/hal/lib/

cactusInc = cactus/lib
cafLib = cactus/lib/stCaf.a
cactusLib = cactus/lib/cactusLib.a

pinchesAndCactiInc = cactus/submodules/sonLib/lib/
pinchesAndCactiLib = cactus/submodules/sonLib/lib/stPinchesAndCacti.a

liblpo = poaV2/liblpo.a

all: cactus poa bin/RepeatScout RepeatMaskerRule halBinaries bin/lastz local


local: bin/neighborJoining bin/denseBundles bin/clusterByAlignmentDistances bin/getThreadPartitions bin/tests bin/getHeaviestBundles bin/minhash bin/poToGraphViz bin/getTECandidates bin/getSequencesFromHAL bin/build_clusters bin/filterNs bin/getElementsFromPinchGraph bin/getCoveredSeeds bin/repeatGraphTests bin/getAlignmentDistances bin/buildClusters

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

bin/lastz:
	cd cPecan/externalTools/lastz-distrib-1.03.54 && make
	cp cPecan/externalTools/lastz-distrib-1.03.54/src/lastz bin/

${objs}: ${sources} ${libSonLib}
	g++ -I ${sonLibInc} -I smhasher/src -c ${sources}

bin/tests: impl/tests.cpp ${libSonLib} ${objs}
	g++ -g -o bin/tests -I ${sonLibInc} impl/tests.cpp ${murmurHashSources} ${objs} ${libSonLib} -lm

bin/minhash: impl/MinhashMain.cpp ${libSonLib} ${objs}
	${cpp} ${cppflags} -I smhasher/src -I ${sonLibInc} -o ./bin/minhash impl/MinhashMain.cpp ${murmurHashSources} ${objs} ${libSonLib}

bin/poToGraphViz: impl/poToGraphViz.c ${liblpo} ${libSonLib}
	${cc} ${cflags} -o bin/poToGraphViz -I poaV2/ -I ${sonLibInc} impl/poToGraphViz.c ${liblpo} ${libSonLib} -lm

bin/getTECandidates: impl/getTECandidates.cpp ${libHal} ${libSonLib}
	PATH=${PWD}/cactus/submodules/hdf5/bin:${PATH} h5c++ ${cppflags} -o bin/getTECandidates -I ${halInc} -I ${sonLibInc} impl/getTECandidates.cpp ${libHal} ${libSonLib} -lm

bin/getSequencesFromHAL: impl/getSequencesFromHAL.cpp ${libHal} ${libSonLib}
	PATH=${PWD}/cactus/submodules/hdf5/bin:${PATH} h5c++ ${cppflags} -o bin/getSequencesFromHAL -I ${halInc} -I ${sonLibInc} impl/getSequencesFromHAL.cpp ${libHal} ${libSonLib} -lm


bin/neighborJoining: impl/neighborJoining.c ${libSonLib}
	${cc} ${cflags} -o bin/neighborJoining -I ${sonLibInc} impl/neighborJoining.c ${libSonLib} -lm

bin/denseBundles: impl/denseBundles.c ${libSonLib}
	${cc} ${cflags} -o bin/denseBundles -I ${sonLibInc} -I poaV2/ impl/denseBundles.c ${libSonLib} ${liblpo} -lm

bin/clusterByAlignmentDistances: impl/clusterByAlignmentDistances.cpp ${libSonLib} ${liblpo} ${objs}
	g++ -g -o bin/clusterByAlignmentDistances -I ${sonLibInc} -I poaV2/ impl/clusterByAlignmentDistances.cpp ${objs} ${libSonLib} ${murmurHashSources} ${liblpo} -lm

bin/getThreadPartitions: impl/getThreadPartitions.c ${liblpo}
	${cc} ${cflags} -o bin/getThreadPartitions -I poaV2/ impl/getThreadPartitions.c ${liblpo} -lm

bin/getHeaviestBundles: impl/getHeaviestBundles.c poaV2/liblpo.a
	${cc} ${cflags} -o bin/getHeaviestBundles -I poaV2/ impl/getHeaviestBundles.c ${liblpo} -lm

bin/filterNs: impl/filterNs.c ${libSonLib}
	${cc} ${cflags} -o bin/filterNs -I ${sonLibInc} impl/filterNs.c ${libSonLib} -lm

impl/repeatGraphs.o: impl/repeatGraphs.c
	${cc} ${cflags} -DDEBUG_ -I ${sonLibInc} -I ${pinchesAndCactiInc} -I ${cactusInc} -c impl/repeatGraphs.c
	mv repeatGraphs.o impl/repeatGraphs.o

bin/repeatGraphTests: impl/repeatGraphs.o impl/repeatGraphTests.c ${libSonLib} ${pinchesAndCactiLib} ${cafLib} ${cactusLib}
	${cc} ${cflags} -o bin/repeatGraphTests -I impl/ -I ${sonLibInc} -I ${sonLibTestInc} -I ${pinchesAndCactiInc} -I ${cactusInc} impl/repeatGraphTests.c impl/repeatGraphs.o ${pinchesAndCactiLib} ${cafLib} ${cactusLib} ${libSonLib} ${libCu} -lm -lz

bin/getElementsFromPinchGraph: impl/repeatGraphs.o impl/getElementsFromPinchGraph.c ${libSonLib} ${pinchesAndCactiLib} ${cafLib} ${cactusLib}
	${cc} ${cflags} -o bin/getElementsFromPinchGraph -I impl/ -I ${sonLibInc} -I ${pinchesAndCactiInc} -I ${cactusInc} impl/getElementsFromPinchGraph.c impl/repeatGraphs.o ${pinchesAndCactiLib} ${cafLib} ${cactusLib} ${libSonLib} -lm -lz

bin/getCoveredSeeds: impl/getCoveredSeeds.c ${libSonLib}
	${cc} ${cflags} -o bin/getCoveredSeeds -I ${sonLibInc} impl/getCoveredSeeds.c ${libSonLib} -lm

bin/getAlignmentDistances: impl/getAlignmentDistances.c ${libSonLib}
	${cc} ${cflags} -o bin/getAlignmentDistances -I ${sonLibInc} impl/getAlignmentDistances.c ${libSonLib} -lm

bin/buildClusters: impl/buildClusters.c ${libSonLib}
	${cc} ${cflags} -o bin/buildClusters -I ${sonLibInc} impl/buildClusters.c ${libSonLib} -lm

bin/build_clusters: scripts/build_clusters
	cp scripts/build_clusters bin/build_clusters
	chmod +x bin/build_clusters

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
	tar -xvf rmblast-2.9.0+-x64-linux.tar.gz
	rm rmblast-2.9.0+-x64-linux.tar.gz

RepeatMaskerRule: RepeatMasker rmblast-2.9.0 bin/trf
	cd RepeatMasker && sed \
    	-e "s/DEFAULT_SEARCH_ENGINE = \"crossmatch\"/DEFAULT_SEARCH_ENGINE = 'ncbi'/g;" \
		-e "s|RMBLAST_DIR\\s\+=\\s\+\"\/usr\/local\/rmblast\"|RMBLAST_DIR = \"${PWD}/rmblast-2.9.0\"|g;" \
		-e "s|TRF_PRGM = \"\"|TRF_PRGM = \"${PWD}/bin/trf\"|g;" \
		RepeatMaskerConfig.tmpl > RepeatMaskerConfig.pm

clean:
	rm -f repeats.a cactusRepeatAnnotator repeatAnnotatorTests *.o
	cd cactus && make clean
