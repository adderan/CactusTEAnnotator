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


all: cactus RepeatMaskerRule halBinaries bin/lastz local


local: bin/tests bin/getTECandidates bin/build_clusters bin/filterNs bin/getConsensusFromPairwiseAlignments bin/getAlignmentDistances bin/buildClusters 

halBinaries:
	cd cactus && make
	cp cactus/submodules/hal/bin/* bin/

${libSonLib}:
	cd cactus && make

bin/lastz:
	cd cPecan/externalTools/lastz-distrib-1.03.54 && make
	cp cPecan/externalTools/lastz-distrib-1.03.54/src/lastz bin/

${objs}: ${sources} ${libSonLib}
	g++ -I ${sonLibInc} -I smhasher/src -c ${sources}

bin/getTECandidates: impl/getTECandidates.cpp ${libHal} ${libSonLib}
	PATH=${PWD}/cactus/submodules/hdf5/bin:${PATH} h5c++ ${cppflags} -o bin/getTECandidates -I ${halInc} -I ${sonLibInc} impl/getTECandidates.cpp ${libHal} ${libSonLib} -lm

bin/filterNs: impl/filterNs.c ${libSonLib}
	${cc} ${cflags} -o bin/filterNs -I ${sonLibInc} impl/filterNs.c ${libSonLib} -lm

impl/Consensus.o: impl/Consensus.c
	${cc} ${cflags} -DDEBUG_ -I ${sonLibInc} -I ${pinchesAndCactiInc} -c impl/Consensus.c
	mv Consensus.o impl/Consensus.o

bin/tests: impl/Consensus.o impl/tests.c ${libSonLib} ${pinchesAndCactiLib} ${cafLib} ${cactusLib}
	${cc} ${cflags} -o bin/tests -I impl/ -I ${sonLibInc} -I ${sonLibTestInc} -I ${pinchesAndCactiInc} -I ${cactusInc} impl/tests.c impl/Consensus.o ${pinchesAndCactiLib} ${cafLib} ${cactusLib} ${libSonLib} ${libCu} -lm -lz

bin/getConsensus: impl/Consensus.o impl/getConsensus.c ${libSonLib} ${pinchesAndCactiLib} ${cafLib} ${cactusLib}
	${cc} ${cflags} -o bin/getConsensus -I impl/ -I ${sonLibInc} -I ${pinchesAndCactiInc} -I ${cactusInc} impl/getConsensus.c impl/Consensus.o ${pinchesAndCactiLib} ${cafLib} ${cactusLib} ${libSonLib} -lm -lz

bin/buildClusters: impl/buildClusters.c ${libSonLib}
	${cc} ${cflags} -o bin/buildClusters -I ${sonLibInc} impl/buildClusters.c ${libSonLib} -lm

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
