rootPath = .

libSourcesAll = $(wildcard impl/*.cpp)
libSources = $(subst impl/cactus_repeats.cpp,,${libSourcesAll})
libHeaders = $(wildcard impl/*.h)

murmurHashSources=${PWD}/smhasher/src/MurmurHash3.cpp

cpp=g++

all: ./bin/pairwise_distances ${PWD}/hal/lib/halLib.a ./bin/neighborJoining ./bin/maximizeOverlaps ./bin/getAlignmentDistances ./bin/getThreadPartitions

./bin/pairwise_distances: impl/pairwise_distances.cpp ${PWD}/sonLib/lib/sonLib.a
	PATH=${PWD}/hdf5/bin:${PATH} ${cpp} ${cppflags} -I smhasher/src -I hal/lib -I sonLib/lib -o ./bin/pairwise_distances impl/pairwise_distances.cpp ${murmurHashSources} sonLib/lib/sonLib.a

./bin/neighborJoining: impl/neighborJoining.c ${PWD}/sonLib/lib/sonLib.a
	gcc -g -o bin/neighborJoining -I sonLib/lib impl/neighborJoining.c ${PWD}/sonLib/lib/sonLib.a -lm

./bin/maximizeOverlaps: impl/maximizeOverlaps.c ${PWD}/sonLib/lib/sonLib.a
	gcc -g -o bin/maximizeOverlaps -I sonLib/lib -I poaV2/ impl/maximizeOverlaps.c ${PWD}/sonLib/lib/sonLib.a ${PWD}/poaV2/liblpo.a -lm

./bin/getAlignmentDistances: impl/getAlignmentDistances.c ${PWD}/sonLib/lib/sonLib.a
	gcc -g -o bin/getAlignmentDistances -I sonLib/lib -I poaV2/ impl/getAlignmentDistances.c ${PWD}/sonLib/lib/sonLib.a ${PWD}/poaV2/liblpo.a -lm

./bin/getThreadPartitions: impl/getThreadPartitions.c poaV2/liblpo.a
	gcc -g -o bin/getThreadPartitions -I poaV2/ impl/getThreadPartitions.c ${PWD}/poaV2/liblpo.a -lm

${PWD}/sonLib/lib/sonLib.a:
	cd ${PWD}/sonLib/ && make

${PWD}/hdf5/bin/h5c++:
	cd ${PWD}/hdf5 && ./configure --enable-cxx --prefix=${PWD}/hdf5 && CFLAGS=-std=c99 make -j4 -e && make install

${PWD}/hal/lib/halLib.a: ${PWD}/hdf5/bin/h5c++ ${PWD}/sonLib/lib/sonLib.a
	cd ${PWD}/hal && PATH=${PWD}/hdf5/bin:${PATH} make
	cp ${PWD}/hal/bin/* ./bin

./bin/poa:
	cd ${PWD}/poaV2 && make poa
	cp ${PWD}/poaV2/poa ./bin
docker:
	docker build -t cactus-te-annotator .

clean:
	rm -f repeats.a cactusRepeatAnnotator repeatAnnotatorTests *.o
