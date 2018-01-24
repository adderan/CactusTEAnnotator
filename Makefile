rootPath = .

libSourcesAll = $(wildcard impl/*.cpp)
libSources = $(subst impl/cactus_repeats.cpp,,${libSourcesAll})
libHeaders = $(wildcard impl/*.h)

murmurHashSources=${PWD}/smhasher/src/MurmurHash3.cpp

cpp=g++

all: pairwise_distances

pairwise_distances: impl/pairwise_distances.cpp
	PATH=${PWD}/hdf5/bin:${PATH} ${cpp} ${cppflags} -I smhasher/src -I hal/lib -I sonLib/lib -o pairwise_distances impl/pairwise_distances.cpp ${murmurHashSources} sonLib/lib/sonLib.a

${PWD}/hdf5/bin/h5c++:
	cd ${PWD}/hdf5 && ./configure --enable-cxx --prefix=${PWD}/hdf5 && CFLAGS=-std=c99 make -j4 -e && make install

${PWD}/hal/lib/halLib.a: ${PWD}/hdf5/bin/h5c++ ${PWD}/sonLib/lib/sonLib.a
	cd ${PWD}/hal && PATH=${PWD}/hdf5/bin:${PATH} make
	
${PWD}/poaV2/liblpo.a:
	cd ${PWD}/poaV2 && make poa

clean:
	rm -f repeats.a cactusRepeatAnnotator repeatAnnotatorTests *.o
