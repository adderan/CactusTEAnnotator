rootPath = .

libSourcesAll = $(wildcard src/*.cpp)
libSources = $(subst src/cactus_repeats.cpp,,${libSourcesAll})
libHeaders = $(wildcard src/*.h)

murmurHashSources=${PWD}/smhasher/src/MurmurHash3.cpp

cpp=h5c++
cppflags=-fPIC

all: cactusRepeatAnnotator repeatAnnotatorTests

repeats.a: ${libSources} ${libHeaders} ${PWD}/hal/lib/halLib.a ${PWD}/sonLib/lib/sonLib.a ${PWD}/poaV2/liblpo.a ${basicLibsDependencies}
	rm -f *.o
	PATH=${PWD}/hdf5/bin:${PATH} ${cpp} ${cppflags} -I src/ -I smhasher/src/ -I poaV2/ -I hal/lib -I sonLib/lib -c ${libSources} ${murmurHashSources}
	ar rc repeats.a *.o
	ranlib repeats.a

cactusRepeatAnnotator : src/cactus_repeats.cpp repeats.a hal/lib/halLib.a
	PATH=${PWD}/hdf5/bin:${PATH} ${cpp} ${cppflags} -I src -I src -I hal/lib -I sonLib/lib -I src -I tests -o cactusRepeatAnnotator src/cactus_repeats.cpp repeats.a hal/lib/halLib.a sonLib/lib/sonLib.a poaV2/liblpo.a ${basicLibs}

repeatAnnotatorTests: src/test.cpp repeats.a hal/lib/halLib.a
	PATH=${PWD}/hdf5/bin:${PATH} ${cpp} ${cppflags} -UNDEBUG -I src -I src -I hal/lib -I sonLib/lib -I src -I tests -o repeatAnnotatorTests src/test.cpp ${PWD}/repeats.a ${PWD}/sonLib/lib/sonLib.a ${PWD}/hal/lib/halLib.a ${basicLibs}
	
${PWD}/sonLib/lib/sonLib.a:
	cd ${PWD}/sonLib && make

${PWD}/hdf5/bin/h5c++:
	cd ${PWD}/hdf5 && ./configure --enable-cxx --prefix=${PWD}/hdf5 && CFLAGS=-std=c99 make -j4 -e && make install

${PWD}/hal/lib/halLib.a: ${PWD}/hdf5/bin/h5c++ ${PWD}/sonLib/lib/sonLib.a
	cd ${PWD}/hal && PATH=${PWD}/hdf5/bin:${PATH} make
	
${PWD}/poaV2/liblpo.a:
	cd ${PWD}/poaV2 && make poa

clean:
	rm -f repeats.a cactusRepeatAnnotator repeatAnnotatorTests *.o
