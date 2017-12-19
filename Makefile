rootPath = .
include hal/include.mk

libSourcesAll = $(wildcard src/*.cpp)
libSources = $(subst src/main.cpp,,${libSourcesAll})
libHeaders = $(wildcard src/*.h)

all: cactusRepeatAnnotator repeatAnnotatorTests

repeats.a: ${libSources} ${libHeaders} hal/lib/halLib.a sonLib/lib/sonLib.a ${basicLibsDependencies}
	rm -f *.o
	${cpp} ${cppflags} -I src/ -I poaV2/ -I hal/lib -I sonLib/lib -c ${libSources} poaV2/liblpo.a
	ar rc repeats.a *.o
	ranlib repeats.a

cactusRepeatAnnotator : src/main.cpp repeats.a hal/lib/halLib.a
	${cpp} ${cppflags} -I src -I src -I hal/lib -I sonLib/lib -I src -I tests -o cactusRepeatAnnotator src/main.cpp repeats.a sonLib/lib/sonLib.a hal/lib/halLib.a poaV2/liblpo.a ${basicLibs}

repeatAnnotatorTests: src/test.cpp repeats.a hal/lib/halLib.a
	${cpp} ${cppflags} -UNDEBUG -I src -I src -I hal/lib -I sonLib/lib -I src -I tests -o repeatAnnotatorTests src/test.cpp repeats.a sonLib/lib/sonLib.a hal/lib/halLib.a ${basicLibs}

clean:
	rm -f repeats.a cactusRepeatAnnotator repeatAnnotatorTests *.o