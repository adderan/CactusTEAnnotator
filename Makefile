rootPath = .
include hal/include.mk

libSourcesAll = $(wildcard src/*.cpp)
libSources = $(subst src/main.cpp,,${libSourcesAll})
libHeaders = $(wildcard src/*.h)

all: cactusRepeatAnnotator repeatAnnotatorTests

repeats.a: ${libSources} ${libHeaders} hal/lib/halLib.a sonLib/lib/sonLib.a ${basicLibsDependencies}
	rm -f *.o
	${cpp} ${cppflags} -I src -I src -I hal/lib -I sonLib/lib -c ${libSources} -fopenmp
	ar rc repeats.a *.o
	ranlib repeats.a
	rm *.o

cactusRepeatAnnotator : src/main.cpp repeats.a hal/lib/halLib.a
	${cpp} ${cppflags} -I src -I src -I hal/lib -I sonLib/lib -I src -I tests -o cactusRepeatAnnotator src/main.cpp repeats.a sonLib/lib/sonLib.a hal/lib/halLib.a ${basicLibs} -fopenmp

repeatAnnotatorTests: src/test.cpp repeats.a hal/lib/halLib.a
	${cpp} ${cppflags} -UNDEBUG -I src -I src -I hal/lib -I sonLib/lib -I src -I tests -o repeatAnnotatorTests src/test.cpp repeats.a sonLib/lib/sonLib.a hal/lib/halLib.a ${basicLibs} -fopenmp