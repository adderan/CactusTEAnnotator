rootPath = .
include hal/include.mk

libSourcesAll = $(wildcard impl/*.cpp)
libSources = $(subst impl/main.cpp,,${libSourcesAll})
libHeaders = $(wildcard inc/*.h)

all: cactusRepeatAnnotator

repeats.a: ${libSources} ${libHeaders} hal/lib/halLib.a sonLib/lib/sonLib.a ${basicLibsDependencies}
	rm -f *.o
	${cpp} ${cppflags} -I inc -I impl -I hal/lib -I sonLib/lib -c ${libSources} -fopenmp
	ar rc repeats.a *.o
	ranlib repeats.a
	rm *.o

cactusRepeatAnnotator : impl/main.cpp repeats.a hal/lib/halLib.a
	${cpp} ${cppflags} -I inc -I impl -I hal/lib -I sonLib/lib -I impl -I tests -o cactusRepeatAnnotator impl/main.cpp repeats.a sonLib/lib/sonLib.a hal/lib/halLib.a ${basicLibs} -fopenmp

