rootPath = ../
include ../include.mk

libSourcesAll = $(wildcard impl/*.cpp)
libSources = $(subst impl/halRepeatsMain.cpp,,${libSourcesAll})
libHeaders = $(wildcard inc/*.h)

all: ${binPath}/halRepeats

${libPath}/halRepeats.a: ${libSources} ${libHeaders} ${libPath}/halLib.a ${basicLibsDependencies}
	cp ${libHeaders} ${libPath}/
	rm -f *.o
	${cpp} ${cppflags} -I inc -I impl -I ${libPath}/ -c ${libSources}
	ar rc halRepeats.a *.o
	ranlib halRepeats.a
	rm *.o
	mv halRepeats.a ${libPath}/

${binPath}/halRepeats : impl/halRepeatsMain.cpp ${libPath}/halRepeats.a ${libPath}/halLib.a
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halRepeats impl/halRepeatsMain.cpp ${libPath}/halRepeats.a ${libPath}/halLib.a ${basicLibs}