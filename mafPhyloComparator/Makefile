include ../inc/common.mk
binPath = ../bin
extraAPI = src/blockTree.o src/coalescences.o ../lib/sharedMaf.o ../external/CuTest.a ../lib/common.o ../mafComparator/src/comparatorAPI.o ../mafComparator/src/comparatorRandom.o ${sonLibPath}/sonLib.a
progs = $(foreach f, mafPhyloComparator, ${binPath}/$f)

# TODO: do this properly
extraIncFlags = -I../mafComparator/src

all: ${progs}

${binPath}/%: src/%.c ${extraAPI}
	@mkdir -p $(dir $@)
	${cxx} -o $@.tmp $^ ${cflags} ${extraIncFlags} ${lm}
	mv $@.tmp $@

%.o: %.c %.h
	${cxx} -c $< -o $@.tmp ${cflags} ${extraIncFlags}
	mv $@.tmp $@

clean:
	rm -f *.o ${progs} src/*.o && rm -rf ./test/ src/buildVersion.c src/buildVersion.h
