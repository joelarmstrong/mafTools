include ../inc/common.mk
binPath = ../bin
extraAPI = src/blockTree.o src/coalescences.o ../lib/sharedMaf.o ../external/CuTest.a ../lib/common.o ../mafComparator/src/comparatorAPI.o ../mafComparator/src/comparatorRandom.o ${sonLibPath}/sonLib.a
testAPI = ${extraAPI} src/test.blockTree.o src/test.coalescences.o ../external/CuTest.a
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

test: test/allTests test/mafPhyloComparator
	./test/allTests && python2.7 src/test.mafPhyloComparator.py --verbose && rm -rf ./test/ && rmdir ./tempTestDir

test/%: src/%.c ${testAPI}
	@mkdir -p $(dir $@)
	${cxx} -o $@.tmp $^ ${testFlags} ${extraIncFlags} ../external/CuTest.a ${lm}
	mv $@.tmp $@

test/%.o: ../lib/%.c ../inc/%.h
	mkdir -p $(dir $@)
	${cxx} -c $< -o $@.tmp ${testFlags} ${extraIncFlags} ${lm}
	mv $@.tmp $@

