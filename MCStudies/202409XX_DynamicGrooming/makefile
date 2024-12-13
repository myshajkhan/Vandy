
default: DetEf

objects: CATree.o DrawRandom.o DrawRandom2.o TauHelperFunctions3.o

CATree.o: CATree.cpp CATree.h
	g++ CATree.cpp -c -o CATree.o -std=c++11

DrawRandom.o: DrawRandom.cpp DrawRandom.h
	g++ DrawRandom.cpp -c -o DrawRandom.o -std=c++11

DrawRandom2.o: DrawRandom2.cpp DrawRandom2.h
	g++ DrawRandom2.cpp -c -o DrawRandom2.o -std=c++11

TauHelperFunctions3.o: TauHelperFunctions3.cpp TauHelperFunctions3.h
	g++ TauHelperFunctions3.cpp -c -o TauHelperFunctions3.o -std=c++11

Execute: ExampleExecutable_latekT.cpp objects
	g++ ExampleExecutable_latekT.cpp -o Execute -std=c++11 \
		CATree.o DrawRandom.o DrawRandom2.o TauHelperFunctions3.o \
		`root-config --cflags --libs` -g

TestRun: Execute
	./Execute





detectoreffect: Example.cpp objects
	g++ Example.cpp -o detectoreffect -std=c++11 \
		CATree.o DrawRandom.o DrawRandom2.o TauHelperFunctions3.o \
		`root-config --cflags --libs` -g
		
		
DetEf: detectoreffect
	./detectoreffect
