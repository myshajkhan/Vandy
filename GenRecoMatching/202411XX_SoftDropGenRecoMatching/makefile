
default:RRec

objects: CATree.o DrawRandom.o DrawRandom2.o TauHelperFunctions3.o

CATree.o: CATree.cpp CATree.h
	g++ CATree.cpp -c -o CATree.o -std=c++11

DrawRandom.o: DrawRandom.cpp DrawRandom.h
	g++ DrawRandom.cpp -c -o DrawRandom.o -std=c++11

DrawRandom2.o: DrawRandom2.cpp DrawRandom2.h
	g++ DrawRandom2.cpp -c -o DrawRandom2.o -std=c++11

TauHelperFunctions3.o: TauHelperFunctions3.cpp TauHelperFunctions3.h
	g++ TauHelperFunctions3.cpp -c -o TauHelperFunctions3.o -std=c++11

GroomingLateKT: Grooming_LatekT.cpp LateKTGrooming.h objects
	g++ Grooming_LatekT.cpp -o GroomingLateKT -std=c++11 \
		CATree.o DrawRandom.o DrawRandom2.o TauHelperFunctions3.o \
		`root-config --cflags --libs` -g

GroomingSoftDrop: Grooming_SD.cpp SoftDropGrooming.h objects
	g++ Grooming_SD.cpp -o GroomingSoftDrop -std=c++11 \
		CATree.o DrawRandom.o DrawRandom2.o TauHelperFunctions3.o \
		`root-config --cflags --libs` -g

GroomingDynamic: Grooming_Dynamic.cpp DynamicGrooming.h objects
	g++ Grooming_Dynamic.cpp -o GroomingDynamic -std=c++11 \
		CATree.o DrawRandom.o DrawRandom2.o TauHelperFunctions3.o \
		`root-config --cflags --libs` -g

RealReco: Realdata.cpp objects
	g++ Realdata.cpp -o RealReco -std=c++11 \
		CATree.o DrawRandom.o DrawRandom2.o TauHelperFunctions3.o \
		`root-config --cflags --libs` -g

detectoreffect: Matching.cpp objects
	g++ Matching.cpp -o detectoreffect -std=c++11 \
		CATree.o DrawRandom.o DrawRandom2.o TauHelperFunctions3.o \
		`root-config --cflags --libs` -g
		
SoftDrop: GroomingSoftDrop
	./GroomingSoftDrop

Dynamic: GroomingDynamic
	./GroomingDynamic

LateKT: GroomingLateKT
	./GroomingLateKT
				
DetEf: detectoreffect
	./detectoreffect

RRec: RealReco
	./RealReco
