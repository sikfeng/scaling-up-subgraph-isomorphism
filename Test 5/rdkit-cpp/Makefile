CXXFLAGS=-Iinclude/ -Iinclude/rdkit -std=c++17 -Ofast
LDFLAGS=-Llib -lRDKitFingerprints -lRDKitRDGeneral -lRDKitFileParsers -lRDKitSmilesParse -lRDKitSubstructMatch -lRDKitGraphMol -lRDKitDataStructs -lstdc++ -lstdc++fs

all: clean direct withgraphsig
clean:
	rm -f direct withgraphsig *.o
direct: direct.o
withgraphsig: withgraphsig.o GraphSig.o
