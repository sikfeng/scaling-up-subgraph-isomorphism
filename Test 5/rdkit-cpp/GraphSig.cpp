#include <iostream>
#include <vector>
#include <GraphMol/GraphMol.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseIntVect.h>
#include <DataStructs/BitOps.h>
#include <boost/array.hpp>

#include "GraphSig.h"
#include "CountFingerprint.cpp"

GraphSig::GraphSig(RDKit::ROMol *mol) {
  //molFingerprint = RDKit::RDKFingerprintMol(*mol, 1, 7, 1024, 1, true, 0.0, 128, false, true);
  molFingerprint = CountFingerprintMol(*mol, 1, 7, 512, 1, true, 0.0, 128, false, true);
  for (RDKit::ROMol::AtomIterator ai = mol->beginAtoms(); ai != mol->endAtoms(); ++ai) {
    int atomNum = (*ai)->getAtomicNum();
    for (int i = 0; i < 5; ++i) {
      if (countedAtoms[i] == atomNum) {
        ++countAtoms[i];
        break;
      }
    }  
  }
}

int GraphSig::compare(GraphSig subgraphSig) {
  for (int i = 0; i < 5; ++i) {
    if (this->countAtoms[i] < subgraphSig.countAtoms[i]) {
      return -1;
    }
  }
  for (int i = 0; i < this->molFingerprint.size(); ++i) {
    if (this->molFingerprint[i] < subgraphSig.molFingerprint[i]) {
      return -1;
    }
  }
  return 0;
}
