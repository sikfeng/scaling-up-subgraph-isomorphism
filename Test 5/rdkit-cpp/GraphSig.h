#ifndef GRAPHSIG_H
#define GRAPHSIG_H

#include <vector>
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseIntVect.h>
#include <DataStructs/BitOps.h>

#include "CountFingerprint.h"

class GraphSig {
  public:
    static constexpr int countedAtoms[5] = {6, 8, 14, 15, 18};
    int countAtoms[5] = {0, 0, 0, 0, 0};
    //ExplicitBitVect *molFingerprint;
    std::vector<uint32_t> molFingerprint;
    GraphSig(RDKit::ROMol *mol);
    int compare(GraphSig subgraphSig);
};

#endif
