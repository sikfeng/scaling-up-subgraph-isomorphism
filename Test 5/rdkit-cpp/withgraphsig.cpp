#include "GraphSig.h"

#include <iostream>
#include <cstdio>
#include <filesystem>
#include <string>
#include <fstream>
#include <vector>
#include <chrono>

#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/ExplicitBitVect.h>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>

int main(int argc ,char **argv) {
  std::chrono::high_resolution_clock::time_point startExtractFile = std::chrono::high_resolution_clock::now();
  std::string curr_dir = std::filesystem::current_path();
  std::string graphsFilename = curr_dir + "/sdf/candidates/" + argv[1] + "/0.sdf";
  std::string subgraphsFilename = curr_dir + "/sdf/substructures_all/0.sdf";
  RDKit::SDMolSupplier graphsFile(graphsFilename, false);
  RDKit::SDMolSupplier subgraphsFile(subgraphsFilename, false);
  std::vector<RDKit::ROMOL_SPTR> graphs;
  std::vector<RDKit::ROMOL_SPTR> subgraphs;
  while (!graphsFile.atEnd()) {
    RDKit::ROMOL_SPTR graph(graphsFile.next());
    graphs.push_back(graph);
  }
  while (!subgraphsFile.atEnd()) {
    RDKit::ROMOL_SPTR subgraph(subgraphsFile.next());
    subgraphs.push_back(subgraph);
  }
  std::chrono::high_resolution_clock::time_point endExtractFile = std::chrono::high_resolution_clock::now();
  double durationExtractFile = std::chrono::duration_cast<std::chrono::duration<double>>(endExtractFile - startExtractFile).count();
  printf("Extracted structures in %f seconds\n", durationExtractFile);
  printf("%d fragments and %d candidates\n", subgraphs.size(), graphs.size());

  std::chrono::high_resolution_clock::time_point startSigTime = std::chrono::high_resolution_clock::now();
  std::vector<GraphSig> graphSigs;
  std::vector<GraphSig> subgraphSigs;
  for (int i = 0; i < graphs.size(); ++i) {
    graphSigs.push_back(GraphSig(&*graphs[i]));
  }
  for (int j = 0; j < subgraphs.size(); ++j) {
    subgraphSigs.push_back(GraphSig(&*subgraphs[j]));
  }
  std::chrono::high_resolution_clock::time_point endSigTime = std::chrono::high_resolution_clock::now();
  double durationSig = std::chrono::duration_cast<std::chrono::duration<double>>(endSigTime - startSigTime).count();
  printf("Generated signatures in %f seconds\n", durationSig);
 
  std::chrono::high_resolution_clock::time_point startSGI = std::chrono::high_resolution_clock::now();
  int numSGIRun = 0;
  int numSubIso = 0;
  RDKit::MatchVectType matchVect;
  for (int i = 0; i < graphs.size(); ++i) {
    GraphSig graphSig = graphSigs[i];
    for (int j = 0; j < subgraphs.size(); ++j) {
      GraphSig subgraphSig = subgraphSigs[j];
      if (graphSig.compare(subgraphSig) == 0) {
        ++numSGIRun;
        RDKit::ROMOL_SPTR graph = graphs[i];
        RDKit::ROMOL_SPTR subgraph = subgraphs[j];
        if (RDKit::SubstructMatch(*graph, *subgraph, matchVect)) {
          ++numSubIso;
        }
      }
    }
  }
  std::chrono::high_resolution_clock::time_point endSGI = std::chrono::high_resolution_clock::now();
  double durationSGI = std::chrono::duration_cast<std::chrono::duration<double>>(endSGI - startSGI).count();
  printf("Finished SGI in %f seconds\n", durationSGI);
  printf("Found %d substructures, SGI run %d times\n", numSubIso, numSGIRun);
  return 0;
}

