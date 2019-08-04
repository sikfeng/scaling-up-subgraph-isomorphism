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

  std::chrono::high_resolution_clock::time_point startSGI = std::chrono::high_resolution_clock::now();
  int numSGIRun = 0;
  int numSubIso = 0;
  RDKit::MatchVectType matchVect;
  for (int i = 0; i < graphs.size(); ++i) {
    for (int j = 0; j < subgraphs.size(); ++j) {
      ++numSGIRun;
      RDKit::ROMOL_SPTR graph = graphs[i];
      RDKit::ROMOL_SPTR subgraph = subgraphs[j];
      if (RDKit::SubstructMatch(*graph, *subgraph, matchVect)) {
        ++numSubIso;
      }
    }
  }
  std::chrono::high_resolution_clock::time_point endSGI = std::chrono::high_resolution_clock::now();
  double durationSGI = std::chrono::duration_cast<std::chrono::duration<double>>(endSGI - startSGI).count();
  printf("Finished SGI in %f seconds\n", durationSGI);
  printf("Found %d substructures, SGI run %d times\n", numSubIso, numSGIRun);
  return 0;
}

