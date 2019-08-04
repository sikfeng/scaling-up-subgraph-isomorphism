#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseIntVect.h>
#include <DataStructs/BitOps.h>
#include "CountFingerprint.h"
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/random.hpp>
#include <cstdint>
#include <RDGeneral/BoostEndInclude.h>
#include <limits.h>
#include <RDGeneral/hash/hash.hpp>
#include <RDGeneral/types.h>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <GraphMol/Fingerprints/FingerprintUtil.h>
#include <boost/array.hpp>

// caller owns the result, it must be deleted
std::vector<uint32_t> CountFingerprintMol(
  const RDKit::ROMol &mol, unsigned int minPath, unsigned int maxPath,
  unsigned int fpSize, unsigned int nBitsPerHash, bool useHs,
  double tgtDensity, unsigned int minSize, bool branchedPaths,
  bool useBondOrder, std::vector<std::uint32_t> *atomInvariants,
  const std::vector<std::uint32_t> *fromAtoms) {

  std::vector<std::uint32_t> lAtomInvariants;
  if (!atomInvariants) {
    RDKit::RDKitFPUtils::buildDefaultRDKitFingerprintAtomInvariants(mol,
                                                             lAtomInvariants);
    atomInvariants = &lAtomInvariants;
  }

  std::vector<uint32_t> res(fpSize, 0);

  // get all paths
  RDKit::INT_PATH_LIST_MAP allPaths;
  RDKit::RDKitFPUtils::enumerateAllPaths(mol, allPaths, fromAtoms, branchedPaths,
                                  useHs, minPath, maxPath);

  // identify query bonds
  std::vector<short> isQueryBond(mol.getNumBonds(), 0);
  std::vector<const RDKit::Bond *> bondCache;
  RDKit::RDKitFPUtils::identifyQueryBonds(mol, bondCache, isQueryBond);

  boost::dynamic_bitset<> atomsInPath(mol.getNumAtoms());
  for (RDKit::INT_PATH_LIST_MAP_CI paths = allPaths.begin(); paths != allPaths.end();
       paths++) {
    BOOST_FOREACH (const RDKit::PATH_TYPE &path, paths->second) {
      // the bond hashes of the path
      std::vector<unsigned int> bondHashes = RDKit::RDKitFPUtils::generateBondHashes(
          mol, atomsInPath, bondCache, isQueryBond, path, useBondOrder,
          atomInvariants);
      if (!bondHashes.size()) {
        continue;
      }

      // hash the path to generate a seed:
      unsigned long seed;
      if (path.size() > 1) {
        std::sort(bondHashes.begin(), bondHashes.end());

        // finally, we will add the number of distinct atoms in the path at the
        // end
        // of the vect. This allows us to distinguish C1CC1 from CC(C)C
        bondHashes.push_back(static_cast<unsigned int>(atomsInPath.count()));
        seed = gboost::hash_range(bondHashes.begin(), bondHashes.end());
      } else {
        seed = bondHashes[0];
      }

      unsigned int bit = seed % fpSize;

      ++res[bit];
    }
  }
  return res;
}
